#!/usr/bin/env python
"""
------------------------------------------------------------------------------------------------------------------------
Copyright 2024 Roman Ufimov under the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later version.

cast_analyze — build per-exon alignments/trees and estimate paralog divergence.

  1. Read 31exonic_contigs/all_hits.tsv; clear/rebuild 40aln_orth_par/.
  2. Per exon (fused, parallel): write FASTA → MAFFT → FastTree → pairwise
     distances / KDE peak means.
  3. Write pairwise_distances.tsv and divergence mixture plots.

Each run rebuilds alignments, trees, and divergence outputs (may overwrite).
------------------------------------------------------------------------------------------------------------------------
"""

from __future__ import annotations

import itertools
import logging
import os
import random
import shutil
import subprocess
import tempfile
from collections import defaultdict
from shutil import which
from typing import Dict, List, Optional, Set, Tuple

import matplotlib

matplotlib.use("Agg")
import numpy as np
import pandas as pd
from Bio import SeqIO, SeqRecord
from matplotlib import pyplot as plt
from scipy.signal import argrelextrema
from sklearn.mixture import BayesianGaussianMixture
from sklearn.neighbors import KernelDensity

from ParalogWizard import log_exceptions, managed_pool
from ParalogWizard.cast_call import (
    allocate_workers_and_threads,
    require_dir,
    require_file,
    require_tools,
)

logger = logging.getLogger("ParalogWizard")

REQUIRED_TOOLS_BASE = ("mafft",)
_GAP = ord("-")
PLOT_DPI = 150
BGM_MAX_ITER = 2000
BGM_N_INIT = 3


def _resolve_fasttree(prefer_mp: bool = True) -> str:
    """Return a FastTree executable path/name (prefer OpenMP builds when available)."""
    names_mp = ("FastTreeMP", "fasttreemp", "VeryFastTree")
    names_serial = ("FastTree", "fasttree")
    search = list(names_mp) + list(names_serial) if prefer_mp else list(names_serial) + list(names_mp)

    for name in search:
        path = which(name)
        if path:
            logger.debug("Using FastTree binary: %s -> %s", name, path)
            return path

    here = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".local", "bin"))
    for name in search:
        candidate = os.path.join(here, name)
        if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            logger.debug("Using FastTree binary: %s", candidate)
            return candidate
    raise EnvironmentError(
        "Required executable not found on PATH: FastTree / FastTreeMP / fasttree"
    )


# -----------------------------------------------------------------------------
# Alignment and tree building
# -----------------------------------------------------------------------------
@log_exceptions
def strip_mafft_rev_prefix(alignment_text: str) -> str:
    """Remove MAFFT '_R_' header prefix added when a sequence was reverse-complemented."""
    out_lines = []
    for line in alignment_text.splitlines(keepends=True):
        if line.startswith(">_R_"):
            line = ">" + line[4:]
        out_lines.append(line)
    return "".join(out_lines)


def mafft_align(file: str, n_threads: int = 1) -> str:
    """
    Align a FASTA with MAFFT --auto --adjustdirectionaccurately.
    Writes <basename>.fasta.mafft next to the input. Strips MAFFT '_R_' marks.

    Runs in a private temp cwd so parallel MAFFT jobs do not collide on
    makedirectionlist scratch files (_direction, infile, …).
    """
    require_file(file, "FASTA to align")
    abs_file = os.path.abspath(file)
    out_file = f"{os.path.splitext(file)[0]}.fasta.mafft"
    threads = max(1, int(n_threads))

    def _run(adjust: bool) -> subprocess.CompletedProcess:
        cmd = ["mafft", "--quiet", "--thread", str(threads)]
        if adjust:
            cmd.append("--adjustdirectionaccurately")
        cmd.extend(["--auto", abs_file])
        logger.debug("MAFFT command: %s", " ".join(cmd))
        with tempfile.TemporaryDirectory(prefix="mafft_") as tmp:
            return subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True,
                cwd=tmp,
            )

    try:
        result = _run(adjust=True)
    except subprocess.CalledProcessError as e:
        logger.warning(
            "MAFFT --adjustdirectionaccurately failed for %s (%s); retrying without it",
            os.path.basename(file),
            (e.stderr or "").strip()[:200] or f"exit {e.returncode}",
        )
        try:
            result = _run(adjust=False)
        except subprocess.CalledProcessError as e2:
            logger.error("MAFFT alignment failed for %s: %s", file, e2.stderr)
            raise
    n_rev = result.stdout.count(">_R_")
    aligned = strip_mafft_rev_prefix(result.stdout)
    if n_rev:
        logger.info(
            "MAFFT reverse-complemented %d sequence(s) in %s",
            n_rev,
            os.path.basename(file),
        )
    with open(out_file, "w") as fh:
        fh.write(aligned)
    logger.debug("MAFFT completed: %s (%d bytes)", out_file, os.path.getsize(out_file))
    return out_file


@log_exceptions
def fast_tree(
    file: str,
    binary: Optional[str] = None,
    n_threads: int = 1,
) -> str:
    """Build a tree with FastTree from an alignment FASTA. Writes <file>.tre."""
    require_file(file, "alignment for FastTree")
    out_file = f"{file}.tre"
    exe = binary or _resolve_fasttree()
    cmd = [exe, "-nt", "-gtr", "-fastest", file]
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(max(1, int(n_threads)))
    logger.debug("FastTree command: %s (OMP_NUM_THREADS=%s)", " ".join(cmd), env["OMP_NUM_THREADS"])
    try:
        with open(out_file, "w") as outf:
            result = subprocess.run(
                cmd,
                stdout=outf,
                stderr=subprocess.PIPE,
                text=True,
                check=True,
                env=env,
            )
        if result.stderr:
            logger.debug("%s stderr (tail): %s", exe, result.stderr[-500:])
    except subprocess.CalledProcessError as e:
        logger.error("%s failed for %s: %s", exe, file, e.stderr)
        raise
    logger.debug("Tree written: %s (%d bytes)", out_file, os.path.getsize(out_file))
    return out_file


def _write_exon_fasta(
    out_path: str,
    exon_str: str,
    records: List[Tuple[str, str, str]],
) -> None:
    """Write unaligned exon FASTA from (saccver, sample, sequence) rows."""
    with open(out_path, "w") as f_out:
        for saccver, sample, sequence in records:
            try:
                contig_part = str(saccver).split("_N_")[1]
            except (IndexError, AttributeError):
                contig_part = "unknown"
            header = f"{exon_str}_N_{contig_part}_{sample}"
            f_out.write(f">{header}\n{sequence}\n")


# -----------------------------------------------------------------------------
# Divergence helpers
# -----------------------------------------------------------------------------
@log_exceptions
def percent_dissimilarity(seq1: str, seq2: str) -> Optional[float]:
    """
    Percent dissimilarity between two aligned sequences.
    Ignores shared gap columns; returns None if overlap is too short.
    """
    n = min(len(seq1), len(seq2))
    if n == 0:
        return None
    a = np.frombuffer(seq1[:n].lower().encode("ascii", "ignore"), dtype=np.uint8)
    b = np.frombuffer(seq2[:n].lower().encode("ascii", "ignore"), dtype=np.uint8)
    if a.size != b.size:
        m = min(a.size, b.size)
        a = a[:m]
        b = b[:m]
    keep = ~((a == _GAP) & (b == _GAP))
    a = a[keep]
    b = b[keep]
    n = int(a.size)
    if n == 0:
        return None

    either_gap = (a == _GAP) | (b == _GAP)
    interior = ~either_gap
    if not interior.any():
        return None
    left = int(np.argmax(interior))
    right = n - int(np.argmax(interior[::-1]))
    overlap = right - left
    if overlap < 100:
        return None
    denom_left = n - left
    if denom_left <= 0 or right <= 0:
        return None
    if (overlap / denom_left < 0.5) or (overlap / right < 0.5):
        return None

    aa = a[left:right]
    bb = b[left:right]
    both = (aa != _GAP) & (bb != _GAP)
    mismatches = int(np.count_nonzero(both & (aa != bb)))
    return (mismatches / overlap) * 100


@log_exceptions
def get_distance_matrix(file_to_process: str, blocklist: Set[str]):
    """
    Pairwise within-sample distances for one alignment.
    Returns (cluster_means, distance DataFrame, optional KDE plot data).
    """
    require_file(file_to_process, "alignment FASTA")
    current_matrix_to_plot: List[float] = []
    current_matrix_to_write: List[list] = []
    sum_list: List[float] = []
    logger.debug("Computing distance matrix for %s", file_to_process)

    with open(file_to_process) as fasta_file:
        sequences: Dict[str, SeqRecord.SeqRecord] = SeqIO.to_dict(
            SeqIO.parse(fasta_file, "fasta")
        )
    if not sequences:
        logger.warning("No sequences in %s", file_to_process)
        return sum_list, pd.DataFrame(columns=["seq1", "dist", "seq2"]), None

    by_sample: Dict[str, List[str]] = defaultdict(list)
    for full_name in sequences:
        sample_name = "_".join(full_name.split("_")[-2:])
        by_sample[sample_name].append(full_name)

    for samp, names in by_sample.items():
        if len(names) >= 2:
            for pair in itertools.combinations(names, 2):
                d = percent_dissimilarity(
                    str(sequences[pair[0]].seq), str(sequences[pair[1]].seq)
                )
                if d is None:
                    continue
                if samp not in blocklist:
                    current_matrix_to_plot.append(d)
                current_matrix_to_write.append([pair[0], d, pair[1]])
        else:
            for name in names:
                seq = str(sequences[name].seq)
                if len(seq) < 100:
                    continue
                current_matrix_to_write.append([name, np.nan, np.nan])

    to_plot = None
    if current_matrix_to_plot:
        sorted_vals = np.sort(np.asarray(current_matrix_to_plot, dtype=float))
        arr = sorted_vals.reshape(-1, 1)
        kde = KernelDensity(kernel="gaussian", bandwidth=1.5)
        kde.fit(arr)
        space = np.linspace(-1, float(np.max(arr)) + 1, 1000)
        e = np.exp(kde.score_samples(space.reshape(-1, 1)))
        mi = argrelextrema(e, np.less)[0]
        minimum = space[mi]
        if len(minimum) == 0:
            sum_list.append(float(np.mean(sorted_vals)))
        else:
            for i in range(len(minimum)):
                if i == 0:
                    idx = np.where((sorted_vals < minimum[i]) & (sorted_vals >= 0))[0]
                else:
                    idx = np.where(
                        (sorted_vals < minimum[i]) & (sorted_vals > minimum[i - 1])
                    )[0]
                if len(idx) == 0:
                    continue
                cluster = sorted_vals[min(idx) : max(idx) + 1]
                sum_list.append(float(np.mean(cluster)))
            idx_last = np.where(sorted_vals > minimum[-1])[0]
            if len(idx_last):
                last_cluster = sorted_vals[min(idx_last) : max(idx_last) + 1]
                sum_list.append(float(np.mean(last_cluster)))
        to_plot = [space, e, sorted_vals]
        logger.debug(
            "%s: %d pairwise distances, %d KDE peak mean(s)",
            os.path.basename(file_to_process),
            len(current_matrix_to_plot),
            len(sum_list),
        )

    matrix_df = pd.DataFrame(current_matrix_to_write, columns=["seq1", "dist", "seq2"])
    return sum_list, matrix_df, to_plot


@log_exceptions
def process_exon(
    exon_str: str,
    records: List[Tuple[str, str, str]],
    aln_folder: str,
    blocklist: Set[str],
    n_threads: int,
    fasttree_bin: str,
) -> Tuple[List[float], pd.DataFrame, Optional[list]]:
    """
    Fused per-exon job: write FASTA → MAFFT → FastTree → distance matrix.
    """
    fasta_path = os.path.join(aln_folder, f"{exon_str}.fasta")
    _write_exon_fasta(fasta_path, exon_str, records)
    aln_path = mafft_align(fasta_path, n_threads=n_threads)
    fast_tree(aln_path, binary=fasttree_bin, n_threads=n_threads)
    return get_distance_matrix(aln_path, blocklist)


@log_exceptions
def get_model(array: np.ndarray, num_comp: int) -> BayesianGaussianMixture:
    """Fit a BayesianGaussianMixture on divergence values."""
    model = BayesianGaussianMixture(
        n_components=num_comp,
        max_iter=BGM_MAX_ITER,
        n_init=BGM_N_INIT,
    )
    model.fit(array)
    logger.debug("BayesianGaussianMixture fitted with %d components", num_comp)
    return model


@log_exceptions
def plot_vertical_line(
    ax: plt.Axes, line_name: str, line_value: float, colors: List[str], idx: int
) -> None:
    """Draw a labeled vertical reference line on a divergence plot."""
    ax.axvline(
        x=line_value,
        label=f"{line_name} - {np.round(line_value, 2)}",
        color=colors[idx],
        lw=0.6,
        ls="--",
    )


@log_exceptions
def get_plot(path: str, name: str, matrix: np.ndarray, comp: int) -> None:
    """Fit a mixture model and save PNG/SVG divergence plots."""
    os.makedirs(path, exist_ok=True)
    fig, ax = plt.subplots(figsize=(15, 10))
    ax.set_xlabel("Divergence (%)")
    ax.set_ylabel("Frequency")
    fig.suptitle(name)
    max_val = np.round(np.max(matrix)) + 1
    ax.xaxis.set_ticks(np.arange(0, max_val, 1))
    model = get_model(matrix, comp)
    means = model.means_.flatten().tolist()
    sigmas = [np.sqrt(x) for x in model.covariances_.flatten().tolist()]
    mu_sigma = sorted(zip(means, sigmas), key=lambda x: x[0])
    colors = random.sample(
        [
            "#FF0000",
            "#008000",
            "#0000FF",
            "#FF00FF",
            "#800000",
            "#808000",
            "#008080",
            "#2196F3",
            "#4A148C",
            "#512DA8",
            "#2ECC71",
        ],
        k=comp * 3,
    )
    if mu_sigma:
        first_peak, first_sigma = mu_sigma[0]
        plot_vertical_line(
            ax, "first_peak_minus_sigma", first_peak - first_sigma, colors, 0
        )
        plot_vertical_line(ax, "first_peak", first_peak, colors, 1)
        plot_vertical_line(
            ax, "first_peak_plus_sigma", first_peak + first_sigma, colors, 2
        )
        if comp > 1 and len(mu_sigma) > 1:
            second_peak, second_sigma = mu_sigma[1]
            plot_vertical_line(
                ax, "second_peak_minus_sigma", second_peak - second_sigma, colors, 3
            )
            plot_vertical_line(ax, "second_peak", second_peak, colors, 4)
            plot_vertical_line(
                ax, "second_peak_plus_sigma", second_peak + second_sigma, colors, 5
            )
            if comp > 2 and len(mu_sigma) > 2:
                third_peak, third_sigma = mu_sigma[2]
                plot_vertical_line(
                    ax, "third_peak_minus_sigma", third_peak - third_sigma, colors, 6
                )
                plot_vertical_line(ax, "third_peak", third_peak, colors, 7)
                plot_vertical_line(
                    ax, "third_peak_plus_sigma", third_peak + third_sigma, colors, 8
                )
    ax.legend(loc="upper right")
    ax.hist(
        matrix,
        bins=np.arange(0, max_val, 1),
        density=True,
        histtype="stepfilled",
        alpha=0.4,
    )
    space = np.linspace(0, max_val, 1000)
    logprob = model.score_samples(space.reshape(-1, 1))
    pdf = np.exp(logprob)
    responsibilities = model.predict_proba(space.reshape(-1, 1))
    ax.plot(space, pdf, "-k")
    ax.plot(space, responsibilities * pdf[:, np.newaxis], "--k")
    png = os.path.join(path, f"{name}.png")
    svg = os.path.join(path, f"{name}.svg")
    fig.savefig(png, dpi=PLOT_DPI, format="png")
    fig.savefig(svg, format="svg")
    plt.close(fig)
    logger.info("Plot saved: %s / %s", png, svg)


# -----------------------------------------------------------------------------
# Pipeline stages
# -----------------------------------------------------------------------------
@log_exceptions
def build_and_score_exons(
    data_folder: str,
    n_cpu: int,
    blocklist: Set[str],
    log_queue,
) -> List[Tuple[List[float], pd.DataFrame, Optional[list]]]:
    """
    Clear 40aln_orth_par, then run fused per-exon jobs in parallel.
    Returns list of (cluster_means, distance_df, kde_plot_data).
    """
    hits_file = require_file(
        os.path.join(data_folder, "31exonic_contigs", "all_hits.tsv"),
        "all_hits.tsv",
    )
    all_hits = pd.read_csv(hits_file, sep="\t")
    if all_hits.empty:
        raise ValueError(f"Hits table is empty: {hits_file}")
    for col in ("exon", "saccver", "sample", "sequence"):
        if col not in all_hits.columns:
            raise ValueError(f"Required column '{col}' missing from {hits_file}")

    aln_folder = os.path.join(data_folder, "40aln_orth_par")
    if os.path.isdir(aln_folder):
        shutil.rmtree(aln_folder)
        logger.debug("Removed previous %s", aln_folder)
    os.makedirs(aln_folder, exist_ok=True)

    jobs: List[Tuple[str, List[Tuple[str, str, str]]]] = []
    for exon, df in all_hits.groupby("exon", sort=True):
        exon_str = str(exon)
        records = list(
            zip(
                df["saccver"].astype(str).tolist(),
                df["sample"].astype(str).tolist(),
                df["sequence"].astype(str).tolist(),
            )
        )
        jobs.append((exon_str, records))

    if not jobs:
        raise RuntimeError(f"No exons found in {hits_file}")

    fasttree_bin = _resolve_fasttree(prefer_mp=True)
    n_workers, threads_per = allocate_workers_and_threads(n_cpu, len(jobs))
    logger.info(
        "Exon jobs: %d exon(s), %d worker(s) x %d thread(s) (FastTree=%s)",
        len(jobs),
        n_workers,
        threads_per,
        fasttree_bin,
    )

    failures: List[Tuple[str, Exception]] = []
    results: List[Optional[Tuple[List[float], pd.DataFrame, Optional[list]]]] = [
        None
    ] * len(jobs)
    with managed_pool(n_workers, log_queue) as pool:
        async_results = [
            (
                i,
                exon_str,
                pool.apply_async(
                    process_exon,
                    (
                        exon_str,
                        records,
                        aln_folder,
                        blocklist,
                        threads_per,
                        fasttree_bin,
                    ),
                ),
            )
            for i, (exon_str, records) in enumerate(jobs)
        ]
        total = len(async_results)
        for done, (i, exon_str, async_result) in enumerate(async_results, start=1):
            try:
                results[i] = async_result.get()
                logger.debug("Exon done %d/%d: %s", done, total, exon_str)
                if done == total or done % 25 == 0:
                    logger.info("Exon progress: %d / %d", done, total)
            except Exception as e:
                logger.error("Exon %s failed: %s", exon_str, e)
                failures.append((exon_str, e))

    if failures:
        preview = ", ".join(e for e, _ in failures[:20])
        raise RuntimeError(
            f"cast_analyze aborted: {len(failures)} exon(s) failed ({preview}). See log."
        )

    logger.info(
        "Alignment, trees, and distances complete for %d exon(s) under %s",
        len(jobs),
        aln_folder,
    )
    return [r for r in results if r is not None]


@log_exceptions
def estimate_divergence_from_results(
    data_folder: str,
    results: List[Tuple[List[float], pd.DataFrame, Optional[list]]],
) -> None:
    """Aggregate fused-job distance results into tables and plots."""
    aln_folder = require_dir(
        os.path.join(data_folder, "40aln_orth_par"),
        "40aln_orth_par directory",
    )
    if not results:
        raise RuntimeError("No exon results to aggregate.")

    distance_frames = [r[1] for r in results if r[1] is not None and not r[1].empty]
    if distance_frames:
        distances_df = pd.concat(distance_frames, ignore_index=True)
    else:
        distances_df = pd.DataFrame(columns=["seq1", "dist", "seq2"])
    dist_path = os.path.join(aln_folder, "pairwise_distances.tsv")
    distances_df.to_csv(dist_path, sep="\t", index=False)
    logger.info("Wrote %s (%d rows)", dist_path, len(distances_df))

    fig, ax = plt.subplots(figsize=(15, 10))
    plotted = 0
    for r in results:
        if r[2] is None:
            continue
        space, e, sorted_arr = r[2]
        ax.plot(space, e, "k-", alpha=0.05)
        ax.plot(
            sorted_arr, np.zeros(sorted_arr.shape[0]), marker=2, color="k", alpha=0.5
        )
        plotted += 1
    logger.debug("Overlayed %d individual KDE curve(s)", plotted)
    end = np.round(ax.get_xlim()[1], 0)
    ax.xaxis.set_ticks(np.arange(0, end, 1))
    ax.set_xlabel("Divergence (%)")
    ax.set_ylabel("Frequency")
    ind_png = os.path.join(aln_folder, "individual_distributions.png")
    ind_svg = os.path.join(aln_folder, "individual_distributions.svg")
    fig.savefig(ind_png, dpi=PLOT_DPI, format="png")
    fig.savefig(ind_svg, format="svg")
    plt.close(fig)
    logger.info("Wrote %s / %s", ind_png, ind_svg)

    divergence_array = np.array(
        [
            [x]
            for sublist in [r[0] if isinstance(r[0], list) else [r[0]] for r in results]
            for x in sublist
        ]
    )
    if divergence_array.size == 0:
        raise RuntimeError(
            "No divergence values available for mixture fitting. "
            "Check alignments and blocklist."
        )
    logger.info(
        "Fitting divergence mixture models (1–3 components) on %d peak mean(s)",
        len(divergence_array),
    )
    for comp in [1, 2, 3]:
        get_plot(
            aln_folder,
            f"pairwise_distances_distribution_{comp}_comp",
            divergence_array,
            comp,
        )
    logger.info("Divergence estimation complete under %s", aln_folder)


# -----------------------------------------------------------------------------
# Top-level entry
# -----------------------------------------------------------------------------
@log_exceptions
def analyze(
    data_folder: str,
    num_cores: int,
    log_queue,
    blocklist: Optional[Set[str]] = None,
) -> None:
    """Run cast_analyze: fused alignments/trees/distances then divergence plots."""
    blocklist = blocklist or set()
    logger.info(
        "Starting cast_analyze (cores=%d, blocklist=%d)",
        num_cores,
        len(blocklist),
    )
    logger.debug(
        "analyze args: data_folder=%s num_cores=%d blocklist=%s",
        data_folder,
        num_cores,
        ", ".join(sorted(blocklist)) if blocklist else "(none)",
    )
    require_tools(REQUIRED_TOOLS_BASE)
    _resolve_fasttree()  # fail early if missing
    data_folder = os.path.abspath(data_folder)
    require_dir(data_folder, "data folder")
    require_dir(
        os.path.join(data_folder, "31exonic_contigs"),
        "31exonic_contigs directory",
    )

    results = build_and_score_exons(data_folder, num_cores, blocklist, log_queue)
    estimate_divergence_from_results(data_folder, results)
    logger.info("cast_analyze completed under %s", data_folder)
