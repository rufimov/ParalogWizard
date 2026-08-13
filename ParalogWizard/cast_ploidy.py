#!/usr/bin/env python
"""
------------------------------------------------------------------------------------------------------------------------
Copyright 2024 Roman Ufimov under the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later version.

cast_ploidy — estimate sample ploidy with nQuire from remapped BAMs.

  1. Clear/rebuild 102ploidy/.
  2. Per sample (parallel): merge 100remapped/exons*/{sample}_filtered_uniq_sorted.bam
     → 102ploidy/{sample}_filtered_uniq_sorted.bam (+ index).
  3. nQuire create + denoise → {sample}.bin / {sample}_denoised.bin.
  4. nQuire lrdmodel on all denoised bins → 102ploidy/lrdmodel.tsv with a
     assigned ploidy column (dip / tri / tet).

Requires cast_remap outputs under 100remapped/. Each run rebuilds 102ploidy/.
------------------------------------------------------------------------------------------------------------------------
"""

from __future__ import annotations

import logging
import os
import shutil
import subprocess
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd

from ParalogWizard import log_exceptions, managed_pool
from ParalogWizard.cast_call import require_dir, require_file, require_tools

logger = logging.getLogger("ParalogWizard")

REQUIRED_TOOLS = ("samtools", "nQuire")
PLOIDY_DIRNAME = "102ploidy"
REMAPPED_DIRNAME = "100remapped"
BAM_SUFFIX = "_filtered_uniq_sorted.bam"


def _read_samples_list(data_folder: str) -> List[str]:
    path = require_file(
        os.path.join(data_folder, "10deduplicated_reads", "samples_list.txt"),
        "samples_list.txt",
    )
    with open(path) as fh:
        samples = [ln.strip() for ln in fh if ln.strip()]
    if not samples:
        raise ValueError(f"No samples in {path}")
    duplicates = sorted({s for s in samples if samples.count(s) > 1})
    if duplicates:
        raise ValueError(f"Duplicate sample(s) in {path}: {', '.join(duplicates)}")
    logger.info("Loaded %d sample(s) from %s", len(samples), path)
    return samples


def _exon_chunks(remapped_dir: str) -> List[str]:
    chunks = sorted(
        name
        for name in os.listdir(remapped_dir)
        if name.startswith("exons")
        and os.path.isdir(os.path.join(remapped_dir, name))
    )
    if not chunks:
        raise FileNotFoundError(
            f"No exons* chunk directories under {remapped_dir}. Run cast_remap first."
        )
    return chunks


def _run(cmd: List[str], *, cwd: Optional[str] = None, log_path: Optional[str] = None) -> None:
    logger.debug("Running: %s", " ".join(cmd))
    try:
        result = subprocess.run(
            cmd,
            cwd=cwd,
            capture_output=True,
            text=True,
            check=True,
        )
    except subprocess.CalledProcessError as e:
        if log_path:
            with open(log_path, "a") as fh:
                fh.write("$ " + " ".join(cmd) + "\n")
                fh.write(e.stdout or "")
                fh.write(e.stderr or "")
                fh.write("\n")
        logger.error(
            "Command failed (%s): %s\n%s",
            e.returncode,
            " ".join(cmd),
            (e.stderr or e.stdout or "")[-2000:],
        )
        raise
    if log_path and (result.stdout or result.stderr):
        with open(log_path, "a") as fh:
            fh.write("$ " + " ".join(cmd) + "\n")
            if result.stdout:
                fh.write(result.stdout)
            if result.stderr:
                fh.write(result.stderr)
            fh.write("\n")


@log_exceptions
def process_sample_ploidy(
    sample: str,
    data_folder: str,
    chunk_names: List[str],
) -> str:
    """
    Merge chunk BAMs for one sample, index, run nQuire create + denoise.
    Returns path to the denoised .bin file.
    """
    ploidy_dir = os.path.join(data_folder, PLOIDY_DIRNAME)
    remapped_dir = os.path.join(data_folder, REMAPPED_DIRNAME)
    chunk_bams: List[str] = []
    for chunk in chunk_names:
        bam = os.path.join(remapped_dir, chunk, f"{sample}{BAM_SUFFIX}")
        require_file(bam, f"remapped BAM for {sample}/{chunk}")
        chunk_bams.append(bam)

    merged_bam = os.path.join(ploidy_dir, f"{sample}{BAM_SUFFIX}")
    merged_bai = merged_bam + ".bai"
    sample_log = os.path.join(ploidy_dir, f"{sample}_ploidy.log")
    bin_prefix = os.path.join(ploidy_dir, sample)
    raw_bin = f"{bin_prefix}.bin"
    denoised_prefix = os.path.join(ploidy_dir, f"{sample}_denoised")
    denoised_bin = f"{denoised_prefix}.bin"

    logger.info(
        "Sample %s: merging %d chunk BAM(s) → nQuire create/denoise",
        sample,
        len(chunk_bams),
    )
    merge_cmd = ["samtools", "merge", "-f", "-o", merged_bam, *chunk_bams]
    _run(merge_cmd, log_path=sample_log)
    require_file(merged_bam, f"merged BAM for {sample}")
    _run(["samtools", "index", merged_bam], log_path=sample_log)
    require_file(merged_bai, f"BAM index for {sample}")

    _run(
        ["nQuire", "create", "-b", merged_bam, "-o", bin_prefix],
        log_path=sample_log,
    )
    require_file(raw_bin, f"nQuire .bin for {sample}")
    _run(
        ["nQuire", "denoise", "-o", denoised_prefix, raw_bin],
        log_path=sample_log,
    )
    require_file(denoised_bin, f"nQuire denoised .bin for {sample}")
    logger.debug("Sample %s denoised bin ready: %s", sample, denoised_bin)
    return denoised_bin


def _assign_ploidy(lrdmodel: pd.DataFrame) -> pd.Series:
    """
    Assign dip/tri/tet from nQuire lrdmodel output.

    Prefer delta columns (d_dip, d_tri, d_tet): smallest ΔlogL = best support
    (nQuire documentation). Fall back to max of dip/tri/tet if deltas absent.
    """
    delta_cols = ["d_dip", "d_tri", "d_tet"]
    if all(c in lrdmodel.columns for c in delta_cols):
        deltas = lrdmodel[delta_cols].to_numpy(dtype=float)
        argmin = np.nanargmin(deltas, axis=1)
        labels = np.array(["dip", "tri", "tet"])
        assigned = labels[argmin]
        all_nan = np.isnan(deltas).all(axis=1)
        assigned = np.where(all_nan, "unknown", assigned)
        logger.info("Ploidy assigned using minimum ΔlogL (d_dip/d_tri/d_tet)")
        return pd.Series(assigned, index=lrdmodel.index)

    score_cols = ["dip", "tri", "tet"]
    missing = [c for c in score_cols if c not in lrdmodel.columns]
    if missing:
        raise ValueError(
            f"lrdmodel.tsv missing expected column(s): {', '.join(missing)}. "
            f"Found: {list(lrdmodel.columns)}"
        )
    scores = lrdmodel[score_cols].to_numpy(dtype=float)
    argmax = np.nanargmax(scores, axis=1)
    labels = np.array(["dip", "tri", "tet"])
    logger.warning(
        "lrdmodel.tsv has no d_dip/d_tri/d_tet columns; "
        "falling back to max of dip/tri/tet likelihoods"
    )
    return pd.Series(labels[argmax], index=lrdmodel.index)


@log_exceptions
def ploidy(data_folder: str, num_cores: int, log_queue) -> None:
    """Run cast_ploidy end-to-end into 102ploidy/."""
    logger.info("Starting cast_ploidy (cores=%d)", num_cores)
    require_tools(REQUIRED_TOOLS)
    data_folder = os.path.abspath(data_folder)
    require_dir(data_folder, "data folder")
    remapped_dir = require_dir(
        os.path.join(data_folder, REMAPPED_DIRNAME),
        "100remapped directory",
    )
    samples = _read_samples_list(data_folder)
    chunks = _exon_chunks(remapped_dir)
    logger.info("Using %d remap chunk(s): %s", len(chunks), ", ".join(chunks))

    ploidy_dir = os.path.join(data_folder, PLOIDY_DIRNAME)
    if os.path.isdir(ploidy_dir):
        shutil.rmtree(ploidy_dir)
        logger.debug("Removed previous %s", ploidy_dir)
    os.makedirs(ploidy_dir, exist_ok=True)

    n_workers = max(1, min(int(num_cores), len(samples)))
    logger.info(
        "Ploidy jobs: %d sample(s), %d worker(s)", len(samples), n_workers
    )

    failures: List[Tuple[str, Exception]] = []
    denoised_bins: List[Optional[str]] = [None] * len(samples)
    with managed_pool(n_workers, log_queue) as pool:
        async_results = [
            (
                i,
                sample,
                pool.apply_async(
                    process_sample_ploidy,
                    (sample, data_folder, chunks),
                ),
            )
            for i, sample in enumerate(samples)
        ]
        total = len(async_results)
        for done, (i, sample, async_result) in enumerate(async_results, start=1):
            try:
                denoised_bins[i] = async_result.get()
                logger.debug("Sample done %d/%d: %s", done, total, sample)
                if done == total or done % 10 == 0:
                    logger.info("Ploidy progress: %d / %d", done, total)
            except Exception as e:
                logger.error("Ploidy failed for %s: %s", sample, e)
                failures.append((sample, e))

    if failures:
        preview = ", ".join(s for s, _ in failures[:20])
        raise RuntimeError(
            f"cast_ploidy aborted: {len(failures)} sample(s) failed ({preview}). See log."
        )

    bin_paths = [p for p in denoised_bins if p]
    if len(bin_paths) != len(samples):
        raise RuntimeError("Internal error: missing denoised bin path(s).")

    lrd_path = os.path.join(ploidy_dir, "lrdmodel.tsv")
    lrd_log = os.path.join(ploidy_dir, "lrdmodel.log")
    lrd_cmd = ["nQuire", "lrdmodel", "-t", str(max(1, int(num_cores))), *bin_paths]
    logger.info("Running nQuire lrdmodel on %d denoised bin(s)", len(bin_paths))
    try:
        result = subprocess.run(
            lrd_cmd, capture_output=True, text=True, check=True
        )
    except subprocess.CalledProcessError as e:
        with open(lrd_log, "w") as fh:
            fh.write(e.stdout or "")
            fh.write(e.stderr or "")
        logger.error("nQuire lrdmodel failed: %s", e.stderr)
        raise
    with open(lrd_path, "w") as fh:
        fh.write(result.stdout)
    if result.stderr:
        with open(lrd_log, "w") as fh:
            fh.write(result.stderr)
    require_file(lrd_path, "lrdmodel.tsv")

    lrdmodel = pd.read_csv(lrd_path, sep=r"\s+")
    if len(lrdmodel) != len(samples):
        raise RuntimeError(
            f"lrdmodel has {len(lrdmodel)} row(s) but {len(samples)} sample(s). "
            f"See {lrd_path}"
        )
    lrdmodel["sample"] = samples
    lrdmodel["ploidy"] = _assign_ploidy(lrdmodel)
    cols = ["sample", "ploidy"] + [
        c for c in lrdmodel.columns if c not in {"sample", "ploidy"}
    ]
    lrdmodel = lrdmodel[cols]
    lrdmodel.to_csv(lrd_path, sep="\t", index=False)
    counts = lrdmodel["ploidy"].value_counts().to_dict()
    logger.info("Wrote %s (ploidy counts: %s)", lrd_path, counts)
    logger.info("cast_ploidy completed under %s", ploidy_dir)
