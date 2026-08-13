#!/usr/bin/env python
"""
------------------------------------------------------------------------------------------------------------------------
Copyright 2024 Roman Ufimov under the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later version.

cast_retrieve — collect assemblies and extract exonic contigs via BLAST.

  1. Collect per-gene SPAdes contigs from 20assemblies/ → 30raw_contigs/
     (samples from samples_list.txt only; parallel over samples).
  2. Prepare headers into 31exonic_contigs/{sample}.fasta.
  3. Index probe exons once (makeblastdb); per sample blastn contigs as query.
  4. Filter hits (locus match, k-mer cover, probe coverage), slice exons,
     write {sample}.fas, all_hits.tsv, statistics.tsv, retrieve_qc.tsv.
  5. Clean temporary BLAST DBs / prepared FASTAs / raw hit tables.
------------------------------------------------------------------------------------------------------------------------
"""

from __future__ import annotations

import logging
import os
import re
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from glob import glob
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import Bio.SeqRecord
import pandas
from Bio import SeqIO

from ParalogWizard import log_exceptions, managed_pool
from ParalogWizard.cast_call import (
    allocate_workers_and_threads,
    require_dir,
    require_file,
    require_tools,
)

logger = logging.getLogger("ParalogWizard")

REQUIRED_TOOLS = ("makeblastdb", "blastn")
BLAST_EVALUE = 1e-10
BLAST_MAX_TARGET_SEQS = 50
BLAST_TASK = "blastn"


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
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
    logger.debug("Sample list: %s", ", ".join(samples))
    return samples


def _sample_stem_from_raw_basename(basename: str) -> str:
    """
    '{sample}_contigs.fasta' → sample (do not split on underscores inside sample).
    """
    name = Path(basename).name
    if name.endswith("_contigs.fasta"):
        return name[: -len("_contigs.fasta")]
    if name.endswith(".fasta"):
        return Path(name).stem
    return Path(name).stem


@log_exceptions
def slicing(
    dictionary: Dict[str, Bio.SeqRecord.SeqRecord],
    contig_id: str,
    start: int,
    end: int,
) -> str:
    """Extract contig subsequence [start, end] (1-based inclusive, BLAST-style)."""
    sequence = dictionary[contig_id].seq
    logger.debug("Slicing contig %s: start=%s end=%s", contig_id, start, end)
    if start > end:
        return str(sequence.reverse_complement())[-end : -start - 1 : -1][::-1]
    return str(sequence)[start - 1 : end]


def ensure_probe_blast_db(probe_exons: str, out_dir: str) -> str:
    """
    Build a nucleotide BLAST DB for the probe FASTA once under out_dir.
    Returns the DB basename prefix (without extension).
    """
    require_file(probe_exons, "probe exons FASTA")
    os.makedirs(out_dir, exist_ok=True)
    db_prefix = os.path.join(out_dir, "probes_exons_blast")
    marker = db_prefix + ".nsq"
    if os.path.isfile(marker):
        for ext in (".nhr", ".nin", ".nsq"):
            try:
                os.remove(db_prefix + ext)
            except OSError:
                pass
        for path in glob(db_prefix + ".*"):
            try:
                os.remove(path)
            except OSError:
                pass

    cmd = [
        "makeblastdb",
        "-dbtype",
        "nucl",
        "-in",
        probe_exons,
        "-out",
        db_prefix,
        "-parse_seqids",
    ]
    logger.info("Indexing probe exons once: %s", probe_exons)
    logger.debug("makeblastdb: %s", " ".join(cmd))
    try:
        subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        logger.error("makeblastdb (probes) failed: %s", e.stderr)
        raise
    logger.info("Probe BLAST DB ready: %s", db_prefix)
    return db_prefix


def _probe_lengths(probe_exons: str) -> Dict[str, int]:
    lengths = {
        rec.id: len(rec.seq) for rec in SeqIO.parse(probe_exons, "fasta")
    }
    logger.debug("Loaded lengths for %d probe record(s)", len(lengths))
    return lengths


# -----------------------------------------------------------------------------
# Collect / prepare
# -----------------------------------------------------------------------------
def _collect_one_sample(sample: str, folder20: str, target_folder: str) -> str:
    """Collect one sample's gene contigs into 30raw_contigs/{sample}_contigs.fasta."""
    sample_path = os.path.join(folder20, sample)
    if not os.path.isdir(sample_path):
        raise FileNotFoundError(
            f"Sample '{sample}' listed in samples_list.txt but missing under "
            f"{folder20}"
        )
    out_file = os.path.join(target_folder, f"{sample}_contigs.fasta")
    n_loci = 0
    with open(out_file, "w") as contigs_out:
        for locus in sorted(os.listdir(sample_path)):
            locus_file = os.path.join(
                sample_path, locus, f"{locus}_contigs.fasta"
            )
            if not os.path.isfile(locus_file):
                continue
            with open(locus_file) as lf:
                contigs_out.write(lf.read().replace(">", f">{locus}_"))
            n_loci += 1
    if os.path.getsize(out_file) == 0:
        os.remove(out_file)
        raise RuntimeError(f"No contig files found for sample {sample} in {sample_path}")
    logger.info(
        "Collected %s: %d locus file(s) -> %s (%d bytes)",
        sample,
        n_loci,
        out_file,
        os.path.getsize(out_file),
    )
    return sample


@log_exceptions
def collect_contigs(data_folder: str, samples: List[str], num_workers: int) -> List[str]:
    """
    Parallel collect from 20assemblies/<sample>/ into 30raw_contigs/.
    Only samples in `samples` (from samples_list.txt).
    """
    target_folder = os.path.join(data_folder, "30raw_contigs")
    os.makedirs(target_folder, exist_ok=True)
    folder20 = require_dir(
        os.path.join(data_folder, "20assemblies"),
        "20assemblies directory",
    )
    n_workers = max(1, min(int(num_workers), len(samples)))
    logger.info(
        "Collecting %d sample(s) from %s with %d worker(s)",
        len(samples),
        folder20,
        n_workers,
    )
    done: List[str] = []
    failures: List[Tuple[str, Exception]] = []
    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        futures = {
            executor.submit(_collect_one_sample, s, folder20, target_folder): s
            for s in samples
        }
        for fut in as_completed(futures):
            sample = futures[fut]
            try:
                done.append(fut.result())
            except Exception as e:
                logger.error("Collect failed for %s: %s", sample, e)
                failures.append((sample, e))
    if failures:
        preview = ", ".join(s for s, _ in failures[:20])
        raise RuntimeError(
            f"cast_retrieve aborted: collect failed for {len(failures)} sample(s) "
            f"({preview}). See log."
        )
    logger.info("Contig collection complete: %d sample(s)", len(done))
    return sorted(done)


@log_exceptions
def prepare_contigs(raw_fasta: str, exonic_dir: str) -> str:
    """
    Copy raw contig FASTA into 31exonic_contigs/{sample}.fasta with shortened headers.
    Sample ID is taken from '{sample}_contigs.fasta' stem (not underscore-splitting).
    """
    require_file(raw_fasta, "raw contig FASTA")
    os.makedirs(exonic_dir, exist_ok=True)
    sample = _sample_stem_from_raw_basename(os.path.basename(raw_fasta))
    final_path = os.path.join(exonic_dir, f"{sample}.fasta")
    logger.info("Preparing contigs %s -> %s", raw_fasta, final_path)

    with open(raw_fasta) as fh:
        text = fh.read()
    text = re.sub("NODE", "N", text)
    text = re.sub(
        r"length_([0-9]+)_cov_([0-9]+\.[0-9]{2}).*",
        r"\1_c_\2",
        text,
    )
    with open(final_path, "w") as fh:
        fh.write(text)
    logger.debug("Prepared %s (%d bytes)", final_path, os.path.getsize(final_path))
    return final_path


# -----------------------------------------------------------------------------
# Per-sample fused BLAST + correction
# -----------------------------------------------------------------------------
@log_exceptions
def process_sample_retrieve(
    fasta_file: str,
    probe_db: str,
    probe_lengths: Dict[str, int],
    n_threads: int,
    length_cover: float,
    spades_cover: float,
) -> Tuple[pandas.DataFrame, dict]:
    """
    blastn (contigs vs probe DB) + filter + slice + rename for one sample.
    Returns (hits DataFrame, qc dict).
    """
    require_file(fasta_file, "sample contig FASTA")
    sample = Path(fasta_file).stem
    folder = os.path.dirname(fasta_file)
    blast_out = os.path.join(folder, f"reference_in_{sample}_contigs.txt")
    fas_out = os.path.join(folder, f"{sample}.fas")

    contigs = list(SeqIO.parse(fasta_file, "fasta"))
    n_contigs = len(contigs)
    logger.info(
        "Sample %s: BLAST %d contig(s) vs probes (%d thread(s))",
        sample,
        n_contigs,
        n_threads,
    )

    blastn_cmd = [
        "blastn",
        "-task",
        BLAST_TASK,
        "-query",
        fasta_file,
        "-db",
        probe_db,
        "-out",
        blast_out,
        "-evalue",
        str(BLAST_EVALUE),
        "-max_target_seqs",
        str(BLAST_MAX_TARGET_SEQS),
        "-outfmt",
        "6 qaccver saccver pident qcovhsp evalue bitscore qstart qend sstart send",
        "-num_threads",
        str(max(1, int(n_threads))),
    ]
    logger.debug("blastn: %s", " ".join(blastn_cmd))
    try:
        subprocess.run(blastn_cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        logger.error("blastn failed for %s: %s", sample, e.stderr)
        raise

    colnames = [
        "qaccver",
        "saccver",
        "pident",
        "qcovhsp",
        "evalue",
        "bitscore",
        "qstart",
        "qend",
        "sstart",
        "send",
    ]
    if os.path.getsize(blast_out) == 0:
        blast_results = pandas.DataFrame(columns=colnames)
    else:
        blast_results = pandas.read_csv(blast_out, sep="\t", names=colnames)

    n_blast = len(blast_results)
    n_locus_mismatch = 0
    n_low_cover = 0
    n_low_probe_cov = 0

    if blast_results.empty:
        logger.warning("No BLAST hits for %s", sample)
        open(fas_out, "w").close()
        qc = {
            "sample": sample,
            "n_contigs": n_contigs,
            "n_blast_hits": 0,
            "n_locus_mismatch": 0,
            "n_low_kmer_cover": 0,
            "n_low_probe_cover": 0,
            "n_hits_kept": 0,
            "n_exons": 0,
            "fas_bytes": 0,
        }
        return pandas.DataFrame(), qc

    # Contig = query, probe = subject
    blast_results["contig_locus"] = blast_results["qaccver"].str.split("_").str[0]
    blast_results["locus"] = (
        blast_results["saccver"].str.split("-").str[1].str.split("_").str[0]
    )
    blast_results["exon"] = blast_results["saccver"].str.split("-").str[1]
    blast_results["k-mer_cover"] = (
        blast_results["qaccver"].str.split("_").str[-1].astype(float)
    )
    blast_results["probe_len"] = blast_results["saccver"].map(probe_lengths)
    blast_results["probe_cov"] = (
        (blast_results["send"] - blast_results["sstart"]).abs() + 1
    ) / blast_results["probe_len"] * 100.0

    locus_ok = blast_results["locus"] == blast_results["contig_locus"]
    n_locus_mismatch = int((~locus_ok).sum())
    cover_ok = blast_results["k-mer_cover"] >= spades_cover
    n_low_cover = int((locus_ok & ~cover_ok).sum())
    probe_ok = blast_results["probe_cov"] >= length_cover
    n_low_probe_cov = int((locus_ok & cover_ok & ~probe_ok).sum())

    hits = blast_results[locus_ok & cover_ok & probe_ok].copy()
    hits_sorted = hits.sort_values(
        ["exon", "probe_cov", "pident", "evalue", "bitscore"],
        ascending=[True, False, False, True, False],
    )
    hits_dedup = hits_sorted.drop_duplicates(subset=["exon", "qaccver"]).reset_index(
        drop=True
    )
    hits_dedup["sample"] = sample

    contigs_dict = {rec.id: rec for rec in contigs}
    sequences: List[str] = []
    for row in hits_dedup.itertuples(index=False):
        sequences.append(
            slicing(contigs_dict, row.qaccver, int(row.qstart), int(row.qend))
        )

    hits_dedup = hits_dedup.copy()
    hits_dedup["sequence"] = sequences

    # Restore classic all_hits orientation: qaccver=probe, saccver=contig
    out_hits = hits_dedup.rename(
        columns={
            "qaccver": "saccver",
            "saccver": "qaccver",
        }
    )
    out_hits["sstart"] = hits_dedup[["qstart", "qend"]].min(axis=1).values
    out_hits["send"] = hits_dedup[["qstart", "qend"]].max(axis=1).values

    with open(fas_out, "w") as out:
        for counter, row in enumerate(out_hits.itertuples(index=False), start=1):
            key = f"{row.exon}_{'_'.join(str(row.saccver).split('_')[1:])}"
            out.write(f">Contig{counter}_{sample}-{key.replace('_', '-')}\n")
            out.write(f"{row.sequence}\n")

    keep_cols = [
        c
        for c in out_hits.columns
        if c
        in {
            "qaccver",
            "saccver",
            "pident",
            "qcovhsp",
            "evalue",
            "bitscore",
            "sstart",
            "send",
            "locus",
            "k-mer_cover",
            "exon",
            "sample",
            "sequence",
            "probe_cov",
        }
    ]
    out_hits = out_hits[keep_cols]

    qc = {
        "sample": sample,
        "n_contigs": n_contigs,
        "n_blast_hits": n_blast,
        "n_locus_mismatch": n_locus_mismatch,
        "n_low_kmer_cover": n_low_cover,
        "n_low_probe_cover": n_low_probe_cov,
        "n_hits_kept": len(out_hits),
        "n_exons": int(out_hits["exon"].nunique()) if len(out_hits) else 0,
        "fas_bytes": os.path.getsize(fas_out),
    }
    logger.info(
        "Sample %s done: %d/%d BLAST hits kept, %d exon(s), fas=%d bytes",
        sample,
        qc["n_hits_kept"],
        n_blast,
        qc["n_exons"],
        qc["fas_bytes"],
    )
    logger.debug("QC %s: %s", sample, qc)
    return out_hits, qc


@log_exceptions
def copies_stats(all_hits: pandas.DataFrame) -> pandas.DataFrame:
    """Max hits per exon (copy estimate) for each locus × sample."""
    statistics = pandas.DataFrame([], columns=["locus"])
    statistics.set_index("locus", drop=True, inplace=True)
    for (sample, locus), group_df in all_hits.groupby(["sample", "locus"]):
        count = group_df.groupby("exon")["sequence"].count().max()
        statistics.loc[locus, sample] = count
    logger.debug("Statistics table shape=%s", statistics.shape)
    return statistics


@log_exceptions
def clean(path: str) -> None:
    """Remove temporary BLAST DBs, prepared FASTAs, and raw hit tables."""
    logger.info("Removing temporary files in %s...", path)
    removed = 0
    for pattern in (
        "*.fasta",
        "*.n*",
        "reference_in*",
        "probes_exons_blast*",
    ):
        for file in glob(os.path.join(path, pattern)):
            try:
                os.remove(file)
                removed += 1
                logger.debug("Removed %s", file)
            except OSError as exc:
                logger.debug("Could not remove %s: %s", file, exc)
    logger.info("Temporary files removed (%d path(s))", removed)


# -----------------------------------------------------------------------------
# Top-level retrieve
# -----------------------------------------------------------------------------
@log_exceptions
def retrieve(
    data_folder: str,
    collect: bool,
    probe_exons: str,
    num_cores: int,
    length_cover: float,
    spades_cover: float,
    log_queue,
) -> None:
    """
    Run cast_retrieve end-to-end.

    :param collect: collect from 20assemblies into 30raw_contigs
    """
    logger.info(
        "Starting cast_retrieve (collect=%s, cores=%d, length_cover=%s, spades_cover=%s)",
        collect,
        num_cores,
        length_cover,
        spades_cover,
    )
    require_tools(REQUIRED_TOOLS)
    data_folder = os.path.abspath(data_folder)
    probe_exons = os.path.abspath(probe_exons)
    require_dir(data_folder, "data folder")
    require_file(probe_exons, "probe exons FASTA")
    samples = _read_samples_list(data_folder)

    raw_dir = os.path.join(data_folder, "30raw_contigs")
    if collect:
        collect_contigs(data_folder, samples, num_cores)
    elif not os.path.isdir(raw_dir):
        raise RuntimeError(
            "No raw contigs found. Run cast_retrieve with -c after cast_assemble."
        )
    require_dir(raw_dir, "30raw_contigs directory")

    # Only process listed samples; fail if a listed sample has no raw contigs.
    raw_files: List[str] = []
    for sample in samples:
        path = os.path.join(raw_dir, f"{sample}_contigs.fasta")
        if not os.path.isfile(path) or os.path.getsize(path) == 0:
            raise FileNotFoundError(
                f"Missing raw contigs for sample '{sample}': {path}. "
                f"Re-run with -c after cast_assemble."
            )
        raw_files.append(path)

    exonic_dir = os.path.join(data_folder, "31exonic_contigs")
    if os.path.isdir(exonic_dir):
        shutil.rmtree(exonic_dir)
        logger.debug("Removed previous %s", exonic_dir)
    os.makedirs(exonic_dir, exist_ok=True)

    logger.info("Preparing %d raw contig file(s)", len(raw_files))
    prepared_fastas = [prepare_contigs(f, exonic_dir) for f in raw_files]

    probe_db = ensure_probe_blast_db(probe_exons, exonic_dir)
    probe_lengths = _probe_lengths(probe_exons)
    if not probe_lengths:
        raise ValueError(f"No sequences in probe FASTA: {probe_exons}")

    n_workers, threads_per = allocate_workers_and_threads(
        num_cores, len(prepared_fastas)
    )
    logger.info(
        "Retrieve jobs: %d sample(s), %d worker(s) x %d blastn thread(s)",
        len(prepared_fastas),
        n_workers,
        threads_per,
    )

    failures: List[Tuple[str, Exception]] = []
    hit_frames: List[pandas.DataFrame] = []
    qc_rows: List[dict] = []
    with managed_pool(n_workers, log_queue) as pool:
        async_results = [
            (
                fasta,
                pool.apply_async(
                    process_sample_retrieve,
                    (
                        fasta,
                        probe_db,
                        probe_lengths,
                        threads_per,
                        length_cover,
                        spades_cover,
                    ),
                ),
            )
            for fasta in prepared_fastas
        ]
        total = len(async_results)
        for i, (fasta, async_result) in enumerate(async_results, start=1):
            sample = Path(fasta).stem
            try:
                hits_df, qc = async_result.get()
                if hits_df is not None and not hits_df.empty:
                    hit_frames.append(hits_df)
                qc_rows.append(qc)
                logger.debug("Sample done %d/%d: %s", i, total, sample)
                if i == total or i % 10 == 0:
                    logger.info("Retrieve progress: %d / %d", i, total)
            except Exception as e:
                logger.error("Sample %s failed: %s", sample, e)
                failures.append((sample, e))

    if failures:
        preview = ", ".join(s for s, _ in failures[:20])
        raise RuntimeError(
            f"cast_retrieve aborted: {len(failures)} sample(s) failed "
            f"({preview}). See log."
        )

    if not hit_frames:
        raise RuntimeError("No exon hits retained for any sample.")

    all_hits = pandas.concat(hit_frames, ignore_index=True)
    hits_path = os.path.join(exonic_dir, "all_hits.tsv")
    all_hits.to_csv(hits_path, sep="\t", index=False)
    logger.info("Wrote %s (%d rows)", hits_path, len(all_hits))

    stats = copies_stats(all_hits)
    stats_path = os.path.join(exonic_dir, "statistics.tsv")
    stats.to_csv(stats_path, sep="\t", na_rep="NaN")
    logger.info("Wrote %s", stats_path)

    qc_df = pandas.DataFrame(qc_rows)
    qc_path = os.path.join(exonic_dir, "retrieve_qc.tsv")
    qc_df.to_csv(qc_path, sep="\t", index=False)
    logger.info("Wrote %s", qc_path)

    clean(exonic_dir)
    logger.info("cast_retrieve completed under %s", exonic_dir)
