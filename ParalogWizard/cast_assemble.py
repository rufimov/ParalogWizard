#!/usr/bin/env python
"""
------------------------------------------------------------------------------------------------------------------------
Copyright 2024 Roman Ufimov under the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later version.

cast_assemble — map reads to baits and assemble gene contigs under 20assemblies/.

  1. Index the bait FASTA once (next to the bait file) and reuse for all samples.
  2. For each sample, map paired reads (plain or .gz) with BWA mem | samtools.
  3. Distribute mapped reads into per-gene interleaved FASTAs.
  4. Run SPAdes per gene (pool), with a reduced-k redo for failures.
  5. Contigs: 20assemblies/<sample>/<gene>/<gene>_contigs.fasta
     Status lists: 20assemblies/<sample>/spades_*.txt

Samples are processed in parallel when -nc allows; each sample gets a share of
cores for BWA threads and its SPAdes gene pool. Each run processes samples
fresh (may overwrite existing outputs).
------------------------------------------------------------------------------------------------------------------------
"""

from __future__ import annotations

import gzip
import logging
import multiprocessing
import os
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from contextlib import contextmanager
from glob import glob
from typing import Dict, List, Optional, TextIO, Tuple

from ParalogWizard import log_exceptions, managed_pool
from ParalogWizard.cast_call import (
    allocate_workers_and_threads,
    require_dir,
    require_file,
    require_tools,
)

logger = logging.getLogger("ParalogWizard")

REQUIRED_TOOLS = ("bwa", "samtools", "spades.py")

try:
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
except ImportError:
    logging.error(
        "BioPython not found. Please install BioPython: pip install biopython"
    )
    raise


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
@contextmanager
def change_dir(dest: str):
    old = os.getcwd()
    os.chdir(dest)
    try:
        yield
    finally:
        os.chdir(old)


def _open_text(path: str) -> TextIO:
    """Open a text file, transparently handling .gz."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def _strip_read_id(x: str) -> str:
    x = x.split()[0]
    return x[:-2] if x.endswith(("/1", "/2")) else x


def ensure_bwa_index(baitfile: str) -> str:
    """
    Ensure a BWA index exists next to baitfile; build it once if missing.
    Returns the absolute bait path for subsequent bwa mem calls.
    """
    baitfile = os.path.abspath(baitfile)
    require_file(baitfile, "bait FASTA")
    index_file = baitfile + ".amb"
    if os.path.isfile(index_file):
        logger.info("Reusing existing BWA index for %s", baitfile)
        logger.debug("Index marker: %s", index_file)
        return baitfile

    cmd = ["bwa", "index", baitfile]
    logger.info("Creating BWA index once for %s", baitfile)
    logger.debug("BWA index command: %s", " ".join(cmd))
    try:
        subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        logger.error("BWA index creation failed for %s: %s", baitfile, e.stderr)
        raise
    logger.info("BWA index created for %s", baitfile)
    return baitfile


def resolve_read_pair(reads_dir: str, sample: str) -> Tuple[str, str]:
    """
    Return absolute (R1, R2) paths, preferring .fastq.gz over plain .fastq.
    """
    candidates = sorted(glob(os.path.join(reads_dir, f"{sample}.*.fastq*")))
    by_name = {os.path.basename(p): os.path.abspath(p) for p in candidates}

    def _pick(mate: str) -> Optional[str]:
        gz = f"{sample}.{mate}.fastq.gz"
        plain = f"{sample}.{mate}.fastq"
        if gz in by_name:
            return by_name[gz]
        if plain in by_name:
            return by_name[plain]
        # Fallback: any file containing .R1. / .R2.
        for name, path in sorted(by_name.items()):
            if f".{mate}." in name:
                return path
        return None

    r1, r2 = _pick("R1"), _pick("R2")
    if r1 is None or r2 is None:
        raise FileNotFoundError(
            f"Could not resolve R1/R2 for {sample} in {reads_dir} "
            f"(found: {', '.join(sorted(by_name)) or 'none'})"
        )
    require_file(r1, f"{sample} R1")
    require_file(r2, f"{sample} R2")
    logger.debug("Resolved %s reads: R1=%s R2=%s", sample, r1, r2)
    return r1, r2


# -----------------------------------------------------------------------------
# Read distribution
# -----------------------------------------------------------------------------
@log_exceptions
def distribute_reads_bwa(
    bamfilename: str, readfiles: list, work_dir: str
) -> None:
    """
    Split BAM-mapped reads into per-target interleaved FASTAs under work_dir.

    Keeps one open write handle per gene (avoids open/close per pair).
    Read files may be plain FASTQ or .gz. BAM hits are streamed from samtools.
    Thread-safe: does not change process CWD.
    """
    require_file(bamfilename, "BAM file")
    for readfile in readfiles:
        require_file(readfile, "read file")
    os.makedirs(work_dir, exist_ok=True)

    logger.info("Parsing BAM file %s to extract read mappings", bamfilename)
    samtools_cmd = ["samtools", "view", "-F", "4", bamfilename]
    logger.debug("Streaming: %s", " ".join(samtools_cmd))

    read_hit_dict: Dict[str, List[str]] = {}
    proc = subprocess.Popen(
        samtools_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    assert proc.stdout is not None
    try:
        for line in proc.stdout:
            fields = line.split()
            if len(fields) < 3:
                logger.warning("Unexpected samtools output line: %s", line.rstrip())
                continue
            read_id = fields[0]
            target = fields[2].split("-")[-1]
            read_hit_dict.setdefault(read_id, []).append(target)
        stderr = proc.stderr.read() if proc.stderr else ""
        rc = proc.wait()
    finally:
        if proc.poll() is None:
            proc.kill()

    if rc != 0:
        logger.error("samtools failed for %s: %s", bamfilename, stderr)
        raise subprocess.CalledProcessError(rc, samtools_cmd, stderr)

    logger.info("Found %d unique reads with hits", len(read_hit_dict))
    logger.debug(
        "Targets hit: %d unique gene id(s)",
        len({t for targets in read_hit_dict.values() for t in targets}),
    )

    logger.info("Processing paired read files: %s, %s", readfiles[0], readfiles[1])
    handles: Dict[str, TextIO] = {}
    reads_written = 0
    try:
        with _open_text(readfiles[0]) as f1, _open_text(readfiles[1]) as f2:
            iterator1 = FastqGeneralIterator(f1)
            iterator2 = FastqGeneralIterator(f2)

            for id1_long, seq1, _ in iterator1:
                try:
                    id2_long, seq2, _ = next(iterator2)
                except StopIteration:
                    logger.error(
                        "Paired read file %s ended prematurely", readfiles[1]
                    )
                    raise ValueError("Paired read file ended prematurely")

                id1 = _strip_read_id(id1_long)
                id2 = _strip_read_id(id2_long)
                chosen = read_hit_dict.get(id1) or read_hit_dict.get(id2)
                if not chosen:
                    continue

                record = f">{id1}\n{seq1}\n>{id2}\n{seq2}\n"
                for target in chosen:
                    out = handles.get(target)
                    if out is None:
                        gene_dir = os.path.join(work_dir, target)
                        os.makedirs(gene_dir, exist_ok=True)
                        out_path = os.path.join(
                            gene_dir, f"{target}_interleaved.fasta"
                        )
                        out = open(out_path, "w")
                        handles[target] = out
                        logger.debug("Opened write handle for gene %s", target)
                    out.write(record)

                reads_written += 1
                if reads_written % 5000 == 0:
                    logger.debug("Processed %d mapped read pairs", reads_written)
    except Exception as e:
        logger.error(
            "Error processing read files %s, %s: %s",
            readfiles[0],
            readfiles[1],
            e,
        )
        raise
    finally:
        for target, fh in handles.items():
            try:
                fh.close()
            except Exception as exc:
                logger.debug("Error closing handle for %s: %s", target, exc)

    logger.info(
        "Read distribution completed. Wrote %d mapped pairs into %d gene(s)",
        reads_written,
        len(handles),
    )


# -----------------------------------------------------------------------------
# BWA
# -----------------------------------------------------------------------------
@log_exceptions
def bwa(
    readfiles: List[str],
    baitfile: str,
    basename: str,
    cpu: int,
    output_dir: str,
) -> str:
    """
    Map paired reads (plain or .gz) to an already-indexed baitfile.
    Writes {output_dir}/{basename}.bam via bwa mem | samtools view.
    Thread-safe: uses absolute paths only (no chdir).
    """
    require_file(baitfile, "indexed bait FASTA")
    if not os.path.isfile(baitfile + ".amb"):
        raise FileNotFoundError(
            f"BWA index missing for {baitfile}; call ensure_bwa_index() first"
        )
    for readfile in readfiles:
        require_file(readfile, "read file")
    os.makedirs(output_dir, exist_ok=True)

    if not cpu:
        cpu = multiprocessing.cpu_count()

    output_bam = os.path.join(output_dir, f"{basename}.bam")
    map_log = os.path.join(output_dir, f"{basename}_bwa_map.log")
    logger.info(
        "Running BWA mem with %d thread(s) for sample %s",
        cpu,
        basename,
    )

    bwa_mem_cmd = ["bwa", "mem", "-t", str(cpu), baitfile] + list(readfiles)
    samtools_cmd = ["samtools", "view", "-h", "-b", "-S"]
    logger.debug("BWA mem: %s", " ".join(bwa_mem_cmd))
    logger.debug("samtools: %s", " ".join(samtools_cmd))
    logger.debug("BAM output: %s", output_bam)

    with open(output_bam, "wb") as bam_out, open(map_log, "w") as log_fh:
        bwa_process = subprocess.Popen(
            bwa_mem_cmd, stdout=subprocess.PIPE, stderr=log_fh
        )
        samtools_process = subprocess.Popen(
            samtools_cmd,
            stdin=bwa_process.stdout,
            stdout=bam_out,
            stderr=log_fh,
        )
        if bwa_process.stdout is not None:
            bwa_process.stdout.close()
        try:
            bwa_rc = bwa_process.wait()
            sam_rc = samtools_process.wait()
        finally:
            for proc in (bwa_process, samtools_process):
                if proc.poll() is None:
                    proc.kill()

    logger.debug("BWA/samtools return codes: bwa=%s samtools=%s", bwa_rc, sam_rc)
    if bwa_rc != 0 or sam_rc != 0:
        tail = ""
        try:
            with open(map_log) as fh:
                tail = "".join(fh.readlines()[-40:])
        except OSError:
            pass
        logger.error(
            "BWA mapping failed for %s (bwa=%s samtools=%s). Log tail:\n%s",
            basename,
            bwa_rc,
            sam_rc,
            tail,
        )
        raise subprocess.CalledProcessError(
            bwa_rc or sam_rc, bwa_mem_cmd if bwa_rc else samtools_cmd
        )

    if not os.path.isfile(output_bam) or os.path.getsize(output_bam) == 0:
        logger.error("Output BAM file %s was not created or is empty", output_bam)
        raise RuntimeError(f"Output BAM file {output_bam} was not created or is empty")

    logger.info(
        "BWA mapping completed: %s (%d bytes)",
        output_bam,
        os.path.getsize(output_bam),
    )
    return output_bam


# -----------------------------------------------------------------------------
# SPAdes workers
# -----------------------------------------------------------------------------
@log_exceptions
def run_spades_for_gene(gene_dir: str) -> None:
    """
    Run SPAdes (--only-assembler, 1 thread) for one gene directory (absolute path).
    Success is determined later from contigs.fasta.
    """
    gene = os.path.basename(gene_dir.rstrip(os.sep))
    input_file = os.path.join(gene_dir, f"{gene}_interleaved.fasta")
    output_dir = os.path.join(gene_dir, f"{gene}_spades")

    if not os.path.isfile(input_file):
        logger.error("Input file %s does not exist for SPAdes assembly", input_file)
        return

    cmd = [
        "spades.py",
        "--only-assembler",
        "--threads",
        "1",
        "--12",
        input_file,
        "-o",
        output_dir,
    ]
    logger.debug("SPAdes %s: %s", gene, " ".join(cmd))
    try:
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        logger.debug("SPAdes finished for gene %s", gene)
    except subprocess.CalledProcessError as e:
        logger.error("SPAdes assembly failed for gene %s: %s", gene, e.stderr)
    except Exception as e:
        logger.error("Unexpected error during SPAdes for gene %s: %s", gene, e)


@log_exceptions
def run_spades_redo_for_gene(gene_dir: str) -> str:
    """Restart a failed SPAdes run with a reduced k-mer set. gene_dir is absolute."""
    gene = os.path.basename(gene_dir.rstrip(os.sep))
    sp_dir = os.path.join(gene_dir, f"{gene}_spades")
    try:
        kmers = sorted(int(x[1:]) for x in os.listdir(sp_dir) if x.startswith("K"))
    except Exception as e:
        logger.error("Cannot read k-mer directories for gene %s: %s", gene, e)
        return gene_dir

    if len(kmers) < 2:
        logger.error("Not enough k-mers for redo on gene %s (kmers=%s)", gene, kmers)
        return gene_dir

    redo_kmers = ",".join(map(str, kmers[:-1]))
    restart_k = f"k{kmers[-2]}"
    cmd = [
        "spades.py",
        "--restart-from",
        restart_k,
        "-k",
        redo_kmers,
        "-o",
        sp_dir,
    ]
    logger.debug("SPAdes redo %s: %s", gene, " ".join(cmd))
    try:
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        logger.debug("SPAdes redo finished for gene %s", gene)
    except subprocess.CalledProcessError as e:
        logger.error("SPAdes redo failed for gene %s: %s", gene, e.stderr)
    except Exception as e:
        logger.error("Unexpected error during SPAdes redo for gene %s: %s", gene, e)

    return gene_dir


def _check_and_copy_contigs(gene_dirs: List[str]) -> Tuple[List[str], List[str]]:
    """Copy non-empty contigs.fasta → <gene>_contigs.fasta. Returns (ok, failed) gene names."""
    spades_successful = []
    spades_failed = []

    for gene_dir in gene_dirs:
        gene = os.path.basename(gene_dir.rstrip(os.sep))
        contigs_path = os.path.join(gene_dir, f"{gene}_spades", "contigs.fasta")
        gene_failed = True

        if os.path.isfile(contigs_path) and os.path.getsize(contigs_path) > 0:
            target_path = os.path.join(gene_dir, f"{gene}_contigs.fasta")
            try:
                shutil.copy(contigs_path, target_path)
                spades_successful.append(gene)
                gene_failed = False
                logger.debug(
                    "SPAdes OK %s (%d bytes) -> %s",
                    gene,
                    os.path.getsize(contigs_path),
                    target_path,
                )
            except Exception as e:
                logger.error("Failed to copy contigs for gene %s: %s", gene, e)
        elif os.path.isfile(contigs_path):
            logger.debug("SPAdes empty contigs for gene %s", gene)
        else:
            logger.debug("SPAdes missing contigs for gene %s", gene)

        if gene_failed:
            spades_failed.append(gene)

    return spades_successful, spades_failed


def _gene_dirs(work_dir: str, genes: List[str]) -> List[str]:
    return [os.path.join(work_dir, g) for g in genes]


# -----------------------------------------------------------------------------
# SPAdes pool controllers
# -----------------------------------------------------------------------------
@log_exceptions
def spades_initial(genelist: str, cpu: int, log_queue, work_dir: str) -> List[str]:
    """Pool-run initial SPAdes for every gene; return gene names that still lack contigs."""
    require_file(genelist, "gene list")
    with open(genelist) as fh:
        genes = [ln.strip() for ln in fh if ln.strip()]

    gene_dirs = _gene_dirs(work_dir, genes)
    n_workers = max(1, min(int(cpu), len(gene_dirs)))
    logger.info(
        "Starting SPAdes for %d gene(s) with %d worker(s) under %s",
        len(genes),
        n_workers,
        work_dir,
    )
    logger.debug(
        "SPAdes gene list: %s",
        ", ".join(genes[:20]) + (" ..." if len(genes) > 20 else ""),
    )

    with managed_pool(n_workers, log_queue) as pool:
        async_result = pool.map_async(run_spades_for_gene, gene_dirs)
        async_result.get()

    logger.info("SPAdes subprocess calls completed; checking contig outputs...")
    spades_successful, spades_failed = _check_and_copy_contigs(gene_dirs)
    logger.info(
        "SPAdes initial results: %d successful, %d failed",
        len(spades_successful),
        len(spades_failed),
    )
    if spades_failed:
        logger.warning("Failed genes: %s", ", ".join(spades_failed))
    return spades_failed


@log_exceptions
def rerun_spades(genelist: str, cpu: int, log_queue, work_dir: str):
    """
    Redo failed genes that have ≥2 k-mer dirs; write spades_failed_final.txt in work_dir.
    Returns (still_failed_after_redo, all_failed) as gene names.
    """
    require_file(genelist, "failed genes list")
    with open(genelist) as fh:
        genes = [ln.strip() for ln in fh if ln.strip()]

    logger.info("Preparing SPAdes redo for %d failed gene(s)", len(genes))
    redo_genes = []
    spades_failed_initial = []

    for gene in genes:
        gene_dir = os.path.join(work_dir, gene)
        sp_dir = os.path.join(gene_dir, f"{gene}_spades")
        try:
            kmers = sorted(int(x[1:]) for x in os.listdir(sp_dir) if x.startswith("K"))
            if len(kmers) >= 2:
                redo_genes.append(gene)
                logger.debug("Redo eligible %s (kmers=%s)", gene, kmers)
            else:
                spades_failed_initial.append(gene)
                logger.debug("Redo not possible %s (kmers=%s)", gene, kmers)
        except Exception as e:
            spades_failed_initial.append(gene)
            logger.debug("Redo skip %s: %s", gene, e)

    logger.info(
        "Can redo %d gene(s); %d permanently failed",
        len(redo_genes),
        len(spades_failed_initial),
    )

    failed_final_path = os.path.join(work_dir, "spades_failed_final.txt")
    if not redo_genes:
        logger.warning("No genes can be redone")
        with open(failed_final_path, "w") as fh:
            fh.write("\n".join(spades_failed_initial))
        return genes, spades_failed_initial

    redo_dirs = _gene_dirs(work_dir, redo_genes)
    n_workers = max(1, min(int(cpu), len(redo_dirs)))
    logger.debug("SPAdes redo pool: workers=%d", n_workers)
    with managed_pool(n_workers, log_queue) as pool:
        async_result = pool.map_async(run_spades_redo_for_gene, redo_dirs)
        async_result.get()

    logger.info("SPAdes redo subprocess calls completed")
    spades_successful, redo_spades_failed = _check_and_copy_contigs(redo_dirs)
    logger.info(
        "SPAdes redo results: %d successful, %d failed",
        len(spades_successful),
        len(redo_spades_failed),
    )

    all_spades_failed = list(set(spades_failed_initial + redo_spades_failed))
    with open(failed_final_path, "w") as fh:
        fh.write("\n".join(all_spades_failed))
    logger.info(
        "After redo: %d still failed from redo set; %d total failed",
        len(redo_spades_failed),
        len(all_spades_failed),
    )
    return redo_spades_failed, all_spades_failed


@log_exceptions
def run_spades_for_genes(
    genes: List[str],
    cpu: int,
    log_queue,
    output_dir: str,
) -> List[str]:
    """
    Drive initial (+ optional redo) SPAdes for genes under output_dir.
    Status lists are written under output_dir.
    """
    if not genes:
        logger.error("No genes provided for SPAdes assembly")
        raise ValueError("No genes provided for SPAdes assembly")

    logger.info("Starting SPAdes assembly pipeline for %d genes", len(genes))
    os.makedirs(output_dir, exist_ok=True)
    logger.debug("SPAdes status directory: %s", output_dir)

    initial_file = os.path.join(output_dir, "spades_genes_initial.txt")
    with open(initial_file, "w") as fh:
        fh.write("\n".join(genes))
    logger.debug("Wrote %s", initial_file)

    spades_failed_initial = spades_initial(
        initial_file, cpu, log_queue, output_dir
    )

    if spades_failed_initial:
        logger.info(
            "Running SPAdes redo for %d failed gene(s)", len(spades_failed_initial)
        )
        failed_initial_file = os.path.join(output_dir, "spades_failed_initial.txt")
        with open(failed_initial_file, "w") as fh:
            fh.write("\n".join(spades_failed_initial))

        spades_failed_after_redo, _ = rerun_spades(
            failed_initial_file, cpu, log_queue, output_dir
        )

        failed_final_file = os.path.join(output_dir, "spades_failed_final.txt")
        if spades_failed_after_redo:
            logger.error(
                "SPAdes failed after redo for %d gene(s): %s",
                len(spades_failed_after_redo),
                ", ".join(spades_failed_after_redo),
            )
            with open(failed_final_file, "w") as fh:
                fh.write("\n".join(spades_failed_after_redo) + "\n")
        elif os.path.exists(failed_final_file):
            os.remove(failed_final_file)
            logger.debug("Removed stale %s", failed_final_file)
    else:
        failed_final_file = os.path.join(output_dir, "spades_failed_final.txt")
        if os.path.exists(failed_final_file):
            os.remove(failed_final_file)
            logger.debug("Removed stale %s", failed_final_file)

    failed_final_file = os.path.join(output_dir, "spades_failed_final.txt")
    if os.path.isfile(failed_final_file):
        with open(failed_final_file) as fh:
            spades_failed_final = [ln.strip() for ln in fh if ln.strip()]
    else:
        spades_failed_final = []

    spades_genelist = [g for g in genes if g not in set(spades_failed_final)]
    successful_file = os.path.join(output_dir, "spades_successful.txt")
    with open(successful_file, "w") as fh:
        fh.write("\n".join(spades_genelist))
    logger.debug("Wrote %s (%d genes)", successful_file, len(spades_genelist))

    if not spades_genelist:
        logger.error(
            "No genes produced assembled contigs! All %d genes failed", len(genes)
        )
        raise RuntimeError("No genes produced assembled contigs!")

    logger.info(
        "SPAdes assembly completed: %d successful, %d failed",
        len(spades_genelist),
        len(genes) - len(spades_genelist),
    )
    logger.info("SPAdes status files saved under: %s", output_dir)
    return spades_genelist


@log_exceptions
def spades(
    readfiles: List[str],
    genes: List[str],
    cpu: int,
    log_queue,
    output_dir: str,
) -> List[str]:
    """Assemble genes for one sample; returns successful gene ids."""
    if len(readfiles) != 2:
        logger.error("Expected exactly 2 paired read files, got %d", len(readfiles))
        raise ValueError("Please specify exactly two paired read files!")

    for readfile in readfiles:
        require_file(readfile, "read file")

    logger.info("Starting SPAdes for %d gene(s)", len(genes))
    return run_spades_for_genes(genes, cpu, log_queue, output_dir)


# -----------------------------------------------------------------------------
# Clean-up
# -----------------------------------------------------------------------------
@log_exceptions
def remove_spades(work_dir: str = "."):
    """Remove leftover *_spades directories under work_dir."""
    for d in os.listdir(work_dir):
        path = os.path.join(work_dir, d)
        if d.endswith("spades") and os.path.isdir(path):
            shutil.rmtree(path)
            logger.info("Removed directory %s", path)


@log_exceptions
def clean_up(sample_dir: str):
    """Remove per-gene *_spades work dirs under the sample folder."""
    logger.debug("Cleaning SPAdes work dirs under %s", sample_dir)
    for gene in os.listdir(sample_dir):
        gene_dir = os.path.join(sample_dir, gene)
        if not os.path.isdir(gene_dir):
            continue
        for d in os.listdir(gene_dir):
            path = os.path.join(gene_dir, d)
            if d.endswith("spades") and os.path.isdir(path):
                shutil.rmtree(path)
                logger.debug("Removed %s", path)


# -----------------------------------------------------------------------------
# Per-sample driver
# -----------------------------------------------------------------------------
@log_exceptions
def process_one_sample(
    sample: str,
    data_folder: str,
    baitfile: str,
    n_threads: int,
    log_queue,
) -> str:
    """
    Map + distribute + SPAdes for one sample under 20assemblies/<sample>/.
    Uses absolute paths only (safe under ThreadPoolExecutor).
    """
    assemblies_root = os.path.join(data_folder, "20assemblies")
    sample_dir = os.path.join(assemblies_root, sample)
    os.makedirs(sample_dir, exist_ok=True)

    reads_dir = os.path.join(data_folder, "10deduplicated_reads")
    readfiles = list(resolve_read_pair(reads_dir, sample))
    logger.info(
        "Processing sample %s (%d thread(s) for BWA/SPAdes)", sample, n_threads
    )
    logger.debug(
        "sample=%s dir=%s R1=%s R2=%s",
        sample,
        sample_dir,
        readfiles[0],
        readfiles[1],
    )

    bamfile = bwa(
        readfiles=readfiles,
        baitfile=baitfile,
        basename=sample,
        cpu=n_threads,
        output_dir=sample_dir,
    )
    distribute_reads_bwa(
        bamfilename=bamfile, readfiles=readfiles, work_dir=sample_dir
    )

    genes = [
        name
        for name in os.listdir(sample_dir)
        if os.path.isfile(
            os.path.join(sample_dir, name, f"{name}_interleaved.fasta")
        )
    ]
    logger.debug(
        "Sample %s: %d gene folder(s) with interleaved reads",
        sample,
        len(genes),
    )

    if not genes:
        raise RuntimeError(
            f"No genes with interleaved FASTA files for sample {sample}"
        )

    assembled = spades(
        readfiles=readfiles,
        genes=genes,
        cpu=n_threads,
        log_queue=log_queue,
        output_dir=sample_dir,
    )
    clean_up(sample_dir)
    logger.info(
        "Completed sample %s: %d gene(s) assembled",
        sample,
        len(assembled),
    )
    return sample


# -----------------------------------------------------------------------------
# Top-level assemble entry
# -----------------------------------------------------------------------------
@log_exceptions
def assemble(
    data_folder: str,
    baitfile: str,
    num_cores: int,
    log_queue,
) -> None:
    """
    Run cast_assemble for every sample in samples_list.txt.

    Indexes bait once, then processes samples in a pool sized by
    allocate_workers_and_threads(num_cores, n_samples).
    """
    logger.info("Starting cast_assemble (cores=%d)", num_cores)
    logger.debug(
        "assemble args: data_folder=%s baitfile=%s num_cores=%d",
        data_folder,
        baitfile,
        num_cores,
    )
    require_tools(REQUIRED_TOOLS)
    data_folder = os.path.abspath(data_folder)
    require_dir(data_folder, "data folder")
    baitfile = ensure_bwa_index(os.path.abspath(baitfile))

    reads_dir = require_dir(
        os.path.join(data_folder, "10deduplicated_reads"),
        "10deduplicated_reads directory",
    )
    samples_list_path = require_file(
        os.path.join(reads_dir, "samples_list.txt"),
        "samples_list.txt",
    )
    with open(samples_list_path) as fh:
        samples = [ln.strip() for ln in fh if ln.strip()]
    if not samples:
        raise ValueError(f"No samples found in {samples_list_path}")
    logger.info("Loaded %d sample(s) from %s", len(samples), samples_list_path)
    logger.debug("Sample list: %s", ", ".join(samples))

    # Preflight reads (prefer .gz) so we fail early.
    for sample in samples:
        resolve_read_pair(reads_dir, sample)

    assemblies_root = os.path.join(data_folder, "20assemblies")
    os.makedirs(assemblies_root, exist_ok=True)

    n_sample_workers, threads_per_sample = allocate_workers_and_threads(
        num_cores, len(samples)
    )
    logger.info(
        "Assembling %d sample(s): %d worker(s) x %d thread(s) per sample",
        len(samples),
        n_sample_workers,
        threads_per_sample,
    )
    logger.debug(
        "Sample pool plan: workers=%d threads=%d first=%s",
        n_sample_workers,
        threads_per_sample,
        ", ".join(samples[:5]),
    )

    # Sample-level parallelism uses threads (BWA/SPAdes are external processes).
    # Gene-level SPAdes still uses managed_pool processes — avoids nested ProcessPools.
    failures = []
    completed = 0
    total = len(samples)
    with ThreadPoolExecutor(max_workers=n_sample_workers) as executor:
        future_map = {
            executor.submit(
                process_one_sample,
                sample,
                data_folder,
                baitfile,
                threads_per_sample,
                log_queue,
            ): sample
            for sample in samples
        }
        logger.debug("Submitted %d sample job(s) on %d thread(s)", total, n_sample_workers)
        for future in as_completed(future_map):
            sample = future_map[future]
            completed += 1
            try:
                future.result()
                logger.debug("Sample job done %d/%d: %s", completed, total, sample)
                if completed == total or completed % 5 == 0:
                    logger.info(
                        "Assemble progress: %d / %d samples done", completed, total
                    )
            except Exception as e:
                logger.error("Sample %s failed: %s", sample, e)
                failures.append((sample, e))

    if failures:
        preview = ", ".join(s for s, _ in failures[:20])
        extra = " ..." if len(failures) > 20 else ""
        raise RuntimeError(
            f"cast_assemble aborted: {len(failures)} sample(s) failed "
            f"({preview}{extra}). See log for details."
        )

    logger.info(
        "cast_assemble completed: %d sample(s) under %s",
        len(samples),
        assemblies_root,
    )
