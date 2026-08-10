#!/usr/bin/env python
"""
------------------------------------------------------------------------------------------------------------------------
Copyright 2024 Roman Ufimov under the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

cast_remap — map reads to a chunked reference under 100remapped/.

  1. Select exons from the customized reference using paralog statistics and a
     minimum length threshold.
  2. Split selected exons into ~10 reference chunks (exons1, exons2, …), pairing
     ortholog/paralog sequences with MAFFT when both are kept.
  3. Index each chunk with BWA.
  4. Map each sample’s deduplicated reads to each chunk (bwa mem | samtools),
     writing {sample}_filtered_uniq_sorted.bam (+ .bai) with RG/SM set to the
     sample ID.

Variant calling is done separately by cast_call.
------------------------------------------------------------------------------------------------------------------------
"""

import logging
import os
import subprocess
from math import ceil

import pandas as pd
from Bio import SeqIO, SeqRecord
from Bio.Seq import Seq

from ParalogWizard import log_exceptions, managed_pool
from ParalogWizard.cast_call import (
    allocate_workers_and_threads,
    require_dir,
    require_file,
    require_tools,
)

logger = logging.getLogger("ParalogWizard")

REQUIRED_TOOLS = ("bwa", "samtools", "mafft")


def _check_sample_reads(data_folder, sample_list):
    """Verify R1/R2 gzipped reads exist and are non-empty for every sample."""
    reads_dir = require_dir(
        os.path.join(data_folder, "10deduplicated_reads"),
        "10deduplicated_reads directory",
    )
    problems = []
    for sample in sample_list:
        for mate in ("R1", "R2"):
            path = os.path.join(reads_dir, f"{sample}.{mate}.fastq.gz")
            if not os.path.isfile(path):
                problems.append(f"{sample}: missing {mate} ({path})")
                logger.debug("Missing reads: %s", path)
            elif os.path.getsize(path) == 0:
                problems.append(f"{sample}: empty {mate} ({path})")
                logger.debug("Empty reads: %s", path)
            else:
                logger.debug(
                    "OK reads %s %s (%d bytes)", sample, mate, os.path.getsize(path)
                )
    if problems:
        preview = "\n  - ".join(problems[:30])
        extra = f"\n  - ... and {len(problems) - 30} more" if len(problems) > 30 else ""
        logger.error("Read-file check failed:\n  - %s%s", preview, extra)
        raise FileNotFoundError(
            f"cast_remap aborted: {len(problems)} read-file problem(s). See log."
        )
    logger.info(
        "Verified R1/R2 reads for %d sample(s) in %s", len(sample_list), reads_dir
    )


# -----------------------------------------------------------------------------
# Alignment helpers (used when building ortholog/paralog chunk references)
# -----------------------------------------------------------------------------
@log_exceptions
def trim_ends(seq1: str, seq2: str):
    """Trim leading/trailing gap columns from a pair of aligned sequences."""
    before = (len(seq1), len(seq2))
    while seq1 and seq2 and (seq1[0] == "-" or seq2[0] == "-"):
        seq1 = seq1[1:]
        seq2 = seq2[1:]
    while seq1 and seq2 and (seq1[-1] == "-" or seq2[-1] == "-"):
        seq1 = seq1[:-1]
        seq2 = seq2[:-1]
    logger.debug(
        "trim_ends: %s -> %s", before, (len(seq1), len(seq2))
    )
    return seq1, seq2


@log_exceptions
def mafft_align(file: str) -> None:
    """
    Align a FASTA with MAFFT --auto.
    Writes <basename>.mafft.fasta next to the input file.
    """
    if not os.path.exists(file):
        logger.error("Input file %s does not exist", file)
        raise FileNotFoundError(f"Input file {file} does not exist")

    out_file = f"{os.path.splitext(file)[0]}.mafft.fasta"
    cmd = ["mafft", "--auto", file]
    logger.info("Running MAFFT on %s", file)
    logger.debug("MAFFT command: %s", " ".join(cmd))
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        with open(out_file, "w") as aligned:
            aligned.write(result.stdout)
    except subprocess.CalledProcessError as e:
        logger.error("MAFFT alignment failed for %s: %s", file, e.stderr)
        raise
    logger.info(
        "MAFFT alignment completed for %s; output written to %s", file, out_file
    )
    logger.debug("MAFFT output size: %d bytes", os.path.getsize(out_file))


@log_exceptions
def bwa_index(reference_to_map_to: str) -> None:
    """Build a BWA index for a reference FASTA (bwa index -a is)."""
    if not os.path.exists(reference_to_map_to):
        logger.error("Reference file %s does not exist", reference_to_map_to)
        raise FileNotFoundError(f"Reference file {reference_to_map_to} does not exist")

    cmd = ["bwa", "index", "-a", "is", reference_to_map_to]
    logger.info("Creating BWA index for %s", reference_to_map_to)
    logger.debug("BWA index command: %s", " ".join(cmd))
    try:
        subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        logger.error("BWA indexing failed for %s: %s", reference_to_map_to, e.stderr)
        raise
    logger.info("BWA index created for %s", reference_to_map_to)


def _run_pipeline(commands, log_file):
    """
    Run argv command lists as a stdin/stdout pipeline.
    Merged stderr from all stages is written to log_file.
    """
    logger.debug(
        "Pipeline (%d stage(s)) -> %s\n%s",
        len(commands),
        log_file,
        " | ".join(" ".join(c) for c in commands),
    )
    with open(log_file, "w") as log_fh:
        processes = []
        prev_stdout = None
        try:
            for i, cmd in enumerate(commands):
                is_last = i == len(commands) - 1
                proc = subprocess.Popen(
                    cmd,
                    stdin=prev_stdout,
                    stdout=subprocess.PIPE if not is_last else subprocess.DEVNULL,
                    stderr=log_fh,
                )
                logger.debug("Started pipeline stage %d pid=%s: %s", i, proc.pid, cmd[0])
                if prev_stdout is not None:
                    prev_stdout.close()
                processes.append(proc)
                prev_stdout = proc.stdout

            return_codes = [proc.wait() for proc in processes]
        finally:
            for proc in processes:
                if proc.poll() is None:
                    proc.kill()

    logger.debug("Pipeline return codes: %s", return_codes)
    if any(code != 0 for code in return_codes):
        tail = ""
        try:
            with open(log_file) as fh:
                tail = "".join(fh.readlines()[-40:])
        except OSError:
            pass
        raise subprocess.CalledProcessError(
            return_codes[-1] if return_codes else 1,
            commands,
            output=tail,
        )


@log_exceptions
def bwa_map(
    exon: str,
    sample: str,
    main_data_folder: str,
    n_threads: int = 1,
    force: bool = False,
) -> None:
    """
    Map one sample to one reference chunk under 100remapped/<chunk>/.

    Pipeline: bwa mem (-R RG/SM=<sample>) | samtools view (-q 3 -F 0x90C) |
    samtools sort → {sample}_filtered_uniq_sorted.bam (+ index).

    Skips if a non-empty BAM and .bai already exist, unless force=True.
    exon here is a chunk id (e.g. exons3), not a single locus name.
    """
    output_path = os.path.join(main_data_folder, "100remapped", exon)
    os.makedirs(output_path, exist_ok=True)
    logger.debug(
        "bwa_map sample=%s chunk=%s threads=%d force=%s out=%s",
        sample,
        exon,
        n_threads,
        force,
        output_path,
    )

    sorted_bam = os.path.join(output_path, f"{sample}_filtered_uniq_sorted.bam")
    sorted_bai = sorted_bam + ".bai"
    if (
        not force
        and os.path.isfile(sorted_bam)
        and os.path.isfile(sorted_bai)
        and os.path.getsize(sorted_bam) > 0
    ):
        logger.info(
            "Skipping sample %s / exon %s: sorted BAM already exists", sample, exon
        )
        logger.debug(
            "Skip existing BAM %s (%d bytes)",
            sorted_bam,
            os.path.getsize(sorted_bam),
        )
        return

    reference_to_map_to = require_file(
        os.path.join(output_path, f"reference_{exon}.fas"),
        f"reference for chunk {exon}",
    )
    r1_path = require_file(
        os.path.join(
            main_data_folder, "10deduplicated_reads", f"{sample}.R1.fastq.gz"
        ),
        f"R1 reads for {sample}",
    )
    r2_path = require_file(
        os.path.join(
            main_data_folder, "10deduplicated_reads", f"{sample}.R2.fastq.gz"
        ),
        f"R2 reads for {sample}",
    )

    # SM matches samples_list.txt so cast_call / bcftools sample names stay clean.
    read_group = f"@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA"
    map_log = os.path.join(output_path, f"{sample}_bwa_map.log")
    logger.debug("Read group: %r", read_group)

    logger.info(
        "Mapping sample %s to chunk %s using BWA (%d threads)",
        sample,
        exon,
        n_threads,
    )
    commands = [
        [
            "bwa",
            "mem",
            "-t",
            str(n_threads),
            "-R",
            read_group,
            reference_to_map_to,
            r1_path,
            r2_path,
        ],
        ["samtools", "view", "-b", "-q", "3", "-F", "0x90C"],
        ["samtools", "sort", "-@", str(n_threads), "-o", sorted_bam, "-"],
    ]
    try:
        _run_pipeline(commands, map_log)
    except subprocess.CalledProcessError as e:
        logger.error(
            "Mapping pipeline failed for sample %s, chunk %s. See %s\n%s",
            sample,
            exon,
            map_log,
            e.output or "",
        )
        raise

    if not os.path.isfile(sorted_bam) or os.path.getsize(sorted_bam) == 0:
        raise RuntimeError(
            f"Mapping produced no BAM for sample {sample}, chunk {exon}: {sorted_bam}"
        )

    try:
        subprocess.run(
            ["samtools", "index", sorted_bam],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True,
        )
    except subprocess.CalledProcessError as e:
        logger.error(
            "Indexing failed for sample %s, chunk %s: %s",
            sample,
            exon,
            e.stderr,
        )
        raise
    if not os.path.isfile(sorted_bai):
        raise RuntimeError(f"Missing BAM index after samtools index: {sorted_bai}")

    # Drop intermediates left by older remap implementations.
    for legacy in (
        f"{sample}.sam",
        f"{sample}.bam",
        f"{sample}_filtered.bam",
        f"{sample}_filtered_uniq.bam",
    ):
        legacy_path = os.path.join(output_path, legacy)
        if os.path.isfile(legacy_path):
            try:
                os.remove(legacy_path)
                logger.debug("Removed legacy intermediate: %s", legacy_path)
            except OSError as exc:
                logger.debug("Could not remove %s: %s", legacy_path, exc)

    logger.info(
        "Mapped %s -> %s (%d bytes)",
        sample,
        sorted_bam,
        os.path.getsize(sorted_bam),
    )


# -----------------------------------------------------------------------------
# Reference chunking
# -----------------------------------------------------------------------------
@log_exceptions
def split_to_chunks(input_list, chunk_size):
    """Yield successive slices of input_list with length chunk_size (last may be shorter)."""
    for i in range(0, len(input_list), chunk_size):
        yield input_list[i : i + chunk_size]


@log_exceptions
def split_reference(
    reference_dict: dict, data_folder: str, exons_allowed: list, chunk_size: int
):
    """
    Build concatenated reference FASTAs for each chunk under 100remapped/exonsN/.

    exons_allowed is a list of (exon_id, has_paralog). When has_paralog is True,
    ortholog and paralog are MAFFT-aligned, gap-trimmed, and joined with N-spacers
    before concatenation into the chunk sequence.
    Returns chunk ids such as ['exons1', 'exons2', ...].
    """
    exon_chunks = list(split_to_chunks(exons_allowed, chunk_size))
    logger.debug(
        "split_reference: %d exon(s) -> %d chunk(s) (chunk_size=%d)",
        len(exons_allowed),
        len(exon_chunks),
        chunk_size,
    )
    count = 1
    chunk_ids = []
    for chunk in exon_chunks:
        concatenated_seq = ""
        chunk_index = 1
        logger.debug(
            "Building chunk exons%d with %d exon entry(ies)",
            count,
            len(chunk),
        )
        for exon, has_paralog in chunk:
            logger.debug("  exon=%s has_paralog=%s", exon, has_paralog)
            if has_paralog:
                exon_para = f"{exon.split('_')[0]}para_{'_'.join(exon.split('_')[1:])}"
                if exon in reference_dict and exon_para in reference_dict:
                    seq_exon = reference_dict[exon].seq
                    seq_para = reference_dict[exon_para].seq
                    temp_fasta = f"{exon}.fasta"
                    with open(temp_fasta, "w") as temp_file:
                        SeqIO.write(
                            SeqRecord.SeqRecord(
                                seq_exon, id=exon, name=exon, description=exon
                            ),
                            temp_file,
                            "fasta-2line",
                        )
                        SeqIO.write(
                            SeqRecord.SeqRecord(
                                seq_para,
                                id=exon_para,
                                name=exon_para,
                                description=exon_para,
                            ),
                            temp_file,
                            "fasta-2line",
                        )
                    mafft_align(temp_fasta)
                    aligned_file = f"{os.path.splitext(temp_fasta)[0]}.mafft.fasta"
                    aligned_records = list(SeqIO.parse(aligned_file, "fasta"))
                    os.remove(temp_fasta)
                    os.remove(aligned_file)
                    seq_exon_aln = aligned_records[0].seq
                    seq_para_aln = aligned_records[1].seq
                    trimmed_exon, trimmed_para = trim_ends(
                        str(seq_exon_aln), str(seq_para_aln)
                    )
                    if chunk_index == 1:
                        concatenated_seq = (
                            trimmed_exon.replace("-", "")
                            + ("N" * 200)
                            + trimmed_para.replace("-", "")
                        )
                    else:
                        concatenated_seq += (
                            ("N" * 200)
                            + trimmed_exon.replace("-", "")
                            + ("N" * 200)
                            + trimmed_para.replace("-", "")
                        )
                    logger.debug(
                        "  paired %s + %s -> concat len %d",
                        exon,
                        exon_para,
                        len(concatenated_seq),
                    )
                else:
                    logger.debug(
                        "  skipping pair %s / %s (missing from reference_dict)",
                        exon,
                        exon_para,
                    )
            else:
                seq_exon = reference_dict[exon].seq
                if chunk_index == 1:
                    concatenated_seq += str(seq_exon)
                else:
                    concatenated_seq += ("N" * 200) + str(seq_exon)
                logger.debug(
                    "  single %s -> concat len %d", exon, len(concatenated_seq)
                )
            chunk_index += 1
        chunk_id = f"exons{count}"
        output_dir = os.path.join(data_folder, "100remapped", chunk_id)
        os.makedirs(output_dir, exist_ok=True)
        output_fasta = os.path.join(output_dir, f"reference_{chunk_id}.fas")
        seq_record_to_write = SeqRecord.SeqRecord(
            Seq(concatenated_seq), id=chunk_id, name=chunk_id, description=chunk_id
        )
        with open(output_fasta, "w") as ref_file:
            SeqIO.write(seq_record_to_write, ref_file, "fasta")
        logger.debug(
            "Wrote %s (%d bp, %d bytes)",
            output_fasta,
            len(concatenated_seq),
            os.path.getsize(output_fasta),
        )
        chunk_ids.append(chunk_id)
        count += 1
    return chunk_ids


# -----------------------------------------------------------------------------
# Entry point
# -----------------------------------------------------------------------------
@log_exceptions
def remap(
    reference_file: str,
    data_folder: str,
    num_cores: int,
    read_length: int,
    log_queue,
    force: bool = False,
) -> None:
    """
    Run cast_remap: select exons, build ~10 reference chunks, BWA-index them,
    and map every sample to every chunk under 100remapped/.

    Parallelism: num_cores is split into process workers × BWA/samtools threads
    (same helper as cast_call). Aborts if any index or mapping job fails.
    Does not call variants — use cast_call afterwards.

    :param reference_file: customized probe/reference FASTA
    :param data_folder: pipeline data root (contains 10deduplicated_reads/, etc.)
    :param num_cores: CPUs to use (from -nc / cluster reservation)
    :param read_length: drop reference exons shorter than this
    :param log_queue: multiprocessing logging queue
    :param force: remap even when sorted BAM+.bai already exist
    """
    logger.info("Starting cast_remap (force=%s, min exon length=%d)", force, read_length)
    logger.debug(
        "remap args: reference=%s data_folder=%s num_cores=%d force=%s",
        reference_file,
        data_folder,
        num_cores,
        force,
    )
    require_tools(REQUIRED_TOOLS)
    require_dir(data_folder, "data folder")
    require_file(reference_file, "reference / probes_customized FASTA")

    samples_list_path = os.path.join(
        data_folder, "10deduplicated_reads", "samples_list.txt"
    )
    require_file(samples_list_path, "samples_list.txt")
    with open(samples_list_path) as samples_f:
        sample_list = [line.strip() for line in samples_f if line.strip()]
    if not sample_list:
        raise ValueError(f"No samples found in {samples_list_path}")
    duplicates = sorted({s for s in sample_list if sample_list.count(s) > 1})
    if duplicates:
        raise ValueError(
            f"Duplicate sample name(s) in {samples_list_path}: {', '.join(duplicates)}"
        )
    logger.info("Loaded %d sample(s) from %s", len(sample_list), samples_list_path)
    logger.debug("Sample list: %s", ", ".join(sample_list))
    _check_sample_reads(data_folder, sample_list)

    paralogs_tsv = require_file(
        os.path.join(data_folder, "41detected_par", "all_paralogs_for_reference.tsv"),
        "all_paralogs_for_reference.tsv",
    )

    logger.info("Parsing reference FASTA: %s", reference_file)
    reference_dict = SeqIO.to_dict(SeqIO.parse(reference_file, "fasta"))
    if not reference_dict:
        raise ValueError(f"No sequences found in reference FASTA: {reference_file}")
    logger.debug("Reference FASTA records: %d", len(reference_dict))

    corrected_reference = {}
    dropped_short = 0
    for exon in reference_dict.keys():
        exon_short = "_".join(exon.split("_")[1:4]).replace("Contig", "exon")
        seq_str = str(reference_dict[exon].seq).replace("-", "")
        sequence = Seq(seq_str)
        record = reference_dict[exon]
        record.seq = sequence
        if len(sequence) >= read_length:
            corrected_reference[exon_short] = record
            logger.debug(
                "Keep exon %s (from %s) len=%d", exon_short, exon, len(sequence)
            )
        else:
            dropped_short += 1
            logger.debug(
                "Drop short exon %s (from %s) len=%d < %d",
                exon_short,
                exon,
                len(sequence),
                read_length,
            )
    if not corrected_reference:
        raise ValueError(
            f"No reference exons retained with length >= {read_length} from {reference_file}"
        )
    logger.info(
        "Retained %d / %d reference exon record(s) with length >= %d",
        len(corrected_reference),
        len(reference_dict),
        read_length,
    )
    logger.debug("Dropped %d short exon record(s)", dropped_short)

    exons_in_ref = set(corrected_reference.keys())
    paralogs_df = pd.read_csv(paralogs_tsv, sep="\t")
    if paralogs_df.empty:
        raise ValueError(f"Paralog table is empty: {paralogs_tsv}")
    logger.debug(
        "Paralog table shape=%s columns=%s",
        paralogs_df.shape,
        list(paralogs_df.columns),
    )
    from ParalogWizard.cast_detect import exon_stats

    exon_stats_df = exon_stats(paralogs_df)
    exon_stats_path = os.path.join(data_folder, "41detected_par", "exon_statistics.tsv")
    exon_stats_df.to_csv(exon_stats_path, sep="\t", na_rep="NaN")
    logger.info("Wrote exon statistics: %s", exon_stats_path)

    # Keep single-copy exons (≥90% "No") or clear paralog pairs (≥50% "Yes").
    allowed_exons = []
    for column in exon_stats_df.columns:
        num_yes = 0
        num_no = 0
        try:
            num_yes = exon_stats_df[column].value_counts()["Yes"]
        except KeyError:
            try:
                num_no = exon_stats_df[column].value_counts()["No"]
            except KeyError:
                pass
        perc_yes = num_yes / len(sample_list)
        perc_no = num_no / len(sample_list)
        if perc_no >= 0.9 and (column in exons_in_ref):
            allowed_exons.append((column, False))
            logger.debug(
                "Select single-copy %s (perc_no=%.2f)", column, perc_no
            )
        elif (
            perc_yes >= 0.5
            and (column in exons_in_ref)
            and (
                (column.split("_")[0] + "para_" + "_".join(column.split("_")[1:]))
                in exons_in_ref
            )
        ):
            allowed_exons.append((column, True))
            logger.debug(
                "Select paralog pair %s (perc_yes=%.2f)", column, perc_yes
            )
        else:
            logger.debug(
                "Reject exon %s (perc_yes=%.2f perc_no=%.2f in_ref=%s)",
                column,
                perc_yes,
                perc_no,
                column in exons_in_ref,
            )
    if not allowed_exons:
        raise RuntimeError(
            "No exons passed selection thresholds for remapping. "
            "Check exon_statistics.tsv and the reference FASTA."
        )
    logger.info(
        "Selected %d exon(s) for remapping (%d with paralog partner)",
        len(allowed_exons),
        sum(1 for _, has_para in allowed_exons if has_para),
    )

    remapped_root = os.path.join(data_folder, "100remapped")
    os.makedirs(remapped_root, exist_ok=True)
    chunk_size = max(1, ceil(len(allowed_exons) / 10))
    logger.info(
        "Building reference chunks (chunk_size=%d) under %s",
        chunk_size,
        remapped_root,
    )
    logger.debug(
        "Chunk plan: %d exon(s) / ~10 groups -> chunk_size=%d",
        len(allowed_exons),
        chunk_size,
    )
    chunk_ids = split_reference(
        corrected_reference, data_folder, allowed_exons, chunk_size
    )
    if not chunk_ids:
        raise RuntimeError("split_reference returned no chunk ids")
    logger.info("Created %d reference chunk(s): %s", len(chunk_ids), ", ".join(chunk_ids))

    reference_chunk_files = [
        require_file(
            os.path.join(
                data_folder, "100remapped", chunk_id, f"reference_{chunk_id}.fas"
            ),
            f"chunk reference {chunk_id}",
        )
        for chunk_id in chunk_ids
    ]
    n_index_workers, _ = allocate_workers_and_threads(
        num_cores, len(reference_chunk_files)
    )
    logger.info(
        "Indexing %d reference chunk(s) with %d worker(s)",
        len(reference_chunk_files),
        n_index_workers,
    )
    logger.debug("Index workers=%d files=%s", n_index_workers, reference_chunk_files)
    index_failures = []
    with managed_pool(n_index_workers, log_queue) as pool_index:
        async_results = [
            (ref_file, pool_index.apply_async(bwa_index, (ref_file,)))
            for ref_file in reference_chunk_files
        ]
        logger.debug("Submitted %d index job(s)", len(async_results))
        for ref_file, async_result in async_results:
            try:
                async_result.get()
                logger.info("Indexed %s", ref_file)
            except Exception as e:
                logger.error("Error indexing reference chunk %s: %s", ref_file, e)
                index_failures.append((ref_file, e))
    if index_failures:
        raise RuntimeError(
            f"cast_remap aborted: failed to index {len(index_failures)} "
            f"reference chunk(s). See log for details."
        )

    mapping_args = [
        (exon_chunk, sample, data_folder)
        for exon_chunk in chunk_ids
        for sample in sample_list
    ]
    n_map_workers, threads_per_map = allocate_workers_and_threads(
        num_cores, len(mapping_args)
    )
    logger.info(
        "Mapping %d sample×chunk job(s): %d worker(s) x %d BWA/samtools thread(s)",
        len(mapping_args),
        n_map_workers,
        threads_per_map,
    )
    logger.debug(
        "Map pool plan: workers=%d threads=%d first_jobs=%s",
        n_map_workers,
        threads_per_map,
        [
            f"{sample}/{exon}"
            for exon, sample, _ in mapping_args[:5]
        ],
    )
    map_failures = []
    with managed_pool(n_map_workers, log_queue) as pool_bwa:
        async_results = [
            (
                args,
                pool_bwa.apply_async(
                    bwa_map,
                    args,
                    {"n_threads": threads_per_map, "force": force},
                ),
            )
            for args in mapping_args
        ]
        total = len(async_results)
        logger.debug("Submitted %d mapping job(s)", total)
        for i, (args, async_result) in enumerate(async_results, start=1):
            exon, sample, _ = args
            try:
                async_result.get()
                logger.debug("Mapping job done %d/%d: %s/%s", i, total, sample, exon)
                if i == total or i % 25 == 0:
                    logger.info("Mapping progress: %d / %d jobs done", i, total)
            except Exception as e:
                logger.error(
                    "Error mapping sample %s to chunk %s: %s", sample, exon, e
                )
                map_failures.append((exon, sample, e))
    if map_failures:
        preview = ", ".join(
            f"{sample}/{exon}" for exon, sample, _ in map_failures[:20]
        )
        extra = " ..." if len(map_failures) > 20 else ""
        raise RuntimeError(
            f"cast_remap aborted: {len(map_failures)} mapping job(s) failed "
            f"({preview}{extra}). See log for details."
        )

    logger.info(
        "cast_remap completed: %d chunk(s), %d sample(s) under %s",
        len(chunk_ids),
        len(sample_list),
        remapped_root,
    )
