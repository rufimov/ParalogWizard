#!/usr/bin/env python
"""
------------------------------------------------------------------------------------------------------------------------
Copyright 2024 Roman Ufimov under the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

This module does the following:
  1. Filter exons and split the reference into chunks of a given size.
  2. Map the reads of each sample to the reference in parallel.
  3. Pile bam-files, call and filter variants for each chunk.
  4. Concatenate the variants into a single file.
  5. Adjust the concatenated file.
------------------------------------------------------------------------------------------------------------------------
"""

import multiprocessing
import os
import subprocess
from math import ceil

import pandas as pd
from Bio import SeqIO, SeqRecord
from Bio.Seq import Seq

from ParalogWizard import setup_logging, worker_initializer

# -----------------------------------------------------------------------------
# Module-level logger
# -----------------------------------------------------------------------------
logger = setup_logging()


# -----------------------------------------------------------------------------
# Logging Decorator (with detailed context)
# -----------------------------------------------------------------------------
def log_exceptions(func):
    """
    Decorator that logs function entry, exit, and any exceptions.
    It logs the function's name along with its positional and keyword arguments.
    Instead of calling sys.exit(), it re-raises the exception so that the main process
    can catch the error and shut down the pool gracefully.
    """
    from functools import wraps

    @wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as e:
            arg_str = ", ".join(str(arg) for arg in args)
            kwarg_str = ", ".join(f"{k}={v}" for k, v in kwargs.items())
            logger.exception(
                f"Exception in {func.__name__} (args: {arg_str}; kwargs: {kwarg_str}): {e}"
            )
            raise

    return wrapper


# -----------------------------------------------------------------------------
# Sequence Processing Functions
# -----------------------------------------------------------------------------
@log_exceptions
def trim_ends(seq1: str, seq2: str):
    """
    Trims the ends of two sequences to remove any gaps.
    :param seq1: first sequence (string)
    :param seq2: second sequence (string)
    :return: tuple of trimmed sequences
    """
    while seq1 and seq2 and (seq1[0] == "-" or seq2[0] == "-"):
        seq1 = seq1[1:]
        seq2 = seq2[1:]
    while seq1 and seq2 and (seq1[-1] == "-" or seq2[-1] == "-"):
        seq1 = seq1[:-1]
        seq2 = seq2[:-1]
    return seq1, seq2


@log_exceptions
def mafft_align(file: str) -> None:
    """
    Aligns a FASTA file using MAFFT.
    The output is written to a file named <original_basename>.mafft.fasta.
    This function calls MAFFT directly via subprocess.
    """
    if not os.path.exists(file):
        logger.error("Input file %s does not exist", file)
        raise FileNotFoundError(f"Input file {file} does not exist")

    out_file = f"{os.path.splitext(file)[0]}.mafft.fasta"
    cmd = ["mafft", "--auto", file]
    logger.info("Running MAFFT on %s", file)
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


@log_exceptions
def bwa_index(reference_to_map_to: str) -> None:
    """
    Indexes a reference using BWA.
    This function calls BWA directly via subprocess.
    :param reference_to_map_to: path to the reference file
    """
    if not os.path.exists(reference_to_map_to):
        logger.error("Reference file %s does not exist", reference_to_map_to)
        raise FileNotFoundError(f"Reference file {reference_to_map_to} does not exist")

    cmd = ["bwa", "index", "-a", "is", reference_to_map_to]
    logger.info("Creating BWA index for %s", reference_to_map_to)
    try:
        subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        logger.error("BWA indexing failed for %s: %s", reference_to_map_to, e.stderr)
        raise
    logger.info("BWA index created for %s", reference_to_map_to)


@log_exceptions
def bwa_map(exon: str, sample: str, main_data_folder: str, n_threads: int = 1) -> None:
    """
    Maps a sample to a reference using BWA and then filters, sorts, and indexes the resulting BAM file.
    Output files are written into the corresponding exon folder.
    """
    output_path = os.path.join(main_data_folder, "100remapped", exon)
    os.makedirs(output_path, exist_ok=True)

    reference_to_map_to = os.path.join(output_path, f"reference_{exon}.fas")
    r1_path = os.path.join(
        main_data_folder, "10deduplicated_reads", f"{sample}.R1.fastq.gz"
    )
    r2_path = os.path.join(
        main_data_folder, "10deduplicated_reads", f"{sample}.R2.fastq.gz"
    )

    logger.info("Mapping sample %s to exon %s using BWA", sample, exon)
    bwa_cmd = [
        "bwa",
        "mem",
        "-t",
        str(n_threads),
        reference_to_map_to,
        r1_path,
        r2_path,
    ]
    try:
        bwa_result = subprocess.run(bwa_cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        logger.error(
            "BWA mem failed for sample %s, exon %s: %s", sample, exon, e.stderr
        )
        raise

    sam_output_file = os.path.join(output_path, f"{sample}.sam")
    with open(sam_output_file, "w") as sam_out:
        sam_out.write(bwa_result.stdout)

    # Convert SAM to BAM
    bam_output_file = os.path.join(output_path, f"{sample}.bam")
    cmd_view_bam = [
        "samtools",
        "view",
        "-h",
        "-b",
        "-o",
        bam_output_file,
        sam_output_file,
    ]
    try:
        result_view_bam = subprocess.run(
            cmd_view_bam, capture_output=True, text=True, check=True
        )
        logger.info(
            "samtools view output for sample %s: %s", sample, result_view_bam.stdout
        )
    except subprocess.CalledProcessError as e:
        logger.error("samtools view failed for sample %s: %s", sample, e.stderr)
        raise

    # Filter BAM
    filtered_bam = os.path.join(output_path, f"{sample}_filtered.bam")
    cmd_filter1 = [
        "samtools",
        "view",
        "-h",
        "-b",
        "-q",
        "3",
        "-o",
        filtered_bam,
        bam_output_file,
    ]
    try:
        result_filter1 = subprocess.run(
            cmd_filter1, capture_output=True, text=True, check=True
        )
        logger.info("Filter1 output for sample %s: %s", sample, result_filter1.stdout)
    except subprocess.CalledProcessError as e:
        logger.error("Filter1 failed for sample %s: %s", sample, e.stderr)
        raise

    filtered_uniq_bam = os.path.join(output_path, f"{sample}_filtered_uniq.bam")
    cmd_filter2 = [
        "samtools",
        "view",
        "-h",
        "-b",
        "-F",
        "0x90C",
        "-o",
        filtered_uniq_bam,
        filtered_bam,
    ]
    try:
        result_filter2 = subprocess.run(
            cmd_filter2, capture_output=True, text=True, check=True
        )
        logger.info("Filter2 output for sample %s: %s", sample, result_filter2.stdout)
    except subprocess.CalledProcessError as e:
        logger.error("Filter2 failed for sample %s: %s", sample, e.stderr)
        raise

    # Sort BAM
    sorted_bam = os.path.join(output_path, f"{sample}_filtered_uniq_sorted.bam")
    cmd_sort = ["samtools", "sort", "-o", sorted_bam, filtered_uniq_bam]
    try:
        result_sort = subprocess.run(
            cmd_sort, capture_output=True, text=True, check=True
        )
        logger.info("Sort output for sample %s: %s", sample, result_sort.stdout)
    except subprocess.CalledProcessError as e:
        logger.error("Sorting failed for sample %s: %s", sample, e.stderr)
        raise

    # Index BAM
    cmd_index = ["samtools", "index", sorted_bam]
    try:
        result_index = subprocess.run(
            cmd_index, capture_output=True, text=True, check=True
        )
        logger.info("Index output for sample %s: %s", sample, result_index.stdout)
    except subprocess.CalledProcessError as e:
        logger.error("Indexing failed for sample %s: %s", sample, e.stderr)
        raise

    logger.info(
        "Sample %s mapped, filtered, sorted, and indexed successfully for exon %s",
        sample,
        exon,
    )


# -----------------------------------------------------------------------------
# Reference Splitting Functions
# -----------------------------------------------------------------------------
@log_exceptions
def split_to_chunks(input_list, chunk_size):
    """
    Splits a list into chunks of a given size.
    :param input_list: the list to split
    :param chunk_size: the size of each chunk
    :return: generator yielding chunks
    """
    for i in range(0, len(input_list), chunk_size):
        yield input_list[i : i + chunk_size]


@log_exceptions
def split_reference(
    reference_dict: dict, data_folder: str, exons_allowed: list, chunk_size: int
):
    """
    Splits the reference (a dictionary of SeqRecords) into chunks based on allowed exons,
    writes each chunk as a FASTA file in a separate folder, and returns a list of chunk identifiers.
    """
    exon_chunks = list(split_to_chunks(exons_allowed, chunk_size))
    count = 1
    chunk_ids = []
    for chunk in exon_chunks:
        concatenated_seq = ""
        chunk_index = 1
        for exon, has_paralog in chunk:
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
                    # Perform MAFFT alignment
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
            else:
                seq_exon = reference_dict[exon].seq
                if chunk_index == 1:
                    concatenated_seq += str(seq_exon)
                else:
                    concatenated_seq += ("N" * 200) + str(seq_exon)
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
        chunk_ids.append(chunk_id)
        count += 1
    return chunk_ids


# -----------------------------------------------------------------------------
# Main remapping function
# -----------------------------------------------------------------------------
@log_exceptions
def remap(
    reference_file: str, data_folder: str, num_cores: int, read_length: int, log_queue
) -> None:
    """
    Main function of the pipeline:
      1. Filter exons and split the reference into chunks.
      2. Map reads of each sample to the reference in parallel.
      (Subsequent steps: pile BAM files, call/filter variants, concatenate and adjust variants)
    :param reference_file: path to the reference FASTA file
    :param data_folder: main data folder
    :param num_cores: number of cores to use
    :param read_length: minimum read length for an exon to be included
    :param log_queue: logging queue
    """
    # Read sample list.
    with open(
        os.path.join(data_folder, "10deduplicated_reads", "samples_list.txt")
    ) as samples_f:
        sample_list = samples_f.read().splitlines()

    # Parse and filter the reference.
    reference_dict = SeqIO.to_dict(SeqIO.parse(reference_file, "fasta"))
    corrected_reference = {}
    for exon in reference_dict.keys():
        exon_short = "_".join(exon.split("_")[1:4]).replace("Contig", "exon")
        seq_str = str(reference_dict[exon].seq).replace("-", "")
        sequence = Seq(seq_str)
        record = reference_dict[exon]
        record.seq = sequence
        if len(sequence) >= read_length:
            corrected_reference[exon_short] = record

    exons_in_ref = set(corrected_reference.keys())
    paralogs_df = pd.read_csv(
        os.path.join(data_folder, "41detected_par", "all_paralogs_for_reference.tsv"),
        sep="\t",
    )
    from ParalogWizard.cast_detect import exon_stats

    exon_stats_df = exon_stats(paralogs_df)
    exon_stats_df.to_csv(
        os.path.join(data_folder, "41detected_par", "exon_statistics.tsv"),
        sep="\t",
        na_rep="NaN",
    )

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
        elif (
            perc_yes >= 0.5
            and (column in exons_in_ref)
            and (
                (column.split("_")[0] + "para_" + "_".join(column.split("_")[1:]))
                in exons_in_ref
            )
        ):
            allowed_exons.append((column, True))
    chunk_size = ceil(len(allowed_exons) / 10)
    chunk_ids = split_reference(
        corrected_reference, data_folder, allowed_exons, chunk_size
    )

    # Index each reference chunk in parallel using asynchronous calls with timeout.
    reference_chunk_files = [
        os.path.join(data_folder, "100remapped", chunk_id, f"reference_{chunk_id}.fas")
        for chunk_id in chunk_ids
    ]
    with multiprocessing.Pool(
        processes=num_cores, initializer=worker_initializer, initargs=(log_queue,)
    ) as pool_index:
        async_results = []
        for ref_file in reference_chunk_files:
            async_results.append(pool_index.apply_async(bwa_index, (ref_file,)))
        for ref_file, async_result in zip(reference_chunk_files, async_results):
            try:
                async_result.get(timeout=600)
            except multiprocessing.TimeoutError:
                logger.error("Timeout indexing reference chunk %s", ref_file)
            except Exception as e:
                logger.error("Error indexing reference chunk %s: %s", ref_file, e)
        pool_index.close()
        pool_index.join()

    # Map each sample to each chunk in parallel using asynchronous calls with timeout.
    mapping_args = []
    for exon_chunk in chunk_ids:
        for sample in sample_list:
            mapping_args.append((exon_chunk, sample, data_folder))
    with multiprocessing.Pool(
        processes=num_cores, initializer=worker_initializer, initargs=(log_queue,)
    ) as pool_bwa:
        async_results = []
        for args in mapping_args:
            async_results.append(pool_bwa.apply_async(bwa_map, args))
        for args, async_result in zip(mapping_args, async_results):
            exon, sample, _ = args
            try:
                async_result.get(timeout=600)
            except multiprocessing.TimeoutError:
                logger.error("Timeout mapping sample %s to exon %s", sample, exon)
            except Exception as e:
                logger.error("Error mapping sample %s to exon %s: %s", sample, exon, e)
        pool_bwa.close()
        pool_bwa.join()

    logger.info("Mapping done.")
