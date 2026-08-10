#!/usr/bin/env python
"""
Module for ParalogWizard cast_retrieve command.

This command retrieves and processes contigs assembled during the assembly step.
It collects raw contigs, prepares (copies/renames) them, creates BLAST hit tables,
corrects contigs, computes statistics, renames files, and performs cleanup.
Multiprocessing is used where appropriate, and detailed logging is provided.
"""

import fileinput
import logging
import multiprocessing
import subprocess
from glob import glob
import os
import re
import shutil
import sys
from pathlib import Path
from typing import Dict

import Bio.SeqRecord
import pandas
from Bio import SeqIO

from ParalogWizard import worker_initializer, log_exceptions

# Get logger by name (will be configured by ParalogWizard.py)
logger = logging.getLogger("ParalogWizard")


# Using unified log_exceptions decorator from ParalogWizard.__init__


# -----------------------------------------------------------------------------
# Utility Functions and Context Managers
# -----------------------------------------------------------------------------
@log_exceptions
def slicing(
    dictionary: Dict[str, Bio.SeqRecord.SeqRecord], entry: pandas.Series
) -> str:
    """
    Extract the sequence slice (exon) from the contig using hit boundaries.
    """
    sequence = dictionary[entry["saccver"]].seq
    start = entry["sstart"]
    end = entry["send"]
    logger.debug(
        "Slicing sequence for contig %s: start=%s, end=%s", entry["saccver"], start, end
    )
    if start > end:
        result = str(sequence.reverse_complement())[-end : -start - 1 : -1][::-1]
    else:
        result = str(sequence)[start - 1 : end]
    return result


@log_exceptions
def collect_contigs(data_folder: str) -> None:
    """
    Collect all contigs from each sample in the '20assemblies' folder
    and write them into individual files in the '30raw_contigs' folder.
    """
    target_folder = os.path.join(data_folder, "30raw_contigs")
    os.makedirs(target_folder, exist_ok=True)
    folder20 = os.path.join(data_folder, "20assemblies")
    logger.info("Collecting raw contigs from folder: %s", folder20)
    for sample in os.listdir(folder20):
        sample_path = os.path.join(folder20, sample)
        if os.path.isdir(sample_path):
            logger.info("Processing sample: %s", sample)
            out_file = os.path.join(target_folder, f"{sample}_contigs.fasta")
            with open(out_file, "w") as contigs_out:
                for locus in os.listdir(sample_path):
                    locus_file = os.path.join(
                        sample_path, locus, f"{locus}_contigs.fasta"
                    )
                    if os.path.exists(locus_file):
                        logger.info("\tProcessing locus: %s", locus)
                        with open(locus_file) as lf:
                            content = lf.read().replace(">", f">{locus}_")
                            contigs_out.write(content)
                        logger.debug("\tLocus %s processed", locus)
            if os.stat(out_file).st_size == 0:
                logger.warning("Output file %s is empty, removing it.", out_file)
                os.remove(out_file)
            logger.info("Sample %s processed.", sample)
    logger.info("Contig collection complete.")


@log_exceptions
def prepare_contigs(file_to_open: str, file_to_write: str) -> None:
    """
    Copy a FASTA file from 30raw_contigs to 31exonic_contigs,
    adjust sequence headers by shortening them, and rename the file.
    """
    logger.info("Preparing contigs from %s to %s", file_to_open, file_to_write)
    shutil.copy(file_to_open, file_to_write)
    with fileinput.FileInput(file_to_write, inplace=True) as fasta:
        for line in fasta:
            line = re.sub("NODE", "N", line)
            line = re.sub(r"length_([0-9]+)_cov_([0-9]+\.[0-9]{2}).*", r"\1_c_\2", line)
            print(line, end="")
    base_parts = os.path.basename(file_to_write).split(".")[0].split("_")
    new_name = "_".join(base_parts[:2]) + ".fasta"
    new_path = os.path.join(os.path.dirname(file_to_write), new_name)
    os.rename(file_to_write, new_path)
    logger.info("Contigs prepared and renamed to %s", new_name)


# -----------------------------------------------------------------------------
# BLAST Hit Table Creation
# -----------------------------------------------------------------------------
@log_exceptions
def create_hit_tables(
    fasta_file: str, probe_exons: str, n_cpu: int, length_cover: int
) -> None:
    """
    Create a BLAST hit table for the given FASTA file using BLAST+ command-line tools.
    This function calls makeblastdb and blastn via subprocess.
    """
    path = os.path.dirname(fasta_file)
    fasta_basename = os.path.basename(fasta_file)
    sample = Path(fasta_basename).stem
    db_name = os.path.join(path, sample)

    # Create BLAST database
    makeblastdb_cmd = (
        f"makeblastdb -dbtype nucl -in {os.path.join(path, fasta_basename)} "
        f"-out {db_name} -parse_seqids"
    )
    logger.info(
        "Running makeblastdb for sample %s with command: %s", sample, makeblastdb_cmd
    )
    ret = subprocess.call(makeblastdb_cmd, shell=True)
    if ret:
        logger.error("makeblastdb failed for sample %s with exit code %d", sample, ret)
        raise RuntimeError(f"makeblastdb failed for sample {sample} with exit code {ret}")

    # Run blastn
    blast_out = os.path.join(path, f"reference_in_{sample}_contigs.txt")
    blastn_cmd = (
        f"blastn -task blastn -query {probe_exons} -db {db_name} "
        f"-out {blast_out} -qcov_hsp_perc {length_cover} "
        f"-outfmt '6 qaccver saccver pident qcovhsp evalue bitscore sstart send' "
        f"-num_threads {n_cpu}"
    )
    logger.info("Running blastn for sample %s with command: %s", sample, blastn_cmd)
    ret = subprocess.call(blastn_cmd, shell=True)
    if ret:
        logger.error("blastn failed for sample %s with exit code %d", sample, ret)
        raise RuntimeError(f"blastn failed for sample {sample} with exit code {ret}")
    logger.info("Hit table for sample %s is ready", sample)


# -----------------------------------------------------------------------------
# Contig Correction and Statistics
# -----------------------------------------------------------------------------
@log_exceptions
def correct_contgis(file: str, spades_cover: float) -> pandas.DataFrame:
    """
    Filter the BLAST hit table, keep one hit per contig per exon,
    slice contigs to exonic contigs using hit boundaries, and write a new FASTA file.
    """
    logger.info("\tCorrecting file %s", os.path.basename(file))
    folder = os.path.dirname(file)
    sample = Path(file).stem
    blast_file = os.path.join(folder, f"reference_in_{sample}_contigs.txt")
    blast_results = pandas.read_csv(
        blast_file,
        sep="\t",
        names=[
            "qaccver",
            "saccver",
            "pident",
            "qcovhsp",
            "evalue",
            "bitscore",
            "sstart",
            "send",
        ],
    )
    blast_results["locus"] = (
        blast_results["qaccver"].str.split("-").str[1].str.split("_").str[0]
    )
    blast_results["contig_locus"] = blast_results["saccver"].str.split("_").str[0]
    blast_results["k-mer_cover"] = (
        blast_results["saccver"].str.split("_").str[-1].astype(float)
    )
    hits = blast_results[
        (blast_results["locus"] == blast_results["contig_locus"])
        & (blast_results["k-mer_cover"] >= spades_cover)
    ].reset_index(drop=True)
    hits["exon"] = blast_results["qaccver"].str.split("-").str[1]
    hits_sorted = hits.sort_values(
        ["exon", "qcovhsp", "pident", "evalue", "bitscore"],
        ascending=[True, False, False, True, False],
    )
    hits_dedup = hits_sorted.drop_duplicates(subset=["exon", "saccver"]).reset_index(
        drop=True
    )
    hits_dedup["sample"] = sample
    out_fasta = os.path.join(folder, f"{sample}.fas")
    with open(out_fasta, "w") as result_fasta, open(file) as contigs:
        contigs_dict = SeqIO.to_dict(SeqIO.parse(contigs, "fasta"))
        for i in range(len(hits_dedup)):
            seq_slice = slicing(contigs_dict, hits_dedup.loc[i])
            result_fasta.write(
                f">{hits_dedup.loc[i, 'exon']}_{'_'.join(hits_dedup.loc[i, 'saccver'].split('_')[1:])}\n"
                f"{seq_slice}\n"
            )
            hits_dedup.loc[i, "sequence"] = seq_slice
    logger.info("\tContig correction for %s completed", sample)
    return hits_dedup


@log_exceptions
def copies_stats(all_hits: pandas.DataFrame) -> pandas.DataFrame:
    """
    Create a statistics table with the maximum number of hits per exon (deduced copies)
    for each locus per sample.
    """
    statistics = pandas.DataFrame([], columns=["locus"])
    statistics.set_index("locus", drop=True, inplace=True)
    for (sample, locus), group_df in all_hits.groupby(["sample", "locus"]):
        count = group_df.groupby("exon")["sequence"].count().max()
        statistics.loc[locus, sample] = count
    return statistics


@log_exceptions
def rename_contigs(file: str) -> None:
    """
    Rename contigs in a FASTA file to a shorter format as required by HybPhyloMaker.
    """
    sample = os.path.basename(os.path.splitext(file)[0])
    logger.info("\tRenaming contigs in %s", sample)
    with open(file) as fasta_in:
        records = SeqIO.to_dict(SeqIO.parse(fasta_in, "fasta"))
        counter = 1
        lines = []
        for key in sorted(records.keys()):
            lines.append(f">Contig{counter}_{sample}-{key.replace('_', '-')}\n")
            lines.append(str(records[key].seq) + "\n")
            counter += 1
    with open(file, "w") as fasta_out:
        fasta_out.writelines(lines)
    logger.info("\tContig renaming for %s completed", sample)


@log_exceptions
def clean(path: str) -> None:
    """
    Remove temporary files in the given path.
    """
    logger.info("Removing temporary files in %s...", path)
    for file in glob(os.path.join(path, "*.fasta")):
        os.remove(file)
    for file in glob(os.path.join(path, "*.n*")):
        os.remove(file)
    for file in glob(os.path.join(path, "reference_in*")):
        os.remove(file)
    logger.info("Temporary files removed.")


# -----------------------------------------------------------------------------
# Main Workflow Function
# -----------------------------------------------------------------------------
@log_exceptions
def retrieve(
    data_folder: str,
    collect: bool,
    probe_exons: str,
    num_cores: int,
    length_cover: int,
    spades_cover: float,
    log_queue,
) -> None:
    """
    Execute the complete retrieval process:
      - Optionally collect raw contigs (if 'collect' is True).
      - Prepare contigs (copy and rename) from 30raw_contigs to 31exonic_contigs.
      - Create BLAST hit tables and correct contigs in parallel.
      - Write statistics and rename contigs.
      - Clean up temporary files.
    """
    logger.info("Retrieving data...")
    raw_dir = os.path.join(data_folder, "30raw_contigs")
    if not collect and not os.path.isdir(raw_dir):
        logger.error("No raw contigs found. Run cast_retrieve with -c specified.")
        raise RuntimeError("No raw contigs found. Run cast_retrieve with -c specified.")
    elif collect:
        collect_contigs(data_folder)
    exonic_dir = os.path.join(data_folder, "31exonic_contigs")
    os.makedirs(exonic_dir, exist_ok=True)
    logger.info("Preparing contigs from raw data...")
    for file in glob(os.path.join(raw_dir, "*contigs.fasta")):
        target_file = os.path.join(exonic_dir, os.path.basename(file))
        prepare_contigs(file, target_file)
    logger.info("Contig preparation complete.\n")

    logger.info("Creating hit tables and correcting contigs...")
    files = glob(os.path.join(exonic_dir, "*.fasta"))
    # Create BLAST hit tables in parallel.
    args_blast = [(f, probe_exons, num_cores, length_cover) for f in files]
    with multiprocessing.Pool(
        processes=num_cores, initializer=worker_initializer, initargs=(log_queue,)
    ) as pool:
        pool.starmap(create_hit_tables, args_blast)
    # Correct contigs in parallel.
    args_slice = [(f, spades_cover) for f in files]
    with multiprocessing.Pool(
        processes=num_cores, initializer=worker_initializer, initargs=(log_queue,)
    ) as pool:
        results = pool.starmap(correct_contgis, args_slice)
        all_hits = pandas.concat(results)
    all_hits.reset_index(drop=True, inplace=True)
    all_hits = all_hits.drop(columns=["contig_locus"], errors="ignore")
    all_hits.to_csv(os.path.join(exonic_dir, "all_hits.tsv"), sep="\t", index=False)
    logger.info("Hit table and contig correction complete.\n")

    logger.info("Writing statistics...")
    stats = copies_stats(all_hits)
    stats.to_csv(os.path.join(exonic_dir, "statistics.tsv"), sep="\t", na_rep="NaN")
    logger.info("Statistics file created.\n")

    logger.info("Renaming contigs...")
    for file in glob(os.path.join(exonic_dir, "*.fas")):
        rename_contigs(file)
    logger.info("Contig renaming complete.\n")

    clean(exonic_dir)
    logger.info("Data was successfully retrieved!")


# -----------------------------------------------------------------------------
# End of Module
# -----------------------------------------------------------------------------
