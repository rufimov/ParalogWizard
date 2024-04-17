import fileinput
import logging
import multiprocessing
import os
import re
import shutil
from glob import glob
from pathlib import Path
from typing import Dict, List

import Bio.SeqRecord
import pandas
from Bio import SeqIO
from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastnCommandline
from pandas.core import series, frame


def slicing(
    dictionary: Dict[str, Bio.SeqRecord.SeqRecord], entry: pandas.core.series.Series
) -> str:
    """
    Cuts out sequence of exon from contig given the hit entry from a dataframe.
    :param dictionary: dictionary of SeqRecords
    :param entry: hit entry from a dataframe
    :return: sequence of exon
    """
    sequence = dictionary[entry["saccver"]].seq
    start = entry["sstart"]
    end = entry["send"]
    if start > end:
        result = str(sequence.reverse_complement())[-end : -start - 1 : -1][::-1]
    else:
        result = str(sequence)[start - 1 : end]
    return result


def collect_contigs(data_folder, logger):
    """
    Collects all contigs assembled during the assembly step for each sample and saves them to 30raw_contigs folder.
    :param data_folder: path to the folder with data
    :param logger: logger
    :return: None
    """
    os.makedirs(os.path.join(data_folder, "30raw_contigs"), exist_ok=True)
    folder20 = os.path.join(data_folder, "20assemblies")
    for folder in os.listdir(folder20):
        if os.path.isdir(os.path.join(folder20, folder)):
            logger.info(f"Processing {folder}")
            with open(
                os.path.join(data_folder, "30raw_contigs", f"{folder}_contigs.fasta"),
                "w",
            ) as contigs:
                for locus in os.listdir(os.path.join(folder20, folder)):
                    if os.path.exists(
                        os.path.join(folder20, folder, locus, f"{locus}_contigs.fasta")
                    ):
                        logger.info(f"\tProcessing {locus}")
                        with open(
                            os.path.join(
                                folder20,
                                folder,
                                locus,
                                f"{locus}_contigs.fasta",
                            )
                        ) as locus_contigs:
                            file_content = locus_contigs.read().replace(
                                ">", f">{locus}_"
                            )
                            contigs.write(file_content)
                        logger.info("\tOK")
            if (
                os.stat(
                    os.path.join(
                        data_folder, "30raw_contigs", f"{folder}_contigs.fasta"
                    )
                ).st_size
                == 0
            ):
                os.remove(
                    os.path.join(
                        data_folder, "30raw_contigs", f"{folder}_contigs.fasta"
                    )
                )
            logger.info("Ok")


def prepare_contigs(file_to_open, file_to_write, logger):
    """
    Copy fasta file from folder 30raw_contigs to 31exonic_contigs, rename it,
    correct names of sequences to make it shorter
    """
    shutil.copy(file_to_open, file_to_write)
    with fileinput.FileInput(file_to_write, inplace=True) as fasta:
        for line in fasta:
            corrected_line = re.sub("NODE", "N", line)
            corrected_line = re.sub(
                r"length_([0-9]+)_cov_([0-9]+\.[0-9][0-9]).*",
                r"\1_c_\2",
                corrected_line,
            )
            print(corrected_line, end="")
    name_of_file: str = "_".join(
        os.path.basename(file_to_write).split(".")[0].split("_")[0:2]
    )
    os.rename(
        file_to_write,
        os.path.join(os.path.dirname(file_to_write), f"{name_of_file}.fasta"),
    )


def create_hit_tables(fasta_file, probe_exons, n_cpu, length_cover, log_file):
    """Running blast on every fasta file with contigs. Probe file is blasted against contigs. Blast results are saved in
    text file."""
    logger = create_logger(log_file)
    path = os.path.dirname(fasta_file)
    fasta_file: str = os.path.basename(fasta_file)
    sample: str = Path(fasta_file).stem
    NcbimakeblastdbCommandline(
        dbtype="nucl",
        input_file=os.path.join(path, fasta_file),
        out=os.path.join(path, sample),
        parse_seqids=True,
    )()
    logger.info(f"\t\tCreating hit table for {sample}. Running BLAST...")
    NcbiblastnCommandline(
        task="blastn",
        num_threads=n_cpu,
        query=probe_exons,
        db=os.path.join(path, sample),
        out=os.path.join(path, f"reference_in_{sample}_contigs.txt"),
        qcov_hsp_perc=length_cover,
        outfmt="6 qaccver saccver pident qcovhsp evalue bitscore sstart send",
    )()
    logger.info(f"\t\tHit table for {sample} is ready")


def correct_contgis(file, spades_cover, log_file):
    """
    Filters hit table, keeps only one hit per contig per exon. Slices contigs to exonic contigs according to hit
    boundaries.
    """
    logger = create_logger(log_file)
    logger.info(f"\t\tCorrecting file {os.path.basename(file)}")
    folder = os.path.dirname(file)
    sample: str = Path(file).stem
    blast_results = pandas.read_csv(
        os.path.join(folder, f"reference_in_{sample}_contigs.txt"),
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
    blast_results["k-mer_cover"] = blast_results["saccver"].str.split("_").str[-1]
    blast_results["k-mer_cover"] = blast_results["k-mer_cover"].astype(float)
    hits = blast_results[
        (blast_results["locus"] == blast_results["contig_locus"])
        & (blast_results["k-mer_cover"] >= spades_cover)
    ].reset_index(drop=True)
    hits["exon"] = hits["qaccver"].str.split("-").str[1]
    hits_sorted = hits.sort_values(
        ["exon", "qcovhsp", "pident", "evalue", "bitscore"],
        ascending=[True, False, False, True, False],
    )
    hits_dedup = hits_sorted.drop_duplicates(subset=["exon", "saccver"]).reset_index(
        drop=True
    )
    hits_dedup["sample"] = sample
    with open(os.path.join(folder, f"{sample}.fas"), "w") as result_fasta, open(
        file
    ) as contigs:
        contigs_fasta_parsed = SeqIO.to_dict(SeqIO.parse(contigs, "fasta"))
        for index in range(len(hits_dedup)):
            sequence: str = slicing(contigs_fasta_parsed, hits_dedup.loc[index])
            result_fasta.write(
                f">{hits_dedup['exon'][index]}_{'_'.join(hits_dedup['saccver'][index].split('_')[1:])}"
                f"\n{sequence}\n"
            )
            hits_dedup.loc[index, "sequence"] = sequence
    logger.info("\t\tOK")
    return hits_dedup


def copies_stats(all_hits: pandas.core.frame.DataFrame) -> pandas.core.frame.DataFrame:
    """
    Takes file with all hits and creates table with statistics. For each sample number of possible copies per locus
    is calculated. Number of copies is deduced from the maximum number of hits among all exons for a particular locus.
    """
    statistics: pandas.core.frame.DataFrame = pandas.DataFrame([], columns=["locus"])
    statistics.set_index("locus", drop=True, inplace=True)
    grouped_hits = all_hits.groupby(["sample", "locus"])[["exon", "sequence"]]
    for group in grouped_hits:
        count = group[1].groupby("exon")["sequence"].count().max()
        sample = tuple(group[0])[0]
        locus = tuple(group[0])[1]
        statistics.loc[locus, sample] = count
    return statistics


def rename_contigs(file, logger):
    """
    Rename contigs in each fasta file in accordance to HybPhyloMaker.
    """
    sample = os.path.basename(os.path.splitext(file)[0])
    logger.info(f"\tProcessing {sample}")
    with open(file) as result_fasta:
        fasta_parsed = SeqIO.to_dict(SeqIO.parse(result_fasta, "fasta"))
        counter: int = 1
        fasta_to_write: List[str] = list()
        for line in sorted(fasta_parsed.keys()):
            fasta_to_write.append(
                f">Contig{str(counter)}_{sample}-{line.replace('_', '-')}\n"
            )
            fasta_to_write.append(str(fasta_parsed[line].seq) + "\n")
            counter += 1
    with open(file, "w") as result_fasta:
        result_fasta.writelines(fasta_to_write)
    logger.info("\tOK")


def clean(path, logger):
    """
    Removes unnecessary files: blast databes files, uncorrected contigs and raw hit tables.
    """
    logger.info("Removing temporary files...")
    for file in glob(os.path.join(path, "*.fasta")):
        os.remove(file)
    for file in glob(os.path.join(path, "*.n*")):
        os.remove(file)
    for file in glob(os.path.join(path, "reference_in*")):
        os.remove(file)
    logger.info("Done\n")


def create_logger(log_file):
    logger = multiprocessing.get_logger()
    logger.setLevel(logging.INFO)
    log_handler_info = logging.FileHandler(log_file)
    log_handler_info.setLevel(logging.INFO)
    log_formatter_info = logging.Formatter(
        fmt="%(asctime)s - %(message)s", datefmt="%d-%b-%y %H:%M:%S"
    )
    log_handler_info.setFormatter(log_formatter_info)
    logger.addHandler(log_handler_info)
    return logger


def retrieve(
    data_folder, collect, probe_exons, num_cores, length_cover, spades_cover, log_file
):
    """
    Executes everything given the arguments and log file.
    """
    logger = create_logger(log_file)
    logger.info("Retrieving data...\n")
    if not collect and not os.path.isdir(os.path.join(data_folder, "30raw_contigs")):
        logger.info(
            "ERROR: No raw contigs found. Run ParalogWizard cast_retrieve with -c specified."
        )
        exit(1)
    elif collect:
        collect_contigs(data_folder, logger)
    os.makedirs(os.path.join(data_folder, "31exonic_contigs"), exist_ok=True)
    folder_30 = os.path.join(data_folder, "30raw_contigs")
    folder_31 = os.path.join(data_folder, "31exonic_contigs")
    logger.info("Preparing congits...")
    for file in glob(os.path.join(folder_30, "*contigs.fasta")):
        prepare_contigs(
            file,
            os.path.join(folder_31, os.path.basename(file)),
            logger,
        )
    logger.info("Done\n")
    logger.info("Creating hit tables and correcting contigs...")
    files_to_process = []
    for file in glob(os.path.join(folder_31, "*.fasta")):
        files_to_process.append(file)
    args_blast = list(
        zip(
            files_to_process,
            [probe_exons] * (len(files_to_process)),
            [num_cores] * (len(files_to_process)),
            [length_cover] * (len(files_to_process)),
            [log_file] * (len(files_to_process)),
        )
    )
    args_slice = list(
        zip(
            files_to_process,
            [spades_cover] * (len(files_to_process)),
            [log_file] * (len(files_to_process)),
        )
    )
    with multiprocessing.Pool(processes=num_cores) as pool_blast:
        pool_blast.starmap(create_hit_tables, args_blast)
    with multiprocessing.Pool(processes=num_cores) as pool_slice_contig:
        results = pool_slice_contig.starmap(correct_contgis, args_slice)
        all_hits_for_reference = pandas.concat(results)
    all_hits_for_reference.reset_index(drop=True, inplace=True)
    all_hits_for_reference: pandas.core.frame.DataFrame = all_hits_for_reference.drop(
        columns=["contig_locus"]
    )
    all_hits_for_reference.to_csv(
        os.path.join(folder_31, "all_hits.tsv"), sep="\t", index=False
    )
    logger.info("Done\n")
    logger.info("Writing statistics...")
    stats = copies_stats(all_hits_for_reference)
    stats.to_csv(
        os.path.join(folder_31, "statistics.tsv"),
        sep="\t",
        na_rep="NaN",
    )
    logger.info("Statistics file created\n")
    logger.info("Renaming contigs...")
    for file in glob(os.path.join(folder_31, "*.fas")):
        rename_contigs(file, logger)
    logger.info("All contigs were successfully renamed!\n")
    clean(folder_31, logger)
    logger.info("Data was successfully retrieved!")
