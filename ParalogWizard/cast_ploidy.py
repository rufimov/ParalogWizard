"""
This module contains the functions to estimate the ploidy of the samples. It uses the nQuire package to calculate
likelihoods of different ploidy levels for each sample. The ploidy level with the highest likelihood is then assigned to
the sample.
"""
import multiprocessing
import os.path
import subprocess

import numpy as np
import pandas
import pandas as pd
from Bio import SeqIO
from Bio import SeqRecord
from Bio.Align.Applications import MafftCommandline
from Bio.Sequencing.Applications import BwaIndexCommandline
from Bio.Sequencing.Applications import BwaMemCommandline
from Bio.Sequencing.Applications import SamtoolsIndexCommandline
from Bio.Sequencing.Applications import SamtoolsVersion1xSortCommandline
from Bio.Sequencing.Applications import SamtoolsViewCommandline


def trim_ends(seq1, seq2):
    """
    Trims the ends of two sequences to remove any gaps
    :param seq1: sequence one
    :param seq2: sequence two
    """
    while seq1[0] == "-" or seq2[0] == "-":
        seq1 = seq1[1:]
        seq2 = seq2[1:]
    while seq1[-1] == "-" or seq2[-1] == "-":
        seq1 = seq1[:-1]
        seq2 = seq2[:-1]
    return seq1, seq2


def mafft_align(file):
    """
    Aligns a fasta file using mafft
    :param file: a fasta file
    :return: a file with the same name as the input file with .mafft.fasta appended
    """
    stdout, stderr = MafftCommandline(
        input=file,
        auto=True,
    )()
    with open(f"{os.path.splitext(file)[0]}.mafft.fasta", "w") as aligned:
        aligned.write(stdout)


def bwa_index(reference_to_map_to):
    """
    Indexes a reference using BWA
    :param reference_to_map_to: the reference to index
    """
    BwaIndexCommandline(infile=reference_to_map_to, a="is")()


def bwa_map(sample: str, main_data_folder: str, n_threads=1):
    output_path = os.path.join(main_data_folder, "ploidy")
    reference_to_map_to = os.path.join(
        main_data_folder,
        "ploidy",
        f"reference_exons_prepared.fas",
    )
    stdout, stderr = BwaMemCommandline(
        t=n_threads,
        reference=reference_to_map_to,
        read_file1=os.path.join(
            main_data_folder, "10deduplicated_reads", sample + ".R1.fastq.gz"
        ),
        read_file2=os.path.join(
            main_data_folder, "10deduplicated_reads", sample + ".R2.fastq.gz"
        ),
    )()
    with open(os.path.join(output_path, f"{sample}.sam"), "w") as mapped_sample:
        mapped_sample.write(stdout)
    SamtoolsViewCommandline(
        input_file=os.path.join(output_path, f"{sample}.sam"),
        h=True,
        b=True,
        o=os.path.join(output_path, f"{sample}.bam"),
    )()
    SamtoolsViewCommandline(
        input_file=os.path.join(output_path, f"{sample}.bam"),
        h=True,
        b=True,
        q=3,
        o=os.path.join(output_path, f"{sample}_filtered.bam"),
    )()
    SamtoolsViewCommandline(
        input_file=os.path.join(output_path, f"{sample}_filtered.bam"),
        h=True,
        b=True,
        # f=3,
        F=0x90C,
        o=os.path.join(output_path, f"{sample}_filtered_uniq.bam"),
    )()
    SamtoolsVersion1xSortCommandline(
        input=os.path.join(output_path, f"{sample}_filtered_uniq.bam"),
        o=os.path.join(output_path, f"{sample}_filtered_uniq_sorted.bam"),
    )()
    SamtoolsIndexCommandline(
        input_bam=os.path.join(output_path, f"{sample}_filtered_uniq_sorted.bam")
    )()


def prepare_reference(reference_dict, data_folder, exons_allowed):
    exon_count = 1
    to_write = ""
    for item in exons_allowed:
        exon = item[0]
        if item[1]:
            exon_para = f"{exon.split('_')[0] + 'para'}_{'_'.join(exon.split('_')[1:])}"
            if (exon in reference_dict.keys()) and (exon_para in reference_dict.keys()):
                seq = reference_dict[exon].seq
                seq_para = reference_dict[exon_para].seq
                with open(f"{exon}.fasta", "w") as file:
                    SeqIO.write(
                        SeqRecord.SeqRecord(seq, exon, exon, exon),
                        file,
                        "fasta-2line",
                    )
                    SeqIO.write(
                        SeqRecord.SeqRecord(seq_para, exon_para, exon_para, exon_para),
                        file,
                        "fasta-2line",
                    )
                mafft_align(f"{exon}.fasta")
                fasta = list(SeqIO.parse(f"{exon}.mafft.fasta", "fasta"))
                os.remove(f"{exon}.fasta")
                os.remove(f"{exon}.mafft.fasta")
                seq = fasta[0].seq
                seq_para = fasta[1].seq
                seq_trimmed, seq_para_trimmed = trim_ends(seq, seq_para)
                if exon_count == 1:
                    to_write = (
                        to_write
                        + seq_trimmed.replace("-", "")
                        + "N" * 200
                        + seq_para_trimmed.replace("-", "")
                    )
                else:
                    to_write = (
                        to_write
                        + "N" * 200
                        + seq_trimmed.replace("-", "")
                        + "N" * 200
                        + seq_para_trimmed.replace("-", "")
                    )
        else:
            seq = reference_dict[exon].seq
            if exon_count == 1:
                to_write = to_write + seq
            else:
                to_write = to_write + "N" * 200 + seq
        exon_count += 1
    seq_record_to_write = SeqRecord.SeqRecord(
        to_write, "exons_prepared", "exons_prepared", "exons_prepared"
    )
    os.makedirs(os.path.join(data_folder, "ploidy"), exist_ok=True)
    with open(
        os.path.join(
            data_folder,
            "ploidy",
            f"reference_exons_prepared.fas",
        ),
        "w",
    ) as exon_ref_file:
        SeqIO.write(seq_record_to_write, exon_ref_file, "fasta")


def nquire_create(sample, data_folder):
    subprocess.run(
        f"nQuire create -b {os.path.join(data_folder, 'ploidy', sample)}_filtered_uniq_sorted.bam -o {sample}\n"
        f"nQuire denoise -o {sample}_denoised {sample}.bin\n"
        f"mv {sample}.bin {os.path.join(data_folder, 'ploidy')}\n"
        f"mv {sample}_denoised.bin {os.path.join(data_folder, 'ploidy')}",
        shell=True,
    )


def ploidy(reference, data_folder, num_cores, read_length):
    with open(
        os.path.join(data_folder, "10deduplicated_reads", "samples_list.txt")
    ) as samples:
        sample_list = samples.read().splitlines()
    reference_dict = SeqIO.to_dict(SeqIO.parse(reference, "fasta"))
    reference_dict_corrected = dict()
    for exon in reference_dict.keys():
        exon_short = "_".join(exon.split("_")[1:4]).replace("Contig", "exon")
        sequence = reference_dict[exon].seq.ungap("-")
        seq_record = reference_dict[exon]
        seq_record.seq = sequence
        exon_length = len(sequence)
        if exon_length >= read_length:
            reference_dict_corrected[exon_short] = seq_record

    exons_in_ref = reference_dict_corrected.keys()
    all_paralogs_for_reference = pandas.read_csv(
        os.path.join(data_folder, "41detected_par", "all_paralogs_for_reference.tsv"),
        sep="\t",
    )
    from ParalogWizard.cast_detect import exon_stats

    all_paralogs_for_reference_exon_stat = exon_stats(all_paralogs_for_reference)
    all_paralogs_for_reference_exon_stat.to_csv(
        os.path.join(
            data_folder,
            "41detected_par",
            "exon_statistics.tsv",
        ),
        na_rep="NaN",
        sep="\t",
    )
    allowed_exons = []
    for column in all_paralogs_for_reference_exon_stat:
        num_yes = 0
        num_no = 0
        try:
            num_yes = all_paralogs_for_reference_exon_stat[column].value_counts()["Yes"]
        except KeyError:
            try:
                num_no = all_paralogs_for_reference_exon_stat[column].value_counts()[
                    "No"
                ]
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
                f"{column.split('_')[0] + 'para'}_{'_'.join(column.split('_')[1:])}"
                in exons_in_ref
            )
        ):
            allowed_exons.append((column, True))

    prepare_reference(reference_dict_corrected, data_folder, allowed_exons)

    bwa_index(
        os.path.join(
            data_folder,
            "ploidy",
            f"reference_exons_prepared.fas",
        )
    )

    args = list(zip(sample_list, [data_folder] * len(sample_list)))
    with multiprocessing.Pool(processes=num_cores) as pool_bwa:
        pool_bwa.starmap(bwa_map, args, chunksize=1)

    with multiprocessing.Pool(processes=num_cores) as pool_nquire_create:
        pool_nquire_create.starmap(nquire_create, args, chunksize=1)

    denoised_bin_list = [
        f"{os.path.join(data_folder, 'ploidy', sample)}_denoised.bin"
        for sample in sample_list
    ]

    nquire_result = os.popen(
        f"nQuire lrdmodel -t {num_cores} {' '.join(denoised_bin_list)}"
    ).read()

    with open(
        os.path.join(data_folder, "ploidy", "lrdmodel.tsv"), "w"
    ) as file_to_write:
        file_to_write.write(nquire_result)

    lrdmodel = pd.read_csv(
        os.path.join(data_folder, "ploidy", "lrdmodel.tsv"), delimiter="\t"
    )
    lrdmodel["sample"] = sample_list
    conditions = [
        (lrdmodel["dip"] > lrdmodel["tri"]) & (lrdmodel["dip"] > lrdmodel["tet"]),
        (lrdmodel["tri"] > lrdmodel["dip"]) & (lrdmodel["tri"] > lrdmodel["tet"]),
        (lrdmodel["tet"] > lrdmodel["dip"]) & (lrdmodel["tet"] > lrdmodel["tri"]),
    ]
    choices = ["dip", "tri", "tet"]
    lrdmodel["ploidy"] = np.select(conditions, choices, default="unknown")
    lrdmodel.to_csv(
        os.path.join(data_folder, "ploidy", "lrdmodel.tsv"), sep="\t", index=False
    )
