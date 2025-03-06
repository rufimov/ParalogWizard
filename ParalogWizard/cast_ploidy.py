"""
This module contains the functions to estimate the ploidy of the samples. It uses the nQuire package to calculate
likelihoods of different ploidy levels for each sample. The ploidy level with the highest likelihood is then assigned to
the sample.
"""
import multiprocessing
import os.path
import shutil
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


def nquire_create(sample, data_folder):
    for folder_name in os.listdir(os.path.join(data_folder, "100remapped")):
        folder_path = os.path.join(
            os.path.join(data_folder, "100remapped"), folder_name
        )
        if os.path.isdir(folder_path) and folder_name.startswith("exons"):
            bam_file_name = f"{sample}_filtered_uniq_sorted.bam"
            bam_file_path = os.path.join(folder_path, bam_file_name)
            new_file_name = f"{bam_file_name.replace('.bam', '')}_{folder_name}.bam"
            new_file_path = os.path.join(
                os.path.join(data_folder, "102ploidy"), new_file_name
            )
            shutil.copy(bam_file_path, new_file_path)
    subprocess.run(
        f"samtools merge -o {os.path.join(data_folder, '102ploidy', sample+'_filtered_uniq_sorted.bam')} "
        f"{os.path.join(data_folder, '102ploidy', sample+'_filtered_uniq_sorted_exons*.bam')} ",
        shell=True,
    )
    subprocess.run(
        f"samtools index {os.path.join(data_folder, '102ploidy', sample+'_filtered_uniq_sorted.bam')}",
        shell=True,
    )
    subprocess.run(
        f"nQuire create -b {os.path.join(data_folder, '102ploidy', sample)}_filtered_uniq_sorted.bam -o {sample}\n"
        f"nQuire denoise -o {sample}_denoised {sample}.bin\n"
        f"mv {sample}.bin {os.path.join(data_folder, '102ploidy')}\n"
        f"mv {sample}_denoised.bin {os.path.join(data_folder, '102ploidy')}",
        shell=True,
    )


def ploidy(data_folder, num_cores):
    with open(
        os.path.join(data_folder, "10deduplicated_reads", "samples_list.txt")
    ) as samples:
        sample_list = samples.read().splitlines()

    args = list(zip(sample_list, [data_folder] * len(sample_list)))

    with multiprocessing.Pool(processes=num_cores) as pool_nquire_create:
        pool_nquire_create.starmap(nquire_create, args, chunksize=1)

    denoised_bin_list = [
        f"{os.path.join(data_folder, '102ploidy', sample)}_denoised.bin"
        for sample in sample_list
    ]

    nquire_result = os.popen(
        f"nQuire lrdmodel -t {num_cores} {' '.join(denoised_bin_list)}"
    ).read()

    with open(
        os.path.join(data_folder, "102ploidy", "lrdmodel.tsv"), "w"
    ) as file_to_write:
        file_to_write.write(nquire_result)

    lrdmodel = pd.read_csv(
        os.path.join(data_folder, "102ploidy", "lrdmodel.tsv"), delimiter="\t"
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
        os.path.join(data_folder, "102ploidy", "lrdmodel.tsv"), sep="\t", index=False
    )
