import glob
import shutil
import subprocess
import sys
from typing import Dict

import Bio
import pandas
import os

from Bio import SeqIO, SeqRecord
from Bio.Align.Applications import MafftCommandline


def mafft_align_adjust_direction(file):
    stdout, stderr = MafftCommandline(
        input=file, auto=True, adjustdirectionaccurately=True
    )()
    with open(f"{os.path.splitext(file)[0]}.fasta.mafft", "w") as aligned:
        aligned.write(stdout)


def exonerate(
    genes,
    basename,
    run_dir,
    exonerate_genefilename,
    cpu=None,
    thresh=55,
    depth_multiplier=0,
    length_pct=100,
    timeout=None,
):
    if len(genes) == 0:
        print(("ERROR: No genes recovered for {}!".format(basename)))
        return 1

    if os.path.isfile("genes_with_seqs.txt"):
        os.remove("genes_with_seqs.txt")

    print(("Running Exonerate to generate sequences for {} genes".format(len(genes))))

    parallel_cmd_list = ["sudo time parallel", "--eta"]
    if cpu:
        parallel_cmd_list.append("-j {}".format(cpu))
    if timeout:
        parallel_cmd_list.append("--timeout {}%".format(timeout))

    exonerate_cmd_list = [
        "python3",
        "{}/exonerate_hits.py".format(run_dir),
        "{}/{}_baits.fasta",
        "{{}}/{{}}_{}".format("contigs.fasta"),
        "--prefix {{}}/{}".format(basename),
        "-t {}".format(thresh),
        "--depth_multiplier {}".format(depth_multiplier),
        "--length_pct {}".format(length_pct),
        "::::",
        exonerate_genefilename,
        "> genes_with_seqs.txt",
    ]

    exonerate_cmd = " ".join(parallel_cmd_list) + " " + " ".join(exonerate_cmd_list)
    print(exonerate_cmd)
    exitcode = subprocess.call(exonerate_cmd, shell=True)
    if exitcode:
        print("ERROR: Something went wrong with Exonerate!")
        return exitcode
    return


def extend(data_folder, baitfile, n_cpu):
    all_paralogs_for_reference_cleaned = pandas.read_csv(
        os.path.join(
            data_folder, "41detected_par", "all_paralogs_for_reference_cleaned.tsv"
        ),
        sep="\t",
    )

    folder_21 = os.path.join(data_folder, "21supercontigs")
    working_directory = os.path.abspath(os.getcwd())
    os.makedirs(folder_21, exist_ok=True)

    loci = (
        (
            all_paralogs_for_reference_cleaned["locus"].astype(str)
            + all_paralogs_for_reference_cleaned["copy"]
        )
        .drop_duplicates()
        .str.replace("main", "")
        .values.tolist()
    )
    samples = (
        all_paralogs_for_reference_cleaned["sample"].drop_duplicates().values.tolist()
    )

    baits = SeqIO.to_dict(SeqIO.parse(baitfile, "fasta"))

    baits_w_para_to_write = []

    for locus in loci:
        if "para" in locus:
            locus_wa_para = locus[:-4]
            seqs_for_locus = [
                value for key, value in baits.items() if locus_wa_para in key
            ]
            for seq_for_locus in seqs_for_locus:
                seq_id = seq_for_locus.id + "para"
                seq_to_write_w_para = SeqRecord.SeqRecord(
                    seq_for_locus.seq, id=seq_id, name=seq_id, description=seq_id
                )
                baits_w_para_to_write.append(seq_to_write_w_para)
        else:
            seqs_for_locus = [value for key, value in baits.items() if locus in key]
            for seq_for_locus in seqs_for_locus:
                seq_id = seq_for_locus.id
                seq_to_write_wo_para = SeqRecord.SeqRecord(
                    seq_for_locus.seq, id=seq_id, name=seq_id, description=seq_id
                )
                baits_w_para_to_write.append(seq_to_write_wo_para)

    with open(
        f"{baitfile.split('.')[0]}_w_para_for_supercontigs.fasta", "w"
    ) as output_handle:
        SeqIO.write(baits_w_para_to_write, output_handle, "fasta")

    for sample in samples:
        os.makedirs(os.path.join(folder_21, sample), exist_ok=True)
        sequences_contigs_file: Dict[str, Bio.SeqRecord.SeqRecord] = SeqIO.to_dict(
            SeqIO.parse(
                os.path.join(data_folder, "30raw_contigs", f"{sample}_contigs.fasta"),
                "fasta",
            )
        )
        dataframe_for_sample = all_paralogs_for_reference_cleaned[
            all_paralogs_for_reference_cleaned["sample"] == sample
        ]
        for locus in loci:
            os.makedirs(os.path.join(folder_21, sample, locus), exist_ok=True)
            if "para" in locus:
                locus_wa_para = locus[:-4]
                contigs = (
                    dataframe_for_sample[
                        (dataframe_for_sample["locus"].astype(str) == locus_wa_para)
                        & (dataframe_for_sample["copy"] == "para")
                    ]["contig"]
                    .drop_duplicates()
                    .str.replace("c", "cov")
                    .values.tolist()
                )
                locus_baits = []
                for key in baits.keys():
                    if locus_wa_para in key:
                        translated_locus_baits = baits[key].translate(
                            id=f"{key}para", description=f"{key}para"
                        )
                        locus_baits.append(translated_locus_baits)
                with open(
                    os.path.join(folder_21, sample, locus, f"{locus}_baits.fasta"), "w"
                ) as output_handle:
                    SeqIO.write(locus_baits, output_handle, "fasta")
            else:
                contigs = (
                    dataframe_for_sample[
                        (dataframe_for_sample["locus"].astype(str) == locus)
                        & (dataframe_for_sample["copy"] == "main")
                    ]["contig"]
                    .drop_duplicates()
                    .str.replace("c", "cov")
                    .values.tolist()
                )
                locus_baits = []
                for key in baits.keys():
                    if locus in key:
                        translated_locus_baits = baits[key].translate(
                            id=key, description=key
                        )
                        locus_baits.append(translated_locus_baits)
                with open(
                    os.path.join(folder_21, sample, locus, f"{locus}_baits.fasta"), "w"
                ) as output_handle:
                    SeqIO.write(locus_baits, output_handle, "fasta")
            sequences_to_write = []
            for contig in contigs:
                contig = (
                    f"{contig.split('_')[0]}_length_{'_'.join(contig.split('_')[1:])}"
                )
                contig_seq = [
                    value
                    for key, value in sequences_contigs_file.items()
                    if contig in key
                ]
                contig_seq_id = contig_seq[0].id
                correct_seq_id = "_".join(contig_seq_id.split("_")[1:])
                contig_seq[0].id = correct_seq_id
                contig_seq[0].description = correct_seq_id
                contig_seq[0].name = correct_seq_id
                sequences_to_write.append(contig_seq[0])
            with open(
                os.path.join(folder_21, sample, locus, f"{locus}_contigs.fasta"), "w"
            ) as output_handle:
                SeqIO.write(sequences_to_write, output_handle, "fasta")
        os.chdir(os.path.join(folder_21, sample))
        exonerate_genefilename = "exonerate_genelist.txt"
        with open(os.path.join(exonerate_genefilename), "w") as file_to_write:
            for gene in loci:
                file_to_write.write(gene + "\n")
        exitcode = exonerate(
            loci,
            sample,
            os.path.join(working_directory, "ParalogWizard"),
            exonerate_genefilename,
            cpu=n_cpu,
        )
        # if exitcode:
        #     return
        sys.stderr.write(
            "Generated sequences from {} genes!\n".format(
                len(open("genes_with_seqs.txt").readlines())
            )
        )
        os.chdir(os.path.join("../"))
        path_to_intronerate = os.path.join(
            working_directory, "ParalogWizard", "intronerate.py"
        )
        subprocess.call(f"python3 {path_to_intronerate} --prefix {sample}", shell=True)
        os.chdir(working_directory)

    os.chdir(folder_21)
    subprocess.call(
        f"python3 {os.path.join(working_directory,'ParalogWizard', 'retrieve_sequences.py')} "
        f"{os.path.join(working_directory, baitfile.split('.')[0])}_w_para_for_supercontigs.fasta . supercontig",
        shell=True,
    )
    os.chdir(working_directory)
    os.makedirs(os.path.join(folder_21, "supercontig_aln"), exist_ok=True)

    for file in glob.glob(os.path.join(folder_21, "*supercontig.fasta")):
        mafft_align_adjust_direction(file)
        shutil.move(
            file, os.path.join(folder_21, "supercontig_aln", os.path.basename(file))
        )
        shutil.move(
            f"{file}.mafft",
            os.path.join(
                folder_21, "supercontig_aln", f"{os.path.basename(file)}.mafft"
            ),
        )
