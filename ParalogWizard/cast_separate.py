from glob import glob
import logging
import multiprocessing
import os
import re
import shutil
import subprocess
import fileinput

import pandas

from ParalogWizard.cast_analyze import mafft_align
from Bio import SeqIO


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


def run_blat(contigfile, probes, minident, log_file):
    """"""
    logger = create_logger(log_file)
    file = os.path.basename(contigfile)
    data_folder = os.path.dirname(os.path.dirname(contigfile))
    if file != probes:
        logger.info(f"Processing {file}...")
        blat_cmd_output = os.popen(
            f"blat -t=DNA -q=DNA -out=pslx -minIdentity={minident} {probes} {contigfile} \
            {os.path.join(data_folder, '50pslx', f'{os.path.basename(contigfile)}.pslx')}"
        ).read()
        logger.info(blat_cmd_output)
        logger.info("Done")


def replace_trailing(seq):
    for id in range(len(seq)):
        symbol = seq[id]
        if symbol == '-':
            seq = seq[:id] + '?' + seq[id + 1:]
        else:
            break
    for id in reversed(range(len(seq))):
        symbol = seq[id]
        if symbol == '-':
            seq = seq[:id] + '?' + seq[id + 1:]
        else:
            break
    return seq


def align(data_folder, probes, n_cpu, log_file):
    logger = create_logger(log_file)
    shutil.rmtree(os.path.join(data_folder, "60mafft"), ignore_errors=True)
    with open(
        os.path.join(data_folder, "50pslx", "corrected", "list_pslx.txt"), "w"
    ) as list_pslx:
        for pslx_file in glob(
            os.path.join(data_folder, "50pslx", "corrected", "*.pslx")
        ):
            list_pslx.write(pslx_file + "\n")
    #     subprocess.call(
    #         f"python3 ParalogWizard/assembled_exons_to_fastas.py \
    # -l {os.path.join(data_folder, '50pslx', 'corrected', 'list_pslx.txt')} -f {probes} \
    # -d {os.path.join(data_folder, '60mafft')}",
    #         shell=True,
    #     )

    exons_to_fastas_output = os.popen(
        f"python3 ParalogWizard/assembled_exons_to_fastas.py \
-l {os.path.join(data_folder, '50pslx', 'corrected', 'list_pslx.txt')} -f {probes} \
-d {os.path.join(data_folder, '60mafft')}"
    ).read()
    logger.info(exons_to_fastas_output)
    all_loci = set()
    files_to_align = []
    for file in glob(os.path.join(data_folder, "60mafft", "*.fasta")):
        with fileinput.FileInput(file, inplace=True) as file_to_correct:
            for line in file_to_correct:
                line = re.sub(r">.+/", ">", line)
                line = re.sub(r"\.fas", "_contigs.fas", line)
                if not line.startswith(">"):
                    line = re.sub(r"[nN]", "-", line)
                print(line, end="")
        locus = os.path.basename(file).split("_")[3]
        all_loci.add(locus)
        files_to_align.append(file)
    with multiprocessing.Pool(processes=n_cpu) as pool_aln:
        pool_aln.map(mafft_align, files_to_align)
    for file in glob(os.path.join(data_folder, "60mafft", "*.mafft")):
        fasta = list(SeqIO.parse(file, 'fasta'))
        SeqIO.write(fasta, file, "fasta-2line")
        with fileinput.FileInput(file, inplace=True) as file_to_correct:
            for line in file_to_correct:
                if not line.startswith(">"):
                    line = replace_trailing(line)
                print(line, end="")
    os.makedirs(
        os.path.join(data_folder, "70concatenated_exon_alignments"), exist_ok=True
    )
    amas_ex = os.path.abspath("ParalogWizard/AMAS.py")
    os.chdir(os.path.join(data_folder, "70concatenated_exon_alignments"))
    for locus in all_loci:
        exons_to_concat = glob(
            os.path.join("..", "60mafft", f"To_align_Assembly_{locus}_*.fasta.mafft")
        )
        exons_to_concat.sort(key=lambda x: int(x.split("_")[5]))
        line_exons_to_concat = " ".join(exons_to_concat)
        # subprocess.call(
        #     f"python3 {amas_ex} concat -i {line_loci_to_concat} -f fasta -d dna -t Assembly_{locus}.fasta "
        #     f"-p Assembly_{locus}.part",
        #     shell=True,
        # )
        amas_concat_output = os.popen(
            f"python3 {amas_ex} concat -i {line_exons_to_concat} -f fasta -d dna -t Assembly_{locus}.fasta "
            f"-p Assembly_{locus}.part"
        ).read()
        logger.info(amas_concat_output)
        # subprocess.call(
        #     f"python3 {amas_ex} convert -i Assembly_{locus}.fasta -f fasta -d dna -u phylip ",
        #     shell=True,
        # )
        amas_convert_output = os.popen(
            f"python3 {amas_ex} convert -i Assembly_{locus}.fasta -f fasta -d dna -u phylip "
        ).read()
        logger.info(amas_convert_output)

        os.rename(f"Assembly_{locus}.fasta-out.phy", f"Assembly_{locus}.phylip")
        with fileinput.FileInput(f"Assembly_{locus}.part", inplace=True) as part_file:
            for line in part_file:
                corrected_line = re.sub(r"^", "DNA, ", line)
                corrected_line = re.sub(r"_To_align", "", corrected_line)                
                print(corrected_line, end="")
    os.chdir(os.path.dirname(os.getcwd()))


def correct_pslx(pslx_file, log_file):
    """"""
    logger = create_logger(log_file)
    foler50 = os.path.dirname(pslx_file)
    file = os.path.basename(pslx_file)
    with open(pslx_file) as original_pslx_file:
        pslx_file_as_list = original_pslx_file.read().splitlines()
    head = pslx_file_as_list[0:5]
    columns = [
        "match",
        "mismatch",
        "rep_match",
        "Ns",
        "Q_gap_count",
        "Q_gap_bases",
        "T_gap_bases",
        "T_gap_count",
        "strand",
        "Q_name",
        "Q_size",
        "Q_start",
        "Q_end",
        "T_name",
        "T_size",
        "T_start",
        "T_end",
        "block_end",
        "blockSizes",
        "qStarts",
        "tStarts",
        "seq1",
        "seq2",
    ]
    data = [line.split() for line in pslx_file_as_list[5:]]
    pslx_file_dataframe = pandas.DataFrame(data, columns=columns)
    pslx_file_dataframe["match"] = pslx_file_dataframe["match"].astype(int)
    pslx_file_dataframe["mismatch"] = pslx_file_dataframe["mismatch"].astype(int)

    pslx_file_dataframe["similarity"] = pslx_file_dataframe["match"] / (
        pslx_file_dataframe["match"] + pslx_file_dataframe["mismatch"]
    )
    pslx_file_dataframe.sort_values(
        ["T_name", "similarity"], ascending=[True, False], inplace=True
    )
    pslx_file_dataframe.drop_duplicates("T_name", inplace=True)
    pslx_file_dataframe.sort_values(
        ["Q_name", "similarity"], ascending=[True, False], inplace=True
    )
    pslx_file_dataframe.drop_duplicates("Q_name", inplace=True)

    pslx_file_dataframe.drop("similarity", axis=1, inplace=True)
    with open(os.path.join(foler50, "corrected", file), "w") as corrected_pslx:
        for line in head:
            corrected_pslx.write(line + "\n")
    pslx_file_dataframe.to_csv(
        os.path.join(foler50, "corrected", file),
        mode="a",
        header=False,
        sep="\t",
        index=False,
    )


def generate_pslx(data_folder, probes, minident, redlist, num_cores, log_file):
    logger = create_logger(log_file)
    logger.info("Generating pslx files using BLAT...\n")
    os.makedirs(os.path.join(data_folder, "50pslx"), exist_ok=True)
    files_for_blat_list = glob(os.path.join(data_folder, "31exonic_contigs", "*.fas"))
    args_blat = list(
        zip(
            files_for_blat_list,
            [probes] * len(files_for_blat_list),
            [minident] * len(files_for_blat_list),
            [log_file] * len(files_for_blat_list),
        )
    )
    with multiprocessing.Pool(processes=num_cores) as pool_blat:
        pool_blat.starmap(run_blat, args_blat)
    os.makedirs(os.path.join(data_folder, "50pslx", "corrected"), exist_ok=True)
    pslx_file_list = []
    for pslx_file in glob(os.path.join(data_folder, "50pslx", "*.pslx")):
        file = os.path.basename(pslx_file)
        sample = file.split(".")[0]
        if sample not in redlist:
            pslx_file_list.append(pslx_file)
        else:
            shutil.copyfile(
                pslx_file, os.path.join(data_folder, "50pslx", "corrected", file)
            )
    args_correct_pslx = list(zip(pslx_file_list, [log_file] * len(pslx_file_list)))
    with multiprocessing.Pool(processes=num_cores) as pool_correct_pslx:
        pool_correct_pslx.starmap(correct_pslx, args_correct_pslx)
