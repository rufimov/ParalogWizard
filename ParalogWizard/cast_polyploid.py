import fileinput
import functools
import multiprocessing
import os
import re
import shutil
from glob import glob
from typing import Dict

import Bio
import pandas
import sklearn
from Bio import SeqIO
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer

from ParalogWizard.cast_analyze import mafft_align, percent_dissimilarity
from ParalogWizard.cast_retrieve import create_logger
from ParalogWizard.cast_separate import replace_trailing


def correct_pslx(data_folder, pslx_file, logger):
    folder50 = os.path.dirname(pslx_file)
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
        ["Q_name", "similarity"], ascending=[True, False], inplace=True
    )
    pslx_file_dataframe.drop_duplicates("Q_name", inplace=True)

    pslx_file_dataframe.drop("similarity", axis=1, inplace=True)
    with open(os.path.join(folder50, "corrected", file), "w") as corrected_pslx:
        for line in head:
            corrected_pslx.write(line + "\n")
    pslx_file_dataframe.to_csv(
        os.path.join(os.path.join(data_folder, "polyploids", "corrected"), file),
        mode="a",
        header=False,
        sep="\t",
        index=False,
    )


def polyploid(data_folder, probes, n_cpu, log_file):
    logger = create_logger(log_file)
    with open(
        os.path.join(data_folder, "polyploids", "polyploids.txt"), "r"
    ) as polyploids:
        list_polyploids = polyploids.read().splitlines()
    with open(
        os.path.join(data_folder, "10deduplicated_reads", "backbone_list.txt")
    ) as backbone_list:
        backbone = backbone_list.read().splitlines()

    logger.info("Correcting pslx files")
    pslx_to_correct = [
        os.path.join(data_folder, "50pslx", f"{polyploid}.fas.pslx")
        for polyploid in list_polyploids
    ]
    os.makedirs(os.path.join(data_folder, "polyploids", "corrected"), exist_ok=True)
    for pslx_file in pslx_to_correct:
        correct_pslx(data_folder, pslx_file, logger)
    with open(
        os.path.join(data_folder, "polyploids", "polyploids_pslx.txt"), "w"
    ) as pslx_polyploids:
        for polyploid in list_polyploids:
            pslx_polyploids.write(
                os.path.join(
                    data_folder, "polyploids", "corrected", f"{polyploid}.fas.pslx\n"
                )
            )
    shutil.rmtree(
        os.path.join(data_folder, "polyploids", "60mafft_all"), ignore_errors=True
    )
    exons_to_fastas_output = os.popen(
        f"python3 ParalogWizard/assembled_exons_to_fastas.py -a \
    -l {os.path.join(data_folder, 'polyploids', 'polyploids_pslx.txt')} -f {probes} \
    -d {os.path.join(data_folder, 'polyploids', '60mafft_all')}"
    ).read()
    logger.info(exons_to_fastas_output)
    all_loci = set()
    files_to_align = []
    for file in glob(os.path.join(data_folder, "polyploids", "60mafft_all", "*.fasta")):
        with fileinput.FileInput(file, inplace=True) as file_to_correct:
            for line in file_to_correct:
                line = re.sub(r">.+/", ">", line)
                line = re.sub(r"\.fas", "_contigs.fas", line)
                if not line.startswith(">"):
                    line = re.sub(r"[nN]", "-", line)
                print(line, end="")
        sequences: Dict[str, Bio.SeqRecord.SeqRecord] = SeqIO.to_dict(
            SeqIO.parse(file, "fasta")
        )
        sequences_ungap = dict()
        for item in sequences.keys():
            sequence = sequences[item]
            sequence = sequence.seq.replace("-", "")
            if len(sequence) != 0:
                sequences_ungap[item] = sequence
        if len(sequences_ungap) < 1:
            logger.info(
                f"File {file} consists of gaps only. Aligning with mafft skipped."
            )
            os.remove(file)
            continue
        locus = os.path.basename(file).split("_")[3]
        all_loci.add(locus)
        files_to_align.append(file)

    os.makedirs(
        os.path.join(data_folder, "polyploids", "60mafft_backbone"), exist_ok=True
    )

    for file in files_to_align:
        with open(
            os.path.join(data_folder, "60mafft", os.path.basename(file)), "r"
        ) as file_to_read, open(
            os.path.join(
                data_folder,
                "polyploids",
                "60mafft_backbone",
                os.path.basename(file),
            ),
            "w",
        ) as file_to_write:
            seq_dict = SeqIO.to_dict(SeqIO.parse(file_to_read, "fasta"))
            seq_dict_keep = {
                k: seq_dict[k] for k in [f"{x}_contigs.fas" for x in backbone]
            }
            for key in seq_dict_keep.keys():
                file_to_write.write(f">{key}\n{seq_dict_keep[key].seq}\n")
        with open(file, "a") as file_to_append, open(
            os.path.join(
                data_folder,
                "polyploids",
                "60mafft_backbone",
                os.path.basename(file),
            ),
            "r",
        ) as file_to_add:
            addition = file_to_add.read()
            file_to_append.write(addition)

    with multiprocessing.Pool(processes=n_cpu) as pool_aln:
        align_function = functools.partial(mafft_align, log_file=log_file)
        pool_aln.map(align_function, files_to_align)

    for file in glob(os.path.join(data_folder, "polyploids", "60mafft_all", "*.mafft")):
        fasta = list(SeqIO.parse(file, "fasta"))
        SeqIO.write(fasta, file, "fasta-2line")
        with fileinput.FileInput(file, inplace=True) as file_to_correct:
            for line in file_to_correct:
                if not line.startswith(">"):
                    line = replace_trailing(line[:-1]) + "\n"
                print(line, end="")

    for poly in list_polyploids:
        matrix_similarity = pandas.DataFrame()
        for file in glob(
            os.path.join(data_folder, "polyploids", "60mafft_all", "*.mafft")
        ):
            all_seq_dict = SeqIO.to_dict(SeqIO.parse(file, "fasta"))
            backbone_dict = {
                key: all_seq_dict[f"{key}_contigs.fas"] for key in backbone
            }
            list_seqs_polyploid = [
                i for i in all_seq_dict.keys() if (poly in i)
            ]  # In case if polyploid sequences present in backbone alignment
            poly_dict = {key: all_seq_dict[key] for key in list_seqs_polyploid}
            for poly_key in poly_dict.keys():
                poly_seq = poly_dict[poly_key].seq
                for backbone_key in backbone_dict.keys():
                    backbone_seq = backbone_dict[backbone_key].seq
                    if len(str(backbone_seq).replace("?", "")) == 0:
                        dis = None
                    else:
                        dis = percent_dissimilarity(str(poly_seq), str(backbone_seq))
                    if dis is not None:
                        similarity = 100 - dis
                        matrix_similarity.loc[poly_key, backbone_key] = similarity

        # Convert DataFrame to matrix
        mat = matrix_similarity.values

        imp = IterativeImputer()
        imputed = imp.fit_transform(mat)
        # Using sklearn
        km = sklearn.cluster.KMeans(n_clusters=2)
        km.fit(imputed)
        # Get cluster assignment labels
        labels = km.labels_
        # Format results as a DataFrame
        results = pandas.DataFrame([matrix_similarity.index, labels]).T

        os.makedirs(
            os.path.join(data_folder, "polyploids", "60mafft_corr"), exist_ok=True
        )

        for file in glob(
            os.path.join(data_folder, "polyploids", "60mafft_all", "*.mafft")
        ):
            all_seq_dict = SeqIO.to_dict(SeqIO.parse(file, "fasta"))
            all_seq_dict_renamed = {}
            backbone_dict = {
                key: all_seq_dict[f"{key}_contigs.fas"] for key in backbone
            }
            poly_contigs = [i for i in all_seq_dict.keys() if (poly in i)]
            clust_poly_contig = results[results[0].isin(poly_contigs)]
            if clust_poly_contig[1].nunique() > 1:
                clust_poly_contig0 = clust_poly_contig[clust_poly_contig[1] == 0]
                clust_poly_contig1 = clust_poly_contig[clust_poly_contig[1] == 1]
                dict_clust0 = {key: all_seq_dict[key] for key in clust_poly_contig0[0]}
                dict_clust0_renamed = {}
                for name, seq in dict_clust0.items():
                    index = name.find("_")
                    new_name = name[:index] + "SubGenA" + name[index:]
                    seq.id = new_name
                    seq.name = new_name  # Update the shorter name
                    seq.description = ""  # Clear description to avoid redundant information in FASTA header
                    all_seq_dict_renamed[new_name] = seq
                dict_clust1 = {key: all_seq_dict[key] for key in clust_poly_contig1[0]}
                dict_clust1_renamed = {}
                for name, seq in dict_clust1.items():
                    index = name.find("_")
                    new_name = name[:index] + "SubGenB" + name[index:]
                    seq.id = new_name
                    seq.name = new_name  # Update the shorter name
                    seq.description = ""  # Clear description to avoid redundant information in FASTA header
                    all_seq_dict_renamed[new_name] = seq
            elif clust_poly_contig[1].unique().item() == 0:
                dict_clust0 = {key: all_seq_dict[key] for key in clust_poly_contig[0]}
                dict_clust0_renamed = {}
                for name, seq in dict_clust0.items():
                    index = name.find("_")
                    new_name = name[:index] + "SubGenA" + name[index:]
                    seq.id = new_name
                    seq.name = new_name  # Update the shorter name
                    seq.description = ""  # Clear description to avoid redundant information in FASTA header
                    all_seq_dict_renamed[new_name] = seq
            elif clust_poly_contig[1].unique().item() == 1:
                dict_clust1 = {key: all_seq_dict[key] for key in clust_poly_contig[0]}
                dict_clust1_renamed = {}
                for name, seq in dict_clust1.items():
                    index = name.find("_")
                    new_name = name[:index] + "SubGenB" + name[index:]
                    seq.id = new_name
                    seq.name = new_name  # Update the shorter name
                    seq.description = ""  # Clear description to avoid redundant information in FASTA header
                    all_seq_dict_renamed[new_name] = seq
            dict_to_write = all_seq_dict_renamed | backbone_dict
            with open(
                os.path.join(
                    data_folder, "polyploids", "60mafft_corr", os.path.basename(file)
                ),
                "w",
            ) as fasta_file:
                SeqIO.write(dict_to_write.values(), fasta_file, "fasta")
            x = 1
            #
            #
            # x = 1

    x = 1
    y = 1
