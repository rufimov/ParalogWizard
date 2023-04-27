import logging
import multiprocessing
import os
from typing import List
import pandas


def find_cluster(array: List[float]):
    array.sort()
    gaps = []
    for id in range(1, len(array)):
        gaps.append(array[id] - array[id - 1])
    biggest_gap_middle = max(gaps) / 2
    biggest_gap_index = gaps.index(max(gaps))
    boundary = array[biggest_gap_index] + biggest_gap_middle
    return boundary


def adjust_orphaned_main(sample_locus_dataframe):
    copies = sample_locus_dataframe["copy"].values.tolist()
    if "para" not in copies:
        return sample_locus_dataframe
    grouped_exons = sample_locus_dataframe.groupby("exon")
    exons_needed_clustering = []
    for group in grouped_exons:
        if "para" not in group[1]["copy"].values.tolist():
            exons_needed_clustering.append(group[0])
    if not exons_needed_clustering:
        return sample_locus_dataframe
    ident_array = sample_locus_dataframe[
        ~sample_locus_dataframe["exon"].isin(exons_needed_clustering)
    ]["pident"].values.tolist()
    boundary = find_cluster(ident_array)
    pidents_main = sample_locus_dataframe[
        (sample_locus_dataframe["copy"] == "main")
        & (~sample_locus_dataframe["exon"].isin(exons_needed_clustering))
    ]["pident"]
    pidents_para = sample_locus_dataframe[(sample_locus_dataframe["copy"] == "para")][
        "pident"
    ]
    if (pidents_main < boundary).any() or (pidents_para > boundary).any():
        sample_locus_dataframe = sample_locus_dataframe[
            ~sample_locus_dataframe["exon"].isin(exons_needed_clustering)
        ]
        return sample_locus_dataframe
    sample_locus_dataframe.loc[
        (sample_locus_dataframe["pident"] > boundary)
        & (sample_locus_dataframe["exon"].isin(exons_needed_clustering)),
        "copy",
    ] = "main"
    sample_locus_dataframe.loc[
        (sample_locus_dataframe["pident"] < boundary)
        & (sample_locus_dataframe["exon"].isin(exons_needed_clustering)),
        "copy",
    ] = "para"
    return sample_locus_dataframe


def delete_ambiguous_contigs(contig_dataframe):
    warning_contig = ""
    if (
        "para" in contig_dataframe["copy"].unique()
        and "main" in contig_dataframe["copy"].unique()
    ):
        warning_contig = contig_dataframe["saccver"].unique()[0]
    return warning_contig


def clean_paralogs(dataframe, n_cpu, log_file):
    logger = create_logger(log_file)
    grouped_samples = dataframe.groupby(["sample", "locus"])
    split_sample_locus_dataframe = []
    for group in grouped_samples:
        sample_locus_dataframe = group[1].reset_index(drop=True)
        split_sample_locus_dataframe.append(sample_locus_dataframe)
    with multiprocessing.Pool(processes=n_cpu) as pool_adjust:
        results = pool_adjust.map(adjust_orphaned_main, split_sample_locus_dataframe)
        dataframe = pandas.concat(results).reset_index(drop=True)
    grouped_contigs = dataframe.groupby("saccver")
    split_contig_dataframe = []
    for group in grouped_contigs:
        contig_dataframe = group[1]
        split_contig_dataframe.append(contig_dataframe)
    with multiprocessing.Pool(processes=n_cpu) as pool_warning_contigs:
        warning_contigs = pool_warning_contigs.map(
            delete_ambiguous_contigs, split_contig_dataframe
        )
    dataframe = dataframe[~dataframe["saccver"].isin(warning_contigs)].reset_index(
        drop=True
    )
    return dataframe


def score_1_2(sample_locus_dataframe):
    avg_ident = sample_locus_dataframe[sample_locus_dataframe["copy"] == "main"][
        "pident"
    ].mean()
    sample_locus_dataframe.loc[:, "score_1"] = avg_ident
    n_ex = len(sample_locus_dataframe[["exon", "copy"]].drop_duplicates())
    sample_locus_dataframe.loc[:, "score_2"] = n_ex
    return sample_locus_dataframe


def score_3(contig_dataframe):
    ex_contigs = contig_dataframe["exon"].unique()
    n_ex_contigs = len(ex_contigs)
    contig_dataframe.loc[:, "score_3"] = n_ex_contigs
    return contig_dataframe


def score_samples(dataframe, n_cpu, log_file):
    logger = create_logger(log_file)
    dataframe = dataframe.sort_values(
        ["exon", "pident", "qcovhsp", "evalue", "bitscore", "k-mer_cover"],
        ascending=(True, False, False, True, False, False),
    ).reset_index(drop=True)
    grouped_samples = dataframe.groupby(["sample", "locus"])
    split_sample_locus_dataframe = []
    for group in grouped_samples:
        sample_locus_dataframe = group[1]
        split_sample_locus_dataframe.append(sample_locus_dataframe)
    with multiprocessing.Pool(processes=n_cpu) as pool_score_1_2:
        results = pool_score_1_2.map(score_1_2, split_sample_locus_dataframe)
        all_hits_for_reference_scored_1_2 = pandas.concat(results).reset_index(
            drop=True
        )
    grouped_contigs = all_hits_for_reference_scored_1_2.groupby("saccver")
    split_contig_dataframe = []
    for group in grouped_contigs:
        contig_dataframe = group[1]
        split_contig_dataframe.append(contig_dataframe)
    with multiprocessing.Pool(processes=n_cpu) as pool_score_3:
        results = pool_score_3.map(score_3, split_contig_dataframe)
        all_hits_for_reference_scored_3 = pandas.concat(results).reset_index(drop=True)
    all_hits_for_reference_scored_3.sort_values(
        ["exon", "copy", "score_2", "score_1", "sample", "score_3"],
        ascending=[True, True, False, False, True, False],
        inplace=True,
    )
    all_hits_for_reference_scored_3.reset_index(drop=True, inplace=True)
    return all_hits_for_reference_scored_3


def phase_wo_paralog(sample_locus_dataframe):
    ident_array = sample_locus_dataframe["pident"].values.tolist()
    if len(ident_array) < 2:
        return sample_locus_dataframe
    boundary = find_cluster(ident_array)
    sample_locus_dataframe = sample_locus_dataframe[
        sample_locus_dataframe["pident"] > boundary
    ]
    return sample_locus_dataframe


def detect_paralogs(
    pairwise_distances,
    grouped_samples,
    paralog_min_divergence,
    paralog_max_divergence,
    log_file,
):
    logger = create_logger(log_file)
    sample = grouped_samples[0]
    grouped_samples_as_list = []
    grouped_samples = (
        grouped_samples[1]
        .sort_values(
            ["exon", "pident", "qcovhsp", "evalue", "bitscore", "k-mer_cover"],
            ascending=(True, False, False, True, False, False),
        )
        .reset_index(drop=True)
    )
    hits_grouped_exon = grouped_samples.groupby(["exon"])
    for group_exon in hits_grouped_exon:
        exon = group_exon[0]
        group_exon_dataframe = group_exon[1].reset_index(drop=True)
        main_copy = f"""{exon}_N_{group_exon_dataframe.loc[0, "contig"]}_{sample}"""
        id = 0
        len_df = len(group_exon_dataframe)
        while main_copy not in pairwise_distances.values:
            if id + 1 == len_df:
                break
            id += 1
            main_copy = (
                f"""{exon}_N_{group_exon_dataframe.loc[id, "contig"]}_{sample}"""
            )
        if main_copy not in pairwise_distances.values:
            continue
        main_copy_entry = group_exon_dataframe.iloc[[0]].reset_index(drop=True)
        main_copy_entry.loc[0, "copy"] = "main"
        grouped_samples_as_list.append(main_copy_entry.loc[0, :].values.tolist())
        paralog_found = False
        for index in range(id + 1, len_df):
            copy_to_compare = (
                f"""{exon}_N_{group_exon_dataframe.loc[index, "contig"]}_{sample}"""
            )
            seq1_main_seq2_para = (
                (pairwise_distances["seq1"] == main_copy)
                & (pairwise_distances["seq2"] == copy_to_compare)
            ).any()
            seq2_main_seq1_para = (
                (pairwise_distances["seq2"] == main_copy)
                & (pairwise_distances["seq1"] == copy_to_compare)
            ).any()
            if not seq1_main_seq2_para and not seq2_main_seq1_para:
                continue
            if seq1_main_seq2_para:
                div = pairwise_distances[
                    (pairwise_distances["seq1"] == main_copy)
                    & (pairwise_distances["seq2"] == copy_to_compare)
                ]["dist"]
            else:
                div = pairwise_distances[
                    (pairwise_distances["seq2"] == main_copy)
                    & (pairwise_distances["seq1"] == copy_to_compare)
                ]["dist"]
            div = div.values[0]
            secondary_copy_entry = group_exon_dataframe.iloc[[index]].reset_index(
                drop=True
            )
            if paralog_min_divergence < div < paralog_max_divergence:
                logger.info(f"Paralog detected in {sample} for {exon}")
                paralog_found = True
                secondary_copy_entry.loc[0, "copy"] = "para"
            elif paralog_min_divergence > div:
                secondary_copy_entry.loc[0, "copy"] = "main"
            grouped_samples_as_list.append(
                secondary_copy_entry.loc[0, :].values.tolist()
            )
        if not paralog_found:
            logger.info(f"Paralog not found in {sample} for {exon}")
    paralog_for_sample = pandas.DataFrame(
        grouped_samples_as_list,
        columns=[
            "qaccver",
            "saccver",
            "pident",
            "qcovhsp",
            "evalue",
            "bitscore",
            "sstart",
            "send",
            "locus",
            "k-mer_cover",
            "exon",
            "sample",
            "sequence",
            "contig",
            "copy",
        ],
    )
    return paralog_for_sample


def locus_stats(all_paralogs_for_reference):
    """"""
    locus_statistics = pandas.DataFrame([], columns=["samples\locus"])
    locus_statistics.set_index("samples\locus", drop=True, inplace=True)
    grouped_sample_locus = all_paralogs_for_reference.groupby(["sample", "locus"])
    for group in grouped_sample_locus:
        sample_locus_dataframe = group[1]
        sample = group[0][0]
        locus = group[0][1]
        if "para" in sample_locus_dataframe["copy"].values.tolist():
            locus_statistics.loc[sample, locus] = "Yes"
        else:
            locus_statistics.loc[sample, locus] = "No"
    return locus_statistics


def exon_stats(all_paralogs_for_reference):
    """"""
    exon_statistics = pandas.DataFrame([], columns=["samples\locus"])
    exon_statistics.set_index("samples\locus", drop=True, inplace=True)
    grouped_sample_exon = all_paralogs_for_reference.groupby(["sample", "exon"])
    for group in grouped_sample_exon:
        sample_exon_dataframe = group[1]
        sample = group[0][0]
        locus = group[0][1]
        if "para" in sample_exon_dataframe["copy"].values.tolist():
            exon_statistics.loc[sample, locus] = "Yes"
        else:
            exon_statistics.loc[sample, locus] = "No"
    return exon_statistics


def paralog_stats(locus_statistics):
    locus_statistics.reset_index(inplace=True)
    paralog_statistics = pandas.DataFrame(
        [], columns=["sample", "number_of_paralogous_loci"]
    )
    paralog_statistics.set_index("sample", drop=True, inplace=True)
    locus_statistics.replace("Yes", 1, inplace=True)
    locus_statistics.replace("No", 0, inplace=True)
    # locus_statistics.fillna("0", inplace=True)
    for id in range(len(locus_statistics)):
        sample = locus_statistics.loc[id, "samples\locus"]
        row = locus_statistics.loc[id:id]
        n_par = row.sum(axis=1, skipna=True, numeric_only=True).iloc[0]
        paralog_statistics.loc[sample] = n_par
    paralog_statistics["number_of_paralogous_loci"] = paralog_statistics[
        "number_of_paralogous_loci"
    ].replace(["0", 0], "0/NaN")
    return paralog_statistics


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


def create_reference_wo_paralogs(
    data_folder,
    all_hits_for_reference,
    blocklist,
    num_cores,
    log_file,
):
    """"""
    logger = create_logger(log_file)
    logger.info("Creating customized reference...")
    all_hits_for_reference = all_hits_for_reference.sort_values(
        ["exon", "pident", "qcovhsp", "evalue", "bitscore", "k-mer_cover"],
        ascending=(True, False, False, True, False, False),
    ).reset_index(drop=True)
    grouped_samples = all_hits_for_reference.groupby(["sample", "locus"])
    split_sample_locus_dataframe = []
    for group in grouped_samples:
        sample_locus_dataframe = group[1].reset_index(drop=True)
        split_sample_locus_dataframe.append(sample_locus_dataframe)
    with multiprocessing.Pool(processes=num_cores) as pool_wo_para:
        reustls = pool_wo_para.map(phase_wo_paralog, split_sample_locus_dataframe)
    all_hits_for_reference = pandas.concat(reustls).reset_index(drop=True)
    all_hits_for_reference["contig"] = (
        all_hits_for_reference["saccver"].str.split("_N_").str[1]
    )
    all_hits_for_reference = all_hits_for_reference[
        ~all_hits_for_reference["sample"].isin(blocklist)
    ]
    all_hits_for_reference["copy"] = "main"
    all_hits_for_reference_scored = score_samples(
        all_hits_for_reference, num_cores, log_file
    )
    all_paralogs_for_reference_to_write = all_hits_for_reference_scored.drop_duplicates(
        subset=["exon"]
    ).reset_index(drop=True)
    all_paralogs_for_reference_to_write["exon_num"] = (
        all_paralogs_for_reference_to_write["exon"].str.split("_").str[-1]
    )
    all_paralogs_for_reference_to_write[
        "exon_num"
    ] = all_paralogs_for_reference_to_write["exon_num"].astype(int)
    all_paralogs_for_reference_to_write = (
        all_paralogs_for_reference_to_write.sort_values(
            ["locus", "exon_num"],
        ).reset_index(drop=True)
    )
    previous_locus = ""
    with open(
        os.path.join(
            data_folder,
            "41without_par",
            f"customized_reference_for_HybPhyloMaker.fas",
        ),
        "w",
    ) as customized_reference_hpm, open(
        os.path.join(
            data_folder,
            "41without_par",
            f"customized_reference_for_ParalogWizard_separated_exons.fas",
        ),
        "w",
    ) as customized_reference_pw_separate, open(
        os.path.join(
            data_folder,
            "41without_par",
            f"customized_reference_for_HybPiper_concatenated_exons.fas",
        ),
        "w",
    ) as customized_reference_hp_concat:
        for id in range(len(all_paralogs_for_reference_to_write)):
            sample = all_paralogs_for_reference_to_write.loc[id, "sample"]
            exon_num = all_paralogs_for_reference_to_write.loc[id, "exon"].split("_")[
                -1
            ]
            locus = all_paralogs_for_reference_to_write.loc[id, "locus"]
            contig = all_paralogs_for_reference_to_write.loc[id, "contig"]
            seq = all_paralogs_for_reference_to_write.loc[id, "sequence"]
            name_seq_hpm = f">Assembly_{locus}_Contig_{exon_num}_{sample}_N_{contig}"
            name_seq_pw_separate = (
                f">{sample.replace('-', '_')}_N_{contig}-{locus}_exon_{exon_num}"
            )
            customized_reference_hpm.write(f"{name_seq_hpm}\n{seq}\n")
            customized_reference_pw_separate.write(f"{name_seq_pw_separate}\n{seq}\n")
            if locus != previous_locus:
                name_seq_hp_concat = f">{locus}-{locus}"
                customized_reference_hp_concat.write(f"\n{name_seq_hp_concat}\n{seq}")
                previous_locus = locus
            else:
                customized_reference_hp_concat.write(f"{seq}")
        with open(
            os.path.join(
                data_folder,
                "41without_par",
                f"customized_reference_for_HybPiper_concatenated_exons.fas",
            ),
            "r",
        ) as customized_reference_hp_concat:
            data = customized_reference_hp_concat.read().splitlines(True)
        with open(
            os.path.join(
                data_folder,
                "41without_par",
                f"customized_reference_for_HybPiper_concatenate_exons.fas",
            ),
            "w",
        ) as customized_reference_hp_concat:
            customized_reference_hp_concat.writelines(data[1:])
        logger.info("Creating customized reference...")


def prepare_to_write(all_paralogs_for_reference_scored):
    grouped_locus = all_paralogs_for_reference_scored.groupby(["locus"])
    prepared_loci = []
    for locus_group in grouped_locus:
        locus_dataframe = locus_group[1]
        grouped_exons = locus_dataframe.groupby("exon")
        prepared_exons = []
        for exon_group in grouped_exons:
            exon_dataframe = exon_group[1]
            if (exon_dataframe["copy"] == "para").any():
                samples_w_para = set(
                    exon_dataframe[exon_dataframe["copy"] == "para"]["sample"].unique()
                )
                prepared_exons.append(
                    exon_dataframe[exon_dataframe["sample"].isin(samples_w_para)]
                )
            else:
                prepared_exons.append(exon_dataframe)
        prepared_loci.append(pandas.concat(prepared_exons))
    all_paralogs_for_reference_to_write = pandas.concat(prepared_loci)
    all_paralogs_for_reference_to_write = (
        all_paralogs_for_reference_to_write.drop_duplicates(subset=["sequence"])
    )
    all_paralogs_for_reference_to_write = (
        all_paralogs_for_reference_to_write.drop_duplicates(
            subset=["exon", "copy"]
        ).reset_index(drop=True)
    )
    return all_paralogs_for_reference_to_write


def create_reference_w_paralogs(
    data_folder,
    all_hits_for_reference,
    paralog_min_divergence,
    paralog_max_divergence,
    blocklist,
    num_cores,
    log_file,
):
    """"""
    logger = create_logger(log_file)
    pairwise_distances = pandas.read_csv(
        os.path.join(data_folder, "40aln_orth_par", "pairwise_distances.tsv"), sep="\t"
    )
    logger.info("Creating customized reference...")
    all_hits_for_reference["contig"] = (
        all_hits_for_reference["saccver"].str.split("_N_").str[1]
    )
    grouped_samples = all_hits_for_reference.groupby("sample")
    list_groupes = []
    for group in grouped_samples:
        list_groupes.append(group)
    args_detect = list(
        zip(
            [pairwise_distances] * len(list_groupes),
            list_groupes,
            [paralog_min_divergence] * len(list_groupes),
            [paralog_max_divergence] * len(list_groupes),
            [log_file] * len(list_groupes),
        )
    )
    with multiprocessing.Pool(processes=num_cores) as pool_detect:
        results = pool_detect.starmap(detect_paralogs, args_detect)
        all_paralogs_for_reference = pandas.concat(results).reset_index(drop=True)
    all_paralogs_for_reference.dropna(subset=["copy"], inplace=True)
    all_paralogs_for_reference.to_csv(
        os.path.join(data_folder, "41detected_par", "all_paralogs_for_reference.tsv"),
        sep="\t",
        index=False,
    )
    all_paralogs_for_reference_cleaned = clean_paralogs(
        all_paralogs_for_reference, num_cores, log_file
    )
    all_paralogs_for_reference_cleaned.to_csv(
        os.path.join(
            data_folder, "41detected_par", "all_paralogs_for_reference_cleaned.tsv"
        ),
        sep="\t",
        index=False,
    )
    all_paralogs_for_reference_scored = score_samples(
        all_paralogs_for_reference_cleaned, num_cores, log_file
    )
    all_paralogs_for_reference_scored.to_csv(
        os.path.join(
            data_folder, "41detected_par", "all_paralogs_for_reference_scored.tsv"
        ),
        sep="\t",
        index=False,
    )
    all_paralogs_for_reference_scored = all_paralogs_for_reference_scored[
        ~all_paralogs_for_reference_scored["sample"].isin(blocklist)
    ].reset_index(drop=True)
    all_paralogs_for_reference_to_write = prepare_to_write(
        all_paralogs_for_reference_scored
    )
    with open(
        os.path.join(
            data_folder,
            "41detected_par",
            f"customized_reference_div_{paralog_min_divergence}_{paralog_max_divergence}.fas",
        ),
        "w",
    ) as customized_reference:
        for id in range(len(all_paralogs_for_reference_to_write)):
            sample = all_paralogs_for_reference_to_write.loc[id, "sample"]
            exon_num = all_paralogs_for_reference_to_write.loc[id, "exon"].split("_")[
                -1
            ]
            locus = all_paralogs_for_reference_to_write.loc[id, "locus"]
            contig = all_paralogs_for_reference_to_write.loc[id, "contig"]
            copy = all_paralogs_for_reference_to_write.loc[id, "copy"]
            seq = all_paralogs_for_reference_to_write.loc[id, "sequence"]
            if copy == "main":
                copy = ""
            name_seq = f">Assembly_{locus}{copy}_Contig_{exon_num}_{sample}_N_{contig}"
            customized_reference.write(f"{name_seq}\n{seq}\n")
    locus_statistics = locus_stats(all_paralogs_for_reference)
    locus_statistics.to_csv(
        os.path.join(
            data_folder,
            "41detected_par",
            f"locus_statistics_div_{paralog_min_divergence}_{paralog_max_divergence}.tsv",
        ),
        sep="\t",
        na_rep="NaN",
    )
    exon_statistics = exon_stats(all_paralogs_for_reference)
    exon_statistics.to_csv(
        os.path.join(
            data_folder,
            "41detected_par",
            f"exon_statistics_div_{paralog_min_divergence}_{paralog_max_divergence}.tsv",
        ),
        sep="\t",
        na_rep="NaN",
    )
    paralog_statistics = paralog_stats(locus_statistics)
    paralog_statistics.to_csv(
        os.path.join(
            data_folder,
            "41detected_par",
            f"paralog_statistics_div_{paralog_min_divergence}_{paralog_max_divergence}.tsv",
        ),
        sep="\t",
        na_rep="NaN",
    )
    logger.info("Customized reference created!\n")
