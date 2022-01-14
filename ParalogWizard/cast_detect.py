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
    pidents = sample_locus_dataframe["copy"].values.tolist()
    if "para" in pidents:
        clustering_needed = False
        grouped_exons = sample_locus_dataframe.groupby("exon")
        exons_needed_clustering = []
        for group in grouped_exons:
            if "para" not in group[1]["copy"].values.tolist():
                exons_needed_clustering.append(group[0])
                clustering_needed = True
        if clustering_needed:
            ident_array = sample_locus_dataframe["pident"].values.tolist()
            boundary = find_cluster(ident_array)
            for id in range(len(sample_locus_dataframe)):
                if sample_locus_dataframe.loc[id, "exon"] in exons_needed_clustering:
                    ident = sample_locus_dataframe.loc[id, "pident"]
                    if ident > boundary:
                        sample_locus_dataframe.loc[id, "copy"] = "main"
                    else:
                        sample_locus_dataframe.loc[id, "copy"] = "para"
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


def score_1(exon_dataframe):
    counters = list(range(1, len(exon_dataframe) + 1))
    exon_dataframe.loc[:, "score_1"] = counters
    return exon_dataframe


def score_2(sample_locus_dataframe):
    ex = set(zip(sample_locus_dataframe["exon"], sample_locus_dataframe["copy"]))
    n_ex = len(ex)
    sample_locus_dataframe.loc[:, "score_2"] = n_ex
    av_rank = sample_locus_dataframe.loc[:, "score_1"].mean()
    sample_locus_dataframe.loc[:, "score_3"] = av_rank
    return sample_locus_dataframe


def score_3(contig_dataframe):
    ex_contigs = contig_dataframe["exon"].unique()
    n_ex_contigs = len(ex_contigs)
    contig_dataframe.loc[:, "score_4"] = n_ex_contigs
    return contig_dataframe


def score_samples(dataframe, n_cpu, log_file):
    logger = create_logger(log_file)
    dataframe = dataframe.sort_values(
        ["exon", "pident", "qcovhsp", "evalue", "bitscore", "k-mer_cover"],
        ascending=(True, False, False, True, False, False),
    ).reset_index(drop=True)
    grouped_exons = dataframe.groupby("exon")
    split_exon_dataframes = []
    for group in grouped_exons:
        exon_dataframe = group[1]
        split_exon_dataframes.append(exon_dataframe)
    with multiprocessing.Pool(processes=n_cpu) as pool_score_1:
        results = pool_score_1.map(score_1, split_exon_dataframes)
        all_hits_for_reference_scored_1 = pandas.concat(results).reset_index(drop=True)
    grouped_samples = all_hits_for_reference_scored_1.groupby(["sample", "locus"])
    split_sample_locus_dataframe = []
    for group in grouped_samples:
        sample_locus_dataframe = group[1]
        split_sample_locus_dataframe.append(sample_locus_dataframe)
    with multiprocessing.Pool(processes=n_cpu) as pool_score_2:
        results = pool_score_2.map(score_2, split_sample_locus_dataframe)
        all_hits_for_reference_scored_2 = pandas.concat(results).reset_index(drop=True)
    grouped_contigs = all_hits_for_reference_scored_2.groupby("saccver")
    split_contig_dataframe = []
    for group in grouped_contigs:
        contig_dataframe = group[1]
        split_contig_dataframe.append(contig_dataframe)
    with multiprocessing.Pool(processes=n_cpu) as pool_score_3:
        results = pool_score_3.map(score_3, split_contig_dataframe)
        all_hits_for_reference_scored_3 = pandas.concat(results).reset_index(drop=True)
    all_hits_for_reference_scored_3.sort_values(
        ["exon", "copy", "score_2", "score_3", "sample", "score_4"],
        ascending=[True, True, False, True, True, False],
        inplace=True,
    )
    all_hits_for_reference_scored_3.reset_index(drop=True, inplace=True)
    return all_hits_for_reference_scored_3


def phase_wo_paralog(sample_locus_dataframe):
    ident_array = sample_locus_dataframe["pident"].values.tolist()
    boundary = find_cluster(ident_array)
    sample_locus_dataframe = sample_locus_dataframe[sample_locus_dataframe['pident'] > boundary]
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
        paralog_found = False
        main_copy = f"""{exon}_N_{f'{group_exon_dataframe.loc[0, "saccver"]}_{sample}'.split('_N_')[1]}"""
        main_copy_saved = False
        for index in range(1, len(group_exon_dataframe)):
            main_copy_not_present = main_copy not in pairwise_distances["seq1"].unique()
            if main_copy_not_present:
                main_copy = f"""{exon}_N_{f'{group_exon_dataframe.loc[index, "saccver"]}_{sample}'.split(
                    '_N_')[1]} """
                continue
            if not main_copy_saved:
                main_copy_entry = group_exon_dataframe.iloc[[0]].reset_index(drop=True)
                main_copy_entry.loc[0, "copy"] = "main"
                grouped_samples_as_list.append(
                    main_copy_entry.loc[0, :].values.tolist()
                )
                main_copy_saved = True
            copy_to_compare = f"""{exon}_N_{f'{group_exon_dataframe.loc[index, "saccver"]}_{sample}'.split(
                '_N_')[1]}"""
            if copy_to_compare not in pairwise_distances["seq1"].unique():
                continue
            div = pairwise_distances[
                (pairwise_distances["seq1"] == main_copy)
                & (pairwise_distances["seq2"] == copy_to_compare)
                ]["dist"]
            if len(div) == 0:
                continue
            div = div.values[0]
            secondary_copy_entry = group_exon_dataframe.iloc[[index]].reset_index(
                drop=True
            )
            if paralog_min_divergence < div < paralog_max_divergence:
                logger.info(f"Paralog detected in {sample} for {exon}")
                paralog_found = True
                secondary_copy_entry.loc[0, "copy"] = "para"
            else:
                secondary_copy_entry.loc[0, "copy"] = "main"
            grouped_samples_as_list.append(
                secondary_copy_entry.loc[0, :].values.tolist()
            )
        if paralog_found is False:
            if main_copy in pairwise_distances["seq1"].unique():
                logger.info(f"Paralog not found in {sample} for {exon}")
                (
                    grouped_samples_as_list.append(
                        main_copy_entry.loc[0, :].values.tolist()
                    )
                )
            else:
                logger.info(
                    f"Not enough data in {sample} for {exon} for paralog detection."
                )
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
    all_hits_for_reference['copy'] = 'main'
    all_hits_for_reference_scored = score_samples(all_hits_for_reference, num_cores, log_file)
    all_paralogs_for_reference_to_write = (
        all_hits_for_reference_scored.drop_duplicates(
            subset=["exon"]
        ).reset_index(drop=True)
    )
    all_paralogs_for_reference_to_write['exon_num'] = all_paralogs_for_reference_to_write["exon"].str.split("_").str[
        -1
    ]
    all_paralogs_for_reference_to_write['exon_num'] = all_paralogs_for_reference_to_write['exon_num'].astype(int)
    all_paralogs_for_reference_to_write = all_paralogs_for_reference_to_write.sort_values(
        ["locus", "exon_num"],
    ).reset_index(drop=True)
    previous_locus = ''
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
            f"customized_reference_for_ParalogWizard_separate_exons.fas",
        ),
        "w",
    ) as customized_reference_pw_separate, open(
        os.path.join(
            data_folder,
            "41without_par",
            f"customized_reference_for_HybPiper_concatenate_exons.fas",
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
            name_seq_pw_separate = f">{sample.replace('-', '_')}_N_{contig}-{locus}_exon_{exon_num}"
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
                    f"customized_reference_for_HybPiper_concatenate_exons.fas",
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
    pairwise_distances_1 = pandas.read_csv(
        os.path.join(data_folder, "40aln_orth_par", "pairwise_distances.tsv"), sep="\t"
    )
    pairwise_distances_2 = pairwise_distances_1[["seq2", "dist", "seq1"]]
    pairwise_distances_2.rename({"seq2": "seq1", "seq1": "seq2"}, axis=1, inplace=True)
    pairwise_distances = pandas.concat(
        [pairwise_distances_1, pairwise_distances_2]
    ).reset_index(drop=True)
    logger.info("Creating customized reference...")
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

    all_paralogs_for_reference["contig"] = (
        all_paralogs_for_reference["saccver"].str.split("_N_").str[1]
    )
    all_paralogs_for_reference_cleaned = clean_paralogs(
        all_paralogs_for_reference, num_cores, log_file
    )
    all_paralogs_for_reference_cleaned.to_csv(
        os.path.join(data_folder, "41detected_par", "all_paralogs_for_reference.tsv"),
        sep="\t",
        index=False,
    )
    all_paralogs_for_reference_cleaned = all_paralogs_for_reference_cleaned[
        ~all_paralogs_for_reference_cleaned["sample"].isin(blocklist)
    ].reset_index(drop=True)

    all_paralogs_for_reference_scored = score_samples(
        all_paralogs_for_reference_cleaned, num_cores, logger
    )
    all_paralogs_for_reference_to_write = (
        all_paralogs_for_reference_scored.drop_duplicates(
            subset=["exon", "copy"]
        ).reset_index(drop=True)
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
