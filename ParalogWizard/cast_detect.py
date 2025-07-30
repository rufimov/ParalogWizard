import multiprocessing
import os
import sys
from typing import List
import pandas as pd
from ParalogWizard import setup_logging, worker_initializer

# -----------------------------------------------------------------------------
# Module-level logger configured once.
# -----------------------------------------------------------------------------
logger = setup_logging()


# =============================================================================
# Logging Decorator
# =============================================================================
def log_exceptions(func):
    """
    Decorator that logs function entry, exit, and any exceptions.
    Exits with code 1 upon an exception.
    """
    from functools import wraps

    @wraps(func)
    def wrapper(*args, **kwargs):
        try:
            result = func(*args, **kwargs)
            return result
        except Exception as e:
            logger.exception(f"Exception in {func.__name__}: {e}")
            sys.exit(1)

    return wrapper


@log_exceptions
def find_cluster(array: List[float]):
    logger.info(f"find_cluster: Received array with {len(array)} elements")
    array.sort()
    gaps = [array[i] - array[i - 1] for i in range(1, len(array))]
    biggest_gap = max(gaps)
    biggest_gap_middle = biggest_gap / 2
    biggest_gap_index = gaps.index(biggest_gap)
    boundary = array[biggest_gap_index] + biggest_gap_middle
    logger.info(f"find_cluster: Computed boundary = {boundary}")
    return boundary


@log_exceptions
def adjust_orphaned_main(sample_locus_dataframe):
    logger.info("adjust_orphaned_main: Processing a sample locus dataframe")
    copies = sample_locus_dataframe["copy"].values.tolist()
    if "para" not in copies:
        logger.info("adjust_orphaned_main: No 'para' copies found, returning original dataframe")
        return sample_locus_dataframe
    grouped_exons = sample_locus_dataframe.groupby("exon")
    exons_needed_clustering = []
    for group in grouped_exons:
        if "para" not in group[1]["copy"].values.tolist():
            exons_needed_clustering.append(group[0])
    if not exons_needed_clustering:
        logger.info("adjust_orphaned_main: No exons require clustering, returning original dataframe")
        return sample_locus_dataframe
    ident_array = sample_locus_dataframe[
        ~sample_locus_dataframe["exon"].isin(exons_needed_clustering)
    ]["pident"].values.tolist()
    if len(ident_array) < 2:
        logger.warning(
            "adjust_orphaned_main: Not enough pident values to compute a cluster boundary. "
            "Skipping clustering for this sample.")
        return sample_locus_dataframe
    boundary = find_cluster(ident_array)
    logger.info(f"adjust_orphaned_main: Using boundary = {boundary}")
    pidents_main = sample_locus_dataframe[
        (sample_locus_dataframe["copy"] == "main")
        & (~sample_locus_dataframe["exon"].isin(exons_needed_clustering))
        ]["pident"]
    pidents_para = sample_locus_dataframe[sample_locus_dataframe["copy"] == "para"]["pident"]
    if (pidents_main < boundary).any() or (pidents_para > boundary).any():
        logger.info("adjust_orphaned_main: Condition met, removing exons needing clustering")
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
    logger.info("adjust_orphaned_main: Clustering adjustment complete")
    return sample_locus_dataframe


@log_exceptions
def delete_ambiguous_contigs(contig_dataframe):
    logger.info("delete_ambiguous_contigs: Processing contig dataframe")
    warning_contig = ""
    if (
            "para" in contig_dataframe["copy"].unique()
            and "main" in contig_dataframe["copy"].unique()
    ):
        warning_contig = contig_dataframe["saccver"].unique()[0]
        logger.info(f"delete_ambiguous_contigs: Ambiguous contig detected: {warning_contig}")
    return warning_contig


@log_exceptions
def clean_paralogs(dataframe, n_cpu, log_queue):
    logger.info("clean_paralogs: Starting cleaning of paralogs")
    grouped_samples = dataframe.groupby(["sample", "locus"])
    split_sample_locus_dataframe = [group[1].reset_index(drop=True) for group in grouped_samples]
    with multiprocessing.Pool(processes=n_cpu, initializer=worker_initializer, initargs=(log_queue,)) as pool_adjust:
        results = pool_adjust.map(adjust_orphaned_main, split_sample_locus_dataframe)
        dataframe = pd.concat(results).reset_index(drop=True)
        pool_adjust.close()
        pool_adjust.join()
    grouped_contigs = dataframe.groupby("saccver")
    split_contig_dataframe = [group[1] for group in grouped_contigs]
    with multiprocessing.Pool(processes=n_cpu, initializer=worker_initializer,
                              initargs=(log_queue,)) as pool_warning_contigs:
        warning_contigs = pool_warning_contigs.map(delete_ambiguous_contigs, split_contig_dataframe)
        pool_warning_contigs.close()
        pool_warning_contigs.join()
    dataframe = dataframe[~dataframe["saccver"].isin(warning_contigs)].reset_index(drop=True)
    logger.info("clean_paralogs: Cleaning complete")
    return dataframe


@log_exceptions
def score_1_2(sample_locus_dataframe):
    logger.info("score_1_2: Scoring a sample locus dataframe")
    avg_ident = sample_locus_dataframe[sample_locus_dataframe["copy"] == "main"]["pident"].mean()
    sample_locus_dataframe.loc[:, "score_1"] = avg_ident
    n_ex = len(sample_locus_dataframe[["exon", "copy"]].drop_duplicates())
    sample_locus_dataframe.loc[:, "score_2"] = n_ex
    logger.info(f"score_1_2: score_1 set to {avg_ident} and score_2 set to {n_ex}")
    return sample_locus_dataframe


@log_exceptions
def score_3(contig_dataframe):
    logger.info("score_3: Scoring contig dataframe")
    ex_contigs = contig_dataframe["exon"].unique()
    n_ex_contigs = len(ex_contigs)
    contig_dataframe.loc[:, "score_3"] = n_ex_contigs
    logger.info(f"score_3: score_3 set to {n_ex_contigs}")
    return contig_dataframe


@log_exceptions
def score_samples(dataframe, n_cpu, log_queue):
    logger.info("score_samples: Starting sample scoring")
    dataframe = dataframe.sort_values(
        ["exon", "pident", "qcovhsp", "evalue", "bitscore", "k-mer_cover"],
        ascending=(True, False, False, True, False, False),
    ).reset_index(drop=True)
    grouped_samples = dataframe.groupby(["sample", "locus"])
    split_sample_locus_dataframe = [group[1] for group in grouped_samples]
    with multiprocessing.Pool(processes=n_cpu, initializer=worker_initializer, initargs=(log_queue,)) as pool_score_1_2:
        results = pool_score_1_2.map(score_1_2, split_sample_locus_dataframe)
        all_hits_for_reference_scored_1_2 = pd.concat(results).reset_index(drop=True)
        pool_score_1_2.close()
        pool_score_1_2.join()
    grouped_contigs = all_hits_for_reference_scored_1_2.groupby("saccver")
    split_contig_dataframe = [group[1] for group in grouped_contigs]
    with multiprocessing.Pool(processes=n_cpu, initializer=worker_initializer, initargs=(log_queue,)) as pool_score_3:
        results = pool_score_3.map(score_3, split_contig_dataframe)
        all_hits_for_reference_scored_3 = pd.concat(results).reset_index(drop=True)
        pool_score_3.close()
        pool_score_3.join()
    all_hits_for_reference_scored_3.sort_values(
        ["exon", "copy", "score_2", "score_1", "sample", "score_3"],
        ascending=[True, True, False, False, True, False],
        inplace=True,
    )
    all_hits_for_reference_scored_3.reset_index(drop=True, inplace=True)
    logger.info(f"score_samples: Scoring complete with {len(all_hits_for_reference_scored_3)} entries")
    return all_hits_for_reference_scored_3


@log_exceptions
def phase_wo_paralog(sample_locus_dataframe):
    logger.info("phase_wo_paralog: Phasing sample without paralogs")
    ident_array = sample_locus_dataframe["pident"].values.tolist()
    if len(ident_array) < 2:
        logger.info("phase_wo_paralog: Not enough pident values; returning original dataframe")
        return sample_locus_dataframe
    boundary = find_cluster(ident_array)
    logger.info(f"phase_wo_paralog: Using boundary = {boundary}")
    sample_locus_dataframe = sample_locus_dataframe[sample_locus_dataframe["pident"] > boundary]
    return sample_locus_dataframe


@log_exceptions
def detect_paralogs(pairwise_distances, grouped_samples, paralog_min_divergence, paralog_max_divergence):
    logger.info("detect_paralogs: Starting paralog detection")
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
    hits_grouped_exon = grouped_samples.groupby("exon")
    for group_exon in hits_grouped_exon:
        exon = group_exon[0]
        group_exon_dataframe = group_exon[1].reset_index(drop=True)
        main_copy = f"{exon}_N_{group_exon_dataframe.loc[0, 'contig']}_{sample}"
        idx = 0
        len_df = len(group_exon_dataframe)
        while main_copy not in pairwise_distances.values:
            if idx + 1 == len_df:
                break
            idx += 1
            main_copy = f"{exon}_N_{group_exon_dataframe.loc[idx, 'contig']}_{sample}"
        if main_copy not in pairwise_distances.values:
            continue
        main_copy_entry = group_exon_dataframe.iloc[[0]].reset_index(drop=True)
        main_copy_entry.loc[0, "copy"] = "main"
        grouped_samples_as_list.append(main_copy_entry.loc[0, :].values.tolist())
        paralog_found = False
        for index in range(idx + 1, len_df):
            copy_to_compare = f"{exon}_N_{group_exon_dataframe.loc[index, 'contig']}_{sample}"
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
            secondary_copy_entry = group_exon_dataframe.iloc[[index]].reset_index(drop=True)
            if paralog_min_divergence < div < paralog_max_divergence:
                logger.info(f"detect_paralogs: Paralog detected in {sample} for {exon}")
                paralog_found = True
                secondary_copy_entry.loc[0, "copy"] = "para"
            elif paralog_min_divergence > div:
                secondary_copy_entry.loc[0, "copy"] = "main"
            grouped_samples_as_list.append(secondary_copy_entry.loc[0, :].values.tolist())
        if not paralog_found:
            logger.info(f"detect_paralogs: Paralog not found in {sample} for {exon}")
    paralog_for_sample = pd.DataFrame(
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
    logger.info(f"detect_paralogs: Detection complete for sample {sample}")
    return paralog_for_sample


@log_exceptions
def locus_stats(all_paralogs_for_reference):
    logger.info("locus_stats: Calculating locus statistics")
    data = {}
    grouped_sample_locus = all_paralogs_for_reference.groupby(["sample", "locus"])
    for (sample, locus), sample_locus_dataframe in grouped_sample_locus:
        if sample not in data:
            data[sample] = {}
        data[sample][locus] = "Yes" if "para" in sample_locus_dataframe["copy"].values.tolist() else "No"
    locus_statistics = pd.DataFrame.from_dict(data, orient="index")
    # Set the index name to "samples\locus"
    locus_statistics.index.name = r"samples\locus"
    logger.info("locus_stats: Completed locus statistics")
    return locus_statistics


def exon_stats(all_paralogs_for_reference):
    logger.info("exon_stats: Calculating exon statistics")
    data = {}
    grouped_sample_exon = all_paralogs_for_reference.groupby(["sample", "exon"])
    for (sample, exon), sample_exon_dataframe in grouped_sample_exon:
        if sample not in data:
            data[sample] = {}
        data[sample][exon] = "Yes" if "para" in sample_exon_dataframe["copy"].values.tolist() else "No"
    exon_statistics = pd.DataFrame.from_dict(data, orient="index")
    logger.info("exon_stats: Completed exon statistics")
    return exon_statistics


@log_exceptions
def paralog_stats(locus_statistics):
    logger.info("paralog_stats: Calculating paralog statistics")
    # Reset the index; the first column will be named "samples\locus"
    locus_statistics = locus_statistics.reset_index()
    paralog_statistics = pd.DataFrame([], columns=["samples\\locus", "number_of_paralogous_loci"])
    paralog_statistics.set_index("samples\\locus", drop=True, inplace=True)
    locus_statistics = locus_statistics.replace("Yes", 1)
    locus_statistics = locus_statistics.replace("No", 0)
    locus_statistics = locus_statistics.infer_objects(copy=False)
    for idx in range(len(locus_statistics)):
        sample = locus_statistics.loc[idx, r"samples\locus"]
        row = locus_statistics.loc[idx:idx]
        n_par = row.sum(axis=1, skipna=True, numeric_only=True).iloc[0]
        paralog_statistics.loc[sample] = n_par
    paralog_statistics["number_of_paralogous_loci"] = paralog_statistics["number_of_paralogous_loci"].replace(["0", 0],
                                                                                                              "0/NaN")
    logger.info("paralog_stats: Paralog statistics computed")
    return paralog_statistics


@log_exceptions
def create_reference_wo_paralogs(data_folder, all_hits_for_reference, blocklist, num_cores, log_queue):
    logger.info("create_reference_wo_paralogs: Creating customized reference without paralogs")
    all_hits_for_reference = all_hits_for_reference.sort_values(
        ["exon", "pident", "qcovhsp", "evalue", "bitscore", "k-mer_cover"],
        ascending=(True, False, False, True, False, False),
    ).reset_index(drop=True)
    grouped_samples = all_hits_for_reference.groupby(["sample", "locus"])
    split_sample_locus_dataframe = [group[1].reset_index(drop=True) for group in grouped_samples]
    with multiprocessing.Pool(processes=num_cores, initializer=worker_initializer,
                              initargs=(log_queue,)) as pool_wo_para:
        results = pool_wo_para.map(phase_wo_paralog, split_sample_locus_dataframe)
        pool_wo_para.close()
        pool_wo_para.join()
    all_hits_for_reference = pd.concat(results).reset_index(drop=True)
    all_hits_for_reference["contig"] = all_hits_for_reference["saccver"].str.split("_N_").str[1]
    all_hits_for_reference = all_hits_for_reference[~all_hits_for_reference["sample"].isin(blocklist)]
    all_hits_for_reference["copy"] = "main"
    all_hits_for_reference_scored = score_samples(all_hits_for_reference, num_cores, log_queue)
    all_paralogs_for_reference_to_write = all_hits_for_reference_scored.drop_duplicates(subset=["exon"]).reset_index(
        drop=True)
    all_paralogs_for_reference_to_write["exon_num"] = all_paralogs_for_reference_to_write["exon"].str.split("_").str[
        -1].astype(int)
    all_paralogs_for_reference_to_write = all_paralogs_for_reference_to_write.sort_values(
        ["locus", "exon_num"]).reset_index(drop=True)
    previous_locus = ""
    with open(os.path.join(data_folder, "41without_par", "customized_reference_for_HybPhyloMaker.fas"),
              "w") as customized_reference_hpm, \
            open(os.path.join(data_folder, "41without_par",
                              "customized_reference_for_ParalogWizard_separated_exons.fas"),
                 "w") as customized_reference_pw_separate, \
            open(os.path.join(data_folder, "41without_par", "customized_reference_for_HybPiper_concatenated_exons.fas"),
                 "w") as customized_reference_hp_concat:
        for idx in range(len(all_paralogs_for_reference_to_write)):
            sample = all_paralogs_for_reference_to_write.loc[idx, "sample"]
            exon_num = all_paralogs_for_reference_to_write.loc[idx, "exon"].split("_")[-1]
            locus = all_paralogs_for_reference_to_write.loc[idx, "locus"]
            contig = all_paralogs_for_reference_to_write.loc[idx, "contig"]
            seq = all_paralogs_for_reference_to_write.loc[idx, "sequence"]
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
                os.path.join(data_folder, "41without_par", "customized_reference_for_HybPiper_concatenated_exons.fas"),
                "r") as customized_reference_hp_concat:
            data = customized_reference_hp_concat.read().splitlines(True)
        with open(os.path.join(data_folder, "41without_par", "customized_reference_for_HybPiper_concatenate_exons.fas"),
                  "w") as customized_reference_hp_concat:
            customized_reference_hp_concat.writelines(data[1:])
    logger.info("create_reference_wo_paralogs: Customized reference without paralogs created")


@log_exceptions
def prepare_to_write(all_paralogs_for_reference_scored):
    logger.info("prepare_to_write: Preparing paralog data to write to reference")
    grouped_locus = all_paralogs_for_reference_scored.groupby("locus")
    prepared_loci = []
    for locus_group in grouped_locus:
        locus_dataframe = locus_group[1]
        grouped_exons = locus_dataframe.groupby("exon")
        prepared_exons = []
        for exon_group in grouped_exons:
            exon_dataframe = exon_group[1]
            if (exon_dataframe["copy"] == "para").any():
                samples_w_para = set(exon_dataframe[exon_dataframe["copy"] == "para"]["sample"].unique())
                prepared_exons.append(exon_dataframe[exon_dataframe["sample"].isin(samples_w_para)])
            else:
                prepared_exons.append(exon_dataframe)
        prepared_loci.append(pd.concat(prepared_exons))
    all_paralogs_for_reference_to_write = pd.concat(prepared_loci)
    all_paralogs_for_reference_to_write = all_paralogs_for_reference_to_write.drop_duplicates(subset=["sequence"])
    all_paralogs_for_reference_to_write = all_paralogs_for_reference_to_write.drop_duplicates(
        subset=["exon", "copy"]).reset_index(drop=True)
    logger.info("prepare_to_write: Data prepared for writing")
    return all_paralogs_for_reference_to_write


@log_exceptions
def create_reference_w_paralogs(data_folder, all_hits_for_reference, paralog_min_divergence, paralog_max_divergence,
                                blocklist, num_cores, log_queue):
    logger.info("create_reference_w_paralogs: Creating customized reference with paralogs")
    pairwise_distances = pd.read_csv(os.path.join(data_folder, "40aln_orth_par", "pairwise_distances.tsv"), sep="\t")
    all_hits_for_reference["contig"] = all_hits_for_reference["saccver"].str.split("_N_").str[1]
    grouped_samples = all_hits_for_reference.groupby("sample")
    list_groupes = [group for group in grouped_samples]
    args_detect = list(
        zip(
            [pairwise_distances] * len(list_groupes),
            list_groupes,
            [paralog_min_divergence] * len(list_groupes),
            [paralog_max_divergence] * len(list_groupes),
        )
    )
    with multiprocessing.Pool(processes=num_cores, initializer=worker_initializer,
                              initargs=(log_queue,)) as pool_detect:
        results = pool_detect.starmap(detect_paralogs, args_detect)
        all_paralogs_for_reference = pd.concat(results).reset_index(drop=True)
        pool_detect.close()
        pool_detect.join()
    all_paralogs_for_reference.dropna(subset=["copy"], inplace=True)
    all_paralogs_for_reference.to_csv(os.path.join(data_folder, "41detected_par", "all_paralogs_for_reference.tsv"),
                                      sep="\t", index=False)
    all_paralogs_for_reference_cleaned = clean_paralogs(all_paralogs_for_reference, num_cores, log_queue)
    all_paralogs_for_reference_cleaned.to_csv(
        os.path.join(data_folder, "41detected_par", "all_paralogs_for_reference_cleaned.tsv"), sep="\t", index=False)
    all_paralogs_for_reference_scored = score_samples(all_paralogs_for_reference_cleaned, num_cores, log_queue)
    all_paralogs_for_reference_scored.to_csv(
        os.path.join(data_folder, "41detected_par", "all_paralogs_for_reference_scored.tsv"), sep="\t", index=False)
    all_paralogs_for_reference_scored = all_paralogs_for_reference_scored[
        ~all_paralogs_for_reference_scored["sample"].isin(blocklist)].reset_index(drop=True)
    all_paralogs_for_reference_to_write = prepare_to_write(all_paralogs_for_reference_scored)
    with open(os.path.join(data_folder, "41detected_par",
                           f"customized_reference_div_{paralog_min_divergence}_{paralog_max_divergence}.fas"),
              "w") as customized_reference:
        for idx in range(len(all_paralogs_for_reference_to_write)):
            sample = all_paralogs_for_reference_to_write.loc[idx, "sample"]
            exon_num = all_paralogs_for_reference_to_write.loc[idx, "exon"].split("_")[-1]
            locus = all_paralogs_for_reference_to_write.loc[idx, "locus"]
            contig = all_paralogs_for_reference_to_write.loc[idx, "contig"]
            copy = all_paralogs_for_reference_to_write.loc[idx, "copy"]
            seq = all_paralogs_for_reference_to_write.loc[idx, "sequence"]
            if copy == "main":
                copy = ""
            name_seq = f">Assembly_{locus}{copy}_Contig_{exon_num}_{sample}_N_{contig}"
            customized_reference.write(f"{name_seq}\n{seq}\n")
    locus_statistics = locus_stats(all_paralogs_for_reference)
    locus_statistics.to_csv(os.path.join(data_folder, "41detected_par",
                                         f"locus_statistics_div_{paralog_min_divergence}_{paralog_max_divergence}.tsv"),
                            sep="\t", na_rep="NaN")
    exon_statistics = exon_stats(all_paralogs_for_reference)
    exon_statistics.to_csv(os.path.join(data_folder, "41detected_par",
                                        f"exon_statistics_div_{paralog_min_divergence}_{paralog_max_divergence}.tsv"),
                           sep="\t", na_rep="NaN")
    paralog_statistics = paralog_stats(locus_statistics)
    paralog_statistics.to_csv(os.path.join(data_folder, "41detected_par",
                                           f"paralog_statistics_div_{paralog_min_divergence}_"
                                           f"{paralog_max_divergence}.tsv"),
                              sep="\t", na_rep="NaN")
    logger.info("create_reference_w_paralogs: Customized reference with paralogs created")
