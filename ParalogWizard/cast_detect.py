import os
import numpy as np
import pandas as pd
import multiprocessing
from typing import List, Tuple
from functools import wraps
from ParalogWizard import setup_logging, worker_initializer


# ---------------------------
# Logging Decorator
# ---------------------------
def log_exceptions(func):
    """
    Decorator that logs function entry, exit, and any exceptions.
    """

    @wraps(func)
    def wrapper(*args, **kwargs):
        logger = setup_logging()
        try:
            result = func(*args, **kwargs)
            return result
        except Exception as e:
            logger.exception(f"Exception in {func.__name__}: {e}")
            raise

    return wrapper


# ---------------------------
# Helper Functions
# ---------------------------
@log_exceptions
def find_cluster(array: List[float]) -> float:
    if not array:
        raise ValueError("Empty array passed to find_cluster")
    arr = np.sort(np.array(array))
    gaps = np.diff(arr)  # This produces a NumPy array of differences.
    max_gap_idx = np.argmax(gaps)
    # Cast the scalar values to Python float
    max_gap = float(gaps[max_gap_idx])
    boundary = float(arr[max_gap_idx]) + max_gap / 2.0
    return boundary


@log_exceptions
def adjust_orphaned_main(df: pd.DataFrame) -> pd.DataFrame:
    """
    Adjusts the DataFrame by checking for exons that need re‐clustering.
    Uses vectorized operations where possible.
    """
    if "para" not in df["copy"].values:
        return df

    # Determine for each row if its exon group has any 'para'
    has_para = df.groupby("exon")["copy"].transform(lambda x: x.eq("para").any())
    exons_needing_clustering = df.loc[~has_para, "exon"].unique()
    if len(exons_needing_clustering) == 0:
        return df

    non_clustered_mask = ~df["exon"].isin(exons_needing_clustering)
    ident_array = df.loc[non_clustered_mask, "pident"].values
    if len(ident_array) == 0:
        return df
    boundary = find_cluster(ident_array.tolist())

    mask_main = (df["copy"] == "main") & df["exon"].isin(exons_needing_clustering)
    mask_para = (df["copy"] == "para") & df["exon"].isin(exons_needing_clustering)
    if (df.loc[mask_main, "pident"] < boundary).any() or (
        df.loc[mask_para, "pident"] > boundary
    ).any():
        # Inconsistent: drop those exons entirely
        return df[~df["exon"].isin(exons_needing_clustering)]
    else:
        df.loc[
            (df["pident"] > boundary) & df["exon"].isin(exons_needing_clustering),
            "copy",
        ] = "main"
        df.loc[
            (df["pident"] < boundary) & df["exon"].isin(exons_needing_clustering),
            "copy",
        ] = "para"
        return df


@log_exceptions
def delete_ambiguous_contigs(df: pd.DataFrame) -> str:
    """
    If both 'para' and 'main' are present for a contig (saccver), return that contig.
    Otherwise, return an empty string.
    """
    copies = df["copy"].unique()
    if "para" in copies and "main" in copies:
        return df["saccver"].iloc[0]
    return ""


@log_exceptions
def clean_paralogs(df: pd.DataFrame, n_cpu: int, log_queue) -> pd.DataFrame:
    """
    Clean the DataFrame by applying adjust_orphaned_main to each
    (sample, locus) group and then removing contigs that are ambiguous.
    """
    # Process (sample, locus) groups in parallel.
    groups = [g.reset_index(drop=True) for _, g in df.groupby(["sample", "locus"])]
    with multiprocessing.Pool(
        processes=n_cpu, initializer=worker_initializer, initargs=(log_queue,)
    ) as pool:
        adjusted = pool.map(adjust_orphaned_main, groups)
    df_adjusted = pd.concat(adjusted, ignore_index=True)

    # Process each contig group to detect ambiguous contigs.
    contig_groups = [g for _, g in df_adjusted.groupby("saccver")]
    with multiprocessing.Pool(
        processes=n_cpu, initializer=worker_initializer, initargs=(log_queue,)
    ) as pool:
        warnings = pool.map(delete_ambiguous_contigs, contig_groups)
    # Remove rows where 'saccver' is in the list of ambiguous contigs.
    df_cleaned = df_adjusted[~df_adjusted["saccver"].isin(warnings)].reset_index(
        drop=True
    )
    # Use logger from the decorator (no need to reassign here)
    setup_logging().debug(
        f"clean_paralogs() cleaned DataFrame shape: {df_cleaned.shape}"
    )
    return df_cleaned


# ---------------------------
# Scoring Functions
# ---------------------------
@log_exceptions
def score_1_2(df: pd.DataFrame) -> pd.DataFrame:
    """
    Computes score_1 (average pident for main copies) and score_2 (number of distinct
    exon-copy combos) per (sample, locus) group and merges these scores back into the DataFrame.
    """
    agg_main = df[df["copy"] == "main"].groupby(["sample", "locus"])["pident"].mean()
    agg_main.name = "score_1"  # Set the name via the .name attribute

    agg_exon = df.groupby(["sample", "locus"]).apply(
        lambda g: g[["exon", "copy"]].drop_duplicates().shape[0]
    )
    agg_exon.name = "score_2"  # Set the name via the .name attribute

    agg = pd.concat([agg_main, agg_exon], axis=1).reset_index()
    df_scored = df.merge(agg, on=["sample", "locus"], how="left")
    return df_scored


@log_exceptions
def score_3(df: pd.DataFrame) -> pd.DataFrame:
    """
    Computes score_3: number of unique exons per contig (saccver).
    """
    agg = df.groupby("saccver")["exon"].nunique()
    agg.name = "score_3"  # Set the name via the .name attribute
    agg = agg.reset_index()
    df_scored = df.merge(agg, on="saccver", how="left")
    return df_scored


@log_exceptions
def score_samples(df: pd.DataFrame, n_cpu: int, log_queue) -> pd.DataFrame:
    """
    Sorts the DataFrame, applies scoring on (sample, locus) groups in parallel, then on
    contig groups in parallel, and finally sorts the combined results.
    """
    df = df.sort_values(
        ["exon", "pident", "qcovhsp", "evalue", "bitscore", "k-mer_cover"],
        ascending=(True, False, False, True, False, False),
    ).reset_index(drop=True)

    groups = [g.reset_index(drop=True) for _, g in df.groupby(["sample", "locus"])]
    with multiprocessing.Pool(
        processes=n_cpu, initializer=worker_initializer, initargs=(log_queue,)
    ) as pool:
        scored1 = pool.map(score_1_2, groups)
    df_scored = pd.concat(scored1, ignore_index=True)

    contig_groups = [g for _, g in df_scored.groupby("saccver")]
    with multiprocessing.Pool(
        processes=n_cpu, initializer=worker_initializer, initargs=(log_queue,)
    ) as pool:
        scored2 = pool.map(score_3, contig_groups)
    df_scored = pd.concat(scored2, ignore_index=True)

    df_scored = df_scored.sort_values(
        ["exon", "copy", "score_2", "score_1", "sample", "score_3"],
        ascending=[True, True, False, False, True, False],
    ).reset_index(drop=True)
    return df_scored


@log_exceptions
def phase_wo_paralog(df: pd.DataFrame) -> pd.DataFrame:
    """
    If there are at least two rows, computes a boundary using pident and returns only rows
    with pident greater than the boundary.
    """
    if df.shape[0] < 2:
        return df
    boundary = find_cluster(df["pident"].tolist())
    return df[df["pident"] > boundary].reset_index(drop=True)


# ---------------------------
# Detect Paralogs Function
# ---------------------------
@log_exceptions
def detect_paralogs(
    pairwise_distances: pd.DataFrame,
    grouped_samples: Tuple[str, pd.DataFrame],
    paralog_min_divergence: float,
    paralog_max_divergence: float,
) -> pd.DataFrame:
    """
    For each exon in the sample group, determines the main copy by checking its presence in
    pairwise_distances. Then, for each subsequent hit in that exon, uses vectorized comparisons
    (with .eq()) to determine a match and labels accordingly.
    """
    logger = setup_logging()
    sample = grouped_samples[0]
    result_rows = []
    df_sample = (
        grouped_samples[1]
        .sort_values(
            ["exon", "pident", "qcovhsp", "evalue", "bitscore", "k-mer_cover"],
            ascending=(True, False, False, True, False, False),
        )
        .reset_index(drop=True)
    )

    for exon, group in df_sample.groupby("exon"):
        group = group.reset_index(drop=True)
        if group.empty:
            continue
        main_copy = f"{exon}_N_{group.loc[0, 'contig']}_{sample}"
        idx = 0
        while main_copy not in pairwise_distances.values and idx < len(group) - 1:
            idx += 1
            main_copy = f"{exon}_N_{group.loc[idx, 'contig']}_{sample}"
        if main_copy not in pairwise_distances.values:
            logger.debug(
                f"Main copy '{main_copy}' not found for exon '{exon}'; skipping."
            )
            continue
        # Mark main copy.
        main_row = group.iloc[[0]].copy()
        main_row.loc[0, "copy"] = "main"
        result_rows.append(main_row.iloc[0].tolist())
        paralog_found = False
        for i in range(idx + 1, len(group)):
            copy_to_compare = f"{exon}_N_{group.loc[i, 'contig']}_{sample}"
            cond1 = pairwise_distances["seq1"].eq(main_copy) & pairwise_distances[
                "seq2"
            ].eq(copy_to_compare)
            cond2 = pairwise_distances["seq2"].eq(main_copy) & pairwise_distances[
                "seq1"
            ].eq(copy_to_compare)
            seq1_main_seq2 = cond1.any()
            seq2_main_seq1 = cond2.any()
            if not seq1_main_seq2 and not seq2_main_seq1:
                continue
            if seq1_main_seq2:
                dist_series = pairwise_distances.loc[cond1, "dist"]
            else:
                dist_series = pairwise_distances.loc[cond2, "dist"]
            if dist_series.empty:
                continue
            div = dist_series.iloc[0]
            row = group.iloc[[i]].copy()
            if paralog_min_divergence < div < paralog_max_divergence:
                row.loc[0, "copy"] = "para"
                paralog_found = True
            elif div < paralog_min_divergence:
                row.loc[0, "copy"] = "main"
            result_rows.append(row.iloc[0].tolist())
        if not paralog_found:
            logger.info(f"Paralog not found in sample '{sample}' for exon '{exon}'.")

    columns = [
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
    ]
    return pd.DataFrame(result_rows, columns=columns)


# ---------------------------
# Statistics Functions
# ---------------------------
@log_exceptions
def locus_stats(df: pd.DataFrame) -> pd.DataFrame:
    stats = df.groupby(["sample", "locus"])["copy"].apply(
        lambda x: "Yes" if x.eq("para").any() else "No"
    )
    return stats.reset_index().pivot(index="sample", columns="locus", values="copy")


@log_exceptions
def exon_stats(df: pd.DataFrame) -> pd.DataFrame:
    stats = df.groupby(["sample", "exon"])["copy"].apply(
        lambda x: "Yes" if x.eq("para").any() else "No"
    )
    return stats.reset_index().pivot(index="sample", columns="exon", values="copy")


@log_exceptions
def paralog_stats(locus_statistics: pd.DataFrame) -> pd.DataFrame:
    counts = locus_statistics.replace({"Yes": 1, "No": 0}).sum(axis=1)
    stats = pd.DataFrame({"number_of_paralogous_loci": counts})
    stats["number_of_paralogous_loci"] = stats["number_of_paralogous_loci"].replace(
        {0: "0/NaN"}
    )
    return stats


# ---------------------------
# Prepare to Write Reference Functions
# ---------------------------
@log_exceptions
def prepare_to_write(df: pd.DataFrame) -> pd.DataFrame:
    """
    For each locus group, if any exon has a 'para' copy then keep only rows corresponding
    to samples with a 'para' copy for that exon; otherwise, keep all.
    """
    logger = setup_logging()
    logger.debug(f"prepare_to_write() DataFrame shape: {df.shape}")
    prepared = []
    for _, locus_df in df.groupby("locus"):
        exons = []
        for _, exon_df in locus_df.groupby("exon"):
            if exon_df["copy"].eq("para").any():
                samples_w_para = set(
                    exon_df.loc[exon_df["copy"].eq("para"), "sample"].unique()
                )
                exons.append(exon_df[exon_df["sample"].isin(samples_w_para)])
            else:
                exons.append(exon_df)
        prepared.append(pd.concat(exons))
    result = pd.concat(prepared)
    result.drop_duplicates(subset=["sequence"], inplace=True)
    result.drop_duplicates(subset=["exon", "copy"], inplace=True)
    result.reset_index(drop=True, inplace=True)
    logger.debug(f"Final shape after dropping duplicates: {result.shape}")
    return result


@log_exceptions
def create_reference_wo_paralogs(
    data_folder: str,
    all_hits: pd.DataFrame,
    blocklist: List[str],
    num_cores: int,
    log_queue,
) -> None:
    setup_logging().info("Creating customized reference without paralogs...")

    all_hits = all_hits.sort_values(
        ["exon", "pident", "qcovhsp", "evalue", "bitscore", "k-mer_cover"],
        ascending=(True, False, False, True, False, False),
    ).reset_index(drop=True)
    groups = [
        g.reset_index(drop=True) for _, g in all_hits.groupby(["sample", "locus"])
    ]
    with multiprocessing.Pool(
        processes=num_cores, initializer=worker_initializer, initargs=(log_queue,)
    ) as pool:
        results = pool.map(phase_wo_paralog, groups)
    all_hits = pd.concat(results, ignore_index=True)
    all_hits["contig"] = all_hits["saccver"].str.split("_N_").str[1]
    all_hits = all_hits[~all_hits["sample"].isin(blocklist)]
    all_hits["copy"] = "main"
    scored = score_samples(all_hits, num_cores, log_queue)
    unique_exons = scored.drop_duplicates(subset=["exon"]).reset_index(drop=True)
    unique_exons["exon_num"] = unique_exons["exon"].str.split("_").str[-1].astype(int)
    unique_exons.sort_values(["locus", "exon_num"], inplace=True)
    previous_locus = ""
    out_dir = os.path.join(data_folder, "41without_par")
    os.makedirs(out_dir, exist_ok=True)
    f_hpm = os.path.join(out_dir, "customized_reference_for_HybPhyloMaker.fas")
    f_pw_sep = os.path.join(
        out_dir, "customized_reference_for_ParalogWizard_separated_exons.fas"
    )
    f_hp_concat = os.path.join(
        out_dir, "customized_reference_for_HybPiper_concatenated_exons.fas"
    )
    with open(f_hpm, "w") as f1, open(f_pw_sep, "w") as f2, open(
        f_hp_concat, "w"
    ) as f3:
        for _, row in unique_exons.iterrows():
            sample = row["sample"]
            exon_num = row["exon"].split("_")[-1]
            # Get the locus value from the row.
            locus_val = row["locus"]
            # If the value is a float, and it is an integer, convert it to int.
            if isinstance(locus_val, float) and locus_val.is_integer():
                locus = str(int(locus_val))
            # If the value is numeric (an int or a non‑integer float), convert normally.
            elif isinstance(locus_val, (int, float)):
                locus = str(locus_val)
            # Otherwise, treat it as text.
            else:
                locus = str(locus_val)
            contig = row["contig"]
            seq = row["sequence"]
            name_hpm = f">Assembly_{locus}_Contig_{exon_num}_{sample}_N_{contig}"
            name_pw = (
                f">{sample.replace('-', '_')}_N_{contig}-" f"{locus}_exon_{exon_num}"
            )
            f1.write(f"{name_hpm}\n{seq}\n")
            f2.write(f"{name_pw}\n{seq}\n")
            if locus != previous_locus:
                f3.write(f"\n>{locus}-{locus}\n{seq}")
                previous_locus = locus
            else:
                f3.write(f"{seq}")
    setup_logging().info("Customized reference without paralogs created.")


@log_exceptions
def create_reference_w_paralogs(
    data_folder: str,
    all_hits: pd.DataFrame,
    paralog_min_divergence: float,
    paralog_max_divergence: float,
    blocklist: List[str],
    num_cores: int,
    log_queue,
) -> None:
    logger = setup_logging()
    pairwise_path = os.path.join(
        data_folder, "40aln_orth_par", "pairwise_distances.tsv"
    )
    pairwise_distances = pd.read_csv(pairwise_path, sep="\t")
    all_hits["contig"] = all_hits["saccver"].str.split("_N_").str[1]
    grouped_samples = list(all_hits.groupby("sample"))
    args = [
        (pairwise_distances, group, paralog_min_divergence, paralog_max_divergence)
        for group in grouped_samples
    ]
    with multiprocessing.Pool(
        processes=num_cores, initializer=worker_initializer, initargs=(log_queue,)
    ) as pool:
        detected = pool.starmap(detect_paralogs, args)
    all_paralogs = pd.concat(detected, ignore_index=True)
    all_paralogs.dropna(subset=["copy"], inplace=True)
    out_dir = os.path.join(data_folder, "41detected_par")
    os.makedirs(out_dir, exist_ok=True)
    all_paralogs.to_csv(
        os.path.join(out_dir, "all_paralogs_for_reference.tsv"), sep="\t", index=False
    )
    cleaned = clean_paralogs(all_paralogs, num_cores, log_queue)
    cleaned.to_csv(
        os.path.join(out_dir, "all_paralogs_for_reference_cleaned.tsv"),
        sep="\t",
        index=False,
    )
    scored = score_samples(cleaned, num_cores, log_queue)
    scored.to_csv(
        os.path.join(out_dir, "all_paralogs_for_reference_scored.tsv"),
        sep="\t",
        index=False,
    )
    scored = scored[~scored["sample"].isin(blocklist)].reset_index(drop=True)
    to_write = prepare_to_write(scored)
    ref_path = os.path.join(
        out_dir,
        f"customized_reference_div_{paralog_min_divergence}_"
        f"{paralog_max_divergence}.fas",
    )
    with open(ref_path, "w") as f:
        for _, row in to_write.iterrows():
            sample = row["sample"]
            exon_num = row["exon"].split("_")[-1]
            # Get the locus value from the row.
            locus_val = row["locus"]
            # If the value is a float, and it is an integer, convert it to int.
            if isinstance(locus_val, float) and locus_val.is_integer():
                locus = str(int(locus_val))
            # If the value is numeric (an int or a non‑integer float), convert normally.
            elif isinstance(locus_val, (int, float)):
                locus = str(locus_val)
            # Otherwise, treat it as text.
            else:
                locus = str(locus_val)
            contig = row["contig"]
            copy = row["copy"]
            seq = row["sequence"]
            copy_str = "" if copy == "main" else copy
            name_seq = (
                f">Assembly_{locus}{copy_str}_Contig_{exon_num}_" f"{sample}_N_{contig}"
            )
            f.write(f"{name_seq}\n{seq}\n")
    locus_statistics = locus_stats(all_paralogs)
    locus_statistics.to_csv(
        os.path.join(
            out_dir,
            f"locus_statistics_div_{paralog_min_divergence}_"
            f"{paralog_max_divergence}.tsv",
        ),
        sep="\t",
        na_rep="NaN",
    )
    exon_statistics = exon_stats(all_paralogs)
    exon_statistics.to_csv(
        os.path.join(
            out_dir,
            f"exon_statistics_div_{paralog_min_divergence}_"
            f"{paralog_max_divergence}.tsv",
        ),
        sep="\t",
        na_rep="NaN",
    )
    paralog_statistics = paralog_stats(locus_statistics)
    paralog_statistics.to_csv(
        os.path.join(
            out_dir,
            f"paralog_statistics_div_{paralog_min_divergence}_"
            f"{paralog_max_divergence}.tsv",
        ),
        sep="\t",
        na_rep="NaN",
    )
    logger.info("Customized reference with paralogs created.")
