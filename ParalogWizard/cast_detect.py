#!/usr/bin/env python
"""
------------------------------------------------------------------------------------------------------------------------
Copyright 2024 Roman Ufimov under the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later version.

cast_detect — assign main/para copies and write customized references.

  * Without -p: cluster by percent identity per sample×locus → 41without_par/.
  * With -p: use pairwise distances from cast_analyze to label paralogs within
    divergence bounds → 41detected_par/.

Each run clears and rebuilds the output directory.
------------------------------------------------------------------------------------------------------------------------
"""

from __future__ import annotations

import logging
import os
import shutil
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd

from ParalogWizard import log_exceptions, managed_pool
from ParalogWizard.cast_call import require_dir, require_file

logger = logging.getLogger("ParalogWizard")

DistanceIndex = Dict[Tuple[str, str], float]


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
def _clear_outdir(path: str) -> str:
    if os.path.isdir(path):
        shutil.rmtree(path)
        logger.debug("Removed previous %s", path)
    os.makedirs(path, exist_ok=True)
    return path


def _contig_from_saccver(saccver: str) -> str:
    parts = str(saccver).split("_N_")
    return parts[1] if len(parts) > 1 else str(saccver)


def find_cluster(array: List[float]) -> float:
    """Boundary at the midpoint of the largest gap in a sorted 1-D array."""
    if len(array) < 2:
        raise ValueError("find_cluster requires at least 2 values")
    ordered = sorted(array)
    gaps = [ordered[i] - ordered[i - 1] for i in range(1, len(ordered))]
    biggest_gap = max(gaps)
    biggest_gap_index = gaps.index(biggest_gap)
    boundary = ordered[biggest_gap_index] + biggest_gap / 2.0
    logger.debug("find_cluster: boundary=%.6f (n=%d)", boundary, len(ordered))
    return boundary


def build_distance_index(
    pairwise_distances: pd.DataFrame,
) -> Tuple[DistanceIndex, Set[str]]:
    """
    Build bidirectional (seq_a, seq_b) → dist lookup and the set of sequence names
    present in the distance table (including single-copy placeholder rows).
    """
    require_columns = {"seq1", "dist", "seq2"}
    missing = require_columns - set(pairwise_distances.columns)
    if missing:
        raise ValueError(
            f"pairwise_distances.tsv missing column(s): {', '.join(sorted(missing))}"
        )
    index: DistanceIndex = {}
    names: Set[str] = set()
    for row in pairwise_distances.itertuples(index=False):
        s1 = row.seq1
        if isinstance(s1, str) and s1:
            names.add(s1)
        s2 = row.seq2
        s2_ok = isinstance(s2, str) and bool(s2) and s2.lower() != "nan"
        if s2_ok:
            names.add(s2)
        if pd.isna(row.dist) or not s2_ok:
            continue
        d = float(row.dist)
        index[(s1, s2)] = d
        index[(s2, s1)] = d
    logger.info(
        "Distance index: %d directed pair(s), %d sequence name(s)",
        len(index),
        len(names),
    )
    return index, names


# -----------------------------------------------------------------------------
# Cleaning / scoring
# -----------------------------------------------------------------------------
@log_exceptions
def adjust_orphaned_main(sample_locus_dataframe: pd.DataFrame) -> pd.DataFrame:
    """Reassign orphaned exons using a pident gap boundary when paras exist."""
    df = sample_locus_dataframe
    if "para" not in df["copy"].values:
        return df

    exons_needed_clustering = [
        exon
        for exon, g in df.groupby("exon", sort=False)
        if "para" not in g["copy"].values
    ]
    if not exons_needed_clustering:
        return df

    ident_array = df.loc[~df["exon"].isin(exons_needed_clustering), "pident"].tolist()
    if len(ident_array) < 2:
        logger.warning(
            "Not enough pident values to cluster orphaned exons for a sample×locus; "
            "leaving labels unchanged"
        )
        return df

    boundary = find_cluster(ident_array)
    pidents_main = df.loc[
        (df["copy"] == "main") & (~df["exon"].isin(exons_needed_clustering)),
        "pident",
    ]
    pidents_para = df.loc[df["copy"] == "para", "pident"]
    if (pidents_main < boundary).any() or (pidents_para > boundary).any():
        return df.loc[~df["exon"].isin(exons_needed_clustering)].copy()

    out = df.copy()
    mask = out["exon"].isin(exons_needed_clustering)
    out.loc[mask & (out["pident"] > boundary), "copy"] = "main"
    out.loc[mask & (out["pident"] < boundary), "copy"] = "para"
    return out


@log_exceptions
def clean_paralogs(dataframe: pd.DataFrame, n_cpu: int, log_queue) -> pd.DataFrame:
    """Adjust orphaned mains in parallel; drop contigs labeled both main and para."""
    logger.info("Cleaning paralog labels")
    groups = [g.reset_index(drop=True) for _, g in dataframe.groupby(["sample", "locus"])]
    n_workers = max(1, min(int(n_cpu), len(groups))) if groups else 1
    if groups:
        with managed_pool(n_workers, log_queue) as pool:
            results = pool.map(adjust_orphaned_main, groups)
        dataframe = pd.concat(results, ignore_index=True)
    else:
        dataframe = dataframe.copy()

    copy_sets = dataframe.groupby("saccver")["copy"].agg(lambda s: set(s.dropna()))
    ambiguous = copy_sets.index[
        copy_sets.map(lambda s: "main" in s and "para" in s)
    ]
    if len(ambiguous):
        logger.info("Dropping %d ambiguous contig(s) labeled both main and para", len(ambiguous))
        dataframe = dataframe.loc[~dataframe["saccver"].isin(ambiguous)].reset_index(
            drop=True
        )
    logger.info("Paralog cleaning complete (%d row(s))", len(dataframe))
    return dataframe


@log_exceptions
def score_samples(dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    Vectorized scores:
      score_1 = mean main-copy pident per sample×locus
      score_2 = n unique (exon, copy) per sample×locus
      score_3 = n unique exons per saccver
    """
    logger.info("Scoring samples (%d row(s))", len(dataframe))
    df = dataframe.sort_values(
        ["exon", "pident", "qcovhsp", "evalue", "bitscore", "k-mer_cover"],
        ascending=(True, False, False, True, False, False),
    ).reset_index(drop=True)

    score_1 = (
        df.loc[df["copy"] == "main"]
        .groupby(["sample", "locus"], sort=False)["pident"]
        .mean()
    )
    score_2 = (
        df.groupby(["sample", "locus"], sort=False)[["exon", "copy"]]
        .apply(lambda g: int(g.drop_duplicates().shape[0]))
    )
    df = df.copy()
    idx = pd.MultiIndex.from_frame(df[["sample", "locus"]])
    df["score_1"] = idx.map(score_1)
    df["score_2"] = idx.map(score_2)
    df["score_3"] = df.groupby("saccver", sort=False)["exon"].transform("nunique")

    df.sort_values(
        ["exon", "copy", "score_2", "score_1", "sample", "score_3"],
        ascending=[True, True, False, False, True, False],
        inplace=True,
    )
    df.reset_index(drop=True, inplace=True)
    logger.info("Scoring complete")
    return df


@log_exceptions
def phase_wo_paralog(sample_locus_dataframe: pd.DataFrame) -> pd.DataFrame:
    """Keep hits above the largest pident gap (no-paralog mode)."""
    ident_array = sample_locus_dataframe["pident"].tolist()
    if len(ident_array) < 2:
        return sample_locus_dataframe
    boundary = find_cluster(ident_array)
    return sample_locus_dataframe.loc[
        sample_locus_dataframe["pident"] > boundary
    ].copy()


# -----------------------------------------------------------------------------
# Paralog detection
# -----------------------------------------------------------------------------
@log_exceptions
def detect_paralogs(
    sample: str,
    sample_hits: pd.DataFrame,
    dist_index: DistanceIndex,
    names_in_distances: Set[str],
    paralog_min_divergence: float,
    paralog_max_divergence: float,
) -> pd.DataFrame:
    """
    Label main/para copies for one sample using the distance index.
    Returns a DataFrame with a 'copy' column (rows without a usable main are omitted).
    """
    grouped_samples = sample_hits.sort_values(
        ["exon", "pident", "qcovhsp", "evalue", "bitscore", "k-mer_cover"],
        ascending=(True, False, False, True, False, False),
    ).reset_index(drop=True)

    rows: List[pd.Series] = []
    n_para = 0
    for exon, group_exon_dataframe in grouped_samples.groupby("exon", sort=False):
        group_exon_dataframe = group_exon_dataframe.reset_index(drop=True)
        len_df = len(group_exon_dataframe)
        idx = 0
        main_copy = f"{exon}_N_{group_exon_dataframe.loc[0, 'contig']}_{sample}"
        while main_copy not in names_in_distances:
            if idx + 1 >= len_df:
                break
            idx += 1
            main_copy = (
                f"{exon}_N_{group_exon_dataframe.loc[idx, 'contig']}_{sample}"
            )
        if main_copy not in names_in_distances:
            continue

        main_entry = group_exon_dataframe.iloc[idx].copy()
        main_entry["copy"] = "main"
        rows.append(main_entry)

        paralog_found = False
        for index in range(idx + 1, len_df):
            copy_to_compare = (
                f"{exon}_N_{group_exon_dataframe.loc[index, 'contig']}_{sample}"
            )
            div = dist_index.get((main_copy, copy_to_compare))
            if div is None:
                continue
            secondary = group_exon_dataframe.iloc[index].copy()
            if paralog_min_divergence < div < paralog_max_divergence:
                secondary["copy"] = "para"
                paralog_found = True
                n_para += 1
            elif div < paralog_min_divergence:
                secondary["copy"] = "main"
            else:
                # above max divergence: leave unlabeled (dropped later via dropna)
                continue
            rows.append(secondary)
        if paralog_found:
            logger.debug("Paralog(s) in %s / %s", sample, exon)

    if not rows:
        logger.debug("No labeled copies for sample %s", sample)
        return pd.DataFrame()

    out = pd.DataFrame(rows).reset_index(drop=True)
    logger.debug(
        "detect_paralogs %s: %d row(s), %d para label(s)", sample, len(out), n_para
    )
    return out


@log_exceptions
def locus_stats(all_paralogs_for_reference: pd.DataFrame) -> pd.DataFrame:
    logger.info("Computing locus statistics")
    data: Dict[str, Dict[str, str]] = {}
    for (sample, locus), g in all_paralogs_for_reference.groupby(
        ["sample", "locus"], sort=False
    ):
        data.setdefault(sample, {})[locus] = (
            "Yes" if "para" in g["copy"].values else "No"
        )
    locus_statistics = pd.DataFrame.from_dict(data, orient="index")
    locus_statistics.index.name = r"samples\locus"
    return locus_statistics


def exon_stats(all_paralogs_for_reference: pd.DataFrame) -> pd.DataFrame:
    logger.info("Computing exon statistics")
    data: Dict[str, Dict[str, str]] = {}
    for (sample, exon), g in all_paralogs_for_reference.groupby(
        ["sample", "exon"], sort=False
    ):
        data.setdefault(sample, {})[exon] = (
            "Yes" if "para" in g["copy"].values else "No"
        )
    return pd.DataFrame.from_dict(data, orient="index")


@log_exceptions
def paralog_stats(locus_statistics: pd.DataFrame) -> pd.DataFrame:
    logger.info("Computing paralog statistics")
    numeric = locus_statistics.copy()
    numeric = numeric.where(numeric != "Yes", 1)
    numeric = numeric.where(numeric != "No", 0)
    numeric = numeric.apply(pd.to_numeric, errors="coerce")
    counts = numeric.sum(axis=1, skipna=True)
    out = counts.to_frame(name="number_of_paralogous_loci")
    out.index.name = r"samples\locus"
    out["number_of_paralogous_loci"] = out["number_of_paralogous_loci"].replace(
        {0: "0/NaN", 0.0: "0/NaN"}
    )
    return out


@log_exceptions
def prepare_to_write(all_paralogs_for_reference_scored: pd.DataFrame) -> pd.DataFrame:
    """Prefer samples that carry a para at exons where paras exist."""
    logger.info("Preparing reference sequences for writing")
    prepared_loci = []
    for _, locus_dataframe in all_paralogs_for_reference_scored.groupby(
        "locus", sort=False
    ):
        prepared_exons = []
        for _, exon_dataframe in locus_dataframe.groupby("exon", sort=False):
            if (exon_dataframe["copy"] == "para").any():
                samples_w_para = set(
                    exon_dataframe.loc[exon_dataframe["copy"] == "para", "sample"]
                )
                prepared_exons.append(
                    exon_dataframe.loc[exon_dataframe["sample"].isin(samples_w_para)]
                )
            else:
                prepared_exons.append(exon_dataframe)
        prepared_loci.append(pd.concat(prepared_exons, ignore_index=True))
    to_write = pd.concat(prepared_loci, ignore_index=True)
    to_write = to_write.drop_duplicates(subset=["sequence"])
    to_write = to_write.drop_duplicates(subset=["exon", "copy"]).reset_index(drop=True)
    logger.info("Reference write set: %d sequence(s)", len(to_write))
    return to_write


# -----------------------------------------------------------------------------
# Reference writers
# -----------------------------------------------------------------------------
@log_exceptions
def create_reference_wo_paralogs(
    data_folder: str,
    all_hits_for_reference: pd.DataFrame,
    blocklist: Set[str],
    num_cores: int,
    log_queue,
) -> None:
    """Build customized references without paralog labeling → 41without_par/."""
    out_dir = _clear_outdir(os.path.join(data_folder, "41without_par"))
    logger.info("Creating reference WITHOUT paralogs → %s", out_dir)

    hits = all_hits_for_reference.sort_values(
        ["exon", "pident", "qcovhsp", "evalue", "bitscore", "k-mer_cover"],
        ascending=(True, False, False, True, False, False),
    ).reset_index(drop=True)
    groups = [g.reset_index(drop=True) for _, g in hits.groupby(["sample", "locus"])]
    n_workers = max(1, min(int(num_cores), len(groups))) if groups else 1
    if not groups:
        raise RuntimeError("No sample×locus groups in all_hits.tsv")

    with managed_pool(n_workers, log_queue) as pool:
        results = pool.map(phase_wo_paralog, groups)
    hits = pd.concat(results, ignore_index=True)
    hits["contig"] = hits["saccver"].map(_contig_from_saccver)
    hits = hits.loc[~hits["sample"].isin(blocklist)].copy()
    hits["copy"] = "main"
    scored = score_samples(hits)
    to_write = scored.drop_duplicates(subset=["exon"]).reset_index(drop=True)
    to_write["exon_num"] = to_write["exon"].str.split("_").str[-1].astype(int)
    to_write = to_write.sort_values(["locus", "exon_num"]).reset_index(drop=True)

    path_hpm = os.path.join(out_dir, "customized_reference_for_HybPhyloMaker.fas")
    path_pw = os.path.join(
        out_dir, "customized_reference_for_ParalogWizard_separated_exons.fas"
    )
    path_hp = os.path.join(
        out_dir, "customized_reference_for_HybPiper_concatenate_exons.fas"
    )
    # Keep legacy filename typo variant as a copy of the corrected name for compatibility
    path_hp_legacy = os.path.join(
        out_dir, "customized_reference_for_HybPiper_concatenated_exons.fas"
    )

    previous_locus = None
    with open(path_hpm, "w") as hpm, open(path_pw, "w") as pw, open(path_hp, "w") as hp:
        for rec in to_write.to_dict(orient="records"):
            sample = rec["sample"]
            exon_num = str(rec["exon"]).split("_")[-1]
            locus = rec["locus"]
            contig = rec["contig"]
            seq = rec["sequence"]
            hpm.write(
                f">Assembly_{locus}_Contig_{exon_num}_{sample}_N_{contig}\n{seq}\n"
            )
            pw.write(
                f">{sample.replace('-', '_')}_N_{contig}-{locus}_exon_{exon_num}\n{seq}\n"
            )
            header = f">{locus}-{locus}"
            if locus != previous_locus:
                if previous_locus is None:
                    hp.write(f"{header}\n{seq}")
                else:
                    hp.write(f"\n{header}\n{seq}")
                previous_locus = locus
            else:
                hp.write(seq)

    shutil.copyfile(path_hp, path_hp_legacy)
    logger.info(
        "Wrote references without paralogs (%d exon(s)) under %s",
        len(to_write),
        out_dir,
    )


@log_exceptions
def create_reference_w_paralogs(
    data_folder: str,
    all_hits_for_reference: pd.DataFrame,
    paralog_min_divergence: float,
    paralog_max_divergence: float,
    blocklist: Set[str],
    num_cores: int,
    log_queue,
) -> None:
    """Detect paralogs from distances and write customized reference → 41detected_par/."""
    out_dir = _clear_outdir(os.path.join(data_folder, "41detected_par"))
    logger.info(
        "Creating reference WITH paralogs (min=%.3f, max=%.3f) → %s",
        paralog_min_divergence,
        paralog_max_divergence,
        out_dir,
    )

    dist_path = require_file(
        os.path.join(data_folder, "40aln_orth_par", "pairwise_distances.tsv"),
        "pairwise_distances.tsv",
    )
    pairwise_distances = pd.read_csv(dist_path, sep="\t")
    dist_index, names_in_distances = build_distance_index(pairwise_distances)

    hits = all_hits_for_reference.copy()
    hits["contig"] = hits["saccver"].map(_contig_from_saccver)
    sample_groups = [(sample, g.copy()) for sample, g in hits.groupby("sample", sort=False)]
    if not sample_groups:
        raise RuntimeError("No samples in all_hits.tsv")

    n_workers = max(1, min(int(num_cores), len(sample_groups)))
    logger.info(
        "Paralog detection: %d sample(s), %d worker(s)",
        len(sample_groups),
        n_workers,
    )
    with managed_pool(n_workers, log_queue) as pool:
        async_results = [
            (
                sample,
                pool.apply_async(
                    detect_paralogs,
                    (
                        sample,
                        sample_df,
                        dist_index,
                        names_in_distances,
                        paralog_min_divergence,
                        paralog_max_divergence,
                    ),
                ),
            )
            for sample, sample_df in sample_groups
        ]
        failures: List[Tuple[str, Exception]] = []
        results: List[pd.DataFrame] = []
        total = len(async_results)
        for i, (sample, async_result) in enumerate(async_results, start=1):
            try:
                results.append(async_result.get())
                if i == total or i % 10 == 0:
                    logger.info("Detect progress: %d / %d", i, total)
            except Exception as e:
                logger.error("detect_paralogs failed for %s: %s", sample, e)
                failures.append((sample, e))
    if failures:
        preview = ", ".join(s for s, _ in failures[:20])
        raise RuntimeError(
            f"cast_detect aborted: {len(failures)} sample(s) failed ({preview}). See log."
        )

    frames = [r for r in results if r is not None and not r.empty]
    if not frames:
        raise RuntimeError("Paralog detection produced no labeled copies.")
    all_paralogs = pd.concat(frames, ignore_index=True)
    all_paralogs.dropna(subset=["copy"], inplace=True)

    all_paralogs.to_csv(
        os.path.join(out_dir, "all_paralogs_for_reference.tsv"), sep="\t", index=False
    )
    cleaned = clean_paralogs(all_paralogs, num_cores, log_queue)
    cleaned.to_csv(
        os.path.join(out_dir, "all_paralogs_for_reference_cleaned.tsv"),
        sep="\t",
        index=False,
    )
    scored = score_samples(cleaned)
    scored.to_csv(
        os.path.join(out_dir, "all_paralogs_for_reference_scored.tsv"),
        sep="\t",
        index=False,
    )
    scored = scored.loc[~scored["sample"].isin(blocklist)].reset_index(drop=True)
    to_write = prepare_to_write(scored)

    ref_path = os.path.join(
        out_dir,
        f"customized_reference_div_{paralog_min_divergence}_{paralog_max_divergence}.fas",
    )
    with open(ref_path, "w") as customized_reference:
        for rec in to_write.to_dict(orient="records"):
            sample = rec["sample"]
            exon_num = str(rec["exon"]).split("_")[-1]
            locus = rec["locus"]
            contig = rec["contig"]
            copy = rec["copy"]
            seq = rec["sequence"]
            copy_tag = "" if copy == "main" else str(copy)
            customized_reference.write(
                f">Assembly_{locus}{copy_tag}_Contig_{exon_num}_{sample}_N_{contig}\n{seq}\n"
            )

    loc_stats = locus_stats(all_paralogs)
    loc_stats.to_csv(
        os.path.join(
            out_dir,
            f"locus_statistics_div_{paralog_min_divergence}_{paralog_max_divergence}.tsv",
        ),
        sep="\t",
        na_rep="NaN",
    )
    ex_stats = exon_stats(all_paralogs)
    ex_stats.to_csv(
        os.path.join(
            out_dir,
            f"exon_statistics_div_{paralog_min_divergence}_{paralog_max_divergence}.tsv",
        ),
        sep="\t",
        na_rep="NaN",
    )
    par_stats = paralog_stats(loc_stats)
    par_stats.to_csv(
        os.path.join(
            out_dir,
            f"paralog_statistics_div_{paralog_min_divergence}_{paralog_max_divergence}.tsv",
        ),
        sep="\t",
        na_rep="NaN",
    )
    logger.info(
        "Wrote paralog reference (%d sequence(s)) and statistics under %s",
        len(to_write),
        out_dir,
    )
