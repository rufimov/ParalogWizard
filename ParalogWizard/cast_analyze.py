#!/usr/bin/env python
"""
Module for ParalogWizard cast_analyze command.

This module performs several analyses:
  - It aligns sequences using MAFFT.
  - It builds phylogenetic trees using FastTree.
  - It calculates pairwise percent dissimilarities,
    estimates divergence using kernel density estimation,
    and produces plots of the distribution.
Multiprocessing is used where appropriate, and detailed logging is provided.
"""

import glob
import itertools
import multiprocessing
import os
import random
import subprocess
from glob import glob
from typing import Dict, List, Set

import matplotlib
import numpy as np
import pandas as pd
from Bio import SeqIO, SeqRecord
from matplotlib import pyplot as plt
from scipy.signal import argrelextrema
from sklearn.mixture import BayesianGaussianMixture
from sklearn.neighbors import KernelDensity

# Use our setup_logging from ParalogWizard.
from ParalogWizard import setup_logging, worker_initializer

# -----------------------------------------------------------------------------
# Module-level logger
# -----------------------------------------------------------------------------
logger = setup_logging()


# -----------------------------------------------------------------------------
# Logging Decorator (with detailed context)
# -----------------------------------------------------------------------------
def log_exceptions(func):
    """
    Decorator that logs function entry, exit, and any exceptions.
    It logs the function's name along with its positional and keyword arguments.
    This way you know which file, sample, locus, etc., caused the error.
    Instead of calling sys.exit(), it re-raises the exception so that the main process
    can catch the error and shut down the pool gracefully.
    """
    from functools import wraps

    @wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as e:
            arg_str = ", ".join(str(arg) for arg in args)
            kwarg_str = ", ".join(f"{k}={v}" for k, v in kwargs.items())
            logger.exception(
                f"Exception in {func.__name__} (args: {arg_str}; kwargs: {kwarg_str}): {e}"
            )
            raise

    return wrapper


# -----------------------------------------------------------------------------
# Alignment and Tree Building Functions (Using subprocess)
# -----------------------------------------------------------------------------
@log_exceptions
def mafft_align(file: str) -> None:
    """
    Aligns a FASTA file using MAFFT via a subprocess call.
    The output is saved as <original_basename>.fasta.mafft.
    """
    if not os.path.exists(file):
        logger.error("Input file %s does not exist", file)
        raise FileNotFoundError(f"Input file {file} does not exist")
    out_file = f"{os.path.splitext(file)[0]}.fasta.mafft"
    cmd = ["mafft", "--auto", file]
    logger.info("Running MAFFT on %s", file)
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        with open(out_file, "w") as aligned:
            aligned.write(result.stdout)
        logger.info(
            "MAFFT alignment completed for %s, output saved to %s", file, out_file
        )
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        logger.error("MAFFT alignment failed for %s: %s", file, e)
        raise


@log_exceptions
def fast_tree(file: str) -> None:
    """
    Builds a tree using FastTree. First tries fasttreemp; if it fails (or isn’t found),
    falls back to fasttree.
    The tree is saved as <file>.tre.
    """
    if not os.path.exists(file):
        logger.error("Input file %s does not exist", file)
        raise FileNotFoundError(f"Input file {file} does not exist")
    out_file = f"{file}.tre"
    # Try fasttreemp first.
    cmd_mp = ["fasttreemp", "-nt", "-gtr", "-fastest", file]
    logger.info("Running FastTreeMP on %s", file)
    try:
        with open(out_file, "w") as outf:
            subprocess.run(cmd_mp, stdout=outf, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        logger.warning("FastTreeMP failed for %s (%s), trying FastTree", file, e)
        cmd = ["fasttree", "-nt", "-gtr", "-fastest", file]
        try:
            with open(out_file, "w") as outf:
                subprocess.run(cmd, stdout=outf, check=True)
        except subprocess.CalledProcessError as e2:
            logger.error("FastTree failed for %s with error: %s", file, e2)
            raise
    logger.info("Tree built for %s, output saved to %s", file, out_file)


@log_exceptions
def build_alignments(data_folder: str, n_cpu: int, log_queue) -> None:
    """
    Builds exon alignments from the 31exonic_contigs/all_hits.tsv file.
    For each exon, creates a FASTA file with all sequences, then aligns with MAFFT,
    and finally builds a tree with FastTree.

    Files are processed individually to avoid a single hanging file from blocking the pool.
    """
    logger.info("Building individual exon alignments...")
    hits_file = os.path.join(data_folder, "31exonic_contigs", "all_hits.tsv")
    if not os.path.exists(hits_file):
        logger.error("Hits file %s does not exist", hits_file)
        raise FileNotFoundError(f"Hits file {hits_file} does not exist")
    all_hits = pd.read_csv(hits_file, sep="\t")
    grouped_exons = all_hits.groupby("exon")
    aln_folder = os.path.join(data_folder, "40aln_orth_par")
    os.makedirs(aln_folder, exist_ok=True)
    # Write FASTA files for each exon.
    for exon, df in grouped_exons:
        exon_str = str(exon)
        out_path = os.path.join(aln_folder, f"{exon_str}.fasta")
        logger.info("Creating alignment FASTA for exon %s", exon_str)
        with open(out_path, "w") as f_out:
            for _, row in df.iterrows():
                # Construct a sequence header from exon, contig and sample.
                try:
                    contig_part = row["saccver"].split("_N_")[1]
                except IndexError:
                    contig_part = "unknown"
                header = f"{exon_str}_N_{contig_part}_{row['sample']}"
                seq = row["sequence"]
                f_out.write(f">{header}\n{seq}\n")
        logger.debug("FASTA for exon %s written", exon_str)

    # Process MAFFT alignments individually.
    files_to_align = glob(os.path.join(aln_folder, "*fasta"))
    logger.info("Running MAFFT alignment on %d files", len(files_to_align))
    mafft_failed_files = []
    with multiprocessing.Pool(
        processes=n_cpu, initializer=worker_initializer, initargs=(log_queue,)
    ) as pool:
        for file in files_to_align:
            async_result = pool.apply_async(mafft_align, (file,))
            try:
                async_result.get(timeout=600)  # timeout in seconds per file
            except multiprocessing.TimeoutError:
                logger.error("Timeout during MAFFT alignment for file: %s", file)
                mafft_failed_files.append(file)
            except Exception as e:
                logger.error(
                    "Error during MAFFT alignment for file: %s; Error: %s", file, e
                )
                mafft_failed_files.append(file)
        pool.close()
        pool.join()
    if mafft_failed_files:
        logger.error(
            "MAFFT alignment failed for these files and they will be skipped: %s",
            mafft_failed_files,
        )
    else:
        logger.info("MAFFT alignment completed for all files.")

    # Process FastTree individually on the MAFFT output files.
    files_to_tree = glob(os.path.join(aln_folder, "*fasta.mafft"))
    logger.info("Running FastTree on %d alignment files", len(files_to_tree))
    fasttree_failed_files = []
    with multiprocessing.Pool(
        processes=n_cpu, initializer=worker_initializer, initargs=(log_queue,)
    ) as pool:
        for file in files_to_tree:
            async_result = pool.apply_async(fast_tree, (file,))
            try:
                async_result.get(timeout=600)
            except multiprocessing.TimeoutError:
                logger.error("Timeout during FastTree tree building for file: %s", file)
                fasttree_failed_files.append(file)
            except Exception as e:
                logger.error(
                    "Error during FastTree tree building for file: %s; Error: %s",
                    file,
                    e,
                )
                fasttree_failed_files.append(file)
        pool.close()
        pool.join()
    if fasttree_failed_files:
        logger.error(
            "FastTree building failed for these files and they will be skipped: %s",
            fasttree_failed_files,
        )
    else:
        logger.info("FastTree building completed for all files.")

    logger.info("Alignment and tree building complete.")


# -----------------------------------------------------------------------------
# Divergence Analysis Functions
# -----------------------------------------------------------------------------
@log_exceptions
def percent_dissimilarity(seq1: str, seq2: str) -> float or None:
    """
    Compute the percent dissimilarity between two sequences.
    Positions with gaps in both sequences are ignored.
    Returns None if the alignment is too short.
    """
    seq1 = seq1.lower()
    seq2 = seq2.lower()
    filtered = [(a, b) for a, b in zip(seq1, seq2) if not (a == "-" and b == "-")]
    if not filtered:
        return None
    seq1_list, seq2_list = zip(*filtered)
    seq1_list = list(seq1_list)
    seq2_list = list(seq2_list)
    n = len(seq1_list)
    left = 0
    while left < n and (seq1_list[left] == "-" or seq2_list[left] == "-"):
        left += 1
    right = n
    while right > left and (seq1_list[right - 1] == "-" or seq2_list[right - 1] == "-"):
        right -= 1
    overlap = right - left
    if (overlap < 100) or (overlap / (n - left) < 0.5) or (overlap / right < 0.5):
        return None
    mismatches = sum(
        1
        for a, b in zip(seq1_list[left:right], seq2_list[left:right])
        if a != b and a != "-" and b != "-"
    )
    return (mismatches / overlap) * 100


@log_exceptions
def get_distance_matrix(file_to_process: str, blocklist: Set[str]):
    """
    For the sequences in the given FASTA file, compute a distance matrix based on percent dissimilarity.
    Also, detect local minima in the divergence distribution.
    Returns a tuple: (list of cluster means, DataFrame of pairwise distances, plotting data).
    """
    if not os.path.exists(file_to_process):
        logger.error("Input file %s does not exist", file_to_process)
        raise FileNotFoundError(f"Input file {file_to_process} does not exist")
    current_matrix_to_plot = []
    current_matrix_to_write = []
    sum_list = []
    logger.info("Computing distance matrix for %s", file_to_process)
    with open(file_to_process) as fasta_file:
        sequences: Dict[str, SeqRecord.SeqRecord] = SeqIO.to_dict(
            SeqIO.parse(fasta_file, "fasta")
        )
    seq_names = pd.DataFrame(list(sequences.keys()), columns=["full_name"])
    seq_names["sample_name"] = (
        seq_names["full_name"].str.split("_").str[-2:].str.join("_")
    )
    duplicated = set(
        seq_names[seq_names.duplicated(subset="sample_name")]["sample_name"].values
    )
    non_duplicated = set(
        seq_names[~seq_names["sample_name"].isin(duplicated)]["sample_name"].values
    )
    for samp in duplicated:
        names = seq_names[seq_names["sample_name"] == samp]["full_name"].tolist()
        for pair in itertools.combinations(names, 2):
            d = percent_dissimilarity(
                str(sequences[pair[0]].seq), str(sequences[pair[1]].seq)
            )
            if d is None:
                continue
            if samp not in blocklist:
                current_matrix_to_plot.append(d)
            current_matrix_to_write.append([pair[0], d, pair[1]])
    for samp in non_duplicated:
        names = seq_names[seq_names["sample_name"] == samp]["full_name"].tolist()
        for name in names:
            seq = str(sequences[name].seq)
            if len(seq) < 100:
                continue
            current_matrix_to_write.append([name, np.nan, np.nan])
    if current_matrix_to_plot:
        sorted_vals = np.sort(np.array(current_matrix_to_plot))
        arr = sorted_vals.reshape(-1, 1)
        kde = KernelDensity(kernel="gaussian", bandwidth=1.5)
        kde.fit(arr)
        space = np.linspace(-1, np.max(arr) + 1, 1000)
        e = np.exp(kde.score_samples(space.reshape(-1, 1)))
        mi = argrelextrema(e, np.less)[0]
        minimum = space[mi]
        if len(minimum) == 0:
            sum_list.append(float(np.mean(sorted_vals)))
        else:
            for i in range(len(minimum)):
                if i == 0:
                    idx = np.where((sorted_vals < minimum[i]) & (sorted_vals >= 0))[0]
                else:
                    idx = np.where(
                        (sorted_vals < minimum[i]) & (sorted_vals > minimum[i - 1])
                    )[0]
                if len(idx) == 0:
                    continue
                cluster = sorted_vals[min(idx) : max(idx) + 1]
                sum_list.append(float(np.mean(cluster)))
            idx_last = np.where(sorted_vals > minimum[-1])[0]
            if len(idx_last):
                last_cluster = sorted_vals[min(idx_last) : max(idx_last) + 1]
                sum_list.append(float(np.mean(last_cluster)))
        to_plot = [space, e, sorted_vals]
    else:
        to_plot = None
    matrix_df = pd.DataFrame(current_matrix_to_write, columns=["seq1", "dist", "seq2"])
    logger.info("Distance matrix computed for %s", file_to_process)
    return sum_list, matrix_df, to_plot


@log_exceptions
def get_model(array: np.ndarray, num_comp: int) -> BayesianGaussianMixture:
    """
    Fit and return a BayesianGaussianMixture model for the given array.
    """
    model = BayesianGaussianMixture(n_components=num_comp, max_iter=10000, n_init=10)
    model.fit(array)
    logger.debug("BayesianGaussianMixture model fitted with %d components", num_comp)
    return model


@log_exceptions
def plot_vertical_line(
    ax: plt.Axes, line_name: str, line_value: float, colors: List[str], idx: int
) -> None:
    """
    Plot a vertical line on the given axis with a label.
    """
    ax.axvline(
        x=line_value,
        label=f"{line_name} - {np.round(line_value, 2)}",
        color=colors[idx],
        lw=0.6,
        ls="--",
    )
    logger.debug("Plotted vertical line %s at %.2f", line_name, line_value)


@log_exceptions
def get_plot(path: str, name: str, matrix: np.ndarray, comp: int) -> None:
    """
    Generate and save plots of divergence distributions using the fitted model.
    """
    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)
        logger.info("Created output directory %s", path)
    matplotlib.use("Agg")
    fig, ax = plt.subplots(figsize=(15, 10))
    ax.set_xlabel("Divergence (%)")
    ax.set_ylabel("Frequency")
    fig.suptitle(name)
    max_val = np.round(np.max(matrix)) + 1
    ax.xaxis.set_ticks(np.arange(0, max_val, 1))
    model = get_model(matrix, comp)
    means = model.means_.flatten().tolist()
    sigmas = [np.sqrt(x) for x in model.covariances_.flatten().tolist()]
    mu_sigma = sorted(zip(means, sigmas), key=lambda x: x[0])
    colors = random.sample(
        [
            "#FF0000",
            "#008000",
            "#0000FF",
            "#FF00FF",
            "#800000",
            "#808000",
            "#008080",
            "#2196F3",
            "#4A148C",
            "#512DA8",
            "#2ECC71",
        ],
        k=comp * 3,
    )
    if mu_sigma:
        first_peak, first_sigma = mu_sigma[0]
        plot_vertical_line(
            ax, "first_peak_minus_sigma", first_peak - first_sigma, colors, 0
        )
        plot_vertical_line(ax, "first_peak", first_peak, colors, 1)
        plot_vertical_line(
            ax, "first_peak_plus_sigma", first_peak + first_sigma, colors, 2
        )
        if comp > 1 and len(mu_sigma) > 1:
            second_peak, second_sigma = mu_sigma[1]
            plot_vertical_line(
                ax, "second_peak_minus_sigma", second_peak - second_sigma, colors, 3
            )
            plot_vertical_line(ax, "second_peak", second_peak, colors, 4)
            plot_vertical_line(
                ax, "second_peak_plus_sigma", second_peak + second_sigma, colors, 5
            )
            if comp > 2 and len(mu_sigma) > 2:
                third_peak, third_sigma = mu_sigma[2]
                plot_vertical_line(
                    ax, "third_peak_minus_sigma", third_peak - third_sigma, colors, 6
                )
                plot_vertical_line(ax, "third_peak", third_peak, colors, 7)
                plot_vertical_line(
                    ax, "third_peak_plus_sigma", third_peak + third_sigma, colors, 8
                )
    ax.legend(loc="upper right")
    ax.hist(
        matrix,
        bins=np.arange(0, max_val, 1),
        density=True,
        histtype="stepfilled",
        alpha=0.4,
    )
    space = np.linspace(0, max_val, 1000)
    logprob = model.score_samples(space.reshape(-1, 1))
    pdf = np.exp(logprob)
    responsibilities = model.predict_proba(space.reshape(-1, 1))
    ax.plot(space, pdf, "-k")
    ax.plot(space, responsibilities * pdf[:, np.newaxis], "--k")
    fig.savefig(os.path.join(path, f"{name}.png"), dpi=300, format="png")
    fig.savefig(os.path.join(path, f"{name}.svg"), dpi=300, format="svg")
    plt.close(fig)
    logger.info("Plot saved as %s.png and %s.svg", name, name)


@log_exceptions
def estimate_divergence(
    data_folder: str, blocklist: Set[str], num_cores: int, log_queue
) -> None:
    """
    Estimate divergence of paralogs by processing all alignment files,
    plotting individual distributions and overall distributions,
    and saving the resulting matrices and plots.
    """
    logger.info("Estimating divergence of paralogs...")
    matplotlib.use("Agg")
    fig, ax = plt.subplots(figsize=(15, 10))
    files = sorted(glob(os.path.join(data_folder, "40aln_orth_par", "*.fasta.mafft")))
    if not files:
        logger.error(
            "No alignment files found in %s",
            os.path.join(data_folder, "40aln_orth_par"),
        )
        raise FileNotFoundError(
            f"No alignment files found in {os.path.join(data_folder, '40aln_orth_par')}"
        )
    args = [(f, blocklist) for f in files]
    with multiprocessing.Pool(
        processes=num_cores, initializer=worker_initializer, initargs=(log_queue,)
    ) as pool:
        async_result = pool.starmap_async(get_distance_matrix, args)
        try:
            results = async_result.get(timeout=600)
        except multiprocessing.TimeoutError as e:
            logger.error(
                f"Timeout during divergence estimation on files: %s\n{e}", files
            )
            pool.terminate()
            pool.join()
            raise
        except Exception as e:
            logger.error(
                "Error during divergence estimation on files: %s; Error: %s", files, e
            )
            pool.terminate()
            pool.join()
            raise
    distances_df = pd.concat([r[1] for r in results]).reset_index(drop=True)
    distances_df.to_csv(
        os.path.join(data_folder, "40aln_orth_par", "pairwise_distances.tsv"),
        sep="\t",
        index=False,
    )
    for r in results:
        if r[2] is None:
            continue
        space, e, sorted_arr = r[2]
        ax.plot(space, e, "k-", alpha=0.05)
        ax.plot(
            sorted_arr, np.zeros(sorted_arr.shape[0]), marker=2, color="k", alpha=0.5
        )
    end = np.round(ax.get_xlim()[1], 0)
    ax.xaxis.set_ticks(np.arange(0, end, 1))
    ax.set_xlabel("Divergence (%)")
    ax.set_ylabel("Frequency")
    fig.savefig(
        os.path.join(data_folder, "40aln_orth_par", "individual_distributions.png"),
        dpi=300,
        format="png",
    )
    fig.savefig(
        os.path.join(data_folder, "40aln_orth_par", "individual_distributions.svg"),
        dpi=300,
        format="svg",
    )
    plt.close(fig)
    divergence_array = np.array(
        [
            [x]
            for sublist in [r[0] if isinstance(r[0], list) else [r[0]] for r in results]
            for x in sublist
        ]
    )
    logger.info("Fitting divergence distribution models with 1, 2, and 3 components")
    for comp in [1, 2, 3]:
        get_plot(
            os.path.join(data_folder, "40aln_orth_par"),
            f"pairwise_distances_distribution_{comp}_comp",
            divergence_array,
            comp,
        )
    logger.info("Divergence estimation complete.")


# -----------------------------------------------------------------------------
# End of Module
# -----------------------------------------------------------------------------
