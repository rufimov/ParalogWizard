import glob
import itertools
import multiprocessing
import os
import random
from glob import glob
from typing import Dict, List
from typing import Set
from typing import Tuple

import Bio.Application
import Bio.Application
import matplotlib
import numpy
import pandas
from Bio import SeqIO, SeqRecord
from Bio.Align.Applications import MafftCommandline
from Bio.Phylo.Applications import FastTreeCommandline
from matplotlib import axes
from matplotlib import pyplot
from pandas.core import frame
from scipy.signal import argrelextrema
from sklearn.mixture import BayesianGaussianMixture
from sklearn.neighbors import KernelDensity


def mafft_align(file):
    """
    Aligns a fasta file using mafft and then adjusts the direction of the sequences in the file
    to match the direction of the first sequence in the file.
    :param file: a fasta file
    :return: a file with the same name as the input file with .fasta.mafft appended
    """
    stdout, stderr = MafftCommandline(
        input=file,
        auto=True,
    )()
    with open(f"{os.path.splitext(file)[0]}.fasta.mafft", "w") as aligned:
        aligned.write(stdout)


def fast_tree(file):
    """
    Builds a tree using FastTreeMP, and if that fails, falls back to FastTree.
    :param file: a fasta file
    :return: a file with the same name as the input file with .tre appended
    """
    try:
        FastTreeCommandline(
            "fasttreemp",
            nt=True,
            gtr=True,
            input=file,
            out=f"{file}.tre",
            fastest=True,
        )()
    except Bio.Application.ApplicationError:
        FastTreeCommandline(
            "fasttree",
            nt=True,
            gtr=True,
            input=file,
            out=f"{file}.tre",
            fastest=True,
        )()


def build_alignments(data_folder, n_cpu, logger):
    """
    Builds exon alignments using mafft and builds trees using FastTreeMP
    :param data_folder: the folder where the data is stored
    :param n_cpu: the number of CPUs to use
    :param logger: the logger
    :return: None
    """
    logger.info("Building individual exon alignments...")
    all_hits_for_reference: pandas.core.frame.DataFrame = pandas.read_csv(
        os.path.join(data_folder, "31exonic_contigs", "all_hits.tsv"), sep="\t"
    )
    grouped_exons = all_hits_for_reference.groupby("exon")
    os.makedirs(os.path.join(data_folder, "40aln_orth_par"), exist_ok=True)
    for exon in grouped_exons:
        exon_dataframe = exon[1].reset_index(drop=True)
        exon_name = exon[0]
        with open(
            os.path.join(data_folder, "40aln_orth_par", f"{exon_name}.fasta"), "w"
        ) as exon_aln_fasta:
            for index in range(len(exon_dataframe)):
                contig_name = exon_dataframe.loc[index, "saccver"]
                sample = exon_dataframe.loc[index, "sample"]
                seq_name = f"{exon_name}_N_{contig_name.split('_N_')[1]}_{sample}"
                seq = exon_dataframe["sequence"][index]
                exon_aln_fasta.write(f">{seq_name}\n{seq}\n")
    files_to_align = glob(os.path.join(data_folder, "40aln_orth_par", "*fasta"))
    with multiprocessing.Pool(processes=n_cpu) as pool_aln:
        pool_aln.map(mafft_align, files_to_align)
    files_to_align = glob(os.path.join(data_folder, "40aln_orth_par", "*fasta.mafft"))
    with multiprocessing.Pool(processes=n_cpu) as pool_tree:
        pool_tree.map(fast_tree, files_to_align)
    logger.info("Done\n")


def percent_dissimilarity(seq1: str, seq2: str) -> float or None:
    """
    :param seq1: sequence one
    :param seq2: sequence two
    :return: the percent dissimilarity between the two sequences
    """

    seq1 = seq1.lower()
    seq2 = seq2.lower()
    # Removing positions with gaps in both sequences
    seq1_wo_mutual_gaps = str()
    seq2_wo_mutual_gaps = str()
    for nucl in zip(seq1, seq2):
        if nucl[0] == "-" and nucl[1] == "-":
            continue
        seq1_wo_mutual_gaps = seq1_wo_mutual_gaps + nucl[0]
        seq2_wo_mutual_gaps = seq2_wo_mutual_gaps + nucl[1]
    # Length of alignment without mutual gaps
    len_aln = len(seq1_wo_mutual_gaps)
    # Removing parts of alignments with gaps at the end and count removed positions from both sides
    count_left = 0
    count_right = 0
    while seq1_wo_mutual_gaps[0] == "-" or seq2_wo_mutual_gaps[0] == "-":
        seq1_wo_mutual_gaps = seq1_wo_mutual_gaps[1:]
        seq2_wo_mutual_gaps = seq2_wo_mutual_gaps[1:]
        if seq1_wo_mutual_gaps == "" or seq2_wo_mutual_gaps == "":
            return None
        count_left += 1
    while seq1_wo_mutual_gaps[-1] == "-" or seq2_wo_mutual_gaps[-1] == "-":
        seq1_wo_mutual_gaps = seq1_wo_mutual_gaps[:-1]
        seq2_wo_mutual_gaps = seq2_wo_mutual_gaps[:-1]
        if seq1_wo_mutual_gaps == "" or seq2_wo_mutual_gaps == "":
            return None
        count_right += 1
    # Checking if resulting alignment is not less than 0.75 than the length of either or 2 initial sequence without end
    # and mutual gaps
    overlap = len_aln - count_left - count_right
    if (
        (overlap / (len_aln - count_left) < 0.5)
        or (overlap / (len_aln - count_right) < 0.5)
        or overlap < 100
    ):
        return None
    # Counting mismatches, ignoring positions with gaps in one of the sequences
    count = 0
    for nucl in zip(seq1_wo_mutual_gaps, seq2_wo_mutual_gaps):
        if (
            nucl[0] != nucl[1] and nucl[0] != "-" and nucl[1] != "-"
        ):  # gaps are not counted
            # if nucl[0] != nucl[1]:  # gaps are counted
            count += 1
    dissimilarity = (count / len(seq1_wo_mutual_gaps)) * 100
    return dissimilarity


def get_distance_matrix(
    file_to_process: str,
    blocklist: Set[str],
):
    """ """
    current_matrix_to_plot: List[float] = list()
    current_matrix_to_write = []
    sum_list = []
    with open(file_to_process) as fasta_file:
        sequences: Dict[str, Bio.SeqRecord.SeqRecord] = SeqIO.to_dict(
            SeqIO.parse(fasta_file, "fasta")
        )
    seq_names = pandas.DataFrame(list(sequences.keys()))
    seq_names[1] = seq_names[0].str.split("_").str[-2:].str.join("_")
    duplicated_samp = set(seq_names[seq_names.duplicated(subset=1)][1].values)
    non_duplicated_samp = set(seq_names[~seq_names[1].isin(duplicated_samp)][1].values)
    for samp in duplicated_samp:
        seqs_to_pairs = seq_names[seq_names[1] == samp][0].values.tolist()
        for pair in list(itertools.combinations(seqs_to_pairs, 2)):
            sequence1: str = str(sequences[pair[0]].seq)
            sequence2: str = str(sequences[pair[1]].seq)
            distance: float = percent_dissimilarity(sequence1, sequence2)
            if distance is None:
                continue
            if samp not in blocklist:
                current_matrix_to_plot.append(distance)
            current_matrix_to_write.append([pair[0], distance, pair[1]])
    for samp in non_duplicated_samp:
        seqs_non_paired = seq_names[seq_names[1] == samp][0].values.tolist()
        for seq in seqs_non_paired:
            sequence = str(sequences[seq].seq)
            if len(sequence) < 100:
                continue
            current_matrix_to_write.append([seq, numpy.nan, numpy.nan])
    if len(current_matrix_to_plot) > 0:
        # Search for local minima
        current_distance_array = numpy.array(current_matrix_to_plot).reshape(-1, 1)
        sorted_current_distance_array = numpy.array(sorted(current_matrix_to_plot))
        kde = KernelDensity(kernel="gaussian", bandwidth=1.5)
        kde.fit(current_distance_array)
        linear_space = numpy.linspace(-1, numpy.max(current_distance_array) + 1, 1000)
        e = numpy.exp(kde.score_samples(linear_space.reshape(-1, 1)))
        mi = argrelextrema(e, numpy.less)[0]
        minimum = linear_space[mi]
        to_plot = [
            linear_space,
            e,
            sorted_current_distance_array,
        ]
        if len(minimum) == 0:
            sum_list.append(sum(current_matrix_to_plot) / len(current_matrix_to_plot))
        else:
            # Calculate mean in every cluster
            for i in range(0, len(minimum)):
                if i == 0:
                    indices = numpy.where(
                        numpy.logical_and(
                            sorted_current_distance_array < minimum[i],
                            sorted_current_distance_array >= 0,
                        )
                    )[0].tolist()
                else:
                    indices = numpy.where(
                        numpy.logical_and(
                            sorted_current_distance_array < minimum[i],
                            sorted_current_distance_array > minimum[i - 1],
                        )
                    )[0].tolist()
                if len(indices) == 1:
                    cluster = [list(sorted_current_distance_array)[indices[0]]]
                else:
                    cluster = sorted_current_distance_array.tolist()[
                        min(indices) : max(indices) + 1
                    ]
                cluster_mean = sum(cluster) / len(cluster)
                sum_list.append(cluster_mean)
            last_cluster_indices = numpy.where(
                sorted_current_distance_array > minimum[-1]
            )[0].tolist()
            last_cluster = sorted_current_distance_array.tolist()[
                min(last_cluster_indices) : max(last_cluster_indices) + 1
            ]
            last_cluster_mean = sum(last_cluster) / len(last_cluster)
            sum_list.append(last_cluster_mean)
    else:
        to_plot = None
    current_matrix_to_write = pandas.DataFrame(
        current_matrix_to_write, columns=["seq1", "dist", "seq2"]
    )
    return sum_list, current_matrix_to_write, to_plot


def get_model(array: numpy.ndarray, num_comp: int) -> BayesianGaussianMixture:
    """
    Get BayesianGaussianMixture object for given array.
    :param array:  to fit
    :param num_comp: number of components
    :return: BayesianGaussianMixture model with fitted data
    """
    from sklearn.mixture import BayesianGaussianMixture

    mix = BayesianGaussianMixture(
        n_components=num_comp,
        max_iter=10000,
        n_init=10,
    )
    mix.fit(array)
    return mix


def plot_vertical_line(
    plot: matplotlib.axes,
    line_name: str,
    line_value: float,
    list_of_colors: List[str],
    num_for_color: int,
):
    """
    Plot vertical line on given plot.
    :param plot: plot to draw on
    :param line_name: name of line
    :param line_value: value of line
    :param list_of_colors: list of colors to use
    :param num_for_color: number of color to use
    :return:
    """

    import numpy

    plot.axvline(
        x=line_value,
        label=f"{line_name} - {numpy.round(line_value, 2)}",
        c=list_of_colors[num_for_color],
        lw=0.6,
        ls="--",
    )


def get_plot(
    path: str,
    name: str,
    matrix: numpy.ndarray,
    comp: int,
):
    """"""

    # General plot settings
    matplotlib.use("Agg")
    max_of_matrix: float = numpy.round(numpy.max(matrix)) + 1
    fig, axis = pyplot.subplots(figsize=(15, 10))
    axis.set_xlabel("Divergence (%)")
    axis.set_ylabel("Frequency")
    axis.xaxis.set_ticks(numpy.arange(-1, max_of_matrix, 1))
    fig.suptitle(f"{name}")
    # Data to plot
    divergency_distribution_mix: BayesianGaussianMixture = get_model(matrix, comp)
    means: List[float] = divergency_distribution_mix.means_.flatten().tolist()
    sigmas: List[float] = [
        numpy.sqrt(x)
        for x in divergency_distribution_mix.covariances_.flatten().tolist()
    ]
    mu_sigma: List[Tuple[float, float]] = sorted(
        list(zip(means, sigmas)), key=lambda x: x[0]
    )
    # weights = divergency_distribution_mix.weights_.flatten().tolist()
    # Set colours for vertical lines to be plotted
    colors: List[str] = random.sample(
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
    # Plotting vertical lines
    first_peak: float = mu_sigma[0][0]
    first_peak_minus_sigma: float = first_peak - mu_sigma[0][1]
    first_peak_plus_sigma: float = first_peak + mu_sigma[0][1]
    plot_vertical_line(
        axis, "first_peak_minus_sigma", first_peak_minus_sigma, colors, 0
    )
    plot_vertical_line(axis, "first_peak", first_peak, colors, 1)
    plot_vertical_line(axis, "first_peak_plus_sigma", first_peak_plus_sigma, colors, 2)
    if comp > 1:
        second_peak: float = mu_sigma[1][0]
        second_peak_minus_sigma: float = second_peak - mu_sigma[1][1]
        second_peak_plus_sigma: float = second_peak + mu_sigma[1][1]
        plot_vertical_line(
            axis, "second_peak_minus_sigma", second_peak_minus_sigma, colors, 3
        )
        plot_vertical_line(axis, "second_peak", second_peak, colors, 4)
        plot_vertical_line(
            axis, "second_peak_plus_sigma", second_peak_plus_sigma, colors, 5
        )
        if comp > 2:
            third_peak: float = mu_sigma[2][0]
            third_peak_minus_sigma: float = third_peak - mu_sigma[2][1]
            third_peak_plus_sigma: float = third_peak + mu_sigma[2][1]
            plot_vertical_line(
                axis, "third_peak_minus_sigma", third_peak_minus_sigma, colors, 6
            )
            plot_vertical_line(axis, "third_peak", third_peak, colors, 7)
            plot_vertical_line(
                axis, "third_peak_plus_sigma", third_peak_plus_sigma, colors, 8
            )
    # Plot graph and other
    axis.legend(loc="upper right")
    axis.plot(matrix, [0] * matrix.shape[0], marker=2, color="k")
    space: numpy.ndarray = numpy.linspace(0, max_of_matrix, 1000)
    logprob: numpy.ndarray = divergency_distribution_mix.score_samples(
        space.reshape(-1, 1)
    )
    pdf: numpy.ndarray = numpy.exp(logprob)
    responsibilities: numpy.ndarray = divergency_distribution_mix.predict_proba(
        space.reshape(-1, 1)
    )
    axis.hist(
        matrix,
        bins=numpy.arange(0, max_of_matrix, 1),
        density=True,
        histtype="stepfilled",
        alpha=0.4,
    )
    axis.plot(space, pdf, "-k")
    pdf_individual: numpy.ndarray = responsibilities * pdf[:, numpy.newaxis]
    axis.plot(space, pdf_individual, "--k")
    # Save figure
    fig.savefig(os.path.join(path, f"{name}.png"), dpi=300, format="png")
    fig.savefig(os.path.join(path, f"{name}.svg"), dpi=300, format="svg")
    pyplot.close(fig)


def estimate_divergence(data_folder, blocklist, num_cores, logger):
    """"""

    logger.info("Estimating divergence of paralogs...")
    matplotlib.use("Agg")
    fig, axis = pyplot.subplots(figsize=(15, 10))
    files = sorted(glob(os.path.join(data_folder, "40aln_orth_par", "*.fasta.mafft")))
    args_est_div = list(
        zip(
            files,
            [blocklist] * len(files),
        )
    )
    with multiprocessing.Pool(processes=num_cores) as pool_est_div:
        results = pool_est_div.starmap(get_distance_matrix, args_est_div)
    divergencies_to_file = pandas.concat([result[1] for result in results]).reset_index(
        drop=True
    )
    divergency_distribution = [
        item for sublist in [result[0] for result in results] for item in sublist
    ]
    to_plot = [result[2] for result in results]
    for el in to_plot:
        if el is None:
            continue
        linear_space = el[0]
        e = el[1]
        sorted_current_distance_array = el[2]
        axis.plot(linear_space, e, "k-", alpha=0.05)
        axis.plot(
            sorted_current_distance_array,
            [0] * sorted_current_distance_array.shape[0],
            marker=2,
            color="k",
            alpha=0.5,
        )
    divergencies_to_file.to_csv(
        os.path.join(data_folder, "40aln_orth_par", "pairwise_distances.tsv"),
        sep="\t",
        index=False,
    )
    end = axis.get_xlim()[1]
    end = numpy.round(end, 0)
    axis.xaxis.set_ticks(numpy.arange(0, end, 1))
    axis.set_xlabel("Divergence (%)")
    axis.set_ylabel("Frequency")
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
    pyplot.close(fig)
    divergency_distribution_array: numpy.ndarray = numpy.array(
        [[x] for x in divergency_distribution]
    )
    get_plot(
        os.path.join(data_folder, "40aln_orth_par"),
        "pairwise_distances_distribution_1_comp",
        divergency_distribution_array,
        1,
    )
    get_plot(
        os.path.join(data_folder, "40aln_orth_par"),
        "pairwise_distances_distribution_2_comp",
        divergency_distribution_array,
        2,
    )
    get_plot(
        os.path.join(data_folder, "40aln_orth_par"),
        "pairwise_distances_distribution_3_comp",
        divergency_distribution_array,
        3,
    )

    logger.info("Done\n")
