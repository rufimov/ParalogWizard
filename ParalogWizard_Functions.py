import itertools
import random
from typing import List, Dict, Tuple, Set

import Bio
import matplotlib
import numpy
from Bio import SeqRecord, SeqIO
from Bio.Alphabet import generic_dna
from matplotlib import axes, pyplot
from sklearn.mixture import BayesianGaussianMixture
from sklearn.neighbors import KernelDensity
from scipy.signal import argrelextrema


def sort_hit_table_cover(hittable_as_list: List[str], primary_field: str):
    """

    :param hittable_as_list:
    :type hittable_as_list:
    :param primary_field:
    :type primary_field:
    """
    hittable_as_list.sort(
        key=lambda x: float(x.split()[1].split("_")[-1]), reverse=True
    )
    hittable_as_list.sort(key=lambda x: float(x.split()[5]), reverse=True)
    hittable_as_list.sort(key=lambda x: float(x.split()[4]))
    hittable_as_list.sort(key=lambda x: float(x.split()[2]), reverse=True)
    hittable_as_list.sort(key=lambda x: float(x.split()[3]), reverse=True)
    hittable_as_list.sort(key=globals()[primary_field])


def sort_hit_table_ident(hittable_as_list: List[str], primary_field: str):
    """

    :param hittable_as_list:
    :type hittable_as_list:
    :param primary_field:
    :type primary_field:
    """
    hittable_as_list.sort(
        key=lambda x: float(x.split()[1].split("_")[-1]), reverse=True
    )
    hittable_as_list.sort(key=lambda x: float(x.split()[5]), reverse=True)
    hittable_as_list.sort(key=lambda x: float(x.split()[4]))
    hittable_as_list.sort(key=lambda x: float(x.split()[3]), reverse=True)
    hittable_as_list.sort(key=lambda x: float(x.split()[2]), reverse=True)
    hittable_as_list.sort(key=globals()[primary_field])


def exon(string: str) -> str:
    """

    :param string:
    :type string:
    :return:
    :rtype:
    """
    return string.split()[0].split("-")[1]


def locus(string: str) -> str:
    """

    :rtype: object
    """
    return string.split()[0].split("-")[1].split("_")[0]


def contig(string: str) -> str:
    """

    :param string:
    :type string:
    :return:
    :rtype:
    """
    return string.split()[1]


def contig_locus(string: str) -> str:
    """

    :param string:
    :type string:
    :return:
    :rtype:
    """
    return string.split()[1].split("_N_")[0]


def slicing(
    dictionary: Dict[str, Bio.SeqRecord.SeqRecord],
    current_string: str,
    rev: bool,
) -> str:
    """

    :param dictionary:
    :type dictionary:
    :param current_string:
    :type current_string:
    :param rev:
    :type rev:
    :return:
    :rtype:
    """
    sequence = dictionary[current_string.split()[1]].seq
    if rev:
        start = int(current_string.split()[7])
        end = int(current_string.split()[6])
        result = str(sequence.reverse_complement())[-start : -end - 1 : -1][::-1]
    else:
        start = int(current_string.split()[6])
        end = int(current_string.split()[7])
        result = str(sequence)[start - 1 : end]
    return result


def percent_dissimilarity(seq1: str, seq2: str) -> float:
    """

    :param seq1:
    :type seq1:
    :param seq2:
    :type seq2:
    :return:
    :rtype:
    """
    seq1 = seq1.lower()
    seq2 = seq2.lower()
    # Removing parts of alignments with gaps at the end
    while seq1[0] == "-" or seq2[0] == "-":
        seq1 = seq1[1:]
        seq2 = seq2[1:]
    while seq1[-1] == "-" or seq2[-1] == "-":
        seq1 = seq1[:-1]
        seq2 = seq2[:-1]
    # Removing positions with gaps in both sequences
    seq1_wo_mutual_gaps = str()
    seq2_wo_mutual_gaps = str()
    for nucl in zip(seq1, seq2):
        if nucl[0] != "-" and nucl[1] != "-":
            seq1_wo_mutual_gaps = seq1_wo_mutual_gaps + nucl[0]
            seq2_wo_mutual_gaps = seq2_wo_mutual_gaps + nucl[1]
    # Counting mismatches, ignoring positions with gaps in one of the sequences
    count = 0
    for nucl in zip(seq1_wo_mutual_gaps, seq2_wo_mutual_gaps):
        if nucl[0] != nucl[1] and nucl[0] != "-" and nucl[1] != "-":
            count += 1
    dissimilarity = (count / len(seq1)) * 100
    return dissimilarity


def get_distance_matrix(
    file_to_process: str,
    sum_list: List[float],
    list_to_write: List[str],
    blacklist: Set[str],
    axis: matplotlib.axes,
):
    """

    :param file_to_process:
    :type file_to_process:
    :param sum_list:
    :type sum_list:
    :param list_to_write:
    :type list_to_write:
    :param blacklist:
    :type blacklist:
    :param axis:
    :type axis:
    """
    current_distance_matrix: List[float] = list()
    current_list_to_write = []
    with open(file_to_process) as fasta_file:
        sequences: Dict[str, Bio.SeqRecord.SeqRecord] = SeqIO.to_dict(
            SeqIO.parse(fasta_file, "fasta", generic_dna)
        )
        for pair in list(itertools.combinations(sequences.keys(), 2)):
            if "_".join(pair[0].split("_")[-2:]) == "_".join(pair[1].split("_")[-2:]):
                sequence1: str = str(sequences[pair[0]].seq)
                sequence2: str = str(sequences[pair[1]].seq)
                distance: float = percent_dissimilarity(sequence1, sequence2)
                if "_".join(pair[0].split("_")[-2:]) not in blacklist:
                    current_distance_matrix.append(distance)
                current_list_to_write.append(f"{pair[0]}\t{str(distance)}\t{pair[1]}")

    if len(current_distance_matrix) > 0:
        # Search for local minima
        current_distance_array = numpy.array(current_distance_matrix).reshape(-1, 1)
        sorted_current_distance_array = numpy.array(sorted(current_distance_matrix))
        kde = KernelDensity(kernel="gaussian", bandwidth=1.5).fit(
            current_distance_array
        )
        linear_space = numpy.linspace(-1, numpy.max(current_distance_array) + 1, 1000)
        e = numpy.exp(kde.score_samples(linear_space.reshape(-1, 1)))
        mi = argrelextrema(e, numpy.less)[0]
        minimum = linear_space[mi]
        axis.plot(linear_space, e, "k-", alpha=0.05)
        axis.plot(
            sorted_current_distance_array,
            [0] * sorted_current_distance_array.shape[0],
            marker=2,
            color="k",
            alpha=0.5,
        )
        if len(minimum) == 0:
            sum_list.append(sum(current_distance_matrix) / len(current_distance_matrix))
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
                    cluster = [sorted_current_distance_array.tolist()[indices[0]]]
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
    list_to_write.extend(current_list_to_write)


def get_model(array: numpy.ndarray, num_comp: int) -> BayesianGaussianMixture:
    """

    :param array:
    :type array:
    :param num_comp:
    :type num_comp:
    :return:
    :rtype:
    """
    return BayesianGaussianMixture(
        n_components=num_comp,
        max_iter=10000,
        n_init=10,
    ).fit(array)


def plot_vertical_line(
    plot: matplotlib.axes,
    line_name: str,
    line_value: float,
    list_of_colors: List[str],
    num_for_color: int,
):
    """

    :param plot:
    :type plot:
    :param line_name:
    :type line_name:
    :param line_value:
    :type line_value:
    :param list_of_colors:
    :type list_of_colors:
    :param num_for_color:
    :type num_for_color:
    """
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
    """

    :param path:
    :type path:
    :param name:
    :type name:
    :param matrix:
    :type matrix:
    :param comp:
    :type comp:
    """
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
    weights = divergency_distribution_mix.weights_.flatten().tolist()
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
        k=comp * 2,
    )
    # Plotting vertical lines
    first_peak: float = mu_sigma[0][0]
    first_peak_minus_sigma: float = first_peak - mu_sigma[0][1]
    plot_vertical_line(
        axis, "first_peak_minus_sigma", first_peak_minus_sigma, colors, 0
    )
    plot_vertical_line(axis, "first_peak", first_peak, colors, 1)
    if comp > 1:
        second_peak: float = mu_sigma[1][0]
        second_peak_minus_sigma: float = second_peak - mu_sigma[1][1]
        plot_vertical_line(
            axis, "second_peak_minus_sigma", second_peak_minus_sigma, colors, 2
        )
        plot_vertical_line(axis, "second_peak", second_peak, colors, 3)
        if comp > 2:
            third_peak: float = mu_sigma[2][0]
            third_peak_minus_sigma: float = third_peak - mu_sigma[2][1]
            plot_vertical_line(
                axis, "third_peak_minus_sigma", third_peak_minus_sigma, colors, 4
            )
            plot_vertical_line(axis, "third_peak", third_peak, colors, 5)
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
    fig.savefig(f"{path}{name}.png", dpi=300, format="png")
    fig.savefig(f"{path}{name}.svg", dpi=300, format="svg")
    pyplot.close(fig)
