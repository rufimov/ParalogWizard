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
from Bio import SeqIO, SeqRecord
from Bio.Align.Applications import MafftCommandline
from Bio.Alphabet import generic_dna
from Bio.Phylo.Applications import FastTreeCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from matplotlib import axes
from matplotlib import pyplot
from scipy.signal import argrelextrema
from sklearn.mixture import BayesianGaussianMixture
from sklearn.neighbors import KernelDensity

from ParalogWizard.cast_retrieve import exon, contig


def sort_hit_table_ident(hittable_as_list: List[str], primary_field: str):
    """"""
    hittable_as_list.sort(
        key=lambda x: float(x.split()[1].split("_")[-1]), reverse=True
    )
    hittable_as_list.sort(key=lambda x: float(x.split()[5]), reverse=True)
    hittable_as_list.sort(key=lambda x: float(x.split()[4]))
    hittable_as_list.sort(key=lambda x: float(x.split()[3]), reverse=True)
    hittable_as_list.sort(key=lambda x: float(x.split()[2]), reverse=True)
    hittable_as_list.sort(key=globals()[primary_field])


def mafft_align(file):
    stdout, stderr = MafftCommandline(
        input=file,
        adjustdirectionaccurately=True,
        auto=True,
    )()
    with open(f"{os.path.splitext(file)[0]}.mafft.fasta", "w") as aligned:
        aligned.write(stdout.replace(">_R_", ">"))
    return f"{os.path.splitext(file)[0]}.mafft.fasta"


def fast_tree(file):
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


def build_alignments(path_to_data, n_cpu, logger):
    """"""

    logger.info("Building individual exon alignments...")
    with open(
        os.path.join(path_to_data, "31exonic_contigs", "all_hits.txt")
    ) as all_hits:
        all_hits_for_reference: List[str] = [x[:-1] for x in all_hits.readlines()]
    sort_hit_table_ident(all_hits_for_reference, "exon")
    exons: Dict[str, List[SeqRecord]] = dict()
    for hit in all_hits_for_reference:
        if exon(hit) not in set(exons.keys()):
            exons[exon(hit)]: List[SeqRecord] = [
                SeqRecord(
                    Seq(hit.split()[-1], generic_dna),
                    id=contig(hit) + "_" + hit.split()[-2],
                    description="",
                )
            ]
        else:
            new_list: List[SeqRecord] = exons[exon(hit)]
            new_list.append(
                SeqRecord(
                    Seq(hit.split()[-1], generic_dna),
                    id=contig(hit) + "_" + hit.split()[-2],
                    description="",
                )
            )
            exons[exon(hit)] = new_list
    os.makedirs(os.path.join(path_to_data, "40aln_orth_par"), exist_ok=True)
    for key in exons.keys():
        for record in exons[key]:
            new_id = f"{key}_N_{record.id.split('_N_')[1]}"
            record.id = new_id
    for key in exons.keys():
        SeqIO.write(
            exons[key],
            os.path.join(path_to_data, "40aln_orth_par", f"{key}.fasta"),
            "fasta-2line",
        )
    pool_aln = multiprocessing.Pool(processes=n_cpu)
    for file in glob(os.path.join(path_to_data, "40aln_orth_par", "*.fasta")):
        pool_aln.apply_async(mafft_align, (file,))
    pool_aln.close()
    pool_aln.join()
    pool_tree = multiprocessing.Pool(processes=n_cpu)
    for file in glob(os.path.join(path_to_data, "40aln_orth_par", "*.mafft.fasta")):
        pool_tree.apply_async(fast_tree, (file,))
    pool_tree.close()
    pool_tree.join()

    logger.info("Done\n")


def percent_dissimilarity(seq1: str, seq2: str) -> float:
    """"""

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
    blocklist: Set[str],
    axis: matplotlib.axes,
):
    """"""

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
                if "_".join(pair[0].split("_")[-2:]) not in blocklist:
                    current_distance_matrix.append(distance)
                current_list_to_write.append(f"{pair[0]}\t{str(distance)}\t{pair[1]}")

    if len(current_distance_matrix) > 0:
        # Search for local minima
        current_distance_array = numpy.array(current_distance_matrix).reshape(-1, 1)
        sorted_current_distance_array = numpy.array(sorted(current_distance_matrix))
        kde = KernelDensity(kernel="gaussian", bandwidth=1.5)
        kde.fit(current_distance_array)
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
    """"""
    from sklearn.mixture import BayesianGaussianMixture

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
    """"""

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


def estimate_divergence(path_to_data, blocklist, logger):
    """"""

    logger.info("Estimating divergence of paralogs...")
    divergency_distribution: List[float] = []
    divergencies_to_write: List[str] = []
    matplotlib.use("Agg")
    fig, axis = pyplot.subplots(figsize=(15, 10))
    for file in sorted(
        glob(os.path.join(path_to_data, "40aln_orth_par", "*.mafft.fasta"))
    ):
        get_distance_matrix(
            file,
            divergency_distribution,
            divergencies_to_write,
            blocklist,
            axis,
        )
    end = axis.get_xlim()[1]
    end = numpy.round(end, 0)
    axis.xaxis.set_ticks(numpy.arange(0, end, 1))
    axis.set_xlabel("Divergence (%)")
    axis.set_ylabel("Frequency")
    fig.savefig(
        os.path.join(path_to_data, "40aln_orth_par", "individual_distributions.png"),
        dpi=300,
        format="png",
    )
    fig.savefig(
        os.path.join(path_to_data, "40aln_orth_par", "individual_distributions.svg"),
        dpi=300,
        format="svg",
    )
    pyplot.close(fig)
    with open(
        os.path.join(path_to_data, "40aln_orth_par", "pairwise_distances.txt"), "w"
    ) as divergency_distribution_to_write:
        for i in divergencies_to_write:
            divergency_distribution_to_write.write(str(i) + "\n")
    divergency_distribution_array: numpy.ndarray = numpy.array(
        [[x] for x in divergency_distribution]
    )
    get_plot(
        os.path.join(path_to_data, "40aln_orth_par"),
        "pairwise_distances_distribution_1_comp",
        divergency_distribution_array,
        1,
    )
    get_plot(
        os.path.join(path_to_data, "40aln_orth_par"),
        "pairwise_distances_distribution_2_comp",
        divergency_distribution_array,
        2,
    )
    get_plot(
        os.path.join(path_to_data, "40aln_orth_par"),
        "pairwise_distances_distribution_3_comp",
        divergency_distribution_array,
        3,
    )

    logger.info("Done\n")
