import argparse
import copy
import glob
import itertools
import os
import random
import re
import shutil
import sys
from typing import Set, Dict, List, Tuple, Union


import Bio
import Bio.Application
import matplotlib
import numpy
from Bio import SeqIO, SeqRecord
from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastnCommandline
from Bio.Phylo.Applications import FastTreeCommandline
from Bio.Align.Applications import MafftCommandline
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from datetime import datetime
from matplotlib import axes, pyplot
from scipy.signal import argrelextrema
from sklearn.mixture import BayesianGaussianMixture
from sklearn.neighbors import KernelDensity
import logging


class ParsedArgs:
    def __init__(self):
        parser = argparse.ArgumentParser(
            usage="""ParalogWizard <command> [<args>]
The ParalogWizard commands are:
    cast_collect
    cast_retrieve
    cast_analyze
    cast_create
    cast_correct         
Use ParalogWizard <command> -h for help with arguments of the command of interest
"""
        )
        parser.add_argument(
            "command",
            help="Subcommand to run",
        )
        self.args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, self.args.command):
            print("Unrecognized command")
            parser.print_help()
            exit(1)

    def common_args(self, parser):
        parser.add_argument("-d", "--data_folder", required=True)
        parser.add_argument(
            "-nc",
            "--num_cores",
            default=1,
            type=int,
            help="Number of cores used. Default: 1",
        )

    def cast_collect(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("-c", "--contig_folder", required=True)
        self.common_args(parser)
        args = parser.parse_args(sys.argv[2:])
        return args

    def cast_retrieve(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("-pe", "--probes_exons", required=True)
        parser.add_argument("-l", "--length_cover", required=True, type=float)
        parser.add_argument("-s", "--spades_cover", required=True, type=float)
        self.common_args(parser)
        args = parser.parse_args(sys.argv[2:])
        return args

    def cast_analyze(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("-b", "--blocklist", nargs="+", required=False)
        self.common_args(parser)
        args = parser.parse_args(sys.argv[2:])
        if args.blocklist is None:
            args.blocklist = set()
        else:
            args.blocklist = set(args.blocklist)
        return args

    def cast_create(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("-b", "--blocklist", nargs="+", required=False)
        parser.add_argument("-pe", "--probes_exons", required=True)
        parser.add_argument(
            "-p",
            "--paralogs",
            default=False,
            action="store_true",
            help="Paralog detection. Default: off",
        )
        parser.add_argument(
            "-mi",
            "--minimum_divergence",
            required="--paralogs" in sys.argv,
            type=float,
        )
        parser.add_argument(
            "-ma",
            "--maximum_divergence",
            required="--paralogs" in sys.argv,
            type=float,
        )
        self.common_args(parser)
        args = parser.parse_args(sys.argv[2:])
        if args.blocklist is None:
            args.blocklist = set()
        else:
            args.blocklist = set(args.blocklist)
        if args.paralogs is True and (
            args.minimum_divergence is None or args.maximum_divergence is None
        ):
            print("Minimum and maximum divergence is not specified")
            parser.print_help()
            exit(1)
        return args

    def cast_correct(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("-r", "--redlist", nargs="+", required=False)
        parser.add_argument("-i", "--min_identity", required=True, type=float)
        parser.add_argument("-pp", "--probes_paralogs", required=True)
        self.common_args(parser)
        args = parser.parse_args(sys.argv[2:])
        if args.redlist is None:
            args.redlist = set()
        else:
            args.redlist = set(args.redlist)
        return args

    def get_args_dict(self):
        command = self.args.__dict__
        arguments = getattr(self, self.args.command)().__dict__
        argument_dictionary = command.copy()
        argument_dictionary.update(arguments)
        return argument_dictionary


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
    blocklist: Set[str],
    axis: matplotlib.axes,
):
    """

    :param file_to_process:
    :type file_to_process:
    :param sum_list:
    :type sum_list:
    :param list_to_write:
    :type list_to_write:
    :param blocklist:
    :type blocklist:
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
    fig.savefig(f"{path}{name}.png", dpi=300, format="png")
    fig.savefig(f"{path}{name}.svg", dpi=300, format="svg")
    pyplot.close(fig)


def collect_contigs(path_to_contigs, path_to_data, logger):
    """

    :param path_to_contigs:
    :type path_to_contigs:
    :param path_to_data:
    :type path_to_data:
    """
    os.makedirs(f"{path_to_data}/HybPiper_contigs", exist_ok=True)
    for folder in os.listdir(path_to_contigs):
        if os.path.isdir(f"{path_to_contigs}/{folder}"):
            logger.info(f"Processing {folder}")
            with open(
                f"{path_to_data}/HybPiper_contigs/{folder}_contigs.fasta", "w"
            ) as contigs:
                for locus in os.listdir(f"{path_to_contigs}/{folder}"):
                    if os.path.isdir(
                        f"{path_to_contigs}/{folder}/{locus}"
                    ) and os.path.exists(
                        f"{path_to_contigs}/{folder}/{locus}/{locus}_contigs.fasta"
                    ):
                        logger.info(f"\tProcessing {locus}")
                        with open(
                            f"{path_to_contigs}/{folder}/{locus}/{locus}_contigs.fasta"
                        ) as locus_contigs:
                            file_content = locus_contigs.read().replace(
                                ">", f">{locus}_"
                            )
                            contigs.write(file_content)
                        logger.info("\tOK")
            if (
                os.stat(
                    f"{path_to_data}/HybPiper_contigs/{folder}_contigs.fasta"
                ).st_size
                == 0
            ):
                os.remove(f"{path_to_data}/HybPiper_contigs/{folder}_contigs.fasta")
            logger.info("Ok")


def prepare_contigs(path_to_data, logger):
    """

    :param path_to_data:
    :type path_to_data:
    """
    main_path: str = f"{path_to_data}/exons/40contigs/"
    os.makedirs(main_path, exist_ok=True)
    logger.info("Preparing congits...")
    for file in glob.glob(f"{path_to_data}/HybPiper_contigs/*contigs.fasta"):
        shutil.copy(file, os.path.join(main_path, os.path.basename(file)))
    for file in glob.glob(f"{main_path}*.fasta"):
        with open(file, "r") as fasta:
            lines: List[str] = fasta.readlines()
        with open(file, "w") as fasta:
            for line in lines:
                line = line.replace("NODE", "N")
                line = re.sub(
                    r"length_([0-9]+)_cov_([0-9]+\.[0-9][0-9]).*", r"\1_c_\2", line
                )
                fasta.write(line)
        name_of_file: str = "_".join(
            os.path.basename(file).split(".")[0].split("_")[0:2]
        )
        os.rename(file, os.path.join(os.path.dirname(file), f"{name_of_file}.fasta"))
    logger.info("Done\n")


def create_hit_tables(path_to_data, probe_exons, length_cover, n_cpu, logger):
    """

    :param path_to_data:
    :type path_to_data:
    :param probe_exons:
    :type probe_exons:
    :param length_cover:
    :type length_cover:
    :param n_cpu:
    :type n_cpu:
    """
    main_path: str = f"{path_to_data}/exons/40contigs/"
    os.makedirs(main_path, exist_ok=True)
    logger.info("Creating hit tables...")
    for file in glob.glob(f"{main_path}*.fasta"):
        file: str = os.path.basename(file)
        sample: str = file[:-6]
        logger.info(f"\tProcessing {sample}")
        NcbimakeblastdbCommandline(
            dbtype="nucl",
            input_file=f"{main_path}{file}",
            out=f"{main_path}{sample}",
            parse_seqids=True,
        )()
        logger.info("\tRunning BLAST...")
        NcbiblastnCommandline(
            task="blastn",
            query=probe_exons,
            db=f"{main_path}{sample}",
            out=f"{main_path}reference_in_{sample}_contigs.txt",
            qcov_hsp_perc=length_cover,
            num_threads=n_cpu,
            outfmt="6 qaccver saccver pident qcovhsp evalue bitscore sstart send",
        )()
        logger.info("\tOK")
    logger.info("Done\n")


def correct_contgis(
    path_to_data, statistics, spades_cover, all_hits_for_reference, logger
):
    """

    :param path_to_data:
    :type path_to_data:
    :param statistics:
    :type statistics:
    :param spades_cover:
    :type spades_cover:
    :param all_hits_for_reference:
    :type all_hits_for_reference:
    """
    main_path: str = f"{path_to_data}/exons/40contigs/"
    os.makedirs(main_path, exist_ok=True)
    logger.info("Correcting contigs...")
    for file in glob.glob(f"{main_path}*.fasta"):
        sample: str = os.path.basename(os.path.splitext(file)[0])
        logger.info(f" Processing {sample}")
        statistics[sample]: Dict[str, Dict[str, List[str]]] = dict()
        hits: List[str] = list()
        with open(f"{main_path}reference_in_{sample}_contigs.txt") as blast_results:
            blast_results_as_list = [x[:-1] for x in blast_results.readlines()]
            for line in blast_results_as_list:
                if (
                    contig_locus(line) == locus(line)
                    and float(line.split()[1].split("_c_")[1]) >= spades_cover
                ):
                    hits.append(line)
        sort_hit_table_cover(hits, "exon")
        hits_exons_contigs: Set[str] = set()
        hits_dedup: List[str] = list()
        for hit in hits:
            if str(exon(hit) + contig(hit)) not in hits_exons_contigs:
                hits_dedup.append(hit)
            else:
                pass
            hits_exons_contigs.add(exon(hit) + contig(hit))
        with open(f"{main_path}{sample}.fas", "w") as result_fasta, open(
            file
        ) as contigs:
            contigs_fasta_parsed = SeqIO.to_dict(
                SeqIO.parse(contigs, "fasta", generic_dna)
            )
            for hit_dedup in hits_dedup:
                if int(hit_dedup.split()[6]) > int(hit_dedup.split()[7]):
                    sequence: str = slicing(contigs_fasta_parsed, hit_dedup, True)
                else:
                    sequence: str = slicing(contigs_fasta_parsed, hit_dedup, False)
                result_fasta.write(
                    f">{exon(hit_dedup)}_{'_'.join(contig(hit_dedup).split('_')[1:])}\n{sequence}\n"
                )
                all_hits_for_reference.append(f"{hit_dedup}\t{sample}\t{sequence}")
        sort_hit_table_cover(hits_dedup, "locus")
        hits_loci: Set[str] = set()
        hits_exons: Set[str] = set()
        for hit_dedup in hits_dedup:
            if locus(hit_dedup) not in hits_loci:
                statistics[sample][locus(hit_dedup)] = dict()
                if exon(hit_dedup) not in hits_exons:
                    statistics[sample][locus(hit_dedup)][exon(hit_dedup)] = [
                        contig(hit_dedup)
                    ]
                    hits_exons.add(exon(hit_dedup))
                else:
                    statistics[sample][locus(hit_dedup)][exon(hit_dedup)].append(
                        contig(hit_dedup)
                    )
                hits_loci.add(locus(hit_dedup))
            else:
                if exon(hit_dedup) not in hits_exons:
                    statistics[sample][locus(hit_dedup)][exon(hit_dedup)] = [
                        contig(hit_dedup)
                    ]
                    hits_exons.add(exon(hit_dedup))
                else:
                    statistics[sample][locus(hit_dedup)][exon(hit_dedup)].append(
                        contig(hit_dedup)
                    )
        for i in statistics[sample].keys():
            stat = []
            for j in statistics[sample][i].keys():
                stat.append(len(statistics[sample][i][j]))
            statistics[sample][i]: int = max(stat)
        with open(
            f"{main_path}reference_against_{sample}_contigs.txt", "w"
        ) as hittable:
            sort_hit_table_cover(hits, "locus")
            for hit in hits:
                hittable.write(f"{hit}\n")
    with open(f"{path_to_data}/exons/all_hits.txt", "w") as all_hits_to_write:
        for hit in all_hits_for_reference:
            all_hits_to_write.write(f"{hit}\n")
    logger.info(" OK")
    logger.info("All contigs were successfully corrected!\n")


def write_stats(path_to_data, probe_exons, statistics, logger):
    """

    :param path_to_data:
    :type path_to_data:
    :param probe_exons:
    :type probe_exons:
    :param statistics:
    :type statistics:
    """
    logger.info("Writing statistics...")
    main_path: str = f"{path_to_data}/exons/40contigs/"
    os.makedirs(main_path, exist_ok=True)
    with open(f"{main_path}statistics.tsv", "w") as stats, open(
        probe_exons
    ) as reference:
        stats_dict: Dict[str, str] = dict([("gene\t", "")])
        loci: Set[str] = set()
        samples: List[str] = list()
        reference_as_list: List[str] = [x[:-1] for x in reference.readlines()]
        for line in reference_as_list:
            if line.startswith(">"):
                loci.add("_".join(line.split("-")[1].split("_")[:-2]))
        for key in statistics.keys():
            samples.append(key)
        samples.sort()
        for sample in samples:
            stats_dict["gene\t"]: str = stats_dict["gene\t"] + sample + "\t"
        for loc in loci:
            stats_dict[loc + "\t"]: str = ""
            for sample in samples:
                loci_in_sample: Set[str] = set(statistics[sample].keys())
                if loc in loci_in_sample:
                    stats_dict[loc + "\t"]: str = (
                        stats_dict[loc + "\t"] + str(statistics[sample][loc]) + "\t"
                    )
                else:
                    stats_dict[loc + "\t"]: str = stats_dict[loc + "\t"] + "NA" + "\t"
        stats.write("gene\t" + stats_dict["gene\t"] + "\n")
        del stats_dict["gene\t"]
        for key in sorted(list(stats_dict.keys())):
            stats.write(f"{key}{stats_dict[key]}\n")
    logger.info("Statistics file created!\n")


def rename_contigs(path_to_data, logger):
    """

    :param path_to_data:
    :type path_to_data:
    """
    logger.info("Renaming contigs...")
    main_path: str = f"{path_to_data}/exons/40contigs/"
    for file in glob.glob(f"{main_path}*.fas"):
        sample = os.path.basename(os.path.splitext(file)[0])
        logger.info(f" Processing {sample}")
        with open(file) as result_fasta:
            fasta_parsed = SeqIO.to_dict(
                SeqIO.parse(result_fasta, "fasta", generic_dna)
            )
            counter: int = 1
            fasta_to_write: List[str] = list()
            for line in sorted(fasta_parsed.keys()):
                fasta_to_write.append(
                    f">Contig{str(counter)}_{sample}-{line.replace('_', '-')}\n"
                )
                fasta_to_write.append(str(fasta_parsed[line].seq) + "\n")
                counter += 1
        with open(file, "w") as result_fasta:
            result_fasta.writelines(fasta_to_write)
        logger.info(" OK")
    logger.info("All contigs were successfully renamed!\n")


def clean(path_to_data, logger):
    """

    :param path_to_data:
    :type path_to_data:
    """
    logger.info("Removing temporary files...")
    main_path: str = f"{path_to_data}/exons/40contigs/"
    for file in glob.glob(f"{main_path}*.fasta"):
        os.remove(file)
    for file in glob.glob(f"{main_path}reference_in*"):
        os.remove(file)
    for file in glob.glob(f"{main_path}*.n*"):
        os.remove(file)
    logger.info("Done\n")


def build_alignments(path_to_data, n_cpu, logger):
    """

    :param path_to_data:
    :type path_to_data:
    :param n_cpu:
    :type n_cpu:
    """
    logger.info("Building individual exon alignments...")
    with open(path_to_data + "/exons/all_hits.txt") as all_hits:
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
    os.makedirs(path_to_data + "/exons/aln_orth_par", exist_ok=True)
    for key in exons.keys():
        for record in exons[key]:
            new_id = f"{key}_N_{record.id.split('_N_')[1]}"
            record.id = new_id
    for key in exons.keys():
        SeqIO.write(
            exons[key], f"{path_to_data}/exons/aln_orth_par/{key}.fasta", "fasta"
        )
        stdout, stderr = MafftCommandline(
            input=f"{path_to_data}/exons/aln_orth_par/{key}.fasta",
            adjustdirectionaccurately=True,
            auto=True,
            thread=n_cpu,
        )()
        with open(
            f"{path_to_data}/exons/aln_orth_par/{key}.mafft.fasta", "w"
        ) as aligned:
            aligned.write(stdout.replace(">_R_", ">"))
        try:
            FastTreeCommandline(
                "fasttreemp",
                nt=True,
                gtr=True,
                input=f"{path_to_data}/exons/aln_orth_par/{key}.mafft.fasta",
                out=f"{path_to_data}/exons/aln_orth_par/{key}.mafft.fasta.tre",
            )()
        except Bio.Application.ApplicationError:
            FastTreeCommandline(
                "fasttree",
                nt=True,
                gtr=True,
                input=f"{path_to_data}/exons/aln_orth_par/{key}.mafft.fasta",
                out=f"{path_to_data}/exons/aln_orth_par/{key}.mafft.fasta.tre",
            )()

    logger.info("Done\n")


def estimate_divergence(path_to_data, blocklist, logger):
    """

    :param path_to_data:
    :type path_to_data:
    :param blocklist:
    :type blocklist:
    """
    logger.info("Estimating divergence of paralogs...")
    divergency_distribution: List[float] = []
    divergencies_to_write: List[str] = []
    matplotlib.use("Agg")
    fig, axis = matplotlib.pyplot.subplots(figsize=(15, 10))
    for file in sorted(glob.glob(f"{path_to_data}/exons/aln_orth_par/*.mafft.fasta")):
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
        f"{path_to_data}/exons/aln_orth_par/individual_distributions.png",
        dpi=300,
        format="png",
    )
    fig.savefig(
        f"{path_to_data}/exons/aln_orth_par/individual_distributions.svg",
        dpi=300,
        format="svg",
    )
    matplotlib.pyplot.close(fig)
    with open(
        path_to_data + "/exons/aln_orth_par/pairwise_distances.txt", "w"
    ) as divergency_distribution_to_write:
        for i in divergencies_to_write:
            divergency_distribution_to_write.write(str(i) + "\n")
    divergency_distribution_array: numpy.ndarray = numpy.array(
        [[x] for x in divergency_distribution]
    )
    get_plot(
        path_to_data + "/exons/aln_orth_par/",
        "pairwise_distances_distribution_1_comp",
        divergency_distribution_array,
        1,
    )
    get_plot(
        path_to_data + "/exons/aln_orth_par/",
        "pairwise_distances_distribution_2_comp",
        divergency_distribution_array,
        2,
    )
    get_plot(
        path_to_data + "/exons/aln_orth_par/",
        "pairwise_distances_distribution_3_comp",
        divergency_distribution_array,
        3,
    )

    logger.info("Done\n")


def score_samples(list_with_hits: List[str]) -> List[str]:
    """

    :param list_with_hits:
    :type list_with_hits:
    :return:
    :rtype:
    """
    sort_hit_table_ident(list_with_hits, "exon")
    # Add counters to hits for each exon
    counter: int = 0
    current_exon: str = ""
    for hit in list_with_hits:
        if exon(hit) != current_exon:
            counter = 1
            current_exon = exon(hit)
        else:
            counter += 1
        list_with_hits[list_with_hits.index(hit)] = f"{hit}\t{counter}"
    # Prepare dictionary for counting scores for contigs in each locus
    for_samples_scores = dict()
    samples: Set[str] = set()
    loci_samples: Set[str] = set()
    for hit in list_with_hits:
        sample: str = hit.split()[8]
        locus: str = exon(hit).split("_")[0]
        locus_sample: str = locus + sample
        spades_contig: str = contig(hit)
        if sample not in samples:
            for_samples_scores[sample] = dict()
            for_samples_scores[sample][locus] = [dict(), set()]
            for_samples_scores[sample][locus][0][spades_contig] = []
            samples.add(sample)
            loci_samples.add(locus_sample)
        else:
            if locus_sample not in loci_samples:
                for_samples_scores[sample][locus] = [dict(), set()]
                for_samples_scores[sample][locus][0][spades_contig] = []
                loci_samples.add(locus_sample)
            else:
                for_samples_scores[sample][locus][0][spades_contig] = []
    # Count scores for contigs in each locus
    for hit in list_with_hits:
        sample: str = hit.split()[8]
        locus: str = exon(hit).split("_")[0]
        spades_contig: str = contig(hit)
        for_samples_scores[sample][locus][0][spades_contig].append(int(hit.split()[-1]))
        for_samples_scores[sample][locus][1].add(exon(hit))
    samples_scores = copy.deepcopy(for_samples_scores)
    for sample in for_samples_scores.keys():
        for locus in for_samples_scores[sample].keys():
            samples_scores[sample][locus][1]: int = len(
                for_samples_scores[sample][locus][1]
            )
            sample_locus_scores: List[float] = []
            for spades_contig in for_samples_scores[sample][locus][0].keys():
                samples_scores[sample][locus][0][spades_contig] = sum(
                    for_samples_scores[sample][locus][0][spades_contig]
                ) / len(for_samples_scores[sample][locus][0][spades_contig])
                sample_locus_scores.append(
                    samples_scores[sample][locus][0][spades_contig]
                )
            sample_locus_score: float = sum(sample_locus_scores) / len(
                sample_locus_scores
            )
            samples_scores[sample][locus].append(sample_locus_score)
    # Adding scores to each hit
    for i in range(0, len(list_with_hits)):
        hit: str = list_with_hits[i]
        locus: str = exon(hit).split("_")[0]
        spades_contig: str = contig(hit)
        sample: str = hit.split()[8]
        list_with_hits[i] = (
            f"{hit}\t{str(samples_scores[sample][locus][1])}\t{str(samples_scores[sample][locus][2])}"
            f"\t{str(samples_scores[sample][locus][0][spades_contig])}"
        )
    list_with_hits.sort(key=lambda x: float(x.split()[-1]))
    list_with_hits.sort(key=lambda x: x.split()[8])
    list_with_hits.sort(key=lambda x: float(x.split()[-2]))
    list_with_hits.sort(key=lambda x: float(x.split()[-3]), reverse=True)
    list_with_hits.sort(key=lambda x: exon(x))
    return ["\t".join(x.split()[:10]) for x in list_with_hits]


def create_reference_wo_paralogs(
    path_to_data, all_hits_for_reference_scored, blocklist, logger
):
    """

    :param path_to_data:
    :type path_to_data:
    :param all_hits_for_reference_scored:
    :type all_hits_for_reference_scored:
    :param blocklist:
    :type blocklist:
    """
    logger.info("Creating new reference...")
    exons: Set[str] = set()
    with open(
        f"{path_to_data}/exons/new_reference_for_HybPhyloMaker.fas", "w"
    ) as new_reference_HPM, open(
        f"{path_to_data}/exons/new_reference_for_HybPiper.fas", "w"
    ) as fasta_to_concatenate:
        for hit in all_hits_for_reference_scored:
            sample: str = hit.split()[8]
            if sample not in blocklist:
                if exon(hit) not in exons:
                    name_of_locus: str = (
                        exon(hit)
                        .replace("exon", "Contig")
                        .replace("Exon", "Contig")
                        .replace("contig", "Contig")
                        .replace("_", "")
                        .replace("Contig", "_Contig_")
                    )
                    new_reference_HPM.write(
                        f">Assembly_{name_of_locus}_{sample}_{contig(hit)}\n{hit.split()[9]}\n"
                    )
                    fasta_to_concatenate.write(
                        f">{sample.replace('-', '_')}_{contig(hit)}-{name_of_locus.replace('Contig', 'exon')}\n"
                        f"{hit.split()[9]}\n"
                    )
                    exons.add(exon(hit))
    with open(
        f"{path_to_data}/exons/new_reference_for_HybPiper.fas"
    ) as fasta_to_concatenate, open(
        f"{path_to_data}/exons/new_reference_for_HybPiper_concatenated.fas", "w"
    ) as concatenated_fasta:
        fasta_parsed = SeqIO.to_dict(
            SeqIO.parse(fasta_to_concatenate, "fasta", generic_dna)
        )
        current_locus = ""
        list_of_keys = list(fasta_parsed.keys())
        list_of_keys.sort(key=lambda x: int(x.split("-")[1].split("_")[2]))
        list_of_keys.sort(key=lambda x: x.split("-")[1].split("_")[0])
        count = 1
        for key in list_of_keys:
            locus = key.split("-")[1].split("_")[0]
            if count == 1:
                if locus != current_locus:
                    concatenated_fasta.write(
                        ">"
                        + "_".join(key.split("-")[0].split("_")[0:3])
                        + "-"
                        + locus
                        + "\n"
                        + str(fasta_parsed[key].seq)
                    )
                else:
                    concatenated_fasta.write(str(fasta_parsed[key].seq))
            else:
                if locus != current_locus:
                    concatenated_fasta.write(
                        "\n>"
                        + "_".join(key.split("-")[0].split("_")[0:3])
                        + "-"
                        + locus
                        + "\n"
                        + str(fasta_parsed[key].seq)
                    )
                else:
                    concatenated_fasta.write(str(fasta_parsed[key].seq))
            current_locus = locus
            count += 1

    logger.info("New reference created!\n")


def create_reference_w_paralogs(
    path_to_data,
    all_hits_for_reference_scored,
    paralog_statistic,
    paralog_min_divergence,
    paralog_max_divergence,
    blocklist,
    logger,
):
    """

    :param path_to_data:
    :type path_to_data:
    :param all_hits_for_reference_scored:
    :type all_hits_for_reference_scored:
    :param paralog_statistic:
    :type paralog_statistic:
    :param paralog_min_divergence:
    :type paralog_min_divergence:
    :param paralog_max_divergence:
    :type paralog_max_divergence:
    :param blocklist:
    :type blocklist:
    """
    with open(f"{path_to_data}/exons/aln_orth_par/pairwise_distances.txt") as distances:
        pairwise_distances: Dict[str, float] = dict()
        for line in distances.read().splitlines():
            pairwise_distances[f"{line.split()[0]}_{line.split()[2]}"] = float(
                line.split()[1]
            )
            pairwise_distances[f"{line.split()[2]}_{line.split()[0]}"] = float(
                line.split()[1]
            )
    logger.info("Creating new reference...")
    all_paralogs_for_reference = []
    exons: Set[str] = set()
    count: int = 0
    paralog_found: bool = bool()
    current_best: str = ""
    current_locus: Dict[str, str] = dict()
    samples_current_locus: Set[str] = set()
    samples_with_paralogs: Set[str] = set()
    set_of_samples = set()
    for hit in all_hits_for_reference_scored:
        sample: str = hit.split()[8]
        if sample not in set_of_samples:
            set_of_samples.add(sample)
            paralog_statistic[sample] = set()
        if exon(hit) not in exons:
            if not paralog_found and count != 0:
                logger.info(
                    f"No paralog found for {current_best.split()[0].split('-')[1]}"
                )
                for sample_wo_paralog in current_locus.keys():
                    all_paralogs_for_reference.append(current_locus[sample_wo_paralog])
            current_locus: Dict[str, str] = dict()
            samples_current_locus: Set[str] = set()
            samples_with_paralogs: Set[str] = set()
            current_best: str = hit
            current_locus[sample] = hit
            samples_current_locus.add(sample)
            paralog_found: bool = False
            exons.add(exon(hit))
        else:
            if sample in samples_current_locus and sample not in samples_with_paralogs:
                current_exonic_contig = (
                    f"{hit.split()[0].split('-')[1]}_"
                    f"{'_'.join(hit.split()[1].split('_')[1:])}_{sample}"
                )
                exonic_contigs_to_compare = (
                    f"{current_locus[sample].split()[0].split('-')[1]}_"
                    f"{'_'.join(current_locus[sample].split()[1].split('_')[1:])}"
                    f"_{sample}"
                )
                if (
                    paralog_max_divergence
                    > pairwise_distances[
                        f"{current_exonic_contig}_{exonic_contigs_to_compare}"
                    ]
                    > paralog_min_divergence
                ):
                    logger.info(
                        f"Paralog detected for {hit.split()[0].split('-')[1]} in {sample}"
                    )
                    paralog_statistic[sample].add(hit.split()[1].split("_")[0])
                    all_paralogs_for_reference.append(current_locus[sample])
                    name_of_locus_para: str = "_".join(
                        [exon(hit).split("_")[0] + "_para"] + exon(hit).split("_")[1:]
                    )
                    name_of_locus_para = (
                        hit.split()[0].split("-")[0] + "-" + name_of_locus_para
                    )
                    hit = "\t".join([name_of_locus_para] + hit.split()[1:])
                    all_paralogs_for_reference.append(hit)
                    samples_with_paralogs.add(sample)
                    paralog_found: bool = True
            else:
                current_locus[sample] = hit
                samples_current_locus.add(sample)
        count += 1
    all_paralogs_for_reference = score_samples(all_paralogs_for_reference)
    exons: Set[str] = set()
    with open(
        f"{path_to_data}/exons/new_reference_for_HybPhyloMaker_"
        f"div_{paralog_min_divergence}_{paralog_max_divergence}.fas",
        "w",
    ) as new_reference:
        for hit in all_paralogs_for_reference:
            sample = hit.split()[8]
            if sample not in blocklist:
                if exon(hit) not in exons:
                    name_of_locus = (
                        exon(hit)
                        .replace("exon", "Contig")
                        .replace("Exon", "Contig")
                        .replace("contig", "Contig")
                        .replace("_", "")
                        .replace("Contig", "_Contig_")
                    )
                    new_reference.write(
                        f">Assembly_{name_of_locus}_{sample}_{contig(hit)}\n{hit.split()[9]}\n"
                    )
                    exons.add(exon(hit))
        logger.info("New reference created!\n")


def refine_phasing(
    path_to_data, paralog_min_divergence, paralog_max_divergence, logger
):
    """

    :param path_to_data:
    :type path_to_data:
    :param paralog_min_divergence:
    :type paralog_min_divergence:
    :param paralog_max_divergence:
    :type paralog_max_divergence:
    """
    with open(
        f"{path_to_data}/exons/new_reference_for_HybPhyloMaker_"
        f"div_{paralog_min_divergence}_{paralog_max_divergence}.fas",
    ) as reference_to_check:
        new_ref_parsed = SeqIO.to_dict(
            SeqIO.parse(reference_to_check, "fasta", generic_dna)
        )
        all_seqs = set(new_ref_parsed.keys())
        initial_dict_to_check = dict()
        for seq in all_seqs:
            locus = seq.split("_")[1]
            exon = "_".join(seq.split("_")[2:4])
            name_contig = "_".join(seq.split("_")[4:])
            contig = (
                name_contig,
                new_ref_parsed[f"Assembly_{locus}_{exon}_{name_contig}"],
            )
            if locus not in initial_dict_to_check.keys():
                initial_dict_to_check[locus] = dict()
                initial_dict_to_check[locus][exon] = contig
            else:
                initial_dict_to_check[locus][exon] = contig
    final_dict_to_check = dict()
    for key1 in initial_dict_to_check.keys():
        if key1[-4:] != "para":
            final_dict_to_check[key1] = dict()
            final_dict_to_check[key1]["non-para"] = dict()
            final_dict_to_check[key1]["non-para"]["all"] = list()
            for key2 in initial_dict_to_check[key1]:
                final_dict_to_check[key1]["non-para"][key2] = initial_dict_to_check[
                    key1
                ][key2]
                final_dict_to_check[key1]["non-para"]["all"].append(
                    initial_dict_to_check[key1][key2][0]
                )
            if f"{key1}para" in initial_dict_to_check.keys():
                final_dict_to_check[key1]["para"] = dict()
                final_dict_to_check[key1]["para"]["all"] = list()

                for key2 in initial_dict_to_check[f"{key1}para"]:
                    final_dict_to_check[key1]["para"][key2] = initial_dict_to_check[
                        f"{key1}para"
                    ][key2]
                    final_dict_to_check[key1]["para"]["all"].append(
                        initial_dict_to_check[f"{key1}para"][key2][0]
                    )
    checked_dict = copy.deepcopy(final_dict_to_check)
    warn = list()
    for key in final_dict_to_check.keys():
        for item in [
            x for x in final_dict_to_check[key]["non-para"].keys() if x != "all"
        ]:
            if "para" in final_dict_to_check[key].keys():
                if (
                    final_dict_to_check[key]["non-para"][item][0]
                    in final_dict_to_check[key]["para"]["all"]
                ):
                    swap1 = final_dict_to_check[key]["non-para"][item]
                    checked_dict[key]["para"][item] = swap1
                    if item not in final_dict_to_check[key]["para"].keys():
                        del checked_dict[key]["non-para"][item]
                        warn.append(
                            f"Moving {swap1[0]} to a paralog for {key} {item}\n"
                        )
                    else:
                        swap2 = final_dict_to_check[key]["para"][item]
                        checked_dict[key]["non-para"][item] = swap2
                        warn.append(
                            f"Swaping {swap1[0]} and {swap2[0]} for {key} {item}\n"
                        )

    if len(warn) > 0:
        with open(f"{path_to_data}/exons/warnings.txt", "w") as warnings:
            warnings.write(
                "Following genes seems to be phased improperly based on similarity to reference. "
                "Refining attempted.\n"
            )
            for x in warn:
                warnings.write(x)
        with open(
            f"{path_to_data}/exons/refined_new_reference_for_HybPhyloMaker_"
            f"div_{paralog_min_divergence}_{paralog_max_divergence}.fas",
            "w",
        ) as reference_to_write:
            for key in sorted(list(checked_dict.keys())):
                for item in [
                    x
                    for x in sorted(list(checked_dict[key]["non-para"].keys()))
                    if x != "all"
                ]:
                    name = f"Assembly_{key}_{item}_{checked_dict[key]['non-para'][item][0]}"
                    sequence = str(checked_dict[key]["non-para"][item][1].seq)
                    reference_to_write.write(f">{name}\n{sequence}\n")
                if "para" in checked_dict[key].keys():
                    for item in [
                        x
                        for x in sorted(list(checked_dict[key]["para"].keys()))
                        if x != "all"
                    ]:
                        name = f"Assembly_{key}para_{item}_{checked_dict[key]['para'][item][0]}"
                        sequence = str(checked_dict[key]["para"][item][1].seq)
                        reference_to_write.write(f">{name}\n{sequence}\n")
        with open(
            f"{path_to_data}/exons/refined_new_reference_for_HybPhyloMaker_"
            f"div_{paralog_min_divergence}_{paralog_max_divergence}.fas",
        ) as reference_to_check:
            new_ref_parsed = SeqIO.to_dict(
                SeqIO.parse(reference_to_check, "fasta", generic_dna)
            )
            all_seqs = set(new_ref_parsed.keys())
            initial_dict_to_check = dict()
            for seq in all_seqs:
                locus = seq.split("_")[1]
                exon = "_".join(seq.split("_")[2:4])
                name_contig = "_".join(seq.split("_")[4:])
                contig = (
                    name_contig,
                    new_ref_parsed[f"Assembly_{locus}_{exon}_{name_contig}"],
                )
                if locus not in initial_dict_to_check.keys():
                    initial_dict_to_check[locus] = dict()
                    initial_dict_to_check[locus][exon] = contig
                else:
                    initial_dict_to_check[locus][exon] = contig
        final_dict_to_check = dict()
        for key1 in initial_dict_to_check.keys():
            if key1[-4:] != "para":
                final_dict_to_check[key1] = dict()
                final_dict_to_check[key1]["non-para"] = dict()
                final_dict_to_check[key1]["non-para"]["all"] = list()
                for key2 in initial_dict_to_check[key1]:
                    final_dict_to_check[key1]["non-para"][key2] = initial_dict_to_check[
                        key1
                    ][key2]
                    final_dict_to_check[key1]["non-para"]["all"].append(
                        initial_dict_to_check[key1][key2][0]
                    )
                if f"{key1}para" in initial_dict_to_check.keys():
                    final_dict_to_check[key1]["para"] = dict()
                    final_dict_to_check[key1]["para"]["all"] = list()

                    for key2 in initial_dict_to_check[f"{key1}para"]:
                        final_dict_to_check[key1]["para"][key2] = initial_dict_to_check[
                            f"{key1}para"
                        ][key2]
                        final_dict_to_check[key1]["para"]["all"].append(
                            initial_dict_to_check[f"{key1}para"][key2][0]
                        )
        checked_dict = copy.deepcopy(final_dict_to_check)
        warn = list()
        for key in final_dict_to_check.keys():
            for item in [
                x for x in final_dict_to_check[key]["non-para"].keys() if x != "all"
            ]:
                if "para" in final_dict_to_check[key].keys():
                    if (
                        final_dict_to_check[key]["non-para"][item][0]
                        in final_dict_to_check[key]["para"]["all"]
                    ):
                        swap1 = final_dict_to_check[key]["non-para"][item]
                        checked_dict[key]["para"][item] = swap1
                        if item not in final_dict_to_check[key]["para"].keys():
                            del checked_dict[key]["non-para"][item]
                            warn.append(f"{key} {item}\n")
                        else:
                            swap2 = final_dict_to_check[key]["para"][item]
                            checked_dict[key]["non-para"][item] = swap2
                            warn.append(f"{key} {item}\n")
        if len(warn) > 0:
            with open(f"{path_to_data}/exons/warnings.txt", "a") as warnings:
                warnings.write(
                    "\nFollowing genes and exons seems to be impossible to phase properly with the current algorithms. "
                    "Consider removing.\n"
                )
                for x in warn:
                    warnings.write(x)
                warnings.write(
                    "Pay closer attention to the loci mentioned above. It is worth having a look at the "
                    "alignments for particular exons in aln_orth_par/"
                )
        with open(
            f"{path_to_data}/exons/refined_new_reference_for_HybPhyloMaker_"
            f"div_{paralog_min_divergence}_{paralog_max_divergence}.fas",
            "w",
        ) as reference_to_write:
            for key in sorted(list(checked_dict.keys())):
                for item in [
                    x
                    for x in sorted(list(checked_dict[key]["non-para"].keys()))
                    if x != "all"
                ]:
                    name = f"Assembly_{key}_{item}_{checked_dict[key]['non-para'][item][0]}"
                    sequence = str(checked_dict[key]["non-para"][item][1].seq)
                    reference_to_write.write(f">{name}\n{sequence}\n")
                if "para" in checked_dict[key].keys():
                    for item in [
                        x
                        for x in sorted(list(checked_dict[key]["para"].keys()))
                        if x != "all"
                    ]:
                        name = f"Assembly_{key}para_{item}_{checked_dict[key]['para'][item][0]}"
                        sequence = str(checked_dict[key]["para"][item][1].seq)
                        reference_to_write.write(f">{name}\n{sequence}\n")
    else:
        with open(f"{path_to_data}/exons/warnings.txt", "w") as warnings:
            warnings.write("No warnings generated.\n")


def write_paralog_stats(
    path_to_data,
    paralog_min_divergence,
    paralog_max_divergence,
    paralog_statistic,
    probes,
    logger,
):
    """

    :param path_to_data:
    :type path_to_data:
    :param paralog_min_divergence:
    :type paralog_min_divergence:
    :param paralog_max_divergence:
    :type paralog_max_divergence:
    :param paralog_statistic:
    :type paralog_statistic:
    :param probes:
    :type probes:
    """
    with open(
        f"{path_to_data}/exons/paralog_statistics_div_{paralog_min_divergence}_{paralog_max_divergence}.tsv",
        "w",
    ) as par_stat:
        for sample in sorted(list(paralog_statistic.keys())):
            par_stat.write(f"{sample}\t{str(len(paralog_statistic[sample]))}\n")
    with open(probes) as probe_loci:
        all_probe_exons = SeqIO.to_dict(
            SeqIO.parse(probe_loci, "fasta", generic_dna)
        ).keys()
        all_loci: Set[str] = set()
        for key in all_probe_exons:
            all_loci.add(key.split("-")[1].split("_")[0])
    with open(
        f"{path_to_data}/exons/locus_statistics_div_{paralog_min_divergence}_{paralog_max_divergence}.tsv",
        "w",
    ) as loci_par_stat:
        loci_par_stat.write(r"sample\locus")
        for locus in sorted(list(all_loci)):
            loci_par_stat.write(f"\t{locus}")
        loci_par_stat.write("\n")
        for sample in sorted(list(paralog_statistic.keys())):
            loci_par_stat.write(sample)
            for locus in sorted(list(all_loci)):
                if locus in paralog_statistic[sample]:
                    loci_par_stat.write("\tYes")
                else:
                    loci_par_stat.write("\tN/A")
            loci_par_stat.write("\n")


def aln_similarity(blat_hit_string: str) -> float:
    """

    :param blat_hit_string:
    :type blat_hit_string:
    :return:
    :rtype:
    """
    return (
        float(blat_hit_string.split()[0])
        / (float(blat_hit_string.split()[0]) + float(blat_hit_string.split()[1]))
    ) * 100


def number_locus_for_sort(blat_hit_string: str) -> str:
    """

    :rtype: object
    """
    return blat_hit_string.split()[13].split("_")[1]


def number_exon_for_sort(blat_hit_string: str) -> int:
    """

    :param blat_hit_string:
    :type blat_hit_string:
    :return:
    :rtype:
    """
    return int(blat_hit_string.split()[13].split("_")[3])


def number_contig_for_sort(blat_hit_string: str) -> int:
    """

    :param blat_hit_string:
    :type blat_hit_string:
    :return:
    :rtype:
    """
    return int(blat_hit_string.split()[9].split("_")[0][6:])


def run_blat(path_to_data, probes, minident, logger):
    """

    :param path_to_data:
    :type path_to_data:
    :param probes:
    :type probes:
    :param minident:
    :type minident:
    """
    logger.info("Generating pslx files using BLAT...\n")
    os.makedirs(f"{path_to_data}/exons/50pslx", exist_ok=True)
    for contigfile in glob.glob(f"{path_to_data}/exons/40contigs/*.fas"):
        file = os.path.basename(contigfile)
        if file != probes:
            logger.info(f"Processing {file}...")
            blat_cmd_output = os.popen(
                f"blat -t=DNA -q=DNA -out=pslx -minIdentity={minident} {probes} {contigfile} {contigfile}.pslx"
            ).read()
            logger.info(blat_cmd_output)
            shutil.move(
                f"{path_to_data}/exons/40contigs/{file}.pslx",
                f"{path_to_data}/exons/50pslx/{file}.pslx",
            )
            logger.info("Done")


def correct(path_to_data, probes, redlist, logger):
    """

    :rtype: object
    """
    os.makedirs(f"{path_to_data}/exons/50pslx/corrected", exist_ok=True)
    for file in glob.glob(f"{path_to_data}/exons/50pslx/*.pslx"):
        if os.path.basename(file) != f"{probes}.pslx":
            with open(file) as pslx_file, open(
                f"{path_to_data}/exons/50pslx/corrected/{os.path.basename(file)}",
                "w",
            ) as corrected_pslx_file:
                file = pslx_file.readlines()
                head = file[0:5]
                list_to_work = file[5:]
                list_to_work.sort(key=aln_similarity, reverse=True)
                list_to_work.sort(key=number_exon_for_sort)
                list_to_work.sort(key=number_locus_for_sort)
                list_to_work_cleaned1 = []
                hits1 = set()
                for line in list_to_work:
                    if line.split()[13] not in hits1:
                        list_to_work_cleaned1.append(line)
                        hits1.add(line.split()[13])
                list_to_work_cleaned1.sort(key=aln_similarity, reverse=True)
                list_to_work_cleaned1.sort(key=number_contig_for_sort)
                list_to_work_cleaned2 = []
                hits2 = set()
                for line in list_to_work_cleaned1:
                    sample = f"{line.split()[9].split('-')[0].split('_')[1]}-{line.split()[9].split('-')[1]}"
                    if sample not in redlist:
                        if line.split()[9] not in hits2:
                            list_to_work_cleaned2.append(line)
                            hits2.add(line.split()[9])
                    else:
                        list_to_work_cleaned2.append(line)
                for line in head:
                    corrected_pslx_file.write(line)
                for line in list_to_work_cleaned2:
                    corrected_pslx_file.write(line)


def main():
    arguments = ParsedArgs().get_args_dict()
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    log_handler = logging.FileHandler(
        f'ParalogWizard_{arguments["command"]}_{datetime.now().strftime("%d.%b.%y_%H.%M")}.log',
        "w",
    )
    log_handler.setLevel(logging.INFO)
    log_formatter = logging.Formatter(
        fmt="%(asctime)s - %(message)s", datefmt="%d-%b-%y %H:%M:%S"
    )
    log_handler.setFormatter(log_formatter)
    logger.addHandler(log_handler)
    if arguments["command"] == "cast_collect":
        logger.info(
            f"""ParalogWizard cast_collect running with the following settings
            main data folder - {arguments["data_folder"]}
            folder with contigs - {arguments["contig_folder"]}"""
        )
        collect_contigs( arguments["contig_folder"], arguments["data_folder"], logger)
    elif arguments["command"] == "cast_retrieve":
        logger.info(
            f"""ParalogWizard cast_collect running with the following settings
            main data folder - {arguments["data_folder"]}
            probe file - {arguments["probes_exons"]}
            filter for blast length cover - {arguments["length_cover"]}
            k-mer cover threshold for spades contigs - {arguments["spades_cover"]}
            number of used cores - {arguments["num_cores"]}"""
        )
        logger.info("Retrieving data...\n")
        statistics: Dict[str, Dict[str, Union[Dict[str, List[str]], int]]] = dict()
        all_hits_for_reference: List[str] = list()
        prepare_contigs(
            arguments["data_folder"],
            logger,
        )
        create_hit_tables(
            arguments["data_folder"],
            arguments["probes_exons"],
            arguments["length_cover"],
            arguments["num_cores"],
            logger,
        )
        correct_contgis(
            arguments["data_folder"],
            statistics,
            arguments["spades_cover"],
            all_hits_for_reference,
            logger,
        )
        write_stats(
            arguments["data_folder"], arguments["probes_exons"], statistics, logger
        )
        rename_contigs(arguments["data_folder"], logger)
        clean(arguments["data_folder"], logger)
        logger.info("Data was successfully retrieved!")
    elif arguments["command"] == "cast_analyze":
        if len(arguments["blocklist"]) > 0:
            blocklist_string = ", ".join(sp for sp in list(arguments["blocklist"]))
            logger.info(
                f"""ParalogWizard cast_collect running with the following settings
            main data folder - {arguments["data_folder"]}
            species not taken to paralogs divergency estimation - {blocklist_string}
            number of used cores - {arguments["num_cores"]}"""
            )
        else:
            logger.info(
                f"""ParalogWizard cast_collect running with the following settings
            main data folder - {arguments["data_folder"]}
            all species taken to paralogs divergency estimation
            number of used cores - {arguments["num_cores"]}"""
            )
        # build_alignments(arguments["data_folder"], arguments["num_cores"], logger)
        estimate_divergence(arguments["data_folder"], arguments["blocklist"], logger)
    elif arguments["command"] == "cast_create":
        with open(f"{arguments['data_folder']}/exons/all_hits.txt") as all_hits:
            all_hits_for_reference: List[str] = [x[:-1] for x in all_hits.readlines()]
        all_hits_for_reference_scored = score_samples(all_hits_for_reference)
        if not arguments["paralogs"]:
            if len(arguments["blocklist"]) > 0:
                blocklist_string = ", ".join(sp for sp in list(arguments["blocklist"]))
                logger.info(
                    f"""ParalogWizard cast_collect running with the following settings
                main data folder - {arguments["data_folder"]}
                paralogs are not being searched
                species not taken to paralogs divergency estimation - {blocklist_string}"""
                )
            else:
                logger.info(
                    f"""ParalogWizard cast_collect running with the following settings
                main data folder - {arguments["data_folder"]}
                paralogs are not being searched
                all species taken to paralogs divergency estimation"""
                )
            create_reference_wo_paralogs(
                arguments["data_folder"],
                all_hits_for_reference_scored,
                arguments["blocklist"],
                logger,
            )
        else:
            if len(arguments["blocklist"]) == 0:
                blocklist_string = ", ".join(sp for sp in list(arguments["blocklist"]))
                logger.info(
                    f"""ParalogWizard cast_collect running with the following settings
                main data folder - {arguments["data_folder"]}
                paralogs are searched with minium {arguments["minimum_divergence"]} and maximum {arguments["maximum_divergence"]} divergence
                species not taken to paralogs divergency estimation - {blocklist_string}"""
                )
            else:
                logger.info(
                    f"""ParalogWizard cast_collect running with the following settings
                main data folder - {arguments["data_folder"]}
                paralogs are searched with minium {arguments["minimum_divergence"]} and maximum {arguments["maximum_divergence"]} divergence
                all species taken to paralogs divergency estimation"""
                )
            paralog_statistic: Dict[str, Set[str]] = dict()
            create_reference_w_paralogs(
                arguments["data_folder"],
                all_hits_for_reference_scored,
                paralog_statistic,
                arguments["minimum_divergence"],
                arguments["maximum_divergence"],
                arguments["blocklist"],
                logger,
            )
            refine_phasing(
                arguments["data_folder"],
                arguments["minimum_divergence"],
                arguments["maximum_divergence"],
                logger,
            )
            write_paralog_stats(
                arguments["data_folder"],
                arguments["minimum_divergence"],
                arguments["maximum_divergence"],
                paralog_statistic,
                arguments["probes_exons"],
                logger,
            )
    elif arguments["command"] == "cast_correct":
        logger.info(
            f"""ParalogWizard cast_collect running with the following settings
            main data folder - {arguments["data_folder"]}
            folder with contigs - {arguments["contig_folder"]}"""
        )
        run_blat(
            arguments["data_folder"],
            arguments["probes_paralogs"],
            arguments["min_identity"],
            logger,
        )
        correct(
            arguments["data_folder"],
            arguments["probes_paralogs"],
            arguments["redlist"],
            logger,
        )


if __name__ == "__main__":
    main()
