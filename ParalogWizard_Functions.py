import itertools
import random
from typing import List, Dict

import Bio
import matplotlib.pyplot
import numpy
from Bio import SeqRecord, SeqIO
from Bio.Alphabet import generic_dna
from sklearn.mixture import BayesianGaussianMixture


def sort_hit_table_cover(hittable_as_list: List[str], primary_field: str):
    hittable_as_list.sort(
        key=lambda x: float(x.split()[1].split("_")[-1]), reverse=True
    )
    hittable_as_list.sort(key=lambda x: float(x.split()[5]), reverse=True)
    hittable_as_list.sort(key=lambda x: float(x.split()[4]))
    hittable_as_list.sort(key=lambda x: float(x.split()[2]), reverse=True)
    hittable_as_list.sort(key=lambda x: float(x.split()[3]), reverse=True)
    hittable_as_list.sort(key=globals()[primary_field])


def sort_hit_table_ident(hittable_as_list: List[str], primary_field: str):
    hittable_as_list.sort(
        key=lambda x: float(x.split()[1].split("_")[-1]), reverse=True
    )
    hittable_as_list.sort(key=lambda x: float(x.split()[5]), reverse=True)
    hittable_as_list.sort(key=lambda x: float(x.split()[4]))
    hittable_as_list.sort(key=lambda x: float(x.split()[3]), reverse=True)
    hittable_as_list.sort(key=lambda x: float(x.split()[2]), reverse=True)
    hittable_as_list.sort(key=globals()[primary_field])


def exon(string: str) -> str:
    return string.split()[0].split("-")[1]


def locus(string: str) -> str:
    return string.split()[0].split("-")[1].split("_")[0]


def contig(string: str) -> str:
    return string.split()[1]


def contig_locus(string: str) -> str:
    return string.split()[1].split("_N_")[0]


def slicing(
    dictionary: Dict[str, Bio.SeqRecord.SeqRecord],
    current_string: str,
    key_column: int,
    start_column: int,
    end_column: int,
    rev: str,
) -> str:
    start = int(current_string.split()[start_column])
    end = int(current_string.split()[end_column])
    sequence = dictionary[current_string.split()[key_column]].seq
    if rev == "yes":
        result = str(sequence.reverse_complement())[-start : -end - 1 : -1][::-1]
    else:
        result = str(sequence)[start - 1 : end]
    return result


def percent_dissimilarity(seq1: str, seq2: str) -> float:
    seq1 = seq1.lower()
    seq2 = seq2.lower()
    while seq1[0] == "-" or seq2[0] == "-":
        seq1 = seq1[1:]
        seq2 = seq2[1:]
    while seq1[-1] == "-" or seq2[-1] == "-":
        seq1 = seq1[:-1]
        seq2 = seq2[:-1]
    count = 0
    seq1_wo_mutual_gaps = str()
    seq2_wo_mutual_gaps = str()
    for i in range(0, len(seq1)):
        if seq1[i] == "-" and seq2[i] == "-":
            continue
        else:
            seq1_wo_mutual_gaps = seq1_wo_mutual_gaps + seq1[i]
            seq2_wo_mutual_gaps = seq2_wo_mutual_gaps + seq2[i]
    for i in range(0, len(seq1_wo_mutual_gaps)):
        if (
            seq1_wo_mutual_gaps[i] != seq2_wo_mutual_gaps[i]
            and seq1_wo_mutual_gaps[i] != "-"
            and seq2_wo_mutual_gaps[i] != "-"
        ):
            count += 1
    return (count / len(seq1)) * 100


def get_distance_matrix(
    file_to_process: str, sum_list: List[float], list_to_write: List[str]
):
    current_distance_matrix: List[float] = []
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
                current_distance_matrix.append(distance)
                current_list_to_write.append(f"{pair[0]}\t{str(distance)}\t{pair[1]}")
    list_to_write.extend(current_list_to_write)
    sum_list.extend(current_distance_matrix)


def get_model(array: numpy.ndarray, num_comp: int) -> BayesianGaussianMixture:
    return BayesianGaussianMixture(
        n_components=num_comp, max_iter=10000, n_init=10
    ).fit(array)


def get_plot(
    path: str,
    name: str,
    matrix: numpy.ndarray,
    mix: BayesianGaussianMixture,
    lines: Dict[str, float],
    num_contig,
    mode,
):
    matplotlib.use("Agg")
    max_of_matrix: float = numpy.round(numpy.max(matrix)) + 1
    fig, axis = matplotlib.pyplot.subplots(figsize=(15, 10))
    axis.plot(matrix, [0] * matrix.shape[0], marker=2, color="k")
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
        k=len(lines),
    )
    count: int = 0
    for i in lines.keys():
        axis.axvline(
            x=lines[i],
            label=f"{i} - {numpy.round(lines[i], 2)}",
            c=colors[count],
            lw=0.6,
            ls="--",
        )
        count += 1
    axis.legend(loc="upper right")
    axis.set_xlabel("Divergence (%)")
    axis.xaxis.set_ticks(numpy.arange(-1, max_of_matrix, 1))
    if num_contig == "":
        fig.suptitle(f"{name}")
    else:
        fig.suptitle(f"{name} ({num_contig} contigs)")
    x: numpy.ndarray = numpy.linspace(0, max_of_matrix, 1000)
    logprob: numpy.ndarray = mix.score_samples(x.reshape(-1, 1))
    responsibilities: numpy.ndarray = mix.predict_proba(x.reshape(-1, 1))
    pdf: numpy.ndarray = numpy.exp(logprob)
    axis.hist(
        matrix,
        bins=numpy.arange(0, max_of_matrix, 1),
        density=True,
        histtype="stepfilled",
        alpha=0.4,
    )
    axis.plot(x, pdf, "-k")
    if mode == "individual":
        pdf_individual: numpy.ndarray = responsibilities * pdf[:, numpy.newaxis]
        axis.plot(x, pdf_individual, "--k")
    elif mode == "mix":
        pass
    fig.savefig(path + name + ".png", dpi=300)
    matplotlib.pyplot.close(fig)
