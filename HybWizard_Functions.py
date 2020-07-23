import itertools
import random
from typing import List, Dict

import Bio
import matplotlib.pyplot
import numpy
from Bio import pairwise2, SeqRecord, SeqIO
from Bio.Alphabet import generic_dna
from sklearn.mixture import GaussianMixture, BayesianGaussianMixture


def sort_hit_table_cover(hittable_as_list: List[str], primary_field: str):
    hittable_as_list.sort(key=lambda x: float(x.split()[5]), reverse=True)
    hittable_as_list.sort(key=lambda x: float(x.split()[4]))
    hittable_as_list.sort(key=lambda x: float(x.split()[2]), reverse=True)
    hittable_as_list.sort(key=lambda x: float(x.split()[3]), reverse=True)
    hittable_as_list.sort(key=globals()[primary_field])


def exon(string: str) -> str:
    return string.split()[0].split('-')[1]


def locus(string: str) -> str:
    return '_'.join(string.split()[0].split('-')[1].split('_')[:-2])


def contig(string: str) -> str:
    return string.split()[1]


def contig_locus(string: str) -> str:
    return string.split()[1].split('_N_')[0]


def slicing(dictionary: Dict[str, Bio.SeqRecord.SeqRecord], current_string: str, key_column: int, start_column: int,
            end_column: int, rev: str) -> str:
    if rev == 'no':
        return str(dictionary[current_string.split()[key_column]].seq)[int(current_string.split()[start_column]) - 1:
                                                                       int(current_string.split()[end_column])]
    elif rev == 'yes':
        return str(dictionary[current_string.split()[key_column]].seq.
                   reverse_complement())[-int(current_string.split()[start_column]):
                                         -int(current_string.split()[end_column]) - 1:-1][::-1]


def percent_dissimilarity(seq1: str, seq2: str) -> float:
    return 100 - (pairwise2.align.globalxx(seq1, seq2)[0][2] / min(len(seq1), len(seq2))) * 100


def get_distance_matrix(file_to_process, sum_list):
    current_distance_matrix = []
    with open(file_to_process) as fasta_file:
        sequences = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta', generic_dna))
        for pair in list(itertools.combinations(sequences.keys(), 2)):
            if '_'.join(pair[0].split('_')[-2:]) == '_'.join(pair[1].split('_')[-2:]):
                sequence1 = str(sequences[pair[0]].seq)
                sequence2 = str(sequences[pair[1]].seq)
                distance = percent_dissimilarity(sequence1, sequence2)
                current_distance_matrix.append(distance)
    sum_list.extend(current_distance_matrix)


def get_best_model(array, num_comp):
    return BayesianGaussianMixture(n_components=num_comp, max_iter=10000, n_init=10).fit(array)


def get_plot(path: str, name: str, matrix: numpy.ndarray, mix: GaussianMixture, lines: Dict[str, float],
             num_contig, mode):
    max_of_matrix = numpy.round(numpy.max(matrix)) + 1
    fig, axis = matplotlib.pyplot.subplots(figsize=(15, 10))
    axis.plot(matrix, [0] * matrix.shape[0], marker=2, color='k')
    colors = random.sample(['#FF0000', '#008000', '#0000FF', '#FF00FF', '#800000', '#808000', '#008080', '#2196F3',
                            '#4A148C', '#512DA8', '#2ECC71'], k=len(lines))
    count = 0
    for i in lines.keys():
        axis.axvline(x=lines[i], label=f"{i} - {numpy.round(lines[i], 2)}", c=colors[count], lw=0.6, ls='--')
        count += 1
    axis.legend(loc='upper right')
    axis.set_xlabel('Divergence (%)')
    axis.xaxis.set_ticks(numpy.arange(-1, max_of_matrix, 1))
    if num_contig == '':
        fig.suptitle(f'{name}')
    else:
        fig.suptitle(f'{name} ({num_contig} contigs)')
    x = numpy.linspace(-1, max_of_matrix, 1000)
    logprob = mix.score_samples(x.reshape(-1, 1))
    responsibilities = mix.predict_proba(x.reshape(-1, 1))
    pdf = numpy.exp(logprob)
    axis.hist(matrix, bins=numpy.arange(-1, max_of_matrix, 1),
              density=True, histtype='stepfilled', alpha=0.4)
    axis.plot(x, pdf, '-k')
    if mode == 'individual':
        pdf_individual = responsibilities * pdf[:, numpy.newaxis]
        axis.plot(x, pdf_individual, '--k')
    elif mode == 'mix':
        pass
    fig.savefig(path + name + '.png', dpi=300)
    matplotlib.pyplot.close(fig)


def sort_hit_table_ident(hittable_as_list: List[str], primary_field: str):
    hittable_as_list.sort(key=lambda x: float(x.split()[5]), reverse=True)
    hittable_as_list.sort(key=lambda x: float(x.split()[4]))
    hittable_as_list.sort(key=lambda x: float(x.split()[3]), reverse=True)
    hittable_as_list.sort(key=lambda x: float(x.split()[2]), reverse=True)
    hittable_as_list.sort(key=globals()[primary_field])
