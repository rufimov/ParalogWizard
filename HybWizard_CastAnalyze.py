import glob
import itertools
import os
import random
import sys
from multiprocessing.pool import ThreadPool
from typing import List, Dict

import Bio
import matplotlib.pyplot
import numpy
from Bio import pairwise2, SeqRecord, SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from sklearn.mixture import GaussianMixture, BayesianGaussianMixture

matplotlib.use('Agg')

path_to_data_HPM = sys.argv[1].strip()
blacklist: set = set([x.strip() for x in sys.argv[2].split(',')])


def exon(string: str) -> str:
    return string.split()[0].split('-')[1]


def contig(string: str) -> str:
    return string.split()[1]


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


def get_plot(name: str, matrix: numpy.ndarray, mix: GaussianMixture, lines: Dict[str, float], num_contig, mode):
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
    fig.savefig(path_to_data_HPM + '/exons/aln_orth_par/' + name + '.png', dpi=300)
    matplotlib.pyplot.close(fig)


def sort_hit_table_ident(hittable_as_list: List[str], primary_field: str):
    hittable_as_list.sort(key=lambda x: float(x.split()[5]), reverse=True)
    hittable_as_list.sort(key=lambda x: float(x.split()[4]))
    hittable_as_list.sort(key=lambda x: float(x.split()[3]), reverse=True)
    hittable_as_list.sort(key=lambda x: float(x.split()[2]), reverse=True)
    hittable_as_list.sort(key=globals()[primary_field])


with open(path_to_data_HPM + '/exons/all_hits.txt') as all_hits:
    all_hits_for_reference = [x[:-1] for x in all_hits.readlines()]
print('Building individual exon alignments and estimating their divergence...')
sort_hit_table_ident(all_hits_for_reference, 'exon')
exons: Dict[str, List[Bio.SeqRecord.SeqRecord]] = dict()
for hit in all_hits_for_reference:
    if hit.split()[-2] not in blacklist:
        if exon(hit) not in set(exons.keys()):
            exons[exon(hit)] = [SeqRecord(Seq(hit.split()[-1], generic_dna), id=contig(hit) + '_' + hit.split()[-2],
                                          description='')]
        else:
            new_list = exons[exon(hit)]
            new_list.append(SeqRecord(Seq(hit.split()[-1], generic_dna), id=contig(hit) + '_' + hit.split()[-2],
                                      description=''))
            exons[exon(hit)] = new_list
os.makedirs(path_to_data_HPM + '/exons/aln_orth_par', exist_ok=True)
for key in exons.keys():
    SeqIO.write(exons[key], path_to_data_HPM + '/exons/aln_orth_par/' + key + '.fasta', 'fasta')
divergency_distribution = []
pool = ThreadPool(12)
for file in sorted(glob.glob(path_to_data_HPM + '/exons/aln_orth_par/*.fasta')):
    pool.apply_async(get_distance_matrix, (file, divergency_distribution))
pool.close()
pool.join()
with open(path_to_data_HPM + '/exons/aln_orth_par/pairwise_distances.txt', 'w') as divergency_distribution_to_write:
    for i in divergency_distribution:
        divergency_distribution_to_write.write(str(i) + '\n')
divergency_distribution_array = numpy.array([[x] for x in divergency_distribution])
divergency_distribution_mix = get_best_model(divergency_distribution_array, 2)
means = divergency_distribution_mix.means_.flatten().tolist()
sigmas = [numpy.sqrt(x) for x in divergency_distribution_mix.covariances_.flatten().tolist()]
mu_sigma = sorted(list(zip(means, sigmas)), key=lambda x: x[0])
first_peak = mu_sigma[0][0]
first_peak_minus_sigma = first_peak - mu_sigma[0][1]
second_peak = mu_sigma[1][0]
second_peak_minus_sigma = second_peak - mu_sigma[1][1]
get_plot('pairwise_distances_distribution_2_comp', divergency_distribution_array, divergency_distribution_mix, {
    'first peak minus sigma ': numpy.round(first_peak_minus_sigma, 2),
    'first peak ': numpy.round(first_peak, 2),
    'second peak minus sigma': numpy.round(second_peak_minus_sigma, 2),
    'second peak ': numpy.round(second_peak, 2)}, '', 'individual')
divergency_distribution_mix = get_best_model(divergency_distribution_array, 3)
means = divergency_distribution_mix.means_.flatten().tolist()
sigmas = [numpy.sqrt(x) for x in divergency_distribution_mix.covariances_.flatten().tolist()]
mu_sigma = sorted(list(zip(means, sigmas)), key=lambda x: x[0])
first_peak = mu_sigma[0][0]
first_peak_minus_sigma = first_peak - mu_sigma[0][1]
second_peak = mu_sigma[1][0]
second_peak_minus_sigma = second_peak - mu_sigma[1][1]
third_peak = mu_sigma[2][0]
third_peak_minus_sigma = third_peak - mu_sigma[2][1]
get_plot('pairwise_distances_distribution_3_comp', divergency_distribution_array, divergency_distribution_mix, {
    'first peak minus sigma ': numpy.round(first_peak_minus_sigma, 2),
    'first peak ': numpy.round(first_peak, 2),
    'second peak minus sigma': numpy.round(second_peak_minus_sigma, 2),
    'second peak ': numpy.round(second_peak, 2),
    'third peak minus sigma ': numpy.round(third_peak_minus_sigma, 2),
    'third peak ': numpy.round(third_peak, 2)}, '', 'individual')
print('Done\n')
