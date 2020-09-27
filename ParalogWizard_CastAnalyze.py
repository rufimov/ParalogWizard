import glob
import os
import sys
from multiprocessing.pool import ThreadPool
from typing import List, Dict, Tuple, Set

import numpy
from Bio import SeqRecord, SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from sklearn.mixture import BayesianGaussianMixture

from ParalogWizard_Functions import sort_hit_table_ident, exon, contig, get_distance_matrix, get_model, get_plot


def build_alignments():
    print('Building individual exon alignments...')
    with open(path_to_data_HPM + '/exons/all_hits.txt') as all_hits:
        all_hits_for_reference: List[str] = [x[:-1] for x in all_hits.readlines()]
    sort_hit_table_ident(all_hits_for_reference, 'exon')
    exons: Dict[str, List[SeqRecord]] = dict()
    for hit in all_hits_for_reference:
        if hit.split()[-2] not in blacklist:
            if exon(hit) not in set(exons.keys()):
                exons[exon(hit)]: List[SeqRecord] = [SeqRecord(Seq(hit.split()[-1], generic_dna),
                                                               id=contig(hit) + '_' + hit.split()[-2],
                                                               description='')]
            else:
                new_list: List[SeqRecord] = exons[exon(hit)]
                new_list.append(SeqRecord(Seq(hit.split()[-1], generic_dna), id=contig(hit) + '_' + hit.split()[-2],
                                          description=''))
                exons[exon(hit)] = new_list
    os.makedirs(path_to_data_HPM + '/exons/aln_orth_par', exist_ok=True)
    for key in exons.keys():
        SeqIO.write(exons[key], path_to_data_HPM + '/exons/aln_orth_par/' + key + '.fasta', 'fasta')
    print('Done\n')


def estimate_divergence():
    print('Estimating divergence of paralogs...')
    divergency_distribution: List[float] = []
    pool: ThreadPool = ThreadPool(200)
    for file in sorted(glob.glob(path_to_data_HPM + '/exons/aln_orth_par/*.fasta')):
        pool.apply_async(get_distance_matrix, (file, divergency_distribution))
    pool.close()
    pool.join()
    with open(path_to_data_HPM + '/exons/aln_orth_par/pairwise_distances.txt', 'w') as divergency_distribution_to_write:
        for i in divergency_distribution:
            divergency_distribution_to_write.write(str(i) + '\n')
    divergency_distribution_array: numpy.ndarray = numpy.array([[x] for x in divergency_distribution])
    divergency_distribution_mix: BayesianGaussianMixture = get_model(divergency_distribution_array, 1)
    means: List[float] = divergency_distribution_mix.means_.flatten().tolist()
    sigmas: List[float] = [numpy.sqrt(x) for x in divergency_distribution_mix.covariances_.flatten().tolist()]
    mu_sigma: List[Tuple[float, float]] = sorted(list(zip(means, sigmas)), key=lambda x: x[0])
    first_peak: float = mu_sigma[0][0]
    first_peak_minus_sigma: float = first_peak - mu_sigma[0][1]
    get_plot(path_to_data_HPM + '/exons/aln_orth_par/', 'pairwise_distances_distribution_1_comp',
             divergency_distribution_array, divergency_distribution_mix, {
                 'first peak minus sigma ': numpy.round(first_peak_minus_sigma, 2),
                 'first peak ': numpy.round(first_peak, 2)}, '', 'individual')
    divergency_distribution_mix: BayesianGaussianMixture = get_model(divergency_distribution_array, 2)
    means: List[float] = divergency_distribution_mix.means_.flatten().tolist()
    sigmas: List[float] = [numpy.sqrt(x) for x in divergency_distribution_mix.covariances_.flatten().tolist()]
    mu_sigma: List[Tuple[float, float]] = sorted(list(zip(means, sigmas)), key=lambda x: x[0])
    first_peak: float = mu_sigma[0][0]
    first_peak_minus_sigma: float = first_peak - mu_sigma[0][1]
    second_peak: float = mu_sigma[1][0]
    second_peak_minus_sigma: float = second_peak - mu_sigma[1][1]
    get_plot(path_to_data_HPM + '/exons/aln_orth_par/', 'pairwise_distances_distribution_2_comp',
             divergency_distribution_array, divergency_distribution_mix, {
                 'first peak minus sigma ': numpy.round(first_peak_minus_sigma, 2),
                 'first peak ': numpy.round(first_peak, 2),
                 'second peak minus sigma': numpy.round(second_peak_minus_sigma, 2),
                 'second peak ': numpy.round(second_peak, 2)}, '', 'individual')
    divergency_distribution_mix: BayesianGaussianMixture = get_model(divergency_distribution_array, 3)
    means: List[float] = divergency_distribution_mix.means_.flatten().tolist()
    sigmas: List[float] = [numpy.sqrt(x) for x in divergency_distribution_mix.covariances_.flatten().tolist()]
    mu_sigma: List[Tuple[float, float]] = sorted(list(zip(means, sigmas)), key=lambda x: x[0])
    first_peak: float = mu_sigma[0][0]
    first_peak_minus_sigma: float = first_peak - mu_sigma[0][1]
    second_peak: float = mu_sigma[1][0]
    second_peak_minus_sigma: float = second_peak - mu_sigma[1][1]
    third_peak: float = mu_sigma[2][0]
    third_peak_minus_sigma: float = third_peak - mu_sigma[2][1]
    get_plot(path_to_data_HPM + '/exons/aln_orth_par/', 'pairwise_distances_distribution_3_comp',
             divergency_distribution_array, divergency_distribution_mix, {
                 'first peak minus sigma ': numpy.round(first_peak_minus_sigma, 2),
                 'first peak ': numpy.round(first_peak, 2),
                 'second peak minus sigma': numpy.round(second_peak_minus_sigma, 2),
                 'second peak ': numpy.round(second_peak, 2),
                 'third peak minus sigma ': numpy.round(third_peak_minus_sigma, 2),
                 'third peak ': numpy.round(third_peak, 2)}, '', 'individual')
    print('Done\n')


if __name__ == '__main__':
    path_to_data_HPM: str = sys.argv[1].strip()
    blacklist: Set[str] = set([x.strip() for x in sys.argv[2].split(',')])
    build_alignments()
    estimate_divergence()
