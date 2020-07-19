import os
import glob
import sys
import shutil
import re
import Bio
from Bio import pairwise2, SeqRecord, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from typing import List, Dict, Set
from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastnCommandline
from Bio.Align.Applications import MafftCommandline
import matplotlib.pyplot
from scipy.stats import norm
import numpy
import itertools
from sklearn.mixture import GaussianMixture, BayesianGaussianMixture
from multiprocessing.pool import ThreadPool
import random

matplotlib.use('Agg')

path_to_data_HP: str = sys.argv[1].strip()
path_to_data_HPM: str = sys.argv[2].strip()
probe_HP_one_repr: str = sys.argv[3].strip()
length_cover: float = float(sys.argv[4].strip())
spades_cover: float = float(sys.argv[5].strip())
new_reference_bool: str = sys.argv[6].strip()
blacklist: set = set([x.strip() for x in sys.argv[7].split(',')])
paralogs_bool: str = sys.argv[8].strip()
paralog_min_divergence = float()
paralog_max_divergence = float()


def sort_hit_table_cover(hittable_as_list: List[str], primary_field: str):
    hittable_as_list.sort(key=lambda x: float(x.split()[5]), reverse=True)
    hittable_as_list.sort(key=lambda x: float(x.split()[4]))
    hittable_as_list.sort(key=lambda x: float(x.split()[2]), reverse=True)
    hittable_as_list.sort(key=lambda x: float(x.split()[3]), reverse=True)
    hittable_as_list.sort(key=globals()[primary_field])


def sort_hit_table_ident(hittable_as_list: List[str], primary_field: str):
    hittable_as_list.sort(key=lambda x: float(x.split()[5]), reverse=True)
    hittable_as_list.sort(key=lambda x: float(x.split()[4]))
    hittable_as_list.sort(key=lambda x: float(x.split()[3]), reverse=True)
    hittable_as_list.sort(key=lambda x: float(x.split()[2]), reverse=True)
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


def intersect_of_equal_length(str1, str2):
    while (str1.startswith('-') or str2.startswith('-')) and str1[0] != str2[0]:
        str1 = str1[1:]
        str2 = str2[1:]
    while (str1.endswith('-') or str2.endswith('-')) and str1[0] != str2[0]:
        str1 = str1[:-1]
        str2 = str2[:-1]
    return len(str1)


def percent_dissimilarity(seq1: str, seq2: str) -> float:
    similarity = float()
    lenmax = max(len(seq1), len(seq2))
    for i in range(lenmax):
        if seq1[i] == seq2[i]:
            similarity += 1
    distance = 100 - ((similarity / intersect_of_equal_length(seq1, seq2)) * 100)
    return distance


def percent_dissimilarity_pairwise2(seq1: str, seq2: str) -> float:
    return 100 - (pairwise2.align.globalxx(seq1, seq2)[0][2] / min(len(seq1), len(seq2))) * 100


def get_best_model(array, num_comp):
    # N = numpy.arange(1, num_comp + 1)
    # models = [None for i in range(len(N))]
    # for i in range(len(N)):
    #     models[i] = GaussianMixture(n_components=N[i], max_iter=10000, n_init=10).fit(array)
    # AIC = [m.aic(array) for m in models]
    # BIC = [m.bic(array) for m in models]
    # M_best = models[numpy.argmin(AIC)]
    # return M_best
    return BayesianGaussianMixture(n_components=num_comp, max_iter=10000, n_init=10).fit(array)


def mix_pdf(x, loc, scale, weights):
    d = numpy.zeros_like(x)
    for mu, sigma, pi in zip(loc, scale, weights):
        d += pi * norm.pdf(x, loc=mu, scale=sigma)
    return d


def individual_pdf(x, loc, scale, weights):
    d = numpy.zeros_like(x)
    d += weights * norm.pdf(x, loc=loc, scale=scale)
    return d


def get_plot(name: str, matrix: numpy.ndarray, mix: GaussianMixture, lines: Dict[str, float], num_contig, mode):
    max_of_matrix = numpy.round(numpy.max(matrix)) + 1
    fig, axis = matplotlib.pyplot.subplots(figsize=(15, 10))
    axis.plot(matrix, [0] * matrix.shape[0], marker=2, color='k')
    colors = random.sample(['r', 'g', 'b', 'c', 'm'], k=len(lines))
    count = 0
    for i in lines.keys():
        axis.axvline(x=lines[i], label=f"{i} - {numpy.round(lines[i], 2)}", c=colors[count], lw=0.5, ls='--')
        count += 1
    axis.legend(loc='upper right')
    axis.set_xlabel('Divergence (%)')
    if num_contig == '':
        fig.suptitle(f'{name}')
    else:
        fig.suptitle(f'{name} ({num_contig} contigs)')
    # pi, mu, sigma = mix.weights_.flatten(), mix.means_.flatten(), numpy.sqrt(mix.covariances_.flatten())
    # grid = numpy.arange(-1, max_of_matrix, 0.01)
    # axis.hist(matrix, bins=numpy.arange(-1, max_of_matrix, 1), density=True, histtype='stepfilled', alpha=0.4)
    # axis.plot(grid, mix_pdf(grid, mu, sigma, pi), '-k')
    # for loc, scale, weights in zip(mu, sigma, pi):
    #     axis.plot(grid, individual_pdf(grid, loc, scale, weights), '--k')
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


def get_distance_matrix(file_to_process, max_num_comp):
    mafft_cline = MafftCommandline(input=file_to_process, adjustdirection=True)
    stdout, stderr = mafft_cline()
    with open(file_to_process + '.mafft.fas', "w") as aligned_to_write:
        aligned_to_write.write(stdout)
    os.system('trimal -automated1 -in ' + file_to_process + '.mafft.fas -out ' + file_to_process +
              '.mafft.trim.fas -fasta')
    file_to_process = file_to_process + '.mafft.trim.fas'
    exon_name = file_to_process.split('/')[-1].split('.')[0]
    current_distance_matrix = []
    with open(file_to_process) as fasta_file:
        sequences = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta', generic_dna))
        num_of_contigs = len(sequences)
        for pair in list(itertools.combinations(sequences.keys(), 2)):
            
            sequence1 = str(sequences[pair[0]].seq)
            sequence2 = str(sequences[pair[1]].seq)
            distance = percent_dissimilarity(sequence1, sequence2)
            current_distance_matrix.append(distance)
    with open(file_to_process + '.dist.txt', 'w') as dist_to_write:
        for i in current_distance_matrix:
            dist_to_write.write(str(i) + '\n')
    if len(current_distance_matrix) < max_num_comp:
        max_num_comp = len(current_distance_matrix)
    if 1 < len(current_distance_matrix) != current_distance_matrix. \
            count(current_distance_matrix[0]):
        current_distance_matrix = list([[x] for x in current_distance_matrix])
        current_distance_matrix = numpy.array(current_distance_matrix)
        distribution_mix = get_best_model(current_distance_matrix, max_num_comp)
        means = zip(distribution_mix.means_.flatten().tolist(), distribution_mix.weights_.flatten().tolist())
        max_weight = distribution_mix.weights_.flatten().tolist()[0]
        means = sorted([x[0] for x in means if max_weight/x[1] < 10])
        first_peak = means[0]
        if len(means) == 1:
            get_plot(exon_name, current_distance_matrix, distribution_mix,
                     {'single peak of pairwise distance distribution': numpy.round(first_peak, 2)},
                     num_of_contigs, 'mix')
        elif len(means) == 2:
            second_peak = means[1]
            if numpy.round(first_peak, 2) == numpy.round(second_peak, 2):
                get_plot(exon_name, current_distance_matrix, distribution_mix,
                         {'single peak of pairwise distance distribution': numpy.round(first_peak, 2)},
                         num_of_contigs, 'mix')
            else:
                second_peaks.extend(means[1:])
                get_plot(exon_name, current_distance_matrix, distribution_mix,
                         {'first peak of pairwise distance distribution': numpy.round(first_peak, 2),
                          'second peak of pairwise distance': numpy.round(second_peak, 2)}, num_of_contigs,
                         'individual')
        else:
            second_peak = means[1]
            third_peak = means[2]
            if numpy.round(first_peak, 2) == numpy.round(second_peak, 2) == numpy.round(third_peak, 2):
                get_plot(exon_name, current_distance_matrix, distribution_mix,
                         {'single peak of pairwise distance distribution': numpy.round(first_peak, 2)},
                         num_of_contigs, 'mix')
            elif numpy.round(first_peak, 2) == numpy.round(second_peak, 2) != numpy.round(third_peak, 2):
                second_peaks.extend(means[2:])
                get_plot(exon_name, current_distance_matrix, distribution_mix,
                         {'first peak of pairwise distance distribution': numpy.round(first_peak, 2),
                          'second peak of pairwise distance': numpy.round(third_peak, 2)}, num_of_contigs,
                         'mix')
            elif numpy.round(first_peak, 2) != numpy.round(second_peak, 2) == numpy.round(third_peak, 2):
                second_peaks.extend(means[1:2])
                get_plot(exon_name, current_distance_matrix, distribution_mix,
                         {'first peak of pairwise distance distribution': numpy.round(first_peak, 2),
                          'second peak of pairwise distance': numpy.round(second_peak, 2)}, num_of_contigs,
                         'mix')
            elif numpy.round(first_peak, 2) != numpy.round(second_peak, 2) != numpy.round(third_peak, 2):
                second_peaks.extend(means[1:])
                get_plot(exon_name, current_distance_matrix, distribution_mix,
                         {'first peak of pairwise distance distribution': numpy.round(first_peak, 2),
                          'second peak of pairwise distance': numpy.round(second_peak, 2),
                          'third peak of pairwise distance': numpy.round(third_peak, 2)}, num_of_contigs,
                         'individual')


# print('Converting data...\n')
# print('**********************************************************************************************************')
# print('\n')
# main_path: str = path_to_data_HPM + '/exons/40contigs/'
# if not os.path.exists(main_path):
#     os.makedirs(main_path, exist_ok=True)
# for file in glob.glob(path_to_data_HP + '/*contigs.fasta'):
#     shutil.copy(file, main_path + file.split('/')[-1])
# for file in glob.glob(main_path + '*.fasta'):
#     with open(file, 'r') as fasta:
#         lines: List[str] = fasta.readlines()
#     with open(file, 'w') as fasta:
#         for line in lines:
#             fasta.write(re.sub(r'length_([0-9]+)_cov_([0-9]+\.[0-9][0-9]).*', r'\1_c_\2', line.replace('NODE', 'N')))
#     name_of_file: str = '_'.join(file.split('/')[-1].split('.')[0].split('_')[0:2])
#     path_to_file: str = '/'.join(file.split('/')[:-1])
#     os.rename(file, path_to_file + '/' + name_of_file + '.fasta')
# print('Creating hit tables...')
# for file in glob.glob(main_path + '*.fasta'):
#     file: str = file.split('/')[-1]
#     sample: str = file[:-6]
#     print('\tProcessing ' + sample)
#     NcbimakeblastdbCommandline(dbtype='nucl', input_file=main_path + file,
#                                out=main_path + sample, parse_seqids=True)()
#     print('\tRunning BLAST...')
#     NcbiblastnCommandline(task='blastn', query=probe_HP_one_repr, db=main_path + sample,
#                           out=main_path + 'reference_in_' + sample + '_contigs.txt', qcov_hsp_perc=length_cover,
#                           num_threads=4, outfmt='6 qaccver saccver pident qcovhsp evalue bitscore sstart send')()
#     print('\tOK')
# print('Done\n')
# print('Correcting contigs...')
# statistics: Dict[str, dict] = dict()
# all_hits_for_reference: List[str] = list()
# for file in glob.glob(main_path + '*.fasta'):
#     sample: str = file.split('/')[-1][:-6]
#     print(' Processing ' + sample)
#     statistics[sample]: Dict[str, int] = dict()
#     hits: List[str] = list()
#     with open(main_path + 'reference_in_' + sample + '_contigs.txt') \
#             as blast_results, \
#             open(main_path + sample + '.fas', 'w') as result_fasta, \
#             open(file) as contigs:
#         contigs_fasta_parsed = SeqIO.to_dict(SeqIO.parse(contigs, 'fasta', generic_dna))
#         for line in blast_results.readlines():
#             if contig_locus(line) == locus(line) and float(line.split()[1].split('_c_')[1]) >= spades_cover:
#                 hits.append(line)
#         sort_hit_table_cover(hits, 'exon')
#         hits_exons_contigs: Set[str] = set()
#         hits_dedup: List[str] = list()
#         for hit in hits:
#             if str(exon(hit) + ' ' + contig(hit)) not in hits_exons_contigs:
#                 hits_dedup.append(hit)
#             else:
#                 pass
#             hits_exons_contigs.add(exon(hit) + ' ' + contig(hit))
#         for hit_dedup in hits_dedup:
#             if int(hit_dedup.split()[6]) > int(hit_dedup.split()[7]):
#                 sequence: str = slicing(contigs_fasta_parsed, hit_dedup, 1, 7, 6, 'yes')
#                 result_fasta.write('>' + exon(hit_dedup) + '_' + '_'.join(contig(hit_dedup).split('_')[1:]) + '\n' +
#                                    sequence + '\n')
#                 all_hits_for_reference.append('{0}\t{1}\t{2}'.format(hit_dedup, sample, sequence))
#             else:
#                 sequence: str = slicing(contigs_fasta_parsed, hit_dedup, 1, 6, 7, 'no')
#                 result_fasta.write('>' + exon(hit_dedup) + '_' + '_'.join(contig(hit_dedup).split('_')[1:]) + '\n' +
#                                    sequence + '\n')
#                 all_hits_for_reference.append('{0}\t{1}\t{2}'.format(hit_dedup, sample, sequence))
#         sort_hit_table_cover(hits_dedup, 'locus')
#         hits_loci: Set[str] = set()
#         for hit_dedup in hits_dedup:
#             if locus(hit_dedup) not in hits_loci:
#                 statistics[sample][locus(hit_dedup)] = {contig(hit_dedup)}
#             else:
#                 new_set = statistics[sample][locus(hit_dedup)]
#                 new_set.add(contig(hit_dedup))
#                 statistics[sample][locus(hit_dedup)] = new_set
#             hits_loci.add(locus(hit_dedup))
#     with open(main_path + 'reference_against_' + sample + '_contigs.txt', 'w') as \
#             hittable:
#         sort_hit_table_cover(hits, 'locus')
#         for hit in hits:
#             hittable.write(hit)
# print(' OK')
# print('All contigs were successfully corrected!\n')
# print('Writing statistics...')
# with open(main_path + 'statistics.csv', 'w') as stats, open(probe_HP_one_repr) as \
#         reference:
#     stats_dict: Dict[str, str] = dict([('gene\t', '')])
#     loci: Set[str] = set()
#     samples: List[str] = list()
#     for line in reference.read().splitlines():
#         if line.startswith('>'):
#             loci.add('_'.join(line.split('-')[1].split('_')[:-2]))
#     for key in statistics.keys():
#         samples.append(key)
#     samples.sort()
#     for sample in samples:
#         stats_dict['gene\t']: str = stats_dict['gene\t'] + sample + '\t'
#     for loc in loci:
#         stats_dict[loc + '\t']: str = ''
#         for sample in samples:
#             loci_in_sample: Set[str] = set(statistics[sample].keys())
#             if loc in loci_in_sample:
#                 stats_dict[loc + '\t']: str = stats_dict[loc + '\t'] + str(len(statistics[sample][loc])) + '\t'
#             else:
#                 stats_dict[loc + '\t']: str = stats_dict[loc + '\t'] + 'NA' + '\t'
#     stats.write('gene\t' + stats_dict['gene\t'] + '\n')
#     del stats_dict['gene\t']
#     for key in sorted(list(stats_dict.keys())):
#         stats.write(key + stats_dict[key] + '\n')
# print('Statistics file created!\n')
if new_reference_bool == 'yes':
    if paralogs_bool == 'yes':
#         print('Building individual exon alignments and estimating their divergence...')
#         sort_hit_table_ident(all_hits_for_reference, 'exon')
#         exons: Dict[str, List[Bio.SeqRecord.SeqRecord]] = dict()
#         for hit in all_hits_for_reference:
#             if hit.split()[-2] not in blacklist:
#                 if exon(hit) not in set(exons.keys()):
#                     exons[exon(hit)] = [
#                         SeqRecord(Seq(hit.split()[-1], generic_dna), id=contig(hit) + '_' + hit.split()[-2],
#                                   description='')]
#                 else:
#                     new_list = exons[exon(hit)]
#                     new_list.append(SeqRecord(Seq(hit.split()[-1], generic_dna), id=contig(hit) + '_' + hit.split()[-2],
#                                               description=''))
#                     exons[exon(hit)] = new_list
#         os.makedirs(path_to_data_HPM + '/exons/aln_orth_par', exist_ok=True)
#         for key in exons.keys():
#             SeqIO.write(exons[key], path_to_data_HPM + '/exons/aln_orth_par/' + key + '.fasta', 'fasta')
        second_peaks = []
        pool = ThreadPool(5)
        for file in sorted(glob.glob(path_to_data_HPM + '/exons/aln_orth_par/*.fasta')):
            pool.apply_async(get_distance_matrix, (file, 3))
        pool.close()
        pool.join()
        second_peaks_array = numpy.array([[x] for x in second_peaks])
        second_peaks_mix = get_best_model(second_peaks_array, 3)
        means = sorted(second_peaks_mix.means_.flatten().tolist())
        first_peak = means[0]
        second_peak = means[1]
        third_peak = means[2]
        get_plot('histogram_second_peaks', second_peaks_array, second_peaks_mix, {
            'first peak of second peaks\' distribution': numpy.round(first_peak, 2),
            'second peak of second peaks\' distribution': numpy.round(second_peak, 2),
            'third peak of second peaks\' distribution': numpy.round(third_peak, 2)}, '', 'individual')
        paralog_min_divergence = second_peak
        paralog_max_divergence = third_peak
        with open(path_to_data_HPM + '/exons/aln_orth_par/paralog_divergence.txt', 'w') as paralog_divergence:
            paralog_divergence.write(str(paralog_min_divergence) + '\n')
            paralog_divergence.write(str(paralog_max_divergence) + '\n')
        print('Done\n')
#     print('Creating new reference...')
#     sort_hit_table_ident(all_hits_for_reference, 'exon')
#     with open(path_to_data_HPM + '/exons/all_hits.txt', 'w') as all_hits:
#         for hit in all_hits_for_reference:
#             all_hits.write(hit + '\n')
#     with open(path_to_data_HPM + '/exons/new_reference_for_HybPhyloMaker.fas', 'w') as new_reference:
#         exons: Set[str] = set()
#         paralog_written: bool = bool()
#         main_copy_in_case_of_no_par: bool = bool()
#         current_best: str = ''
#         count = 0
#         for hit in all_hits_for_reference:
#             if hit.split()[-2] not in blacklist and float(hit.split()[3]) >= 75:
#                 if paralogs_bool != 'yes':
#                     if exon(hit) not in exons:
#                         name_of_locus = exon(hit).replace('exon', 'Contig').replace('Exon', 'Contig') \
#                             .replace('contig', 'Contig').replace('_', '').replace('Contig', '_Contig_')
#                         new_reference.write(
#                             '>Assembly_' + name_of_locus + '_' + hit.split()[-2] + '\n' + hit.split()[-1] +
#                             '\n')
#                         exons.add(exon(hit))
#                 else:
#                     if exon(hit) not in exons:
#                         if not paralog_written and count != 0:
#                             name_of_locus = exon(current_best).replace('exon', 'Contig').replace('Exon', 'Contig') \
#                                 .replace('contig', 'Contig').replace('_', '').replace('Contig', '_Contig_')
#                             new_reference.write('>Assembly_' + name_of_locus + '_' + current_best.split()[-2] + '\n'
#                                                 + current_best.split()[-1] + '\n')
#                         current_locus = dict()
#                         current_best = hit
#                         current_locus[hit.split()[-2]] = hit
#                         exons.add(exon(hit))
#                         paralog_written: bool = False
#                     else:
#                         if hit.split()[-2] in set(current_locus.keys()):
#                             if not paralog_written:
#                                 current_seq: str = hit.split()[-1].replace('-', '')
#                                 seq_to_compare: str = current_locus[hit.split()[-2]].split()[-1].replace('-',
#                                                                                                          '')
#                                 if paralog_max_divergence > percent_dissimilarity_pairwise2(current_seq,
#                                                                                             seq_to_compare) > \
#                                         paralog_min_divergence:
#                                     print('Paralog detected for ' + hit.split()[0].split('-')[1])
#                                     name_of_locus_main: str = exon(current_locus[hit.split()[-2]]).replace('exon',
#                                                                                                            'Contig'). \
#                                         replace('Exon', 'Contig').replace('contig', 'Contig').replace('_', ''). \
#                                         replace('Contig', '_Contig_')
#                                     new_reference.write(
#                                         '>Assembly_' + name_of_locus_main + '_' + hit.split()[-2] + '\n' +
#                                         seq_to_compare + '\n')
#                                     name_of_locus_para = exon(hit).replace('exon', 'Contig').replace('Exon', 'Contig') \
#                                         .replace('contig', 'Contig').replace('_', '').replace('Contig', '_Contig_')
#                                     name_of_locus_para: str = '_'.join([name_of_locus_para.split('_')[0] + 'par2'] +
#                                                                        name_of_locus_para.split('_')[1:])
#                                     new_reference.write(
#                                         '>Assembly_' + name_of_locus_para + '_' + hit.split()[-2] + '\n' +
#                                         hit.split()[-1] + '\n')
#                                     paralog_written: bool = True
#                         else:
#                             current_locus[hit.split()[-2]] = hit
#             count += 1
#     print('New reference created!\n')
# print('Renaming contigs...')
# for file in glob.glob(main_path + '*.fas'):
#     sample = file.split('/')[-1][:-4]
#     print(' Processing ' + sample)
#     with open(file) as result_fasta:
#         fasta_parsed = SeqIO.to_dict(SeqIO.parse(result_fasta, 'fasta', generic_dna))
#         counter: int = 1
#         fasta_to_write: List[str] = list()
#         for line in sorted(fasta_parsed.keys()):
#             fasta_to_write.append('>Contig' + str(counter) + '_' + sample + '-' + line.replace('_', '-')
#                                   + '\n')
#             fasta_to_write.append(str(fasta_parsed[line].seq) + '\n')
#             counter += 1
#     with open(file, 'w') as result_fasta:
#         result_fasta.writelines(fasta_to_write)
#     print(' OK')
# print('All contigs were successfully renamed!\n')
# print('Removing temporary files...')
# for file in glob.glob(main_path + '*.fasta'):
#     os.remove(file)
# for file in glob.glob(main_path + 'reference_in*'):
#     os.remove(file)
# for file in glob.glob(main_path + '*.n*'):
#     os.remove(file)
# print('Done\n')
# print('**********************************************************************************************************')
# print('\nData was successfully converted!')
