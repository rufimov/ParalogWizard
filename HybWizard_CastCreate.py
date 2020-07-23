import sys
from typing import List, Set

from Bio import pairwise2


def sort_hit_table_ident(hittable_as_list: List[str], primary_field: str):
    hittable_as_list.sort(key=lambda x: float(x.split()[5]), reverse=True)
    hittable_as_list.sort(key=lambda x: float(x.split()[4]))
    hittable_as_list.sort(key=lambda x: float(x.split()[3]), reverse=True)
    hittable_as_list.sort(key=lambda x: float(x.split()[2]), reverse=True)
    hittable_as_list.sort(key=globals()[primary_field])


def exon(string: str) -> str:
    return string.split()[0].split('-')[1]


def percent_dissimilarity(seq1: str, seq2: str) -> float:
    return 100 - (pairwise2.align.globalxx(seq1, seq2)[0][2] / min(len(seq1), len(seq2))) * 100


path_to_data_HPM = sys.argv[1]
blacklist: set = set([x.strip() for x in sys.argv[2].split(',')])
paralogs_bool = sys.argv[3]
paralog_min_divergence = float(sys.argv[4].strip())

with open(path_to_data_HPM + '/exons/all_hits.txt') as all_hits:
    all_hits_for_reference = [x[:-1] for x in all_hits.readlines()]
sort_hit_table_ident(all_hits_for_reference, 'exon')
print('Creating new reference...')
with open(path_to_data_HPM + '/exons/new_reference_for_HybPhyloMaker.fas', 'w') as new_reference:
    exons: Set[str] = set()
    paralog_written: bool = bool()
    main_copy_in_case_of_no_par: bool = bool()
    current_best: str = ''
    count = 0
    for hit in all_hits_for_reference:
        if hit.split()[-2] not in blacklist and float(hit.split()[3]) >= 75:
            if paralogs_bool != 'yes':
                if exon(hit) not in exons:
                    name_of_locus = exon(hit).replace('exon', 'Contig').replace('Exon', 'Contig') \
                        .replace('contig', 'Contig').replace('_', '').replace('Contig', '_Contig_')
                    new_reference.write(
                        '>Assembly_' + name_of_locus + '_' + hit.split()[-2] + '\n' + hit.split()[-1] +
                        '\n')
                    exons.add(exon(hit))
            else:
                if exon(hit) not in exons:
                    if not paralog_written and count != 0:
                        name_of_locus = exon(current_best).replace('exon', 'Contig').replace('Exon', 'Contig') \
                            .replace('contig', 'Contig').replace('_', '').replace('Contig', '_Contig_')
                        new_reference.write('>Assembly_' + name_of_locus + '_' + current_best.split()[-2] + '\n'
                                            + current_best.split()[-1] + '\n')
                    current_locus = dict()
                    current_best = hit
                    current_locus[hit.split()[-2]] = hit
                    exons.add(exon(hit))
                    paralog_written: bool = False
                else:
                    if hit.split()[-2] in set(current_locus.keys()):
                        if not paralog_written:
                            current_seq: str = hit.split()[-1].replace('-', '')
                            seq_to_compare: str = current_locus[hit.split()[-2]].split()[-1].replace('-', '')
                            if percent_dissimilarity(current_seq, seq_to_compare) > \
                                    paralog_min_divergence:
                                print('Paralog detected for ' + hit.split()[0].split('-')[1])
                                name_of_locus_main: str = exon(current_locus[hit.split()[-2]]).replace('exon',
                                                                                                       'Contig'). \
                                    replace('Exon', 'Contig').replace('contig', 'Contig').replace('_', ''). \
                                    replace('Contig', '_Contig_')
                                new_reference.write(
                                    '>Assembly_' + name_of_locus_main + '_' + hit.split()[-2] + '\n' +
                                    seq_to_compare + '\n')
                                name_of_locus_para = exon(hit).replace('exon', 'Contig').replace('Exon', 'Contig') \
                                    .replace('contig', 'Contig').replace('_', '').replace('Contig', '_Contig_')
                                name_of_locus_para: str = '_'.join([name_of_locus_para.split('_')[0] + 'par2'] +
                                                                   name_of_locus_para.split('_')[1:])
                                new_reference.write(
                                    '>Assembly_' + name_of_locus_para + '_' + hit.split()[-2] + '\n' +
                                    hit.split()[-1] + '\n')
                                paralog_written: bool = True
                    else:
                        current_locus[hit.split()[-2]] = hit
        count += 1
print('New reference created!\n')
