import sys
from typing import Set

from HybWizard_Functions import sort_hit_table_ident, exon, percent_dissimilarity


def create_reference():
    print('Creating new reference...')
    with open(path_to_data_HPM + '/exons/all_hits.txt') as all_hits:
        all_hits_for_reference = [x[:-1] for x in all_hits.readlines()]
    sort_hit_table_ident(all_hits_for_reference, 'exon')
    with open(path_to_data_HPM + '/exons/new_reference_for_HybPhyloMaker.fas', 'w') as new_reference:
        exons: Set[str] = set()
        paralog_written: bool = bool()
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


if __name__ == '__main__':
    path_to_data_HPM = sys.argv[1]
    blacklist: set = set([x.strip() for x in sys.argv[2].split(',')])
    paralogs_bool = sys.argv[3]
    paralog_min_divergence = float(sys.argv[4].strip())
    create_reference()
