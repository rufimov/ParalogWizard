import os
import glob
import sys
import shutil
import re
import Bio
from Bio import pairwise2, SeqRecord, SeqIO
from Bio.Alphabet import generic_dna
from typing import List, Dict, Set
from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastnCommandline

path_to_data_HP: str = sys.argv[1].strip()
path_to_data_HPM: str = sys.argv[2].strip()
probe_HP_one_repr: str = sys.argv[3].strip()
length_cover: int = int(sys.argv[4].strip())
spades_cover: float = float(sys.argv[5].strip())
new_reference_bool: str = sys.argv[6].strip()
blacklist: set = set([x.strip() for x in sys.argv[7].split(',')])
paralogs_bool: str = sys.argv[8].strip()
paralog_min_divergence: int = int(sys.argv[9].strip())
paralog_max_divergence: int = int(sys.argv[10].strip())


def sort_hit_table_cover(hittable_as_list: List[str], primary_field: str):
    hittable_as_list.sort(key=lambda x: float(x.split()[5]), reverse=True)
    hittable_as_list.sort(key=lambda x: float(x.split()[4]))
    hittable_as_list.sort(key=lambda x: float(x.split()[2]), reverse=True)
    hittable_as_list.sort(key=lambda x: float(x.split()[3]), reverse=True)
    hittable_as_list.sort(key=globals()[primary_field])


def locus(string: str) -> str:
    return string.split()[0].split('-')[1]


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


print('Converting data...\n')
print('**********************************************************************************************************')
print('\n')
main_path: str = path_to_data_HPM + '/exons/40contigs/'
if not os.path.exists(main_path):
    os.makedirs(main_path)
for file in glob.glob(path_to_data_HP + '/*contigs.fasta'):
    shutil.copy(file, main_path + file.split('/')[-1])
for file in glob.glob(main_path + '*.fasta'):
    with open(file, 'r') as fasta:
        lines: List[str] = fasta.readlines()
    with open(file, 'w') as fasta:
        for line in lines:
            fasta.write(re.sub(r'length_([0-9]+)_cov_([0-9]+\.[0-9][0-9]).*', r'\1_c_\2', line.replace('NODE', 'N')))
    name_of_file: str = '_'.join(file.split('/')[-1].split('.')[0].split('_')[0:2])
    path_to_file: str = '/'.join(file.split('/')[:-1])
    os.rename(file, path_to_file + '/' + name_of_file + '.fasta')
print('Creating hit tables...')
for file in glob.glob(main_path + '*.fasta'):
    file: str = file.split('/')[-1]
    sample: str = file[:-6]
    print('\tProcessing ' + sample)
    NcbimakeblastdbCommandline(dbtype='nucl', input_file=main_path + file,
                               out=main_path + sample, parse_seqids=True)()
    print('\tRunning BLAST...')
    NcbiblastnCommandline(task='blastn', query=probe_HP_one_repr, db=main_path + sample,
                          out=main_path + 'reference_in_' + sample + '_contigs.txt',
                          outfmt='6 qaccver saccver pident qcovhsp evalue bitscore sstart send')()
    print('\tOK')
print('Done\n')
print('Correcting contigs...')
statistics: Dict[str, dict] = dict()
all_hits_for_reference: List[str] = list()
for file in glob.glob(main_path + '*.fasta'):
    sample: str = file.split('/')[-1][:-6]
    print(' Processing ' + sample)
    statistics[sample]: Dict[str, int] = dict()
    hits: List[str] = list()
    with open(main_path + 'reference_in_' + sample + '_contigs.txt') \
            as blast_results, \
            open(main_path + sample + '.fas', 'w') as result_fasta, \
            open(file) as contigs:
        contigs_fasta_parsed = SeqIO.to_dict(SeqIO.parse(contigs, 'fasta', generic_dna))
        for line in blast_results.readlines():
            if contig_locus(line) == locus(line) and int(line.split()[3]) >= \
                    length_cover and float(line.split()[1].split('_c_')[1]) >= spades_cover:
                hits.append(line)
        sort_hit_table_cover(hits, 'locus')
        hits_loci_contigs: Set[str] = set()
        hits_dedup: List[str] = list()
        for hit in hits:
            if str(locus(hit) + ' ' + contig(hit)) not in hits_loci_contigs:
                hits_dedup.append(hit)
            else:
                pass
            hits_loci_contigs.add(locus(hit) + ' ' + contig(hit))
        sort_hit_table_cover(hits_dedup, 'contig')
        hits_contigs: Set[str] = set()
        for hit_dedup in hits_dedup:
            if contig(hit_dedup) not in hits_contigs:
                if int(hit_dedup.split()[6]) > int(hit_dedup.split()[7]):
                    sequence: str = slicing(contigs_fasta_parsed, hit_dedup, 1, 7, 6, 'yes')
                    result_fasta.write('>' + contig(hit_dedup) + '\n' + sequence + '\n')
                    all_hits_for_reference.append('{0}\t{1}\t{2}'.format(hit_dedup, sample, sequence))
                else:
                    sequence: str = slicing(contigs_fasta_parsed, hit_dedup, 1, 6, 7, 'no')
                    result_fasta.write('>' + hit_dedup.split()[1] + '\n' + sequence + '\n')
                    all_hits_for_reference.append('{0}\t{1}\t{2}'.format(hit_dedup, sample, sequence))
                hits_contigs.add(contig(hit_dedup))
            else:
                pass
        sort_hit_table_cover(hits_dedup, 'locus')
        hits_loci: Set[str] = set()
        for hit_dedup in hits_dedup:
            if locus(hit_dedup) not in hits_loci:
                statistics[sample][locus(hit_dedup)] = 1
            else:
                statistics[sample][locus(hit_dedup)] += 1
            hits_loci.add(locus(hit_dedup))
    with open(main_path + 'reference_against_' + sample + '_contigs.txt', 'w') as \
            hittable:
        for hit in hits:
            hittable.write(hit)
    # with open(main_path + 'reference_against_' + sample + '_contigs_dedup.txt',
    #           'w') as \
    #         hittable_dedup:
    #     for hit_dedup in hits_dedup:
    #         hittable_dedup.write(hit_dedup)
print(' OK')
print('All contigs were successfully corrected!\n')
print('Writing statistics...')
with open(main_path + 'statistics.csv', 'w') as stats, open(probe_HP_one_repr) as \
        reference:
    stats_dict: Dict[str, str] = dict([('gene\t', '')])
    loci: Set[str] = set()
    samples: List[str] = list()
    for line in reference.read().splitlines():
        if line.startswith('>'):
            loci.add(line.split('-')[1])
    for key in statistics.keys():
        samples.append(key)
    samples.sort()
    for sample in samples:
        stats_dict['gene\t']: str = stats_dict['gene\t'] + sample + '\t'
    for loc in loci:
        stats_dict[loc + '\t']: str = ''
        for sample in samples:
            loci_in_sample: Set[str] = set(statistics[sample].keys())
            if loc in loci_in_sample:
                stats_dict[loc + '\t']: str = stats_dict[loc + '\t'] + str(statistics[sample][loc]) + '\t'
            else:
                stats_dict[loc + '\t']: str = stats_dict[loc + '\t'] + 'NA' + '\t'
    stats.write('gene\t' + stats_dict['gene\t'] + '\n')
    del stats_dict['gene\t']
    for key in sorted(list(stats_dict.keys())):
        stats.write(key + stats_dict[key] + '\n')
print('Statistics file created!\n')
if new_reference_bool == 'yes':
    print('Creating new reference...')
    sort_hit_table_cover(all_hits_for_reference, 'locus')
    with open(path_to_data_HPM + '/exons/new_reference_for_HybPhyloMaker.fas', 'w') as new_reference:
        num_paralog: int = 0
        cover_best_seq: int = 0
        best_seq: str = ''
        exons: Set[str] = set()
        for hit in all_hits_for_reference:
            if hit.split()[-2] not in blacklist:
                if locus(hit) not in exons:
                    name_of_locus = locus(hit).replace('exon', 'Contig').replace('Exon', 'Contig') \
                        .replace('contig', 'Contig').replace('_', '').replace('Contig', '_Contig_')
                    new_reference.write('>Assembly_' + name_of_locus + '_' + hit.split()[-2] + '\n' + hit.split()[-1] +
                                        '\n')
                    best_seq: str = hit.split()[-1]
                    exons.add(locus(hit))
                    num_paralog: int = 1
                else:
                    if paralogs_bool == 'yes' and num_paralog == 1:
                        current_seq: str = hit.split()[-1]
                        if paralog_max_divergence > percent_dissimilarity(current_seq, best_seq) > \
                                paralog_min_divergence:
                            print('Paralog detected for ' + hit.split()[0].split('-')[1])
                            name_of_locus: str = locus(hit).replace('exon', 'Contig'). \
                                replace('Exon', 'Contig').replace('contig', 'Contig').replace('_', ''). \
                                replace('Contig', '_Contig_')
                            name_of_locus: str = '_'.join([name_of_locus.split('_')[0] + 'par2'] +
                                                          name_of_locus.split('_')[1:])
                            new_reference.write('>Assembly_' + name_of_locus + '_' + hit.split()[-2] + '\n' +
                                                hit.split()[-1] + '\n')
                            num_paralog += 1
    print('New reference created!\n')
print('Renaming contigs...')
for file in glob.glob(main_path + '*.fas'):
    sample = file.split('/')[-1][:-4]
    print(' Processing ' + sample)
    with open(file) as result_fasta:
        fasta_as_list: List[str] = result_fasta.read().splitlines()
        fasta_parsed: Dict[str, str] = dict()
        for i in range(0, len(fasta_as_list), 2):
            fasta_parsed[fasta_as_list[i]] = fasta_as_list[i + 1]
        counter: int = 1
        fasta_to_write: List[str] = list()
        fasta_parsed_as_list: List[str] = list(fasta_parsed.keys())
        fasta_as_list.sort()
        for line in fasta_parsed_as_list:
            fasta_to_write.append('>Contig' + str(counter) + '_' + sample + '-' + line[1:].replace('_', '-')
                                  + '\n')
            fasta_to_write.append(fasta_parsed[line] + '\n')
            counter += 1
    with open(file, 'w') as result_fasta:
        result_fasta.writelines(fasta_to_write)
    print(' OK')
print('All contigs were successfully renamed!\n')
print('Removing temporary files...')
for file in glob.glob(main_path + '*.fasta'):
    os.remove(file)
for file in glob.glob(main_path + 'reference_in*'):
    os.remove(file)
for file in glob.glob(main_path + '*.n*'):
    os.remove(file)
print('Done\n')
print('**********************************************************************************************************')
print('\nData was successfully converted!')
