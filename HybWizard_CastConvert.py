import glob
import os
import re
import shutil
import sys
from typing import List, Dict, Set, Union

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastnCommandline

from HybWizard_Functions import sort_hit_table_cover, exon, locus, contig, contig_locus, slicing


def prepare_contigs():
    print('Preparing congits...')
    for file in glob.glob(path_to_data_HP + '/*contigs.fasta'):
        shutil.copy(file, main_path + file.split('/')[-1])
    for file in glob.glob(main_path + '*.fasta'):
        with open(file, 'r') as fasta:
            lines: List[str] = fasta.readlines()
        with open(file, 'w') as fasta:
            for line in lines:
                fasta.write(
                    re.sub(r'length_([0-9]+)_cov_([0-9]+\.[0-9][0-9]).*', r'\1_c_\2', line.replace('NODE', 'N')))
        name_of_file: str = '_'.join(file.split('/')[-1].split('.')[0].split('_')[0:2])
        path_to_file: str = '/'.join(file.split('/')[:-1])
        os.rename(file, path_to_file + '/' + name_of_file + '.fasta')
    print('Done\n')


def create_hit_tables():
    print('Creating hit tables...')
    for file in glob.glob(main_path + '*.fasta'):
        file: str = file.split('/')[-1]
        sample: str = file[:-6]
        print('\tProcessing ' + sample)
        NcbimakeblastdbCommandline(dbtype='nucl', input_file=main_path + file,
                                   out=main_path + sample, parse_seqids=True)()
        print('\tRunning BLAST...')
        NcbiblastnCommandline(task='blastn', query=probe_HP_one_repr, db=main_path + sample,
                              out=main_path + 'reference_in_' + sample + '_contigs.txt', qcov_hsp_perc=length_cover,
                              num_threads=8, outfmt='6 qaccver saccver pident qcovhsp evalue bitscore sstart send')()
        print('\tOK')
    print('Done\n')


def correct_contgis():
    print('Correcting contigs...')
    for file in glob.glob(main_path + '*.fasta'):
        sample: str = file.split('/')[-1][:-6]
        print(' Processing ' + sample)
        statistics[sample]: Dict[str, Dict[str, List[str]]] = dict()
        hits: List[str] = list()
        with open(main_path + 'reference_in_' + sample + '_contigs.txt') as blast_results:
            blast_results_as_list = [x[:-1] for x in blast_results.readlines()]
            for line in blast_results_as_list:
                if contig_locus(line) == locus(line) and float(line.split()[1].split('_c_')[1]) >= spades_cover:
                    hits.append(line)
        sort_hit_table_cover(hits, 'exon')
        hits_exons_contigs: Set[str] = set()
        hits_dedup: List[str] = list()
        for hit in hits:
            if str(exon(hit) + ' ' + contig(hit)) not in hits_exons_contigs:
                hits_dedup.append(hit)
            else:
                pass
            hits_exons_contigs.add(exon(hit) + ' ' + contig(hit))
        with open(main_path + sample + '.fas', 'w') as result_fasta, open(file) as contigs:
            contigs_fasta_parsed = SeqIO.to_dict(SeqIO.parse(contigs, 'fasta', generic_dna))
            for hit_dedup in hits_dedup:
                if int(hit_dedup.split()[6]) > int(hit_dedup.split()[7]):
                    sequence: str = slicing(contigs_fasta_parsed, hit_dedup, 1, 7, 6, 'yes')
                    result_fasta.write('>' + exon(hit_dedup) + '_' + '_'.join(contig(hit_dedup).split('_')[1:]) + '\n' +
                                       sequence + '\n')
                    all_hits_for_reference.append('{0}\t{1}\t{2}'.format(hit_dedup, sample, sequence))
                else:
                    sequence: str = slicing(contigs_fasta_parsed, hit_dedup, 1, 6, 7, 'no')
                    result_fasta.write('>' + exon(hit_dedup) + '_' + '_'.join(contig(hit_dedup).split('_')[1:]) + '\n' +
                                       sequence + '\n')
                    all_hits_for_reference.append('{0}\t{1}\t{2}'.format(hit_dedup, sample, sequence))
        sort_hit_table_cover(hits_dedup, 'locus')
        hits_loci: Set[str] = set()
        hits_exons: Set[str] = set()
        for hit_dedup in hits_dedup:
            if locus(hit_dedup) not in hits_loci:
                statistics[sample][locus(hit_dedup)] = dict()
                if exon(hit_dedup) not in hits_exons:
                    statistics[sample][locus(hit_dedup)][exon(hit_dedup)] = [contig(hit_dedup)]
                    hits_exons.add(exon(hit_dedup))
                else:
                    statistics[sample][locus(hit_dedup)][exon(hit_dedup)].append(contig(hit_dedup))
                hits_loci.add(locus(hit_dedup))
            else:
                if exon(hit_dedup) not in hits_exons:
                    statistics[sample][locus(hit_dedup)][exon(hit_dedup)] = [contig(hit_dedup)]
                    hits_exons.add(exon(hit_dedup))
                else:
                    statistics[sample][locus(hit_dedup)][exon(hit_dedup)].append(contig(hit_dedup))
        for i in statistics[sample].keys():
            stat = []
            for j in statistics[sample][i].keys():
                stat.append(len(statistics[sample][i][j]))
            statistics[sample][i]: int = max(stat)
        with open(main_path + 'reference_against_' + sample + '_contigs.txt', 'w') as \
                hittable:
            sort_hit_table_cover(hits, 'locus')
            for hit in hits:
                hittable.write(hit + '\n')
    with open(path_to_data_HPM + '/exons/all_hits.txt', 'w') as all_hits_to_write:
        for hit in all_hits_for_reference:
            all_hits_to_write.write(hit + '\n')
    print(' OK')
    print('All contigs were successfully corrected!\n')


def write_stats():
    print('Writing statistics...')
    with open(main_path + 'statistics.csv', 'w') as stats, open(probe_HP_one_repr) as \
            reference:
        stats_dict: Dict[str, str] = dict([('gene\t', '')])
        loci: Set[str] = set()
        samples: List[str] = list()
        reference_as_list: List[str] = [x[:-1] for x in reference.readlines()]
        for line in reference_as_list:
            if line.startswith('>'):
                loci.add('_'.join(line.split('-')[1].split('_')[:-2]))
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


def rename_contigs():
    print('Renaming contigs...')
    for file in glob.glob(main_path + '*.fas'):
        sample = file.split('/')[-1][:-4]
        print(' Processing ' + sample)
        with open(file) as result_fasta:
            fasta_parsed = SeqIO.to_dict(SeqIO.parse(result_fasta, 'fasta', generic_dna))
            counter: int = 1
            fasta_to_write: List[str] = list()
            for line in sorted(fasta_parsed.keys()):
                fasta_to_write.append('>Contig' + str(counter) + '_' + sample + '-' + line.replace('_', '-')
                                      + '\n')
                fasta_to_write.append(str(fasta_parsed[line].seq) + '\n')
                counter += 1
        with open(file, 'w') as result_fasta:
            result_fasta.writelines(fasta_to_write)
        print(' OK')
    print('All contigs were successfully renamed!\n')


def clean():
    print('Removing temporary files...')
    for file in glob.glob(main_path + '*.fasta'):
        os.remove(file)
    # for file in glob.glob(main_path + 'reference_in*'):
    #     os.remove(file)
    for file in glob.glob(main_path + '*.n*'):
        os.remove(file)
    print('Done\n')


if __name__ == "__main__":
    print('Converting data...\n')
    print('**********************************************************************************************************')
    print('\n')
    path_to_data_HP: str = sys.argv[1].strip()
    path_to_data_HPM: str = sys.argv[2].strip()
    probe_HP_one_repr: str = sys.argv[3].strip()
    length_cover: float = float(sys.argv[4].strip())
    spades_cover: float = float(sys.argv[5].strip())
    main_path: str = path_to_data_HPM + '/exons/40contigs/'
    os.makedirs(main_path, exist_ok=True)
    statistics: Dict[str, Dict[str, Union[Dict[str, List[str]], int]]] = dict()
    all_hits_for_reference: List[str] = list()
    prepare_contigs()
    create_hit_tables()
    correct_contgis()
    write_stats()
    rename_contigs()
    clean()
    print('**********************************************************************************************************')
    print('\nData was successfully converted!')
