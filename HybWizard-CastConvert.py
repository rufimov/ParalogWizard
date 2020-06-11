import os
import glob
import sys
import shutil
import re

path_to_data_HP = sys.argv[1]
path_to_data_HPM = sys.argv[2]
probe_HP_one_repr = sys.argv[3]
length_cover = int(sys.argv[4])
spades_cover = float(sys.argv[5])
new_reference_bool = sys.argv[6]

os.makedirs(path_to_data_HPM + '/exons/40contigs')
for file in glob.glob(path_to_data_HP + '/*contigs.fasta'):
    shutil.move(file, path_to_data_HPM + '/exons/40contigs/' + file.split('/')[-1])
for file in glob.glob(path_to_data_HPM + '/exons/40contigs/*.fasta'):
    with open(file, 'r') as fasta:
        lines = fasta.readlines()
    with open(file, 'w') as fasta:
        for line in lines:
            fasta.write(re.sub(r'length_([0-9]+)_cov_([0-9]+\.[0-9][0-9]).*', r'\1_c_\2', line.replace('NODE', 'N')))
    name_of_file = file.split('/')[-1]
    path_to_file = file.split('/')[:-1]
    os.rename(file, '/'.join(path_to_file) + '/' + name_of_file.split('.')[0] + '.fasta')
for file in glob.glob(path_to_data_HPM + '/exons/40contigs/*.fasta'):
    file = file.split('/')[-1]
    sample = file[:-6]
    os.system('echo -e "\n\tProcessing %(sample)s"\n'
              'makeblastdb -in %(path_to_data_HPM)s/exons/40contigs/%(file)s -parse_seqids -dbtype nucl '
              '-out %(path_to_data_HPM)s/exons/40contigs/%(sample)s || exit 1\n'
              'echo -e "\tRunning BLAST..."\n'
              'blastn -task blastn '
              '-db %(path_to_data_HPM)s/exons/40contigs/%(sample)s '
              '-query %(probe_HP_one_repr)s '
              '-out %(path_to_data_HPM)s/exons/40contigs/reference_in_%(sample)s_contigs.txt '
              '-outfmt "6 qaccver saccver pident qcovhsp evalue bitscore sstart send" || exit 1\n'
              'echo -e "\tOK"' % {'file': file,
                                  'sample': sample,
                                  'path_to_data_HPM': path_to_data_HPM,
                                  'probe_HP_one_repr': probe_HP_one_repr})

print('Done\n\nCorrecting contigs..')
statistics = {}
all_hits_for_reference = []
for sample in glob.glob(path_to_data_HPM + '/exons/40contigs/*.fasta'):
    sample = sample.split('/')[-1]
    print(' Processing ' + sample)
    statistics[sample[:-6]] = {}
    hits = []
    with open(path_to_data_HPM + '/exons/40contigs/' + 'reference_in_' + sample[:-6] + '_contigs.txt') \
            as blast_results, \
            open(path_to_data_HPM + '/exons/40contigs/' + sample[:-6] + '.fas', 'w') as result_fasta, \
            open(path_to_data_HPM + '/exons/40contigs/' + sample) as contigs:
        corrected_contigs_fasta = []
        for line in contigs.read().splitlines():
            if line.startswith('>'):
                corrected_contigs_fasta.append(line)
                corrected_contigs_fasta.append('')
            else:
                corrected_contigs_fasta[-1] = corrected_contigs_fasta[-1] + line
        contigs_fasta_parsed = dict()
        for i in range(0, len(corrected_contigs_fasta), 2):
            contigs_fasta_parsed[corrected_contigs_fasta[i]] = corrected_contigs_fasta[i + 1]
        for line in blast_results.read().splitlines():
            if line.split()[1].split('_N_')[0] == line.split()[0].split('-')[1] and int(line.split()[3]) >= \
                    length_cover and float(line.split()[1].split('_c_')[1]) >= spades_cover:
                hits.append(line)
        hits.sort(key=lambda x: float(x.split()[5]), reverse=True)
        hits.sort(key=lambda x: float(x.split()[4]))
        hits.sort(key=lambda x: float(x.split()[2]), reverse=True)
        hits.sort(key=lambda x: float(x.split()[3]), reverse=True)
        hits.sort(key=lambda x: x.split()[1])
        contig_hits = set()
        for hit in hits:
            if hit.split()[1] not in contig_hits:
                if int(hit.split()[6]) > int(hit.split()[7]):
                    result_fasta.write('>' + hit.split()[1] + '\n' + contigs_fasta_parsed['>' + hit.split()[1]][
                                                                     int(hit.split()[7]) - 1:int(
                                                                         hit.split()[6])] + '\n')
                    all_hits_for_reference.append(
                        '{0}\t{1}\t{2}'.format(hit, sample[:-6], contigs_fasta_parsed['>' + hit.split()[1]][
                                                                 int(hit.split()[7]) - 1:int(
                                                                     hit.split()[6])]))

                else:
                    result_fasta.write('>' + hit.split()[1] + '\n' + contigs_fasta_parsed['>' + hit.split()[1]][
                                                                     int(hit.split()[6]) - 1:int(
                                                                         hit.split()[7])] + '\n')
                    all_hits_for_reference.append(
                        '{0}\t{1}\t{2}'.format(hit, sample[:-6], contigs_fasta_parsed['>' + hit.split()[1]][
                                                                 int(hit.split()[6]) - 1:int(
                                                                     hit.split()[7])]))
                contig_hits.add(hit.split()[1])
            else:
                pass
        hits.sort(key=lambda x: float(x.split()[5]), reverse=True)
        hits.sort(key=lambda x: float(x.split()[4]))
        hits.sort(key=lambda x: float(x.split()[2]), reverse=True)
        hits.sort(key=lambda x: float(x.split()[3]), reverse=True)
        hits.sort(key=lambda x: x.split()[0].split('-')[1])
        hits_loci_contigs = set()
        hits_dedup = []
        for hit in hits:
            if str(hit.split()[0].split('-')[1] + ' ' + hit.split()[1]) not in hits_loci_contigs:
                hits_dedup.append(hit)
            else:
                pass
            hits_loci_contigs.add(hit.split()[0].split('-')[1] + ' ' + hit.split()[1])
        hits_dedup.sort(key=lambda x: float(x.split()[5]), reverse=True)
        hits_dedup.sort(key=lambda x: float(x.split()[4]))
        hits_dedup.sort(key=lambda x: float(x.split()[3]), reverse=True)
        hits_dedup.sort(key=lambda x: float(x.split()[2]), reverse=True)
        hits_dedup.sort(key=lambda x: x.split()[0].split('-')[1])
        hits_loci = set()
        for hit_dedup in hits_dedup:
            if hit_dedup.split()[0].split('-')[1] not in hits_loci:
                statistics[sample[:-6]][hit_dedup.split()[0].split('-')[1]] = 1
            else:
                statistics[sample[:-6]][hit_dedup.split()[0].split('-')[1]] += 1
            hits_loci.add(hit_dedup.split()[0].split('-')[1])
    with open(path_to_data_HPM + '/exons/40contigs/' + 'reference_against_' + sample[:-6] + '_contigs.txt', 'w') as \
            hittable:
        for hit in hits:
            hittable.write(hit + '\n')
    print(' OK')
print('All contigs were successfully corrected!\n')
print('Writing statistics...')
with open(path_to_data_HPM + '/exons/40contigs/statistics.csv', 'w') as stats, open(probe_HP_one_repr) as \
        reference:
    stats_dict = dict([('gene\t', '')])
    loci = set()
    samples = []
    for line in reference.read().splitlines():
        if line.startswith('>'):
            loci.add(line[1:].split('-')[1])
    for key in statistics.keys():
        samples.append(key)
    samples.sort()
    for sample in samples:
        stats_dict['gene\t'] = stats_dict['gene\t'] + sample + '\t'
    for locus in loci:
        stats_dict[locus + '\t'] = ''
        for sample in samples:
            loci_in_sample = set(statistics[sample].keys())
            if locus in loci_in_sample:
                stats_dict[locus + '\t'] = stats_dict[locus + '\t'] + str(statistics[sample][locus]) + '\t'
            else:
                stats_dict[locus + '\t'] = stats_dict[locus + '\t'] + 'NA' + '\t'
    stats.write('gene\t' + stats_dict['gene\t'] + '\n')
    del stats_dict['gene\t']
    for key in sorted(list(stats_dict.keys())):
        stats.write(key + stats_dict[key] + '\n')
print('Statistics file created!\n')
if new_reference_bool == 'yes':
    print('Creating new reference...')
    all_hits_for_reference.sort(key=lambda x: float(x.split()[5]), reverse=True)
    all_hits_for_reference.sort(key=lambda x: float(x.split()[4]))
    all_hits_for_reference.sort(key=lambda x: float(x.split()[2]), reverse=True)
    all_hits_for_reference.sort(key=lambda x: float(x.split()[3]), reverse=True)
    all_hits_for_reference.sort(key=lambda x: x.split()[0])
    exons = set()
    with open(path_to_data_HPM + '/exons/new_reference_for_HybPhyloMaker.fas', 'w') as new_reference:
        for hit in all_hits_for_reference:
            if hit.split()[0] not in exons:
                name_of_locus = hit.split()[0].split('-')[1].replace('exon', 'Contig').replace('Exon', 'Contig') \
                    .replace('contig', 'Contig').replace('_', '').replace('Contig', '_Contig_')
                new_reference.write('>Assembly_' + name_of_locus + '_' + hit.split()[-2] + '\n' + hit.split()[-1] +
                                    '\n')
            else:
                pass
            exons.add(hit.split()[0])
    print('New reference created!\n')
print('Renaming contigs...')
for sample in glob.glob(path_to_data_HPM + '/exons/40contigs/*.fasta'):
    sample = sample.split('/')[-1]
    print(' Processing ' + sample)
    with open(path_to_data_HPM + '/exons/40contigs/' + sample[:-6] + '.fas') as result_fasta:
        fasta_as_list = result_fasta.read().splitlines()
        fasta_parsed = dict()
        for i in range(0, len(fasta_as_list), 2):
            fasta_parsed[fasta_as_list[i]] = fasta_as_list[i + 1]
        counter = 1
        fasta_to_write = []
        fasta_parsed_as_list = list(fasta_parsed.keys())
        fasta_parsed_as_list = sorted(fasta_parsed_as_list)
        for line in fasta_parsed_as_list:
            fasta_to_write.append('>Contig' + str(counter) + '_' + sample[:-6] + '\n')
            fasta_to_write.append(fasta_parsed[line] + '\n')
            counter += 1
    with open(path_to_data_HPM + '/exons/40contigs/' + sample[:-6] + '.fas', 'w') as result_fasta:
        result_fasta.writelines(fasta_to_write)
    print(' OK')
print('All contigs were successfully renamed!\n')
print('Removing temporary files...')
for file in glob.glob(path_to_data_HPM + '/exons/40contigs/*.fasta'):
    os.remove(file)
for file in glob.glob(path_to_data_HPM + '/exons/40contigs/reference_in*'):
    os.remove(file)
for file in glob.glob(path_to_data_HPM + '/exons/40contigs/*.n*'):
    os.remove(file)
print('Done\n')
print('**********************************************************************************************************')
print('\nData was successfully converted!')
