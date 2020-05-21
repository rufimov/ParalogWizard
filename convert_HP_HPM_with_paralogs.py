import os

import sys

path_to_data_HP = sys.argv[1]
path_to_data_HPM = sys.argv[2]
probe_HP_one_repr = sys.argv[3]

os.system(
    'echo "**********************************************************************************************************"\n'
    'echo "\nCollecting raw contigs..."\n'
    'cd %s\n'
    'for folder in $(find . -maxdepth 1 -type d | sed \'s/.\///\' | tail -n +2); do\n'
    'cd $folder\n'
    'echo "\tProcessing ${folder}"\n'
    'for gene in $(find . -maxdepth 1 -type d | sed \'s/.\///\' | tail -n +2); do\n'
    'locus=$gene\n'
    'if test -f "$locus/${locus}_contigs.fasta"; then\n'
    'cat $locus/${locus}_contigs.fasta | sed "s/>/>${locus}_/g" >> ../${folder}_contigs.fasta\n'
    'fi\n'
    'done\n'
    'cd ..\n'
    'echo "\tOK"\n'
    'done\n'
    'echo "Done"\n'
    'echo \n' % path_to_data_HP)
os.system('mkdir -p %(path_to_data_HPM)s/exons/40contigs\n'
          'mv %(path_to_data_HP)s/*contigs.fasta %(path_to_data_HPM)s/exons/40contigs\n'
          'cd %(path_to_data_HPM)s/exons/40contigs\n'
          'for file in $(ls *.fasta); do\n'
          'sed -i "s/_length.\+$//g" "$file"\n'
          'mv $file $(echo $file | sed \'s/dedup_contigs.//g\')\n'
          'ls *.fasta > list_of_files.txt\n'
          'done\n'
          'echo "Creating hit table for each sample..."' % {'path_to_data_HP': path_to_data_HP,
                                                            'path_to_data_HPM': path_to_data_HPM})
with open('%s/exons/40contigs/list_of_files.txt' % path_to_data_HPM) as list_of_files:
    for file in list_of_files.read().splitlines():
        sample = file[:-6]
        os.system('echo "\n\tProcessing %(sample)s"\n'
                  'makeblastdb -in %(path_to_data_HPM)s/exons/40contigs/%(file)s -parse_seqids -dbtype nucl '
                  '-out %(path_to_data_HPM)s/exons/40contigs/%(blast_database)s\n'
                  'echo "\tRunning BLAST..."\n'
                  'blastn -task blastn -evalue 1e-50 -db '
                  '%(path_to_data_HPM)s/exons/40contigs/%(blast_database)s -query '
                  '%(probe_HP_one_repr)s '
                  '-out %(path_to_data_HPM)s/exons/40contigs/reference_in_%(sample)s_contigs.txt  -outfmt "6 '
                  'qaccver saccver pident qcovs '
                  'evalue bitscore sstart send"\n'
                  'echo "\tOK"' % {'file': file, 'sample': sample, 'blast_database': file[:-6],
                                   'path_to_data_HPM': path_to_data_HPM,
                                   'probe_HP_one_repr': probe_HP_one_repr})

print('Done\n\nCorrecting contigs..')
for sample in open('%s/exons/40contigs/list_of_files.txt' % path_to_data_HPM).read().splitlines():
    print(' Processing ' + sample)
    with open('%s/exons/40contigs/' % path_to_data_HPM + 'reference_in_' + sample[:-6] + '_contigs'
                                                                                         '.txt') as \
            blast_results, open('%s/exons/40contigs/' % path_to_data_HPM + sample[:-6] + '.fas',
                                'w') as result_fasta, open(
        '%s/exons/40contigs/' % path_to_data_HPM + sample) \
            as contigs:
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
        hits = []
        for line in blast_results.read().splitlines():
            if line.split()[1].split('_NODE_')[0] == line.split()[0].split('-')[1] and int(line.split()[3]) >= 70:
                hits.append(line)
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

                else:
                    result_fasta.write('>' + hit.split()[1] + '\n' + contigs_fasta_parsed['>' + hit.split()[1]][
                                                                     int(hit.split()[6]) - 1:int(
                                                                         hit.split()[7])] + '\n')
                contig_hits.add(hit.split()[1])
            else:
                pass
    print(' OK')
print('All contigs were successfully corrected!\n')
print('Renaming contigs...')
for sample in open('%s/exons/40contigs/list_of_files.txt' % path_to_data_HPM).read().splitlines():
    print(' Processing ' + sample)
    with open('%s/exons/40contigs/' % path_to_data_HPM + sample[:-6] + '.fas') as result_fasta:
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
    with open('%s/exons/40contigs/' % path_to_data_HPM + sample[:-6] + '.fas', 'w') as result_fasta:
        result_fasta.truncate(0)
        result_fasta.writelines(fasta_to_write)
    print(' OK')
print('All contigs were successfully renamed!\n')
print('Removing temporary files...')
os.system('cd %s/exons/40contigs\n'
          'rm *.fasta *.txt *.n*\n' % path_to_data_HPM)
print('Done\n')
print('**********************************************************************************************************')
print('\nData was successfully converted!')
