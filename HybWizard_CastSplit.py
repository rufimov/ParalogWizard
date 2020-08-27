from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastnCommandline
from Bio import SeqIO
from Bio.Alphabet import generic_dna
import os
import glob
import sys
concat_exons = sys.argv[1]  # concatenated exons
separat_exons = sys.argv[2]  # separated exons
probes = sys.argv[3]  # probe file
best_separate_exons = 'best_hits_as_exons.fasta'
result_file = sys.argv[4]  # final output file
blast_task = sys.argv[5]  # blastn task ('megablast', 'dc-megablast', 'blastn')

print('Building database for %s...' % concat_exons)
NcbimakeblastdbCommandline(dbtype='nucl', input_file=concat_exons,
                           out=concat_exons, parse_seqids=True)()
print('Done')
print('Blasting %(probes)s against %(database)s' % {'probes': probes,
                                                    'database': concat_exons})
NcbiblastnCommandline(task=blast_task, query=probes, db=concat_exons,
                      out=probes + '_against_' + concat_exons + '.txt',
                      outfmt="6 qaccver saccver pident qcovhsp evalue bitscore",
                      num_threads=4)()
print('Done')
print('Choosing the best hit...')
with open(probes + '_against_' + concat_exons + '.txt') as blast_results:
    hits = blast_results.readlines()
    hits.sort(key=lambda x: float(x.split()[5]), reverse=True)
    hits.sort(key=lambda x: float(x.split()[4]))
    hits.sort(key=lambda x: float(x.split()[2]), reverse=True)
    hits.sort(key=lambda x: float(x.split()[3]), reverse=True)
    hits.sort(key=lambda x: x.split()[0].split('-')[1])
    uniqe_hits = set()
    best_hits = []
    for hit in hits:
        if hit.split()[0].split('-')[1] not in uniqe_hits:
            best_hits.append(hit)
            uniqe_hits.add(hit.split()[0].split('-')[1])
print('Done')
print('Splitting best hits to exons...')
with open(separat_exons) as all_exons:
    all_exons_parsed = SeqIO.to_dict(SeqIO.parse(all_exons, 'fasta', generic_dna))
with open(best_separate_exons, 'w') as best_exons:
    for best_hit in best_hits:
        locus = best_hit.split()[1]
        probe = best_hit.split()[0]
        exons = [val for key, val in all_exons_parsed.items() if locus in key]
        for exon in exons:
            name = str(exon.id)
            sequence = str(exon.seq)
            best_exons.write('>' + probe + '_' + name + '\n' + sequence + '\n')
NcbimakeblastdbCommandline(dbtype='nucl', input_file=probes,
                           out=probes, parse_seqids=True)()
NcbiblastnCommandline(task=blast_task, query=best_separate_exons, db=probes,
                      out=best_separate_exons + '_against_' + probes + '.txt', num_threads=4,
                      outfmt='6 qaccver saccver pident qcovhsp evalue bitscore sstart send qstart qend')()
with open(probes) as probes_to_parse:
    probes_as_dict = SeqIO.to_dict(SeqIO.parse(probes_to_parse, 'fasta', generic_dna))
with open(best_separate_exons + '_against_' + probes + '.txt') as new_blast_results, \
        open(result_file, 'w') as result_file:
    hits = new_blast_results.readlines()
    cleaned_hits = []
    for hit in hits:
        if hit.split()[0].split('_exon-')[0] == hit.split()[1]:
            cleaned_hits.append(hit)
    cleaned_hits.sort(key=lambda x: float(x.split()[5]), reverse=True)
    cleaned_hits.sort(key=lambda x: float(x.split()[4]))
    cleaned_hits.sort(key=lambda x: float(x.split()[2]), reverse=True)
    cleaned_hits.sort(key=lambda x: float(x.split()[3]), reverse=True)
    cleaned_hits.sort(key=lambda x: x.split()[0])
    hits_exons = set()
    cleaned_dedup_hits = []
    for cleaned_hit in cleaned_hits:
        if cleaned_hit.split()[0] not in hits_exons:
            cleaned_dedup_hits.append(cleaned_hit)
            hits_exons.add(cleaned_hit.split()[0])
    cleaned_dedup_hits.sort(key=lambda x: x.split()[0])
    cleaned_dedup_hits.sort(key=lambda x: x.split()[1].split('-')[1])
    for cleaned_dedup_hit in cleaned_dedup_hits:
        name_of_locus = cleaned_dedup_hit.split()[1]
        num_exon = cleaned_dedup_hit.split()[0].split('_exon-')[1].split('-')[1]
        if num_exon.isdigit():
            pass
        else:
            num_exon = '1'
        if int(cleaned_dedup_hit.split()[6]) > int(cleaned_dedup_hit.split()[7]):
            start = int(cleaned_dedup_hit.split()[7])
            end = int(cleaned_dedup_hit.split()[6])
            sequence = str(probes_as_dict[name_of_locus][start - 1:end].seq.reverse_complement())
        else:
            start = int(cleaned_dedup_hit.split()[6])
            end = int(cleaned_dedup_hit.split()[7])
            sequence = str(probes_as_dict[name_of_locus][start - 1:end].seq)
        result_file.write('>' + name_of_locus + '_exon_' + num_exon + '\n' + sequence + '\n')
for file in glob.glob('*.n*'):
    os.remove(file)
os.('best_hits_as_exons.fasta')
print('Done')
