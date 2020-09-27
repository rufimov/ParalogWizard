import sys
import gffutils
from Bio import SeqRecord, SeqIO
from Bio.Alphabet import generic_dna

gff = sys.argv[1]
fasta_file = sys.argv[2]
result_fasta = sys.argv[3]
concatenated_fasta = sys.argv[4]
# print('Creating gff database...')
# db = gffutils.create_db(gff, dbfn=gff[:-4], force=True, keep_order=True,
#                         merge_strategy='merge', sort_attribute_values=True)
# print('Done')
print('Creating feature database...')
db = gffutils.FeatureDB(gff[:-4], keep_order=True)
genes = list(db.features_of_type('gene'))
pseudogenes = list(db.features_of_type('pseudogene'))
genes.extend(pseudogenes)
print('Done')
print('Reading genome fasta...')
with open(fasta_file) as fasta:
    fasta_parsed = SeqIO.to_dict(SeqIO.parse(fasta, 'fasta', generic_dna))
print('Done')
print('Writing fasta with exons...')
correct_fasta_parsed = {}
for key in fasta_parsed:
    correct_fasta_parsed[key.split()[0]] = fasta_parsed[key]
with open(result_fasta, 'w') as all_exons:
    for gene in genes:
        # print(gene)
        type_of_feature = gene.attributes['gene_biotype'][0]
        if type_of_feature == 'protein_coding' or type_of_feature == 'pseudogene' or type_of_feature == 'lncRNA' \
                or type_of_feature == 'transcribed_pseudogene':
            if type_of_feature == 'pseudogene':
                strand = gene.strand
                name = gene.attributes['Name'][0].replace('-', '_')
                chromosome = gene.seqid
                exons = []
                for i in db.children(gene, featuretype='exon', order_by='start'):
                    print(i)
                    exons.append(i)
                exons.sort(key=lambda x: x.start)
                exons.sort(key=lambda x: x.attributes['gene'][0])
                if strand == '-':
                    exons = exons[::-1]
                    count = 1
                    transcripts = set()
                    for exon in exons:
                        transcript = exon.attributes['gene'][0].replace('-', '_')
                        if transcript not in transcripts:
                            count = 1
                            transcripts.add(transcript)
                        start = exon.start
                        end = exon.end
                        sequence = str(correct_fasta_parsed[chromosome][start - 1:end].seq.reverse_complement())
                        all_exons.write('>' + name + '-' + transcript + '-exon_' + str(count) +
                                        '-' + type_of_feature + '-' + chromosome + '\n' + sequence + '\n')
                        count += 1
                elif strand == '+':
                    count = 1
                    transcripts = set()
                    for exon in exons:
                        transcript = exon.attributes['gene'][0].replace('-', '_')
                        if transcript not in transcripts:
                            count = 1
                            transcripts.add(transcript)
                        start = exon.start
                        end = exon.end
                        sequence = str(correct_fasta_parsed[chromosome][start - 1:end].seq)
                        all_exons.write('>' + name + '-' + transcript + '-exon_' + str(count) +
                                        '-' + type_of_feature + '-' + chromosome + '\n' + sequence + '\n')
                        count += 1
            else:
                strand = gene.strand
                name = gene.attributes['Name'][0].replace('-', '_')
                chromosome = gene.seqid
                exons = []
                for i in db.children(gene, featuretype='exon', order_by='start'):
                     # print(i)
                    exons.append(i)
                exons.sort(key=lambda x: x.start)
                try:
                    exons.sort(key=lambda x: x.attributes['transcript_id'][0])
                except KeyError:
                    exons.sort(key=lambda x: x.attributes['orig_transcript_id'][0])
                if strand == '-':
                    exons = exons[::-1]
                    count = 1
                    transcripts = set()
                    for exon in exons:
                        try:
                            transcript = exon.attributes['transcript_id'][0].replace('-', '_')
                        except KeyError:
                            transcript = exon.attributes['orig_transcript_id'][0].replace('-', '_')
                        if transcript not in transcripts:
                            count = 1
                            transcripts.add(transcript)
                        start = exon.start
                        end = exon.end
                        sequence = str(correct_fasta_parsed[chromosome][start - 1:end].seq.reverse_complement())
                        all_exons.write('>' + name + '-' + transcript + '-exon_' + str(count) +
                                        '-' + type_of_feature + '-' + chromosome + '\n' + sequence + '\n')
                        count += 1
                elif strand == '+':
                    count = 1
                    transcripts = set()
                    for exon in exons:
                        try:
                            transcript = exon.attributes['transcript_id'][0].replace('-', '_')
                        except KeyError:
                            transcript = exon.attributes['orig_transcript_id'][0].replace('-', '_')
                        if transcript not in transcripts:
                            count = 1
                            transcripts.add(transcript)
                        start = exon.start
                        end = exon.end
                        sequence = str(correct_fasta_parsed[chromosome][start - 1:end].seq)
                        all_exons.write('>' + name + '-' + transcript + '-exon_' + str(count) +
                                        '-' + type_of_feature + '-' + chromosome + '\n' + sequence + '\n')
                        count += 1
print('Done')
print('Concatenating exons...')
with open(result_fasta) as fasta_to_concatenate, open(concatenated_fasta, 'w') as concat_fasta:
    fasta_parsed = SeqIO.to_dict(SeqIO.parse(fasta_to_concatenate, 'fasta', generic_dna))
    current_transcript = ''
    count = 1
    list_of_keys = list(fasta_parsed.keys())
    list_of_keys.sort(key=lambda x: int(x.split('-')[2].split('_')[1]))
    list_of_keys.sort(key=lambda x: x.split('-')[1])
    list_of_keys.sort(key=lambda x: x.split('-')[0])
    list_of_keys.sort(key=lambda x: x.split('-')[-1])
    for key in list_of_keys:
        transcript = key.split('-')[1]
        locus = key.split('-')[0]
        type_of_feature = key.split('-')[3]
        chromosome = key.split('-')[4]
        if count == 1:
            if transcript != current_transcript:
                concat_fasta.write('>' + locus + '-' + transcript +
                                   '-' + type_of_feature + '-' + chromosome + '\n' + str(fasta_parsed[key].seq))
            else:
                concat_fasta.write(str(fasta_parsed[key].seq))
        else:
            if transcript != current_transcript:
                concat_fasta.write('\n>' + locus + '-' + transcript +
                                   '-' + type_of_feature + '-' + chromosome + '\n' + str(fasta_parsed[key].seq))
            else:
                concat_fasta.write(str(fasta_parsed[key].seq))
        current_transcript = transcript
        count += 1
print('Done')


