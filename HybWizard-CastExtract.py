import sys
import gffutils
from Bio import SeqRecord, SeqIO
from Bio.Alphabet import generic_dna

gff = sys.argv[1]
fasta_file = sys.argv[2]
result_fasta = sys.argv[3]
concatenated_fasta = sys.argv[4]
print('Creating gff database...')
db = gffutils.create_db(gff, dbfn=gff[:-4], force=True, keep_order=True,
                        merge_strategy='merge', sort_attribute_values=True)
print('Done')
print('Creating feature database...')
db = gffutils.FeatureDB(gff[:-4], keep_order=True)
genes = list(db.features_of_type('gene'))
pseudogenes = list(db.features_of_type('pseudogene'))
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
        exons = []
        for i in db.children(gene, featuretype='exon', order_by='start'):
            exons.append(i.astuple())
        exons.sort(key=lambda x: int(x[0].split('-')[-1]))
        exons.sort(key=lambda x: '-'.join(x[0].split('-')[:-1]))
        for exon in exons:
            name = exon[0]
            chromosome = exon[1]
            start = int(exon[4])
            end = int(exon[5])
            strand = exon[7]
            if strand == '-':
                sequence = str(correct_fasta_parsed[chromosome][start - 1:end].seq.reverse_complement())
            if strand == '+':
                sequence = str(correct_fasta_parsed[chromosome][start - 1:end].seq)
            all_exons.write('>' + name + '-' + chromosome + '\n' + sequence + '\n')
    for pseudogene in pseudogenes:
        pseudogene_exons = []
        for i in db.children(pseudogene, featuretype='exon', order_by='start'):
            pseudogene_exon = i.astuple()
            if pseudogene_exon[0].startswith('exon-XR'):
                pseudogene_exons.append(pseudogene_exon)
            else:
                name = pseudogene_exon[0]
                chromosome = pseudogene_exon[1]
                start = int(pseudogene_exon[4])
                end = int(pseudogene_exon[5])
                strand = pseudogene_exon[7]
                sequence = str(correct_fasta_parsed[chromosome][start - 1:end].seq)
                all_exons.write('>' + name + '-' + chromosome + '\n' + sequence + '\n')
        pseudogene_exons.sort(key=lambda x: int(x[0].split('-')[-1]))
        pseudogene_exons.sort(key=lambda x: '-'.join(x[0].split('-')[:-1]))
        for pseudogene_exon in pseudogene_exons:
            name = pseudogene_exon[0]
            chromosome = pseudogene_exon[1]
            start = int(pseudogene_exon[4])
            end = int(pseudogene_exon[5])
            strand = pseudogene_exon[7]
            if strand == '-':
                sequence = str(correct_fasta_parsed[chromosome][start - 1:end].seq.reverse_complement())
            if strand == '+':
                sequence = str(correct_fasta_parsed[chromosome][start - 1:end].seq)
            all_exons.write('>' + name + '-' + chromosome + '\n' + sequence + '\n')
print('Done')
print('Concatenating exons...')
with open(result_fasta) as fasta_to_concatenate, open(concatenated_fasta, 'w') as concat_fasta:
    fasta_parsed = SeqIO.to_dict(SeqIO.parse(fasta_to_concatenate, 'fasta', generic_dna))
    current_transcript = ''
    for key in sorted(list(fasta_parsed.keys())):
        transcript = key.split('-')[1]
        if transcript != current_transcript:
            concat_fasta.write('\n>' + transcript + '\n' + str(fasta_parsed[key].seq))

        else:
            concat_fasta.write(str(fasta_parsed[key].seq))
        current_transcript = transcript
print('Done')


