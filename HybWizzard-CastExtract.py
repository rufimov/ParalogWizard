import gffutils
from Bio import SeqRecord, SeqIO
from Bio.Alphabet import generic_dna
db = gffutils.create_db('GCF_002114115.1_ASM211411v1_genomic.gff', dbfn='Malus_2017', force=True, keep_order=True,
                        merge_strategy='merge', sort_attribute_values=True)
db = gffutils.FeatureDB('Malus_2017', keep_order=True)
genes = list(db.features_of_type('gene'))
pseudogenes = list(db.features_of_type('pseudogene'))
with open('GCF_002114115.1_ASM211411v1_genomic.fas') as fasta:
    fasta_parsed = SeqIO.to_dict(SeqIO.parse(fasta, 'fasta', generic_dna))
for key in fasta_parsed:
    fasta_parsed[key.split()[0]] = fasta_parsed.pop(key)
    del fasta_parsed[key]
with open('all_exons_Malus_2017.fasta', 'w') as all_exons:
    for gene in genes:
        exons = []
        for i in db.children(gene, featuretype='exon', order_by='start'):
            exons.append(i.astuple())
        exons.sort(key=lambda x: int(x[0].split('-')[-1]))
        for exon in exons:
            name = exon[0]
            chromosome = exon[1]
            start = int(exon[4])
            end = int(exon[5])
            strand = exon[7]
            if strand == '-':
                sequence = str(fasta_parsed[chromosome][start - 1:end].seq.reverse_complement())
            if strand == '+':
                sequence = str(fasta_parsed[chromosome][start - 1:end].seq)
            all_exons.write('>' + name + '_' + chromosome + '_' + '\n' + sequence + '\n')


