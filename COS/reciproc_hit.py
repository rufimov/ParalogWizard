from typing import Set, Any

from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastnCommandline
from Bio import SeqIO
from Bio.Alphabet import generic_dna
import glob
import os

NcbimakeblastdbCommandline(dbtype='nucl', input_file='cos_ref_HanXRQr1.0_exons.fasta',
                           out='cos_ref_HanXRQr1.0_exons', parse_seqids=True)()
NcbimakeblastdbCommandline(dbtype='nucl', input_file='cos_ref_Lsat_Salinas_v7_exons.fasta',
                           out='cos_ref_Lsat_Salinas_v7_exons', parse_seqids=True)()
NcbiblastnCommandline(task='blastn', query='cos_ref_Lsat_Salinas_v7_exons.fasta', db='cos_ref_HanXRQr1.0_exons',
                      out='lett_in_sunf.txt', qcov_hsp_perc=80,
                      outfmt="6 qaccver saccver pident qcovhsp evalue bitscore",
                      num_threads=4)()
NcbiblastnCommandline(task='blastn', query='cos_ref_HanXRQr1.0_exons.fasta', db='cos_ref_Lsat_Salinas_v7_exons',
                      out='sunf_in_lett.txt', qcov_hsp_perc=80,
                      outfmt="6 qaccver saccver pident qcovhsp evalue bitscore",
                      num_threads=4)()
for file in glob.glob('*.n*'):
    os.remove(file)
with open('lett_in_sunf.txt') as lett_in_sunf, open('sunf_in_lett.txt') as sunf_in_lett:
    lett_in_sunf_list = lett_in_sunf.readlines()
    sunf_in_lett_list = sunf_in_lett.readlines()
    sunf1 = set()
    for line1 in lett_in_sunf_list:
        sunf1.add(line1.split()[1])
    sunf2 = set()
    for line2 in sunf_in_lett_list:
        sunf2.add(line2.split()[0])
    lett1 = set()
    for line3 in lett_in_sunf_list:
        lett1.add(line3.split()[0])
    lett2 = set()
    for line4 in sunf_in_lett_list:
        lett2.add(line4.split()[1])
    sunf = sunf1.intersection(sunf2)
    lett: Set[str] = lett1.intersection(lett2)
with open('cos_ref_HanXRQr1.0_exons.fasta') as helianthus:
    sunf_parsed = SeqIO.to_dict(SeqIO.parse(helianthus, 'fasta', generic_dna))
with open('cos_ref_Lsat_Salinas_v7_exons.fasta') as lactuca:
    lett_parsed = SeqIO.to_dict(SeqIO.parse(lactuca, 'fasta', generic_dna))
with open('cos_ref_HanXRQr1.0_exons_wo_conflicts.fasta', 'w') as sunf_correct, open(
        'cos_ref_Lsat_Salinas_v7_exons_wo_conflicts.fasta', 'w') as lett_correct:
    sunf_list = list(sunf)
    sunf_list.sort(key=lambda x: int(x.split('-')[1].split('_')[2]))
    sunf_list.sort(key=lambda x: x.split('-')[1].split('_')[0])
    for line in sunf_list:
        SeqIO.write(sunf_parsed[line], sunf_correct, "fasta")
    lett_list = list(lett)
    lett_list.sort(key=lambda x: int(x.split('-')[1].split('_')[2]))
    lett_list.sort(key=lambda x: x.split('-')[1].split('_')[0])
    for line in lett_list:
        SeqIO.write(lett_parsed[line], lett_correct, "fasta")






