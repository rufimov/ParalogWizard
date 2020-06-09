import sys
probes = sys.argv[1]
import os


def aln_similarity(y):
    return (float(y.split('\t')[0]) / (float(y.split('\t')[0]) + float(y.split('\t')[1]))) * 100


def number_locus_for_sort(z):
    return z.split()[13].split('_')[1]


def number_exon_for_sort(a):
    return int(a.split()[13].split('_')[3])


def number_contig_for_sort(b):
    return int(b.split()[9].split('_')[0][6:])


os.system('ls *.pslx > list_of_pslx.txt\n'
          'mkdir corrected')

with open("list_of_pslx.txt") as list:
    list_to_process = list.read().splitlines()
for file in list_to_process:
    if file != '%s.pslx' % probes:
        with open(file) as pslx_file, open('corrected/' + file, 'w') as corrected_pslx_file:
            file = pslx_file.read().splitlines()
            head = file[0:5]
            list_to_work = file[5:]
            list_to_work.sort(key=aln_similarity, reverse=True)
            list_to_work.sort(key=number_exon_for_sort)
            list_to_work.sort(key=number_locus_for_sort)
            list_to_work_cleaned1 = []
            hits1 = set()
            for line in list_to_work:
                if line.split()[13] not in hits1:
                    list_to_work_cleaned1.append(line)
                    hits1.add(line.split()[13])
            list_to_work_cleaned1.sort(key=aln_similarity, reverse=True)
            list_to_work_cleaned1.sort(key=number_contig_for_sort)
            list_to_work_cleaned2 = []
            hits2 = set()
            for line in list_to_work_cleaned1:
                if line.split()[9] not in hits2:
                    list_to_work_cleaned2.append(line)
                    hits2.add(line.split()[9])
            for line in head:
                corrected_pslx_file.write(line + '\n')
            for line in list_to_work_cleaned2:
                corrected_pslx_file.write(line + '\n')
