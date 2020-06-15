from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastnCommandline

NcbimakeblastdbCommandline(dbtype='nucl', input_file='all_exons_Malus_2017.fasta',
                           out='all_exons_Malus_2017', parse_seqids=True)()
NcbiblastnCommandline(task='blastn', query=, db='all_exons_Malus_2017',
                      out=main_path + 'reference_in_' + sample + '_contigs.txt',
                      outfmt='6 qaccver saccver pident qcovhsp evalue bitscore sstart send')()
