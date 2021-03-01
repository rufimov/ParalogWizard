# !/usr/bin/env python
import errno
import os
import shutil
import subprocess
import sys

from Bio.SeqIO.QualityIO import FastqGeneralIterator

"""
This script is part of a pipeline to extract phylogenetically-useful sequences from 
Illumina data using the targeted (liquid-phase) sequence enrichment approach.

After a BWA search of the raw reads against the target sequences, the reads need to be 
sorted according to the successful hits. This script takes the BWA output (BAM format)
and the raw read files, and distributes the reads into FASTA files ready for assembly.

If there are multiple results (for example, one for each read direction),
concatenate them prior to sorting.
"""


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def read_sorting(bamfilename):
    samtools_cmd = "samtools view -F 4 {}".format(bamfilename)
    child = subprocess.Popen(
        samtools_cmd, shell=True, stdout=subprocess.PIPE, universal_newlines=True
    )
    bwa_results = child.stdout.readlines()

    read_hit_dict = {}
    for line in bwa_results:
        line = line.split()
        read_id = line[0]
        target = line[2].split("-")[-1]
        if read_id in read_hit_dict:
            if target not in read_hit_dict[read_id]:
                read_hit_dict[read_id].append(target)
        else:
            read_hit_dict[read_id] = [target]
    return read_hit_dict


def write_paired_seqs(target, id1, seq1, id2, seq2):
    mkdir_p(target)
    outfile = open(os.path.join(target, "{}_interleaved.fasta".format(target)), "a")
    outfile.write(">{}\n{}\n".format(id1, seq1))
    outfile.write(">{}\n{}\n".format(id2, seq2))
    outfile.close()


def distribute_reads(readfiles, read_hit_dict):
    num_reads_to_write = len(read_hit_dict)
    iterator1 = FastqGeneralIterator(open(readfiles[0]))
    reads_written = 0
    sys.stderr.write("Read distributing progress:\n")

    iterator2 = FastqGeneralIterator(open(readfiles[1]))

    for id1_long, seq1, qual1 in iterator1:

        id2_long, seq2, qual2 = next(iterator2)

        id1 = id1_long.split()[0]
        if id1.endswith("/1") or id1.endswith("/2"):
            id1 = id1[:-2]

        id2 = id2_long.split()[0]
        if id2.endswith("/1") or id2.endswith("/2"):
            id2 = id2[:-2]

        if id1 in read_hit_dict:
            for target in read_hit_dict[id1]:
                write_paired_seqs(target, id1, seq1, id2, seq2)
            reads_written += 1
        elif id2 in read_hit_dict:
            for target in read_hit_dict[id2]:
                write_paired_seqs(target, id1, seq1, id2, seq2)
            reads_written += 1
        j = (reads_written + 1) / num_reads_to_write
        if int(100 * j) % 5 == 0:
            sys.stderr.write("\r")
            sys.stderr.write("[%-20s] %d%%" % ("=" * int(20 * j), 100 * j))
            sys.stderr.flush()
    sys.stderr.write("\n")


def distribute_reads_to_targets_bwa(bamfilename, readfiles):
    read_hit_dict = read_sorting(bamfilename)
    print("Unique reads with hits: {}".format(len(read_hit_dict)))
    distribute_reads(readfiles, read_hit_dict)


def bwa(readfiles, baitfile, basename, cpu):
    """Conduct BWA search of reads against the baitfile.
    Returns an error if the second line of the baitfile contains characters other than ACTGN"""
    dna = set("ATCGN")
    if os.path.isfile(baitfile):
        # Quick detection of whether baitfile is DNA.
        with open(baitfile) as bf:
            header = bf.readline()
            seqline = bf.readline().rstrip().upper()
            if set(seqline) - dna:
                print(
                    "ERROR: characters other than ACTGN found in first line. You need a nucleotide bait file for BWA!"
                )
                return None

        if os.path.isfile(os.path.split(baitfile)[0] + ".amb"):
            db_file = baitfile
        else:
            print("Making nucleotide bwa index in current directory.")
            baitfile_dir = os.path.split(baitfile)[0]
            if baitfile_dir:
                if os.path.realpath(baitfile_dir) != os.path.realpath("."):
                    shutil.copy(baitfile, ".")
            db_file = os.path.split(baitfile)[1]
            make_bwa_index_cmd = "bwa index {}".format(db_file)
            print(("[CMD]: {}".format(make_bwa_index_cmd)))
            exitcode = subprocess.call(make_bwa_index_cmd, shell=True)
            if exitcode:
                return None
    else:
        print(("ERROR: Cannot find baitfile at: {}".format(baitfile)))
        return None

    if not cpu:
        import multiprocessing

        cpu = multiprocessing.cpu_count()

    if len(readfiles) < 3:
        bwa_fastq = " ".join(readfiles)
    else:
        bwa_fastq = readfiles

    bwa_commands = [
        "time bwa mem",
        "-t",
        str(cpu),
        db_file,
        bwa_fastq,
        " | samtools view -h -b -S - > ",
        basename + ".bam",
    ]
    full_command = " ".join(bwa_commands)
    print(("[CMD]: {}".format(full_command)))
    exitcode = subprocess.call(full_command, shell=True)
    if exitcode:
        return None

    return basename + ".bam"


def distribute_bwa(bamfile, readfiles):
    """"""
    try:
        distribute_reads_to_targets_bwa(bamfile, readfiles)
    except:
        print(
            "ERROR: Something went wrong with distributing reads to gene directories."
        )
        exitcode = True
        return exitcode
