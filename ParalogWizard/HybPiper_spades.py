#!/usr/bin/env python

import os
import shutil
import subprocess
import sys
from typing import List

helptext = """Run the assembler SPAdes with re-dos if any of the k-mers are unsuccessful.
The re-runs are attempted by removing the largest k-mer and re-running spades. If a final
contigs.fasta file is generated, a 'spades.ok' file is saved."""


def make_spades_cmd(
    genelist,
    cpu,
):

    parallel_cmd_list = ["time", "parallel", "--eta"]
    if cpu:
        parallel_cmd_list.append("-j {}".format(cpu))

    spades_cmd_list = (
        "spades.py --only-assembler --threads 1 --12 {{}}/{{}}_interleaved.fasta -o {{}}/{{}}_spades "
        ":::: {} > spades.log".format(genelist)
    )

    spades_cmd = " ".join(parallel_cmd_list) + " " + spades_cmd_list
    return spades_cmd


def spades_initial(
    genelist,
    cpu,
):
    """Run SPAdes on each gene separately using GNU paralell."""
    if os.path.isfile("spades.log"):
        os.remove("spades.log")

    genes = [x.rstrip() for x in open(genelist)]
    spades_cmd = make_spades_cmd(
        genelist,
        cpu,
    )

    sys.stderr.write("Running SPAdes on {} genes\n".format(len(genes)))
    sys.stderr.write(spades_cmd + "\n")
    exitcode = subprocess.call(spades_cmd, shell=True)

    if exitcode:
        sys.stderr.write(
            "ERROR: One or more genes had an error with SPAdes assembly. This may be due to low coverage. No contigs "
            "found for the following genes:\n "
        )

    spades_successful = []
    spades_failed = []

    for gene in genes:
        gene_failed = False
        if os.path.isfile("{}/{}_spades/contigs.fasta".format(gene, gene)):
            contig_file_size = os.stat(
                "{}/{}_spades/contigs.fasta".format(gene, gene)
            ).st_size
            if contig_file_size > 0:
                shutil.copy(
                    "{}/{}_spades/contigs.fasta".format(gene, gene),
                    "{}/{}_contigs.fasta".format(gene, gene),
                )
                spades_successful.append(gene)
            else:
                gene_failed = True
        else:
            gene_failed = True

        if gene_failed:
            sys.stderr.write("{}\n".format(gene))
            spades_failed.append(gene)
    return spades_failed


def rerun_spades(genelist, cpu):
    genes = [x.rstrip() for x in open(genelist)]
    spades_duds = []
    genes_redos = []
    with open("redo_spades_commands.txt", "w") as redo_cmds_file:
        for gene in genes:
            all_kmers = [
                int(x[1:])
                for x in os.listdir(os.path.join(gene, "{}_spades".format(gene)))
                if x.startswith("K")
            ]
            all_kmers.sort()

            if len(all_kmers) < 2:
                sys.stderr.write("WARNING: All Kmers failed for {}!\n".format(gene))
                spades_duds.append(gene)
                continue
            else:
                genes_redos.append(gene)
            redo_kmers = [str(x) for x in all_kmers[:-1]]
            restart_k = "k{}".format(redo_kmers[-1])
            kvals = ",".join(redo_kmers)
            spades_cmd = "spades.py --restart-from {} -k {} -o {}/{}_spades".format(
                restart_k, kvals, gene, gene
            )
            redo_cmds_file.write(spades_cmd + "\n")
    if cpu:
        redo_spades_cmd = "parallel -j {} --eta --timeout 400% :::: redo_spades_commands.txt > spades_redo.log".format(
            cpu
        )
    else:
        redo_spades_cmd = "parallel --eta --timeout 400% :::: redo_spades_commands.txt > spades_redo.log"

    sys.stderr.write("Re-running SPAdes for {} genes\n".format(len(genes_redos)))
    sys.stderr.write(redo_spades_cmd + "\n")
    exitcode = subprocess.call(redo_spades_cmd, shell=True)

    if exitcode:
        sys.stderr.write(
            "ERROR: One or more genes had an error with SPAdes assembly. This may be due to low coverage. No contigs "
            "found for the following genes:\n "
        )

    spades_successful = []
    spades_failed = []
    for gene in genes_redos:
        gene_failed = False
        if os.path.isfile("{}/{}_spades/contigs.fasta".format(gene, gene)):
            if os.stat("{}/{}_spades/contigs.fasta".format(gene, gene)).st_size > 0:
                shutil.copy(
                    "{}/{}_spades/contigs.fasta".format(gene, gene),
                    "{}/{}_contigs.fasta".format(gene, gene),
                )
                spades_successful.append(gene)
            else:
                gene_failed = True
        else:
            gene_failed = True

        if gene_failed:
            sys.stderr.write("{}\n".format(gene))
            spades_duds.append(gene)
    with open("spades_duds.txt", "w") as spades_duds_file:
        spades_duds_file.write("\n".join(spades_duds))

    return spades_failed, spades_duds


def spades_runner(
    genelist,
    cpu,
):

    spades_failed = spades_initial(
        genelist,
        cpu=cpu,
    )

    if len(spades_failed) > 0:
        with open("failed_spades.txt", "w") as failed_spadefile:
            failed_spadefile.write("\n".join(spades_failed))

        spades_failed, spades_duds = rerun_spades(
            "failed_spades.txt",
            cpu=cpu,
        )
        if len(spades_failed) == 0:
            sys.stderr.write("All redos completed successfully!\n")
        else:
            sys.exit(1)


def execute_spades(
    genes,
    cpu,
):
    """Run SPAdes on each gene separately using GNU parallel."""
    with open("spades_genelist.txt", "w") as spadesfile:
        spadesfile.write("\n".join(genes) + "\n")

    if os.path.isfile("spades.log"):
        os.remove("spades.log")
    if os.path.isfile("spades_redo.log"):
        os.remove("spades_redo.log")

    try:
        spades_runner(
            "spades_genelist.txt",
            cpu=cpu,
        )
        if os.path.isfile("spades_duds.txt"):
            spades_duds = [x.rstrip() for x in open("spades_duds.txt")]
        else:
            spades_duds = []
        spades_genelist = []
        for gene in genes:
            if gene not in set(spades_duds):
                spades_genelist.append(gene)
        return spades_genelist
    except:
        sys.stderr.write(
            "WARNING: Something went wrong with the assemblies! Check for failed assemblies and re-run! \n"
        )
        return None


def spades(
    readfiles: List[str],
    genes: List[str],
    cpu,
):
    if len(readfiles) == 2:
        spades_genelist = execute_spades(
            genes,
            cpu=cpu,
        )

    else:
        print("ERROR: Please specify paired read files! Exiting!")
        return
    if not spades_genelist:
        print("ERROR: No genes had assembled contigs! Exiting!")
        return


def list_sub_dirs(parentdir):
    """Given a parent directory return a list of all subdirectories"""
    return next(os.walk(parentdir))[1]


def remove_spades():
    """In the current directory, remove the spades directory."""
    spades_dirs = [
        s for s in os.listdir(".") if s.endswith("spades") and os.path.isdir(s)
    ]
    for s in spades_dirs:
        shutil.rmtree(s)


def clean_up(sample):
    try:
        os.chdir(sample)
        gene_dirs = list_sub_dirs(".")
        print("Found {} gene directories".format(len(gene_dirs)))
    except OSError:
        print("Directory '{}' does not exist!".format(sample))
        sys.exit(1)
    for gene in gene_dirs:
        os.chdir(gene)
        remove_spades()
        os.chdir("..")
