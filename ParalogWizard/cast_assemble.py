#!/usr/bin/env python
"""
Module for ParalogWizard cast_assemble command.

This command conducts a BWA search against a bait file, distributes the reads,
and runs SPAdes assembly. This version uses multiprocessing (parallel run only)
to run SPAdes. It uses Pool as a context manager and passes a worker initializer
and a log_queue.
"""

import errno
import os
import sys
import subprocess
import shutil
from contextlib import contextmanager
from functools import partial, wraps
import multiprocessing

from ParalogWizard import setup_logging, worker_initializer

# -----------------------------------------------------------------------------
# Module-level logger configured once.
# -----------------------------------------------------------------------------
logger = setup_logging()


# =============================================================================
# Logging Decorator
# =============================================================================
def log_exceptions(func):
    """
    Decorator that logs function entry, exit, and any exceptions.
    Exits with code 1 upon an exception.
    """
    from functools import wraps

    @wraps(func)
    def wrapper(*args, **kwargs):
        # logger.debug(f"Entering {func.__name__}")
        try:
            result = func(*args, **kwargs)
            # logger.debug(f"Exiting {func.__name__}")
            return result
        except Exception as e:
            logger.exception(f"Exception in {func.__name__}: {e}")
            sys.exit(1)

    return wrapper


# =============================================================================
# Context Manager and Utilities
# =============================================================================
@contextmanager
def change_dir(destination: str):
    """Temporarily change the working directory to 'destination'."""
    old_dir = os.getcwd()
    os.chdir(destination)
    try:
        yield
    finally:
        os.chdir(old_dir)


@log_exceptions
def mkdir_p(path: str) -> None:
    """Create a directory (and all intermediate directories) like 'mkdir -p'."""
    try:
        os.makedirs(path, exist_ok=True)
    except Exception as exc:
        logger.error("Error creating directory '%s': %s", path, exc)
        raise


# -----------------------------------------------------------------------------
# BAM Parsing & Read Distribution
# -----------------------------------------------------------------------------
@log_exceptions
def distribute_reads_bwa(bamfilename: str, readfiles: list) -> None:
    """
    Run samtools on the BAM file to map read IDs to target names, then distribute
    paired reads from two FASTQ files accordingly.
    """
    samtools_cmd = f"samtools view -F 4 {bamfilename}"
    try:
        result = subprocess.run(
            samtools_cmd, shell=True, capture_output=True, text=True, check=True
        )
    except subprocess.CalledProcessError as e:
        logger.error("samtools command failed: %s", e)
        sys.exit(1)
    bwa_results = result.stdout.splitlines()
    read_hit_dict = {}
    for line in bwa_results:
        fields = line.split()
        if len(fields) < 3:
            logger.warning("Unexpected samtools output line: %s", line)
            continue
        read_id = fields[0]
        target = fields[2].split("-")[-1]
        read_hit_dict.setdefault(read_id, [])
        if target not in read_hit_dict[read_id]:
            read_hit_dict[read_id].append(target)
    try:
        f1 = open(readfiles[0])
        f2 = open(readfiles[1])
    except Exception as e:
        logger.error("Error opening read files: %s", e)
        sys.exit(1)
    from Bio.SeqIO.QualityIO import FastqGeneralIterator

    iterator1 = FastqGeneralIterator(f1)
    iterator2 = FastqGeneralIterator(f2)
    reads_written = 0
    total_reads = len(read_hit_dict)
    logger.info("Starting read distribution for %d unique reads.", total_reads)
    progress_logged = set()
    for id1_long, seq1, _ in iterator1:
        try:
            id2_long, seq2, _ = next(iterator2)
        except StopIteration:
            logger.error("Paired read file ended prematurely.")
            break
        id1 = id1_long.split()[0]
        if id1.endswith("/1") or id1.endswith("/2"):
            id1 = id1[:-2]
        id2 = id2_long.split()[0]
        if id2.endswith("/1") or id2.endswith("/2"):
            id2 = id2[:-2]
        if id1 in read_hit_dict:
            for target in read_hit_dict[id1]:
                # Write paired sequences to target's interleaved FASTA file.
                mkdir_p(target)
                out_path = os.path.join(target, f"{target}_interleaved.fasta")
                try:
                    with open(out_path, "a") as outfile:
                        outfile.write(f">{id1}\n{seq1}\n")
                        outfile.write(f">{id2}\n{seq2}\n")
                except Exception as e:
                    logger.error(
                        "Error writing paired sequences to %s: %s", out_path, e
                    )
                    sys.exit(1)
            reads_written += 1
        elif id2 in read_hit_dict:
            for target in read_hit_dict[id2]:
                mkdir_p(target)
                out_path = os.path.join(target, f"{target}_interleaved.fasta")
                try:
                    with open(out_path, "a") as outfile:
                        outfile.write(f">{id1}\n{seq1}\n")
                        outfile.write(f">{id2}\n{seq2}\n")
                except Exception as e:
                    logger.error(
                        "Error writing paired sequences to %s: %s", out_path, e
                    )
                    sys.exit(1)
            reads_written += 1
        progress = int(((reads_written + 1) / total_reads) * 100)
        if progress % 5 == 0 and progress not in progress_logged:
            logger.info("Distribution progress: %d%%", progress)
            progress_logged.add(progress)
    f1.close()
    f2.close()
    logger.info(
        "Read distribution completed. Total reads distributed: %d", reads_written
    )


# =============================================================================
# BWA and SPAdes Functions
# =============================================================================
@log_exceptions
def bwa(readfiles: list, baitfile: str, basename: str, cpu: int) -> str:
    """
    Conduct a BWA search of reads against the baitfile.
    Returns the path to the resulting BAM file on success, or an empty string on error.
    """
    dna = set("ATCGN")
    if not os.path.isfile(baitfile):
        logger.error("Baitfile '%s' not found.", baitfile)
        return ""
    try:
        with open(baitfile) as bf:
            _ = bf.readline()  # header
            seqline = bf.readline().rstrip().upper()
    except Exception as e:
        logger.error("Error reading baitfile '%s': %s", baitfile, e)
        return ""
    if set(seqline) - dna:
        logger.error(
            "Baitfile '%s' contains invalid characters. Only ACTGN allowed.", baitfile
        )
        return ""
    baitfile_dir = os.path.dirname(baitfile)
    index_file = os.path.join(baitfile_dir, os.path.basename(baitfile) + ".amb")
    if os.path.isfile(index_file):
        db_file = baitfile
    else:
        logger.info("Making nucleotide BWA index in current directory.")
        if baitfile_dir and os.path.realpath(baitfile_dir) != os.path.realpath("."):
            try:
                shutil.copy(baitfile, ".")
            except Exception as e:
                logger.error("Error copying baitfile: %s", e)
                return ""
        db_file = os.path.basename(baitfile)
        make_bwa_index_cmd = f"bwa index {db_file}"
        logger.info("[CMD]: %s", make_bwa_index_cmd)
        exitcode = subprocess.call(make_bwa_index_cmd, shell=True)
        if exitcode:
            logger.error("BWA index creation failed with exit code %d", exitcode)
            return ""
    if not cpu:
        try:
            cpu = multiprocessing.cpu_count()
        except Exception as e:
            logger.error("Error determining CPU count: %s", e)
            return ""
    bwa_fastq = " ".join(readfiles)
    bwa_cmd = (
        f"time bwa mem -t {cpu} {db_file} {bwa_fastq} | "
        f"samtools view -h -b -S > {basename}.bam"
    )
    logger.info("[CMD]: %s", bwa_cmd)
    exitcode = subprocess.call(bwa_cmd, shell=True)
    if exitcode:
        logger.error("BWA command failed with exit code %d", exitcode)
        return ""
    return f"{basename}.bam"


@log_exceptions
def distribute_bwa(bamfile: str, readfiles: list) -> int:
    """
    Distribute reads from the given BAM file to target directories.
    Returns 0 on success or 1 on error.
    """
    try:
        distribute_reads_bwa(bamfile, readfiles)
    except Exception as e:
        logger.error("Error distributing reads to gene directories: %s", e)
        return 1
    return 0


# -----------------------------------------------------------------------------
# SPAdes Assembly Functions
# -----------------------------------------------------------------------------
@log_exceptions
def run_spades_for_gene(gene: str) -> tuple:
    """
    Build and run the SPAdes command for a single gene.
    Returns a tuple (gene, True) on success and (gene, False) on failure.
    """
    cmd = (
        f"spades.py --only-assembler --threads 1 --12 {gene}/{gene}_interleaved.fasta "
        f"-o {gene}/{gene}_spades > spades.log"
    )
    logger.info("Running SPAdes for gene %s with command: %s", gene, cmd)
    exitcode = subprocess.call(cmd, shell=True)
    if exitcode:
        logger.error(
            "SPAdes assembly for gene %s failed with exit code %d", gene, exitcode
        )
        return gene, False
    return gene, True


@log_exceptions
def spades_initial(genelist: str, cpu: int, log_queue) -> list:
    """
    Run SPAdes on each gene in parallel using multiprocessing.
    Returns a list of genes that failed assembly.
    """
    if os.path.isfile("spades.log"):
        try:
            os.remove("spades.log")
        except Exception as e:
            logger.error("Error removing spades.log: %s", e)
    try:
        with open(genelist) as f:
            genes = [line.rstrip() for line in f if line.strip()]
    except Exception as e:
        logger.error("Error reading genelist file '%s': %s", genelist, e)
        sys.exit(1)
    with multiprocessing.Pool(
        processes=cpu, initializer=worker_initializer, initargs=(log_queue,)
    ) as pool:
        results = pool.map(run_spades_for_gene, genes)
    spades_failed = [gene for gene, success in results if not success]
    for gene in genes:
        gene_failed = False
        contigs_path = os.path.join(gene, f"{gene}_spades", "contigs.fasta")
        if os.path.isfile(contigs_path):
            try:
                contig_file_size = os.stat(contigs_path).st_size
            except Exception as e:
                logger.error("Error getting file size for %s: %s", contigs_path, e)
                contig_file_size = 0
            if contig_file_size > 0:
                dest = os.path.join(gene, f"{gene}_contigs.fasta")
                try:
                    shutil.copy(contigs_path, dest)
                except Exception as e:
                    logger.error("Error copying contigs for gene %s: %s", gene, e)
                    gene_failed = True
            else:
                gene_failed = True
        else:
            gene_failed = True
        if gene_failed and gene not in spades_failed:
            logger.error("SPAdes assembly failed for gene: %s", gene)
            spades_failed.append(gene)
    return spades_failed


@log_exceptions
def rerun_spades(genelist: str, cpu: int, log_queue) -> tuple:
    """
    Attempt to re-run SPAdes on the failed genes using alternative k-mer values.
    Returns a tuple: (spades_failed, spades_duds).
    """
    try:
        with open(genelist) as f:
            genes = [line.rstrip() for line in f if line.strip()]
    except Exception as e:
        logger.error("Error reading genelist file '%s': %s", genelist, e)
        sys.exit(1)
    spades_duds = []
    genes_redos = []
    redo_commands = []
    try:
        with open("redo_spades_commands.txt", "w") as redo_cmds_file:
            for gene in genes:
                spades_dir = os.path.join(gene, f"{gene}_spades")
                try:
                    kmers = [
                        int(x[1:]) for x in os.listdir(spades_dir) if x.startswith("K")
                    ]
                except Exception as e:
                    logger.error("Error listing kmers in %s: %s", spades_dir, e)
                    spades_duds.append(gene)
                    continue
                kmers.sort()
                if len(kmers) < 2:
                    logger.warning("All kmers failed for gene %s!", gene)
                    spades_duds.append(gene)
                    continue
                genes_redos.append(gene)
                redo_kmers = [str(x) for x in kmers[:-1]]
                restart_k = f"k{redo_kmers[-1]}"
                kvals = ",".join(redo_kmers)
                spades_cmd = f"spades.py --restart-from {restart_k} -k {kvals} -o {gene}/{gene}_spades"
                redo_cmds_file.write(spades_cmd + "\n")
                redo_commands.append(spades_cmd)
    except Exception as e:
        logger.error("Error writing redo_spades_commands.txt: %s", e)
        sys.exit(1)
    with multiprocessing.Pool(
        processes=cpu, initializer=worker_initializer, initargs=(log_queue,)
    ) as pool:
        results = pool.map(
            lambda cmd: subprocess.call(cmd, shell=True) == 0, redo_commands
        )
    if not all(results):
        logger.error("One or more SPAdes redo commands failed.")
    spades_failed = []
    for gene in genes_redos:
        gene_failed = False
        contigs_path = os.path.join(gene, f"{gene}_spades", "contigs.fasta")
        if os.path.isfile(contigs_path):
            try:
                if os.stat(contigs_path).st_size > 0:
                    dest = os.path.join(gene, f"{gene}_contigs.fasta")
                    shutil.copy(contigs_path, dest)
                else:
                    gene_failed = True
            except Exception as e:
                logger.error("Error processing contigs for gene %s: %s", gene, e)
                gene_failed = True
        else:
            gene_failed = True
        if gene_failed:
            logger.error("SPAdes assembly failed for gene after redo: %s", gene)
            spades_failed.append(gene)
    try:
        with open("spades_duds.txt", "w") as spades_duds_file:
            spades_duds_file.write("\n".join(spades_duds))
    except Exception as e:
        logger.error("Error writing spades_duds.txt: %s", e)
    return spades_failed, spades_duds


@log_exceptions
def run_spades_for_genes(genes: list, cpu: int, log_queue) -> list:
    """
    Write the gene list to file, run SPAdes on all genes (with a redo for failures),
    and return a list of genes that assembled successfully.
    """
    try:
        with open("spades_genelist.txt", "w") as f:
            f.write("\n".join(genes) + "\n")
    except Exception as e:
        logger.error("Error writing spades_genelist.txt: %s", e)
        sys.exit(1)
    # Remove old log files.
    for logfile in ["spades.log", "spades_redo.log"]:
        if os.path.isfile(logfile):
            try:
                os.remove(logfile)
            except Exception as e:
                logger.warning("Could not remove %s: %s", logfile, e)
    spades_failed = spades_initial("spades_genelist.txt", cpu=cpu, log_queue=log_queue)
    if spades_failed:
        try:
            with open("failed_spades.txt", "w") as f:
                f.write("\n".join(spades_failed))
        except Exception as e:
            logger.error("Error writing failed_spades.txt: %s", e)
        spades_failed, spades_duds = rerun_spades(
            "failed_spades.txt", cpu=cpu, log_queue=log_queue
        )
        if spades_failed:
            logger.error(
                "SPAdes failed for the following genes after redo: %s", spades_failed
            )
            sys.exit(1)
        else:
            logger.info("All SPAdes re-run assemblies completed successfully!")
    spades_duds = []
    if os.path.isfile("spades_duds.txt"):
        try:
            with open("spades_duds.txt") as f:
                spades_duds = [line.rstrip() for line in f if line.strip()]
        except Exception as e:
            logger.error("Error reading spades_duds.txt: %s", e)
    successful_genes = [gene for gene in genes if gene not in set(spades_duds)]
    if not successful_genes:
        logger.error("No genes had assembled contigs!")
        sys.exit(1)
    return successful_genes


@log_exceptions
def spades(readfiles: list, genes: list, cpu: int, log_queue) -> None:
    """
    Run SPAdes assembly on genes.
    Expects exactly two paired read files.
    Exits if the file count is not 2 or if no assemblies were successful.
    """
    if len(readfiles) != 2:
        logger.error("Please specify exactly two paired read files! Exiting!")
        sys.exit(1)
    spades_genelist = run_spades_for_genes(genes, cpu=cpu, log_queue=log_queue)
    logger.info("SPAdes completed for %d genes.", len(spades_genelist))


# =============================================================================
# Cleanup and Directory Functions
# =============================================================================
@log_exceptions
def remove_spades() -> None:
    """Remove directories in the current directory that end with 'spades'."""
    for s in os.listdir("."):
        if s.endswith("spades") and os.path.isdir(s):
            try:
                shutil.rmtree(s)
                logger.info("Removed directory: %s", s)
            except Exception as e:
                logger.error("Error removing directory %s: %s", s, e)


@log_exceptions
def clean_up(sample: str) -> None:
    """
    For each gene directory (in the current directory), remove any subdirectories ending with 'spades'.
    """
    # List subdirectories in current directory.
    gene_dirs = [d for d in os.listdir(".") if os.path.isdir(d)]
    logger.info("Found %d gene directories in '%s'.", len(gene_dirs), sample)
    for gene in gene_dirs:
        try:
            with change_dir(gene):
                for d in os.listdir("."):
                    if d.endswith("spades") and os.path.isdir(d):
                        shutil.rmtree(d)
                        logger.info("Removed spades directory: %s in gene %s", d, gene)
        except Exception as e:
            logger.error("Error cleaning up gene directory '%s': %s", gene, e)
