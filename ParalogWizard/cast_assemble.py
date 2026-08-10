#!/usr/bin/env python
"""
Module for ParalogWizard cast_assemble command.

This command conducts a BWA search against a bait file, distributes the reads,
and runs SPAdes assembly.  Parallel SPAdes runs use multiprocessing.Pool.
"""

from __future__ import annotations

import logging
import multiprocessing
import os
import shutil
import subprocess
from contextlib import contextmanager
from functools import wraps
from typing import List

from ParalogWizard import worker_initializer, log_exceptions

# Get logger by name (will be configured by ParalogWizard.py)
logger = logging.getLogger("ParalogWizard")

try:
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
except ImportError:
    import logging
    logging.error("BioPython not found. Please install BioPython: pip install biopython")
    raise

# --------------------------------------------------------------------------- #
#  Helpers (using unified log_exceptions from ParalogWizard.__init__)        #
# --------------------------------------------------------------------------- #


@contextmanager
def change_dir(dest: str):
    old = os.getcwd()
    os.chdir(dest)
    try:
        yield
    finally:
        os.chdir(old)


# --------------------------------------------------------------------------- #
#  Read distribution helpers                                                  #
# --------------------------------------------------------------------------- #
@log_exceptions
def distribute_reads_bwa(bamfilename: str, readfiles: list[str]) -> None:
    """
    Distribute reads from BAM file to target directories based on BWA alignments.
    Creates interleaved FASTA files for each target gene.
    """
    if not os.path.exists(bamfilename):
        logger.error("BAM file %s does not exist", bamfilename)
        raise FileNotFoundError(f"BAM file {bamfilename} does not exist")
    
    for readfile in readfiles:
        if not os.path.exists(readfile):
            logger.error("Read file %s does not exist", readfile)
            raise FileNotFoundError(f"Read file {readfile} does not exist")

    logger.info("Parsing BAM file %s to extract read mappings", bamfilename)
    samtools_cmd = ["samtools", "view", "-F", "4", bamfilename]
    try:
        result = subprocess.run(
            samtools_cmd, capture_output=True, text=True, check=True
        )
    except subprocess.CalledProcessError as e:
        logger.error("samtools command failed for %s: %s", bamfilename, e.stderr)
        raise

    read_hit_dict = {}
    for line in result.stdout.splitlines():
        fields = line.split()
        if len(fields) < 3:
            logger.warning("Unexpected samtools output line: %s", line)
            continue
        read_id = fields[0]
        target = fields[2].split("-")[-1]
        read_hit_dict.setdefault(read_id, []).append(target)

    logger.info("Found %d unique reads with hits", len(read_hit_dict))

    logger.info("Processing paired read files: %s, %s", readfiles[0], readfiles[1])
    try:
        with open(readfiles[0]) as f1, open(readfiles[1]) as f2:
            iterator1 = FastqGeneralIterator(f1)
            iterator2 = FastqGeneralIterator(f2)

            reads_written = 0
            for id1_long, seq1, _ in iterator1:
                try:
                    id2_long, seq2, _ = next(iterator2)
                except StopIteration:
                    logger.error("Paired read file %s ended prematurely", readfiles[1])
                    raise ValueError("Paired read file ended prematurely")

                def _strip(x: str) -> str:
                    x = x.split()[0]
                    return x[:-2] if x.endswith(("/1", "/2")) else x

                id1 = _strip(id1_long)
                id2 = _strip(id2_long)

                chosen = read_hit_dict.get(id1) or read_hit_dict.get(id2)
                if not chosen:
                    continue

                for target in chosen:
                    os.makedirs(target, exist_ok=True)
                    out_path = os.path.join(target, f"{target}_interleaved.fasta")
                    try:
                        with open(out_path, "a") as out:
                            out.write(f">{id1}\n{seq1}\n>{id2}\n{seq2}\n")
                    except Exception as e:
                        logger.error("Error writing paired sequences to %s: %s", out_path, e)
                        raise

                reads_written += 1
                if reads_written % 1000 == 0:
                    logger.debug("Processed %d read pairs", reads_written)

    except Exception as e:
        logger.error("Error processing read files %s, %s: %s", readfiles[0], readfiles[1], e)
        raise

    logger.info("Read distribution completed. Processed %d read pairs", reads_written)



# --------------------------------------------------------------------------- #
#  BWA                                                                        #
# --------------------------------------------------------------------------- #
@log_exceptions
def bwa(readfiles: List[str], baitfile: str, basename: str, cpu: int) -> str:
    """
    Conduct BWA search of reads against the baitfile.
    Creates a BWA index if needed, then runs BWA mem to create a BAM file.
    """
    if not os.path.isfile(baitfile):
        logger.error("Baitfile '%s' not found", baitfile)
        raise FileNotFoundError(f"Baitfile '{baitfile}' not found")

    for readfile in readfiles:
        if not os.path.exists(readfile):
            logger.error("Read file %s does not exist", readfile)
            raise FileNotFoundError(f"Read file {readfile} does not exist")

    # Check for BWA index files
    baitfile_dir = os.path.dirname(baitfile)
    index_file = os.path.join(baitfile_dir, os.path.basename(baitfile) + ".amb")
    
    if not os.path.isfile(index_file):
        logger.info("Creating BWA index for %s", baitfile)
        if baitfile_dir and os.path.realpath(baitfile_dir) != os.path.realpath("."):
            shutil.copy(baitfile, ".")
            baitfile = os.path.basename(baitfile)
        
        bwa_index_cmd = ["bwa", "index", baitfile]
        logger.info("Running BWA index command: %s", " ".join(bwa_index_cmd))
        try:
            subprocess.run(bwa_index_cmd, capture_output=True, text=True, check=True)
        except subprocess.CalledProcessError as e:
            logger.error("BWA index creation failed for %s: %s", baitfile, e.stderr)
            raise
        logger.info("BWA index created successfully for %s", baitfile)

    if not cpu:
        cpu = multiprocessing.cpu_count()

    output_bam = f"{basename}.bam"
    logger.info("Running BWA mem with %d threads for %d read files", cpu, len(readfiles))
    
    # Build BWA mem command
    bwa_mem_cmd = ["bwa", "mem", "-t", str(cpu), baitfile] + readfiles
    samtools_cmd = ["samtools", "view", "-h", "-b", "-S"]
    
    logger.info("BWA mem command: %s", " ".join(bwa_mem_cmd))
    logger.info("Samtools command: %s", " ".join(samtools_cmd))
    
    try:
        # Run BWA mem and pipe to samtools
        with open(output_bam, "wb") as bam_out:
            bwa_process = subprocess.Popen(
                bwa_mem_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=False
            )
            samtools_process = subprocess.Popen(
                samtools_cmd, stdin=bwa_process.stdout, stdout=bam_out, 
                stderr=subprocess.PIPE, text=False
            )
            bwa_process.stdout.close()  # Allow bwa_process to receive a SIGPIPE if samtools_process exits
            
            # Wait for both processes to complete
            bwa_stderr = bwa_process.communicate()[1]
            samtools_stderr = samtools_process.communicate()[1]
            
            if bwa_process.returncode != 0:
                logger.error("BWA mem failed: %s", bwa_stderr.decode())
                raise subprocess.CalledProcessError(bwa_process.returncode, bwa_mem_cmd)
            
            if samtools_process.returncode != 0:
                logger.error("Samtools view failed: %s", samtools_stderr.decode())
                raise subprocess.CalledProcessError(samtools_process.returncode, samtools_cmd)
                
    except Exception as e:
        logger.error("BWA mapping failed for %s: %s", basename, e)
        raise

    if not os.path.exists(output_bam) or os.path.getsize(output_bam) == 0:
        logger.error("Output BAM file %s was not created or is empty", output_bam)
        raise RuntimeError(f"Output BAM file {output_bam} was not created or is empty")

    logger.info("BWA mapping completed successfully. Output: %s", output_bam)
    return output_bam


# --------------------------------------------------------------------------- #
#  SPAdes helpers                                                             #
# --------------------------------------------------------------------------- #
@log_exceptions
def run_spades_for_gene(gene: str) -> None:
    """
    Run SPAdes assembly for a single gene.
    Returns gene name (success/failure determined later by checking output files).
    """
    input_file = os.path.join(gene, f"{gene}_interleaved.fasta")
    output_dir = os.path.join(gene, f"{gene}_spades")
    
    if not os.path.exists(input_file):
        logger.error("Input file %s does not exist for SPAdes assembly", input_file)
    
    cmd = [
        "spades.py", "--only-assembler", "--threads", "1", 
        "--12", input_file, "-o", output_dir
    ]
    
    logger.debug("Running SPAdes for gene %s: %s", gene, " ".join(cmd))
    try:
        subprocess.run(
            cmd, capture_output=True, text=True, check=True
        )
        logger.debug("SPAdes subprocess completed for gene %s", gene)
            
    except subprocess.CalledProcessError as e:
        logger.error("SPAdes assembly failed for gene %s: %s", gene, e.stderr)
    except Exception as e:
        logger.error("Unexpected error during SPAdes assembly for gene %s: %s", gene, e)
    


@log_exceptions
def run_spades_redo_for_gene(gene: str) -> str:
    """
    Run SPAdes redo assembly for a single gene with reduced k-mer set.
    Returns gene name (success/failure determined later by checking output files).
    """
    sp_dir = os.path.join(gene, f"{gene}_spades")
    
    # Get available kmers from previous run
    try:
        kmers = sorted(int(x[1:]) for x in os.listdir(sp_dir) if x.startswith("K"))
    except Exception as e:
        logger.error("Cannot read k-mer directories for gene %s: %s", gene, e)
        return gene
    
    if len(kmers) < 2:
        logger.error("Not enough k-mers for redo on gene %s", gene)
        return gene
    
    # Build redo command with reduced k-mer set
    redo_kmers = ",".join(map(str, kmers[:-1]))
    restart_k = f"k{kmers[-2]}"
    
    cmd = [
        "spades.py", "--restart-from", restart_k, 
        "-k", redo_kmers, "-o", sp_dir
    ]
    
    logger.debug("Running SPAdes redo for gene %s: %s", gene, " ".join(cmd))
    try:
        subprocess.run(
            cmd, capture_output=True, text=True, check=True
        )
        logger.debug("SPAdes redo subprocess completed for gene %s", gene)
    except subprocess.CalledProcessError as e:
        logger.error("SPAdes redo failed for gene %s: %s", gene, e.stderr)
    except Exception as e:
        logger.error("Unexpected error during SPAdes redo for gene %s: %s", gene, e)
    
    return gene


@log_exceptions
def spades_initial(genelist: str, cpu: int, log_queue) -> List[str]:
    """
    Run initial SPAdes assembly on all genes using multiprocessing.
    Returns list of genes that failed assembly.
    """
    if not os.path.exists(genelist):
        logger.error("Gene list file %s does not exist", genelist)
        raise FileNotFoundError(f"Gene list file {genelist} does not exist")
    
    # Read gene list
    with open(genelist) as fh:
        genes = [ln.strip() for ln in fh if ln.strip()]
    
    logger.info("Starting SPAdes assembly for %d genes using %d processes", len(genes), cpu)
    
    pool = None
    try:
        pool = multiprocessing.Pool(
        processes=cpu, initializer=worker_initializer, initargs=(log_queue,)
        )
        
        # Use map_async for better error handling
        async_result = pool.map_async(run_spades_for_gene, genes)
        async_result.get()  # Wait for completion without timeout
        
        logger.info("SPAdes subprocess calls completed, checking results...")
            
    except Exception as e:
        logger.error("Error during SPAdes initial run: %s", e)
        raise
    finally:
        if pool:
            pool.close()
            pool.join()
    
    # Check results following HybPiper logic: examine actual output files
    spades_successful, spades_failed = _check_and_copy_contigs(genes)
    
    logger.info("SPAdes initial results: %d successful, %d failed", 
               len(spades_successful), len(spades_failed))
    
    if spades_failed:
        logger.warning("Failed genes: %s", ", ".join(spades_failed))
    
    return spades_failed


@log_exceptions
def rerun_spades(genelist: str, cpu: int, log_queue):
    """
    Rerun failed SPAdes assemblies with reduced k-mer sets.
    Returns tuple of (still_failed_genes, dud_genes).
    """
    if not os.path.exists(genelist):
        logger.error("Failed genes list file %s does not exist", genelist)
        raise FileNotFoundError(f"Failed genes list file {genelist} does not exist")
    
    # Read gene list
    with open(genelist) as fh:
        genes = [ln.strip() for ln in fh if ln.strip()]
    
    logger.info("Preparing to rerun SPAdes for %d failed genes", len(genes))
    
    # Filter genes that can be redone (have enough k-mers)
    redo_genes = []
    spades_failed_initial = []
    
    for gene in genes:
        sp_dir = os.path.join(gene, f"{gene}_spades")
        try:
            kmers = sorted(int(x[1:]) for x in os.listdir(sp_dir) if x.startswith("K"))
            if len(kmers) >= 2:
                redo_genes.append(gene)
            else:
                spades_failed_initial.append(gene)
        except Exception:
            spades_failed_initial.append(gene)
    
    logger.info("Can redo %d genes, %d genes failed initially", len(redo_genes), len(spades_failed_initial))
    
    if not redo_genes:
        logger.warning("No genes can be redone - all genes failed initially")
        # Write failed genes list
        with open("spades_failed_final.txt", "w") as fh:
            fh.write("\n".join(spades_failed_initial))
        return genes, spades_failed_initial
    
    pool = None
    try:
        pool = multiprocessing.Pool(
        processes=cpu, initializer=worker_initializer, initargs=(log_queue,)
        )
        
        # Use map_async for better error handling
        async_result = pool.map_async(run_spades_redo_for_gene, redo_genes)
        async_result.get()  # Wait for completion without timeout
        
        logger.info("SPAdes redo subprocess calls completed")
        
    except Exception as e:
        logger.error("Error during SPAdes redo run: %s", e)
        raise
    finally:
        if pool:
            pool.close()
            pool.join()

    # Check redo results following HybPiper logic: examine actual output files
    spades_successful, redo_spades_failed = _check_and_copy_contigs(redo_genes)

    logger.info("SPAdes redo results: %d successful, %d failed", 
               len(spades_successful), len(redo_spades_failed))

    # Combine failed genes from initial failure analysis and redo failures
    all_spades_failed = list(set(spades_failed_initial + redo_spades_failed))
    
    # Write final failed genes list
    with open("spades_failed_final.txt", "w") as fh:
        fh.write("\n".join(all_spades_failed))
    logger.info("After redo: %d genes still failed, %d total failed", 
               len(redo_spades_failed), len(all_spades_failed))
    
    return redo_spades_failed, all_spades_failed


# === utility mini-helpers (keep top-level for pickling) ======================
def _check_and_copy_contigs(genes: List[str]) -> tuple[List[str], List[str]]:
    """Check SPAdes results and copy contigs files. Returns (successful, failed) gene lists."""
    spades_successful = []
    spades_failed = []
    
    for gene in genes:
        gene_failed = False
        contigs_path = os.path.join(gene, f"{gene}_spades", "contigs.fasta")
        
        if os.path.isfile(contigs_path):
            contig_file_size = os.path.getsize(contigs_path)
            if contig_file_size > 0:
                # Copy contigs file as in HybPiper
                target_path = os.path.join(gene, f"{gene}_contigs.fasta")
                try:
                    shutil.copy(contigs_path, target_path)
                    spades_successful.append(gene)
                    logger.debug("SPAdes successful for gene %s, copied contigs", gene)
                except Exception as e:
                    logger.error("Failed to copy contigs for gene %s: %s", gene, e)
                    gene_failed = True
            else:
                gene_failed = True
                logger.debug("SPAdes produced empty contigs file for gene %s", gene)
        else:
            gene_failed = True
            logger.debug("SPAdes did not produce contigs file for gene %s", gene)
        
        if gene_failed:
            spades_failed.append(gene)
    
    return spades_successful, spades_failed





# --------------------------------------------------------------------------- #
#  Main SPAdes controller called by top-level script                          #
# --------------------------------------------------------------------------- #
@log_exceptions
def run_spades_for_genes(genes: List[str], cpu: int, log_queue, output_dir: str = None) -> List[str]:
    """
    Main controller for running SPAdes assembly on multiple genes.
    Returns list of successfully assembled genes.
    """
    if not genes:
        logger.error("No genes provided for SPAdes assembly")
        raise ValueError("No genes provided for SPAdes assembly")
    
    logger.info("Starting SPAdes assembly pipeline for %d genes", len(genes))
    
    # Determine output directory for SPAdes files
    if output_dir is None:
        output_dir = "."
    
    # Write initial gene list
    initial_file = os.path.join(output_dir, "spades_genes_initial.txt")
    with open(initial_file, "w") as fh:
        fh.write("\n".join(genes))

    # Run initial SPAdes attempt
    spades_failed_initial = spades_initial(initial_file, cpu, log_queue)
    
    # Handle failures with redo attempt
    if spades_failed_initial:
        logger.info("Running SPAdes redo for %d failed genes", len(spades_failed_initial))
        # Write failed genes list for redo
        failed_initial_file = os.path.join(output_dir, "spades_failed_initial.txt")
        with open(failed_initial_file, "w") as fh:
            fh.write("\n".join(spades_failed_initial))
        
        spades_failed_after_redo, spades_failed_final = rerun_spades(failed_initial_file, cpu, log_queue)
        
        if spades_failed_after_redo:
            logger.error("SPAdes failed after redo for %d genes: %s", 
                        len(spades_failed_after_redo), ", ".join(spades_failed_after_redo))
            # Write final failures after redo
            failed_final_file = os.path.join(output_dir, "spades_failed_final.txt")
            with open(failed_final_file, "w") as fh:
                fh.write("\n".join(spades_failed_after_redo) + "\n")
        else:
            # Remove old file if no failures after redo
            failed_final_file = os.path.join(output_dir, "spades_failed_final.txt")
            if os.path.exists(failed_final_file):
                os.remove(failed_final_file)
    else:
        # No initial failures, remove any existing failed files
        failed_final_file = os.path.join(output_dir, "spades_failed_final.txt")
        if os.path.exists(failed_final_file):
            os.remove(failed_final_file)

    # Determine successfully assembled genes
    failed_final_file = os.path.join(output_dir, "spades_failed_final.txt")
    if os.path.isfile(failed_final_file):
        with open(failed_final_file) as fh:
            spades_failed_final = [ln.strip() for ln in fh if ln.strip()]
    else:
        spades_failed_final = []
    
    spades_genelist = [g for g in genes if g not in set(spades_failed_final)]
    
    # Write successful genes list
    successful_file = os.path.join(output_dir, "spades_successful.txt")
    with open(successful_file, "w") as fh:
        fh.write("\n".join(spades_genelist))
    
    if not spades_genelist:
        logger.error("No genes produced assembled contigs! All %d genes failed", len(genes))
        raise RuntimeError("No genes produced assembled contigs!")
    
    logger.info("SPAdes assembly completed: %d successful, %d failed", 
               len(spades_genelist), len(genes) - len(spades_genelist))
    logger.info("SPAdes output files saved to: %s", output_dir)
    return spades_genelist


@log_exceptions
def spades(readfiles: List[str], genes: List[str], cpu: int, log_queue, output_dir: str = None) -> List[str]:
    """
    Main SPAdes function for assembling genes from paired read files.
    Returns list of successfully assembled genes.
    """
    if len(readfiles) != 2:
        logger.error("Expected exactly 2 paired read files, got %d", len(readfiles))
        raise ValueError("Please specify exactly two paired read files!")
    
    for readfile in readfiles:
        if not os.path.exists(readfile):
            logger.error("Read file %s does not exist", readfile)
            raise FileNotFoundError(f"Read file {readfile} does not exist")

    logger.info("Starting SPAdes assembly with read files: %s", ", ".join(readfiles))
    return run_spades_for_genes(genes, cpu, log_queue, output_dir)


# --------------------------------------------------------------------------- #
#  Clean-up helpers                                                           #
# --------------------------------------------------------------------------- #
@log_exceptions
def remove_spades():
    for d in os.listdir("."):
        if d.endswith("spades") and os.path.isdir(d):
            shutil.rmtree(d)
            logger.info("Removed directory %s", d)


@log_exceptions
def clean_up(sample: str):
    for gene in [g for g in os.listdir(".") if os.path.isdir(g)]:
        with change_dir(gene):
            for d in os.listdir("."):
                if d.endswith("spades") and os.path.isdir(d):
                    shutil.rmtree(d)
                    logger.info("Removed %s in gene %s", d, gene)
