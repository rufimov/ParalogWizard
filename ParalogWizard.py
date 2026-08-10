#!/usr/bin/env python3
"""
ParalogWizard – paralog-detection pipeline for HybSeq data.

Sub-commands:
  cast_assemble, cast_retrieve, cast_analyze, cast_detect, cast_separate,
  cast_remap, cast_call, cast_ploidy, cast_polyploid, cast_phase
"""

from __future__ import annotations

# -----------------------------------------------------------------------------#
#  Standard library                                                             #
# -----------------------------------------------------------------------------#
import argparse
import functools
import logging
import logging.handlers
import multiprocessing
import os
import signal
import shutil
import sys
import time
from contextlib import contextmanager
from datetime import datetime
from glob import glob
from typing import Callable

# -----------------------------------------------------------------------------#
#  Third-party                                                                 #
# -----------------------------------------------------------------------------#
import pandas as pd

# -----------------------------------------------------------------------------#
#  Own package                                                                 #
# -----------------------------------------------------------------------------#
from ParalogWizard import (  # noqa: E402  (keep local import style)
    compress_fastq_files,
    decompress_fastq_files,
    listener_process,
    log_paths_from_base,
    setup_logging,
)

# =============================================================================
#  Small helpers
# =============================================================================
@contextmanager
def change_dir(destination: str):
    """Temporarily change the working directory."""
    previous = os.getcwd()
    os.chdir(destination)
    try:
        yield
    finally:
        os.chdir(previous)


def fatal(msg: str):
    """
    Log an error and raise RuntimeError so the decorator/main can handle cleanup.

    Replaces all previous `sys.exit(1)` calls inside worker functions.
    """
    logger = logging.getLogger("ParalogWizard")
    logger.error(msg)
    raise RuntimeError(msg)


# =============================================================================
#  Decorator for sub-commands
# =============================================================================
def log_command(func: Callable):
    """
    Log entry/exit around each sub-command and ensure **every** exception
    (including SystemExit, KeyboardInterrupt) is logged then re-raised.
    """

    @functools.wraps(func)
    def wrapper(args, log_queue):
        logger = logging.getLogger("ParalogWizard")
        logger.info("Starting %s", func.__name__)
        try:
            result = func(args, log_queue)
            return result
        except BaseException as exc:  # catches everything
            logger.exception("Error in %s: %s", func.__name__, exc)
            raise

    return wrapper


# =============================================================================
#  Argument parsing
# =============================================================================
def setup_argparser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="ParalogWizard: a pipeline for paralog detection in HybSeq data."
    )

    # Shared logging flags for every subcommand.
    log_parent = argparse.ArgumentParser(add_help=False)
    log_parent.add_argument(
        "--debug",
        action="store_true",
        help=(
            "Enable debug logging: console shows DEBUG, and a *.debug.log file "
            "is written with full detail (INFO/WARNING still go to *.log; "
            "ERROR/CRITICAL to *.errors.log)"
        ),
    )

    # `required=True` for add_subparsers() exists only from Python 3.7 upward
    if sys.version_info >= (3, 7):
        subparsers = parser.add_subparsers(
            dest="command", required=True, help="Subcommand to run"
        )
    else:
        subparsers = parser.add_subparsers(dest="command", help="Subcommand to run")

    # ---- cast_assemble -------------------------------------------------------
    pa = subparsers.add_parser(
        "cast_assemble", parents=[log_parent], help="Run the cast_assemble step"
    )
    pa.add_argument(
        "-pr",
        "--baitfile",
        required=True,
        help=(
            "FASTA with bait sequences (if multiple baits per gene, IDs must look "
            "like >Taxon-geneName)"
        ),
    )
    pa.add_argument("-d", "--data_folder", required=True, help="Main data folder")
    pa.add_argument(
        "-nc", "--num_cores", default=1, type=int, help="Number of cores (default 1)"
    )

    # ---- cast_retrieve -------------------------------------------------------
    pr = subparsers.add_parser(
        "cast_retrieve", parents=[log_parent], help="Run the cast_retrieve step"
    )
    pr.add_argument("-pe", "--probes_exons", required=True, help="probes_exons file")
    pr.add_argument(
        "-l", "--length_cover", default=75, type=float, help="Length cover (default 75)"
    )
    pr.add_argument(
        "-s", "--spades_cover", default=5, type=float, help="Spades cover (default 5)"
    )
    pr.add_argument("-c", "--collect_contigs", action="store_true", help="Collect contigs")
    pr.add_argument("-d", "--data_folder", required=True, help="Main data folder")
    pr.add_argument(
        "-nc", "--num_cores", default=1, type=int, help="Number of cores (default 1)"
    )

    # ---- cast_analyze --------------------------------------------------------
    pan = subparsers.add_parser(
        "cast_analyze", parents=[log_parent], help="Run the cast_analyze step"
    )
    pan.add_argument(
        "-b",
        "--blocklist",
        nargs="+",
        help="Species to exclude from paralog divergence estimation",
    )
    pan.add_argument("-d", "--data_folder", required=True, help="Main data folder")
    pan.add_argument(
        "-nc", "--num_cores", default=1, type=int, help="Number of cores (default 1)"
    )

    # ---- cast_detect ---------------------------------------------------------
    pdct = subparsers.add_parser(
        "cast_detect", parents=[log_parent], help="Run the cast_detect step"
    )
    pdct.add_argument(
        "-b", "--blocklist", nargs="+", help="Species to exclude from divergence calc."
    )
    pdct.add_argument("-pe", "--probes_exons", required=True, help="probes_exons file")
    pdct.add_argument("-p", "--paralogs", action="store_true", help="Enable paralog detection")
    pdct.add_argument("-mi", "--minimum_divergence", type=float, help="Min divergence")
    pdct.add_argument("-ma", "--maximum_divergence", type=float, help="Max divergence")
    pdct.add_argument("-d", "--data_folder", required=True, help="Main data folder")
    pdct.add_argument(
        "-nc", "--num_cores", default=1, type=int, help="Number of cores (default 1)"
    )

    # ---- cast_separate -------------------------------------------------------
    ps = subparsers.add_parser(
        "cast_separate", parents=[log_parent], help="Run the cast_separate step"
    )
    ps.add_argument("-r", "--redlist", nargs="+", help="Taxa to exclude")
    ps.add_argument("-i", "--min_identity", required=True, type=float, help="BLAT id")
    ps.add_argument("-pc", "--probes_customized", required=True, help="probes_customized file")
    ps.add_argument("-d", "--data_folder", required=True, help="Main data folder")
    ps.add_argument("-nc", "--num_cores", default=1, type=int, help="Cores (default 1)")

    # ---- cast_extend (optional; module may be local-only / gitignored) --------
    pe = subparsers.add_parser(
        "cast_extend", parents=[log_parent], help="Run the cast_extend step (local)"
    )
    pe.add_argument("-pr", "--baitfile", required=True, help="bait FASTA")
    pe.add_argument("-d", "--data_folder", required=True, help="Main data folder")
    pe.add_argument("-nc", "--num_cores", default=1, type=int, help="Cores (default 1)")

    # ---- cast_remap ----------------------------------------------------------
    prm = subparsers.add_parser(
        "cast_remap", parents=[log_parent], help="Run the cast_remap step"
    )
    prm.add_argument("-pc", "--probes_customized", required=True, help="probes_customized file")
    prm.add_argument("-e", "--exon_length", type=int, default=150, help="Exon length")
    prm.add_argument("-d", "--data_folder", required=True, help="Main data folder")
    prm.add_argument("-nc", "--num_cores", default=1, type=int, help="Cores (default 1)")
    prm.add_argument(
        "--force",
        action="store_true",
        help="Redo mapping even when sorted BAMs already exist",
    )

    # ---- cast_call -----------------------------------------------------------
    pcl = subparsers.add_parser(
        "cast_call", parents=[log_parent], help="Run the cast_call step"
    )
    pcl.add_argument("-d", "--data_folder", required=True, help="Main data folder")
    pcl.add_argument("-nc", "--num_cores", default=1, type=int, help="Cores (default 1)")
    pcl.add_argument(
        "--force",
        action="store_true",
        help="Redo variant calling even when complete VCFs already exist",
    )

    # ---- cast_phase ----------------------------------------------------------
    pph = subparsers.add_parser(
        "cast_phase", parents=[log_parent], help="Run the cast_phase step"
    )
    pph.add_argument("-d", "--data_folder", required=True, help="Main data folder")
    pph.add_argument("-nc", "--num_cores", default=1, type=int, help="Cores (default 1)")

    # ---- cast_ploidy ---------------------------------------------------------
    ppl = subparsers.add_parser(
        "cast_ploidy", parents=[log_parent], help="Run the cast_ploidy step"
    )
    ppl.add_argument("-d", "--data_folder", required=True, help="Main data folder")
    ppl.add_argument("-nc", "--num_cores", default=1, type=int, help="Cores (default 1)")

    # ---- cast_polyploid ------------------------------------------------------
    ppp = subparsers.add_parser(
        "cast_polyploid", parents=[log_parent], help="Run the cast_polyploid step"
    )
    ppp.add_argument("-pc", "--probes_customized", required=True, help="probes_customized file")
    ppp.add_argument("-d", "--data_folder", required=True, help="Main data folder")
    ppp.add_argument("-nc", "--num_cores", default=1, type=int, help="Cores (default 1)")

    return parser


def parse_args():
    parser = setup_argparser()
    args = parser.parse_args()

    if args.command is None:  # Python < 3.7
        parser.error("A subcommand is required (cast_assemble, cast_call …).")

    return args


# =============================================================================
#  Sub-command implementations
# =============================================================================
@log_command
def run_cast_assemble(args, log_queue):
    from ParalogWizard.cast_assemble import bwa, clean_up, distribute_reads_bwa, spades

    logger = logging.getLogger("ParalogWizard")
    logger.info("Running cast_assemble with data_folder: %s", args.data_folder)

    data_folder = os.path.abspath(args.data_folder)
    os.makedirs(data_folder, exist_ok=True)
    baitfile = os.path.abspath(args.baitfile)
    if not os.path.isfile(baitfile):
        fatal(f"Bait file '{baitfile}' not found.")
    
    with change_dir(data_folder):
        samples_list = os.path.join("10deduplicated_reads", "samples_list.txt")
        if not os.path.isfile(samples_list):
            fatal(f"Samples list file '{samples_list}' not found.")

        with open(samples_list) as fh:
            samples = [ln.strip() for ln in fh if ln.strip()]
        if not samples:
            fatal(f"No samples found in '{samples_list}'.")

        os.makedirs("20assemblies", exist_ok=True)

        with change_dir("20assemblies"):
            for sample in samples:
                sample_dir = os.path.join(os.getcwd(), sample)
                os.makedirs(sample_dir, exist_ok=True)

                with change_dir(sample_dir):
                    pattern = os.path.join("..", "..", "10deduplicated_reads", sample) + ".*.fastq*"
                    readfiles = glob(pattern)
                    if len(readfiles) < 2:
                        logger.error("Not enough readfiles for sample %s (pattern %s)", sample, pattern)
                        continue

                    # two first files, absolute
                    readfiles = [
                        os.path.abspath(os.path.join("..", "..", "10deduplicated_reads", os.path.basename(readfiles[0]))),
                        os.path.abspath(os.path.join("..", "..", "10deduplicated_reads", os.path.basename(readfiles[1]))),
                    ]

                    # Check readfiles exist before proceeding
                    missing_files = [rf for rf in readfiles if not os.path.isfile(rf)]
                    if missing_files:
                        logger.error("Missing readfiles for sample %s: %s", sample, ", ".join(missing_files))
                        continue

                    logger.info("Processing sample %s", sample)
                    readfiles = decompress_fastq_files(readfiles)
                    
                    # These functions already handle their own errors and logging via @log_exceptions
                    bamfile = bwa(readfiles=readfiles, baitfile=baitfile, basename=sample, cpu=args.num_cores)
                    distribute_reads_bwa(bamfilename=bamfile, readfiles=readfiles)

                    genes = [
                        x for x in os.listdir(".")
                        if os.path.isfile(os.path.join(x, f"{x}_interleaved.fasta"))
                    ]
                    
                    if not genes:
                        logger.warning("No genes found with interleaved FASTA files for sample %s", sample)
                    else:
                        assembled_genes = spades(readfiles=readfiles, genes=genes, cpu=args.num_cores, log_queue=log_queue, output_dir=data_folder)
                        logger.info("SPAdes completed for sample %s: %d genes assembled", sample, len(assembled_genes) if assembled_genes else 0)

                    compress_fastq_files(readfiles)
                    clean_up(sample)
                    logger.info("Completed processing sample %s", sample)

    logger.info("ParalogWizard cast_assemble completed.")


@log_command
def run_cast_retrieve(args, log_queue):
    from ParalogWizard.cast_retrieve import retrieve

    logger = logging.getLogger("ParalogWizard")
    
    # Validate inputs
    data_folder = os.path.abspath(args.data_folder)
    if not os.path.isdir(data_folder):
        fatal(f"Data folder '{data_folder}' does not exist or is not a directory.")
    
    probes_exons = os.path.abspath(args.probes_exons)
    if not os.path.isfile(probes_exons):
        fatal(f"Probes exons file '{probes_exons}' not found.")

    logger.info("Running cast_retrieve with data_folder=%s, probes_exons=%s, cores=%d, length_cover=%.1f, spades_cover=%.1f", 
               data_folder, probes_exons, args.num_cores, args.length_cover, args.spades_cover)
    logger.info("Collect contigs: %s", "enabled" if args.collect_contigs else "disabled")
    
    retrieve(
        data_folder,
        args.collect_contigs,
        probes_exons,
        args.num_cores,
        args.length_cover,
        args.spades_cover,
        log_queue,
    )
    logger.info("ParalogWizard cast_retrieve completed successfully.")


@log_command
def run_cast_analyze(args, log_queue):
    from ParalogWizard.cast_analyze import build_alignments, estimate_divergence

    logger = logging.getLogger("ParalogWizard")
    
    # Validate inputs
    data_folder = os.path.abspath(args.data_folder)
    if not os.path.isdir(data_folder):
        fatal(f"Data folder '{data_folder}' does not exist or is not a directory.")
    
    # Check for required input directory
    required_dir = os.path.join(data_folder, "31exonic_contigs")
    if not os.path.isdir(required_dir):
        fatal(f"Required directory '{required_dir}' not found. Run cast_retrieve first.")
    
    blocklist = set(args.blocklist) if args.blocklist else set()

    if blocklist:
        logger.info(
            "Running cast_analyze: data_folder=%s, blocklist=%s, cores=%d",
            data_folder, ", ".join(sorted(blocklist)), args.num_cores
        )
    else:
        logger.info("Running cast_analyze: data_folder=%s, all species included, cores=%d", 
                   data_folder, args.num_cores)

    logger.info("Building alignments...")
    build_alignments(data_folder, args.num_cores, log_queue)
    logger.info("Alignment building completed successfully.")
    
    logger.info("Estimating divergence...")
    estimate_divergence(data_folder, blocklist, args.num_cores, log_queue)
    logger.info("Divergence estimation completed successfully.")
    
    logger.info("ParalogWizard cast_analyze completed successfully.")


@log_command
def run_cast_detect(args, log_queue):
    from ParalogWizard.cast_detect import create_reference_w_paralogs, create_reference_wo_paralogs

    logger = logging.getLogger("ParalogWizard")
    
    # Validate inputs
    data_folder = os.path.abspath(args.data_folder)
    if not os.path.isdir(data_folder):
        fatal(f"Data folder '{data_folder}' does not exist or is not a directory.")
        
    logger.info("Starting cast_detect with data_folder=%s, cores=%d", data_folder, args.num_cores)

    folder_31 = os.path.join(data_folder, "31exonic_contigs")
    if not os.path.isdir(folder_31):
        fatal(f"Required directory '{folder_31}' not found. Run cast_retrieve first.")

    all_hits_path = os.path.join(folder_31, "all_hits.tsv")
    if not os.path.isfile(all_hits_path):
        fatal(f"Required file '{all_hits_path}' not found.")

    try:
        logger.info("Loading hits data from %s", all_hits_path)
        all_hits_for_reference = pd.read_csv(all_hits_path, sep="\t")
        logger.info("Loaded %d hits from all_hits.tsv", len(all_hits_for_reference))
    except Exception as e:
        fatal(f"Failed to read hits file '{all_hits_path}': {e}")
        
    blocklist = set(args.blocklist) if args.blocklist else set()
    if blocklist:
        logger.info("Blocklist contains %d species: %s", len(blocklist), ", ".join(sorted(blocklist)))

    if not args.paralogs:
        logger.info("Running WITHOUT paralog detection")
        folder_41 = os.path.join(data_folder, "41without_par")
        os.makedirs(folder_41, exist_ok=True)
        logger.info("Output directory: %s", folder_41)
        
        create_reference_wo_paralogs(
            data_folder, all_hits_for_reference, blocklist, args.num_cores, log_queue
        )
        logger.info("Reference without paralogs created successfully")
    else:
        if args.minimum_divergence is None or args.maximum_divergence is None:
            fatal("Minimum and maximum divergence must be specified when enabling paralog detection.")

        logger.info("Running WITH paralog detection: min_div=%.3f, max_div=%.3f", 
                   args.minimum_divergence, args.maximum_divergence)
        folder_41 = os.path.join(data_folder, "41detected_par")
        os.makedirs(folder_41, exist_ok=True)
        logger.info("Output directory: %s", folder_41)
        
        create_reference_w_paralogs(
            data_folder,
            all_hits_for_reference,
            args.minimum_divergence,
            args.maximum_divergence,
            blocklist,
            args.num_cores,
            log_queue,
        )
        logger.info("Reference with paralogs created successfully")

    logger.info("ParalogWizard cast_detect completed successfully.")


@log_command
def run_cast_separate(args, log_queue):
    from ParalogWizard.cast_separate import align, generate_pslx

    logger = logging.getLogger("ParalogWizard")
    
    # Validate inputs
    data_folder = os.path.abspath(args.data_folder)
    if not os.path.isdir(data_folder):
        fatal(f"Data folder '{data_folder}' does not exist or is not a directory.")
        
    probes_customized = os.path.abspath(args.probes_customized)
    if not os.path.isfile(probes_customized):
        fatal(f"Probes customized file '{probes_customized}' not found.")

    if not (0 <= args.min_identity <= 100):
        fatal(f"Minimum identity must be between 0 and 100, got {args.min_identity}")

    redlist = set(args.redlist) if args.redlist else set()
    logger.info(
        "Running cast_separate: data_folder=%s, probes_customized=%s, min_identity=%.1f%%, cores=%d",
        data_folder, probes_customized, args.min_identity, args.num_cores
    )
    if redlist:
        logger.info("Redlist contains %d taxa: %s", len(redlist), ", ".join(sorted(redlist)))
    else:
        logger.info("No taxa in redlist - all taxa will be processed")

    logger.info("Generating PSLX files...")
    generate_pslx(
        data_folder, probes_customized, args.min_identity, redlist, args.num_cores, log_queue
    )
    logger.info("PSLX generation completed successfully")
    
    logger.info("Running alignment...")
    align(data_folder, probes_customized, args.num_cores, log_queue)
    logger.info("Alignment completed successfully")
    
    logger.info("ParalogWizard cast_separate completed successfully.")


@log_command
def run_cast_extend(args, log_queue):
    """Local-only step: requires untracked ParalogWizard/cast_extend.py."""
    try:
        from ParalogWizard.cast_extend import extend
    except ImportError:
        fatal(
            "cast_extend is not available. Keep ParalogWizard/cast_extend.py "
            "(and its helper scripts) locally; they are gitignored."
        )

    logger = setup_logging()

    data_folder = os.path.abspath(args.data_folder)
    if not os.path.isdir(data_folder):
        fatal(f"Data folder '{data_folder}' does not exist or is not a directory.")

    baitfile = os.path.abspath(args.baitfile)
    if not os.path.isfile(baitfile):
        fatal(f"Bait file '{baitfile}' not found.")

    logger.info(
        "Running cast_extend: data_folder=%s, baitfile=%s, cores=%d",
        data_folder,
        baitfile,
        args.num_cores,
    )
    extend(data_folder, baitfile, args.num_cores)
    logger.info("ParalogWizard cast_extend completed successfully.")


@log_command
def run_cast_remap(args, log_queue):
    from ParalogWizard.cast_remap import remap

    logger = setup_logging()
    
    # Validate inputs
    data_folder = os.path.abspath(args.data_folder)
    if not os.path.isdir(data_folder):
        fatal(f"Data folder '{data_folder}' does not exist or is not a directory.")
        
    probes_customized = os.path.abspath(args.probes_customized)
    if not os.path.isfile(probes_customized):
        fatal(f"Probes customized file '{probes_customized}' not found.")

    if args.exon_length <= 0:
        fatal(f"Exon length must be positive, got {args.exon_length}")

    logger.info("Running cast_remap: data_folder=%s, probes_customized=%s, exon_length=%d, cores=%d", 
               data_folder, probes_customized, args.exon_length, args.num_cores)
    
    remap(probes_customized, data_folder, args.num_cores, args.exon_length, log_queue, force=args.force)
    logger.info("ParalogWizard cast_remap completed successfully.")


@log_command
def run_cast_call(args, log_queue):
    from ParalogWizard.cast_call import call_variants

    logger = setup_logging()
    
    # Validate inputs
    data_folder = os.path.abspath(args.data_folder)
    if not os.path.isdir(data_folder):
        fatal(f"Data folder '{data_folder}' does not exist or is not a directory.")

    logger.info("Running cast_call: data_folder=%s, cores=%d", data_folder, args.num_cores)
    
    call_variants(data_folder, args.num_cores, log_queue, force=args.force)
    logger.info("ParalogWizard cast_call completed successfully.")


# @log_command
# def run_cast_polyploid(args, log_queue):
#     from ParalogWizard.cast_polyploid import polyploid

#     logger = setup_logging()
    
#     # Validate inputs
#     data_folder = os.path.abspath(args.data_folder)
#     if not os.path.isdir(data_folder):
#         fatal(f"Data folder '{data_folder}' does not exist or is not a directory.")
        
#     probes_customized = os.path.abspath(args.probes_customized)
#     if not os.path.isfile(probes_customized):
#         fatal(f"Probes customized file '{probes_customized}' not found.")

#     logger.info("Running cast_polyploid: data_folder=%s, probes_customized=%s, cores=%d", 
#                data_folder, probes_customized, args.num_cores)
    
#     try:
#         polyploids_dir = os.path.join(data_folder, "polyploids")
#         os.makedirs(polyploids_dir, exist_ok=True)
#         logger.info("Created polyploids directory: %s", polyploids_dir)

#         source_polyploids = os.path.join(data_folder, "10deduplicated_reads", "polyploids.txt")
#         if not os.path.isfile(source_polyploids):
#             fatal(f"Polyploids file '{source_polyploids}' not found.")

#         dest_polyploids = os.path.join(polyploids_dir, "polyploids.txt")
#         try:
#             shutil.copyfile(source_polyploids, dest_polyploids)
#             logger.info("Copied polyploids file from %s to %s", source_polyploids, dest_polyploids)
#         except Exception as exc:
#             fatal(f"Failed to copy polyploids file from '{source_polyploids}' to '{dest_polyploids}': {exc}")

#         polyploid(data_folder, probes_customized, args.num_cores)
#         logger.info("ParalogWizard cast_polyploid completed successfully.")
#     except Exception as e:
#         logger.error("cast_polyploid failed: %s", e)
#         raise


@log_command
def run_cast_phase(args, log_queue):
    from ParalogWizard.cast_phase import phase

    logger = setup_logging()
    
    # Validate inputs
    data_folder = os.path.abspath(args.data_folder)
    if not os.path.isdir(data_folder):
        fatal(f"Data folder '{data_folder}' does not exist or is not a directory.")

    logger.info("Running cast_phase: data_folder=%s, cores=%d", data_folder, args.num_cores)
    
    phased_dir = os.path.join(data_folder, "101phased")
    os.makedirs(phased_dir, exist_ok=True)
    logger.info("Created phased directory: %s", phased_dir)
    
    phase(data_folder, args.num_cores, logger)
    logger.info("ParalogWizard cast_phase completed successfully.")


@log_command
def run_cast_ploidy(args, log_queue):
    from ParalogWizard.cast_ploidy import ploidy

    logger = setup_logging()
    
    # Validate inputs
    data_folder = os.path.abspath(args.data_folder)
    if not os.path.isdir(data_folder):
        fatal(f"Data folder '{data_folder}' does not exist or is not a directory.")

    logger.info("Running cast_ploidy: data_folder=%s, cores=%d", data_folder, args.num_cores)
    
    ploidy_dir = os.path.join(data_folder, "102ploidy")
    os.makedirs(ploidy_dir, exist_ok=True)
    logger.info("Created ploidy directory: %s", ploidy_dir)
    
    ploidy(data_folder, args.num_cores)
    logger.info("ParalogWizard cast_ploidy completed successfully.")


# =============================================================================
#  Input validation helpers
# =============================================================================
def validate_common_args(args):
    """Validate common arguments shared across commands"""
    # Validate core count
    if args.num_cores < 1:
        fatal(f"Number of cores must be at least 1, got {args.num_cores}")
    
    max_cores = multiprocessing.cpu_count()
    if args.num_cores > max_cores:
        logger = logging.getLogger("ParalogWizard")
        logger.warning("Requested %d cores but only %d available. Using %d cores.", 
                      args.num_cores, max_cores, max_cores)
        args.num_cores = max_cores
    
    # Validate data folder
    if hasattr(args, 'data_folder') and args.data_folder:
        data_folder = os.path.abspath(args.data_folder)
        if not os.path.exists(data_folder):
            logger = logging.getLogger("ParalogWizard")
            logger.info("Data folder '%s' does not exist, creating it...", data_folder)
            try:
                os.makedirs(data_folder, exist_ok=True)
            except Exception as e:
                fatal(f"Failed to create data folder '{data_folder}': {e}")


# =============================================================================
#  Signal handling for graceful shutdown
# =============================================================================
def signal_handler(signum, frame, log_queue=None, listener=None):
    """Handle interruption signals gracefully"""
    logger = logging.getLogger("ParalogWizard")
    logger.warning("Received signal %s, shutting down gracefully...", signum)
    
    if log_queue and listener:
        try:
            log_queue.put_nowait(None)
        except Exception:
            pass
        try:
            listener.join(timeout=15)
        except Exception:
            pass
        if listener.is_alive():
            try:
                listener.terminate()
                listener.join(timeout=5)
            except Exception:
                pass
        try:
            log_queue.close()
        except Exception:
            pass
    
    sys.exit(1)


# =============================================================================
#  Dispatcher
# =============================================================================
def main():
    args = parse_args()

    # Validate common arguments early
    validate_common_args(args)

    debug = bool(getattr(args, "debug", False))
    os.environ["PARALOGWIZARD_DEBUG"] = "1" if debug else "0"

    # logfile base name visible to all children
    log_file = f"ParalogWizard_{args.command}_{datetime.now():%d.%b.%y_%H:%M}.log"
    os.environ["PARALOGWIZARD_LOGFILE"] = log_file
    paths = log_paths_from_base(log_file, debug=debug)

    # Queue + listener (writes split *.log / *.errors.log [/ *.debug.log])
    log_queue: multiprocessing.Queue = multiprocessing.Queue(-1)
    listener = multiprocessing.Process(
        target=listener_process, args=(log_queue, log_file, debug)
    )
    listener.start()
    
    # Set up signal handlers for graceful shutdown
    signal.signal(signal.SIGINT, lambda s, f: signal_handler(s, f, log_queue, listener))
    signal.signal(signal.SIGTERM, lambda s, f: signal_handler(s, f, log_queue, listener))

    level = logging.DEBUG if debug else logging.INFO

    # main logger
    main_logger = logging.getLogger("ParalogWizard")
    main_logger.setLevel(level)
    main_logger.handlers.clear()
    main_logger.addHandler(logging.handlers.QueueHandler(log_queue))

    # also echo to console (DEBUG only when --debug)
    console = logging.StreamHandler(sys.stdout)
    console.setLevel(level)
    console.setFormatter(
        logging.Formatter(
            "[%(asctime)s] [%(processName)s:%(process)d] [%(levelname)s] %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
    )
    main_logger.addHandler(console)

    main_logger.info("Starting ParalogWizard command: %s", args.command)
    main_logger.info("Normal log: %s", paths["info"])
    main_logger.info("Error log:  %s", paths["errors"])
    if debug:
        main_logger.info("Debug log:  %s", paths["debug"])
        main_logger.debug("Debug mode enabled")
    main_logger.info("Using %d CPU cores (max available: %d)", args.num_cores, multiprocessing.cpu_count())
    main_logger.info("Python version: %s", sys.version.split()[0])
    main_logger.info("Working directory: %s", os.getcwd())

    # map sub-commands
    commands: dict[str, Callable] = {
        "cast_assemble": run_cast_assemble,
        "cast_retrieve": run_cast_retrieve,
        "cast_analyze": run_cast_analyze,
        "cast_detect": run_cast_detect,
        "cast_separate": run_cast_separate,
        "cast_extend": run_cast_extend,
        "cast_remap": run_cast_remap,
        "cast_call": run_cast_call,
        # "cast_polyploid": run_cast_polyploid,  # Commented out
        "cast_phase": run_cast_phase,
        "cast_ploidy": run_cast_ploidy,
    }

    exit_code = 0
    try:
        if args.command not in commands:
            fatal(f"Unrecognized command: {args.command}")

        commands[args.command](args, log_queue)
        main_logger.info("ParalogWizard command '%s' completed successfully.", args.command)
    except BaseException as exc:  # noqa: BLE001
        exit_code = 1
        main_logger.exception("Error executing command %s: %s", args.command, exc)
    finally:
        # Drain logging before closing the queue so the listener can exit.
        try:
            time.sleep(0.2)
            log_queue.put_nowait(None)
        except Exception:
            pass

        try:
            listener.join(timeout=60)
        except Exception:
            pass
        if listener.is_alive():
            try:
                main_logger.warning(
                    "Log listener still running after 60s; terminating it"
                )
            except Exception:
                pass
            try:
                listener.terminate()
                listener.join(timeout=10)
            except Exception:
                pass

        try:
            log_queue.close()
            log_queue.join_thread()
        except Exception:
            pass

    sys.exit(exit_code)


if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()
