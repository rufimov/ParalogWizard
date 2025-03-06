#!/usr/bin/python3
"""
ParalogWizard is a pipeline for paralog detection in HybSeq data.
It consists of several steps:
    cast_assemble, cast_retrieve, cast_analyze, cast_detect, cast_separate,
    cast_extend, cast_remap, cast_call, cast_ploidy, cast_polyploid, and cast_phase.
"""

import argparse
import logging
import logging.handlers
import multiprocessing
import os
import sys
import time
from datetime import datetime

import shutil
from contextlib import contextmanager
import functools
from glob import glob

import pandas as pd

# Import shared functions from your ParalogWizard package.
from ParalogWizard import (
    decompress_fastq_files,
    compress_fastq_files,
    setup_logging,
    listener_process,
)


# -----------------------------------------------------------------------------
# Helper: Context Manager for Changing Directories
# -----------------------------------------------------------------------------
@contextmanager
def change_dir(destination: str):
    """Temporarily change the working directory."""
    prev_dir = os.getcwd()
    os.chdir(destination)
    try:
        yield
    finally:
        os.chdir(prev_dir)


# -----------------------------------------------------------------------------
# Helper: Decorator for Logging Command Function Entry/Exit
# -----------------------------------------------------------------------------
def log_command(func):
    @functools.wraps(func)
    def wrapper(args, log_queue):
        logger = logging.getLogger("ParalogWizard")
        logger.info("Starting %s", func.__name__)
        try:
            result = func(args, log_queue)
            logger.info("Completed %s successfully.", func.__name__)
            return result
        except Exception as e:
            logger.exception("Error in %s: %s", func.__name__, e)
            sys.exit(1)

    return wrapper


# -----------------------------------------------------------------------------
# Setup Argparse with Subparsers for Each Command
# -----------------------------------------------------------------------------
def setup_argparser():
    parser = argparse.ArgumentParser(
        description="ParalogWizard: A pipeline for paralog detection in HybSeq data."
    )
    subparsers = parser.add_subparsers(
        dest="command", required=True, help="Subcommand to run"
    )

    # --- cast_assemble ---
    parser_assemble = subparsers.add_parser(
        "cast_assemble", help="Run the cast_assemble step"
    )
    parser_assemble.add_argument(
        "-pr",
        "--baitfile",
        required=True,
        help=(
            "FASTA file containing bait sequences for each gene. "
            "If multiple baits exist for a gene, the id must be of the form: >Taxon-geneName"
        ),
    )
    parser_assemble.add_argument(
        "-d", "--data_folder", required=True, help="Main data folder"
    )
    parser_assemble.add_argument(
        "-nc",
        "--num_cores",
        default=1,
        type=int,
        help="Number of cores to use (default: 1)",
    )

    # --- cast_retrieve ---
    parser_retrieve = subparsers.add_parser(
        "cast_retrieve", help="Run the cast_retrieve step"
    )
    parser_retrieve.add_argument(
        "-pe", "--probes_exons", required=True, help="Path to the probes_exons file"
    )
    parser_retrieve.add_argument(
        "-l",
        "--length_cover",
        default=75,
        type=float,
        help="Length cover (default: 75)",
    )
    parser_retrieve.add_argument(
        "-s", "--spades_cover", default=5, type=float, help="Spades cover (default: 5)"
    )
    parser_retrieve.add_argument(
        "-c",
        "--collect_contigs",
        action="store_true",
        default=False,
        help="Collect contigs flag",
    )
    parser_retrieve.add_argument(
        "-d", "--data_folder", required=True, help="Main data folder"
    )
    parser_retrieve.add_argument(
        "-nc",
        "--num_cores",
        default=1,
        type=int,
        help="Number of cores to use (default: 1)",
    )

    # --- cast_analyze ---
    parser_analyze = subparsers.add_parser(
        "cast_analyze", help="Run the cast_analyze step"
    )
    parser_analyze.add_argument(
        "-b",
        "--blocklist",
        nargs="+",
        help="List of species to exclude from paralog divergence estimation",
    )
    parser_analyze.add_argument(
        "-d", "--data_folder", required=True, help="Main data folder"
    )
    parser_analyze.add_argument(
        "-nc",
        "--num_cores",
        default=1,
        type=int,
        help="Number of cores to use (default: 1)",
    )

    # --- cast_detect ---
    parser_detect = subparsers.add_parser(
        "cast_detect", help="Run the cast_detect step"
    )
    parser_detect.add_argument(
        "-b",
        "--blocklist",
        nargs="+",
        help="List of species to exclude from paralog divergence estimation",
    )
    parser_detect.add_argument(
        "-pe", "--probes_exons", required=True, help="Path to the probes_exons file"
    )
    parser_detect.add_argument(
        "-p",
        "--paralogs",
        action="store_true",
        default=False,
        help="Enable paralog detection (default: off)",
    )
    parser_detect.add_argument(
        "-mi",
        "--minimum_divergence",
        type=float,
        help="Minimum divergence (required if --paralogs is set)",
    )
    parser_detect.add_argument(
        "-ma",
        "--maximum_divergence",
        type=float,
        help="Maximum divergence (required if --paralogs is set)",
    )
    parser_detect.add_argument(
        "-d", "--data_folder", required=True, help="Main data folder"
    )
    parser_detect.add_argument(
        "-nc",
        "--num_cores",
        default=1,
        type=int,
        help="Number of cores to use (default: 1)",
    )

    # --- cast_separate ---
    parser_separate = subparsers.add_parser(
        "cast_separate", help="Run the cast_separate step"
    )
    parser_separate.add_argument(
        "-r",
        "--redlist",
        nargs="+",
        help="List of taxa to exclude from paralog separation",
    )
    parser_separate.add_argument(
        "-i",
        "--min_identity",
        required=True,
        type=float,
        help="Minimum identity for BLAT",
    )
    parser_separate.add_argument(
        "-pc",
        "--probes_customized",
        required=True,
        help="Path to the probes_customized file",
    )
    parser_separate.add_argument(
        "-d", "--data_folder", required=True, help="Main data folder"
    )
    parser_separate.add_argument(
        "-nc",
        "--num_cores",
        default=1,
        type=int,
        help="Number of cores to use (default: 1)",
    )

    # --- cast_extend ---
    parser_extend = subparsers.add_parser(
        "cast_extend", help="Run the cast_extend step"
    )
    parser_extend.add_argument(
        "-pr", "--baitfile", required=True, help="FASTA file containing bait sequences"
    )
    parser_extend.add_argument(
        "-d", "--data_folder", required=True, help="Main data folder"
    )
    parser_extend.add_argument(
        "-nc",
        "--num_cores",
        default=1,
        type=int,
        help="Number of cores to use (default: 1)",
    )

    # --- cast_remap ---
    parser_remap = subparsers.add_parser("cast_remap", help="Run the cast_remap step")
    parser_remap.add_argument(
        "-pc",
        "--probes_customized",
        required=True,
        help="Path to the probes_customized file",
    )
    parser_remap.add_argument(
        "-e", "--exon_length", type=int, default=150, help="Exon length (default: 150)"
    )
    parser_remap.add_argument(
        "-d", "--data_folder", required=True, help="Main data folder"
    )
    parser_remap.add_argument(
        "-nc",
        "--num_cores",
        default=1,
        type=int,
        help="Number of cores to use (default: 1)",
    )

    # --- cast_call ---
    parser_call = subparsers.add_parser("cast_call", help="Run the cast_call step")
    parser_call.add_argument(
        "-d", "--data_folder", required=True, help="Main data folder"
    )
    parser_call.add_argument(
        "-nc",
        "--num_cores",
        default=1,
        type=int,
        help="Number of cores to use (default: 1)",
    )

    # --- cast_phase ---
    parser_phase = subparsers.add_parser("cast_phase", help="Run the cast_phase step")
    parser_phase.add_argument(
        "-d", "--data_folder", required=True, help="Main data folder"
    )
    parser_phase.add_argument(
        "-nc",
        "--num_cores",
        default=1,
        type=int,
        help="Number of cores to use (default: 1)",
    )

    # --- cast_ploidy ---
    parser_ploidy = subparsers.add_parser(
        "cast_ploidy", help="Run the cast_ploidy step"
    )
    parser_ploidy.add_argument(
        "-d", "--data_folder", required=True, help="Main data folder"
    )
    parser_ploidy.add_argument(
        "-nc",
        "--num_cores",
        default=1,
        type=int,
        help="Number of cores to use (default: 1)",
    )

    # --- cast_polyploid ---
    parser_polyploid = subparsers.add_parser(
        "cast_polyploid", help="Run the cast_polyploid step"
    )
    parser_polyploid.add_argument(
        "-pc",
        "--probes_customized",
        required=True,
        help="Path to the probes_customized file",
    )
    parser_polyploid.add_argument(
        "-d", "--data_folder", required=True, help="Main data folder"
    )
    parser_polyploid.add_argument(
        "-nc",
        "--num_cores",
        default=1,
        type=int,
        help="Number of cores to use (default: 1)",
    )

    return parser


def parse_args():
    parser = setup_argparser()
    return parser.parse_args()


# -----------------------------------------------------------------------------
# Command Functions
# -----------------------------------------------------------------------------


@log_command
def run_cast_assemble(args, log_queue: multiprocessing.Queue):
    """
    Execute the cast_assemble step.
    This involves:
      - Reading a samples list,
      - Locating and decompressing FASTQ files,
      - Running BWA and SPAdes (or a variant with parallel disabled),
      - Cleaning up and compressing FASTQ files back.
    """
    from ParalogWizard.cast_assemble import (
        run_spades_assembly,
        cleanup_gene_directories,
    )
    from ParalogWizard.cast_assemble import run_bwa_mem, distribute_bwa_reads

    logger = setup_logging()
    logger.info("Running cast_assemble with data_folder: %s", args.data_folder)
    data_folder = os.path.abspath(args.data_folder)
    os.makedirs(data_folder, exist_ok=True)

    with change_dir(data_folder):
        samples_list = os.path.join("10deduplicated_reads", "samples_list.txt")
        if not os.path.isfile(samples_list):
            logger.error("Samples list file '%s' not found.", samples_list)
            sys.exit(1)
        with open(samples_list) as f:
            samples = [line.strip() for line in f if line.strip()]

        if not samples:
            logger.error("No samples found in '%s'.", samples_list)
            sys.exit(1)

        os.makedirs("20assemblies", exist_ok=True)
        with change_dir("20assemblies"):
            for sample in samples:
                sample_dir = os.path.join(os.getcwd(), sample)
                os.makedirs(sample_dir, exist_ok=True)
                with change_dir(sample_dir):
                    pattern = (
                        os.path.join("..", "..", "10deduplicated_reads", sample)
                        + ".*.fastq*"
                    )
                    readfiles = glob(pattern)
                    if len(readfiles) < 2:
                        logger.error(
                            "Not enough readfiles for sample: %s (pattern: %s)",
                            sample,
                            pattern,
                        )
                        continue

                    # Use the first two files and compute absolute paths.
                    readfiles = [
                        os.path.abspath(
                            os.path.join(
                                "..",
                                "..",
                                "10deduplicated_reads",
                                os.path.basename(readfiles[0]),
                            )
                        ),
                        os.path.abspath(
                            os.path.join(
                                "..",
                                "..",
                                "10deduplicated_reads",
                                os.path.basename(readfiles[1]),
                            )
                        ),
                    ]

                    for rf in readfiles:
                        if not os.path.isfile(rf):
                            logger.error("Readfile '%s' does not exist.", rf)
                            sys.exit(1)

                    # Decompress FASTQ files if necessary.
                    readfiles = decompress_fastq_files(readfiles)
                    baitfile = os.path.abspath(
                        os.path.join(data_folder, "..", args.baitfile)
                    )
                    if not os.path.isfile(baitfile):
                        logger.error("Bait file '%s' not found.", baitfile)
                        sys.exit(1)

                    bamfile = run_bwa_mem(
                        read_fastq_files=readfiles,
                        bait_filepath=baitfile,
                        output_basename=sample,
                        num_threads=args.num_cores,
                    )
                    if not bamfile:
                        logger.error("BWA step failed for sample: %s", sample)
                        compress_fastq_files(readfiles)
                        sys.exit(1)

                    for fn in glob("./*/*_interleaved.fasta"):
                        try:
                            os.remove(fn)
                        except Exception as e:
                            logger.error("Failed to remove file '%s': %s", fn, e)
                            sys.exit(1)

                    exitcode = distribute_bwa_reads(
                        bam_filepath=bamfile,
                        fastq_files=readfiles,
                    )
                    if exitcode:
                        logger.error("BWA distribution failed for sample: %s", sample)
                        sys.exit(1)

                    genes = [
                        x
                        for x in os.listdir(".")
                        if os.path.isfile(os.path.join(x, x + "_interleaved.fasta"))
                    ]
                    run_spades_assembly(
                        paired_fastq_files=readfiles,
                        gene_list=genes,
                        num_threads=args.num_cores,
                        log_queue=log_queue,
                    )

                    compress_fastq_files(readfiles)
                    cleanup_gene_directories(sample)
    logger.info("ParalogWizard cast_assemble completed.")


@log_command
def run_cast_retrieve(args, log_queue):
    """
    Execute the cast_retrieve step.
    This step retrieves exonic regions based on the probes file.
    """
    from ParalogWizard.cast_retrieve import retrieve

    logger = setup_logging()

    # Check that probes_exons file exists.
    probes_exons = os.path.abspath(args.probes_exons)
    if not os.path.isfile(probes_exons):
        logger.error("Probes exons file '%s' not found.", probes_exons)
        sys.exit(1)

    logger.info(
        "Running cast_retrieve with data_folder=%s, probes_exons=%s",
        args.data_folder,
        probes_exons,
    )
    retrieve(
        args.data_folder,
        args.collect_contigs,
        probes_exons,
        args.num_cores,
        args.length_cover,
        args.spades_cover,
        log_queue,
    )
    logger.info("ParalogWizard cast_retrieve completed.")


@log_command
def run_cast_analyze(args, log_queue):
    """
    Execute the cast_analyze step.
    This step builds alignments and estimates divergence.
    """
    from ParalogWizard.cast_analyze import estimate_divergence, build_alignments

    logger = setup_logging()
    blocklist = set(args.blocklist) if args.blocklist else set()
    if blocklist:
        logger.info(
            "ParalogWizard cast_analyze running with data_folder: %s, blocklist: %s, num_cores: %d",
            args.data_folder,
            ", ".join(blocklist),
            args.num_cores,
        )
    else:
        logger.info(
            "ParalogWizard cast_analyze running with data_folder: %s, all species included, num_cores: %d",
            args.data_folder,
            args.num_cores,
        )

    build_alignments(args.data_folder, args.num_cores, log_queue)
    estimate_divergence(args.data_folder, blocklist, args.num_cores, log_queue)
    logger.info("ParalogWizard cast_analyze completed.")


@log_command
def run_cast_detect(args, log_queue: multiprocessing.Queue):
    """
    Execute the cast_detect step.
    Depending on whether paralog detection is enabled,
    create reference files with or without paralogs.
    """
    from ParalogWizard.cast_detect import (
        create_reference_w_paralogs,
        create_reference_wo_paralogs,
    )

    logger = setup_logging()
    logger.info("Starting cast_detect")

    folder_31 = os.path.join(args.data_folder, "31exonic_contigs")
    if not os.path.isdir(folder_31):
        logger.error("Required directory '%s' not found.", folder_31)
        sys.exit(1)
    all_hits_path = os.path.join(folder_31, "all_hits.tsv")
    if not os.path.isfile(all_hits_path):
        logger.error("Required file '%s' not found.", all_hits_path)
        sys.exit(1)
    all_hits_for_reference = pd.read_csv(all_hits_path, sep="\t")
    blocklist = set(args.blocklist) if args.blocklist else set()

    if not args.paralogs:
        if blocklist:
            logger.info(
                "Running cast_detect with data_folder: %s, no paralogs, blocklist: %s, num_cores: %d",
                args.data_folder,
                ", ".join(blocklist),
                args.num_cores,
            )
        else:
            logger.info(
                "Running cast_detect with data_folder: %s, no paralogs, all species included, num_cores: %d",
                args.data_folder,
                args.num_cores,
            )
        folder_41 = os.path.join(args.data_folder, "41without_par")
        os.makedirs(folder_41, exist_ok=True)
        create_reference_wo_paralogs(
            args.data_folder,
            all_hits_for_reference,
            blocklist,
            args.num_cores,
            log_queue,
        )
    else:
        if args.minimum_divergence is None or args.maximum_divergence is None:
            logger.error(
                "Minimum and maximum divergence must be specified when enabling paralog detection."
            )
            sys.exit(1)
        if blocklist:
            logger.info(
                "Running cast_detect with data_folder: %s, paralogs enabled (min: %s, max: %s), "
                "blocklist: %s, num_cores: %d",
                args.data_folder,
                args.minimum_divergence,
                args.maximum_divergence,
                ", ".join(blocklist),
                args.num_cores,
            )
        else:
            logger.info(
                "Running cast_detect with data_folder: %s, paralogs enabled (min: %s, max: %s), "
                "all species included, num_cores: %d",
                args.data_folder,
                args.minimum_divergence,
                args.maximum_divergence,
                args.num_cores,
            )
        folder_41 = os.path.join(args.data_folder, "41detected_par")
        os.makedirs(folder_41, exist_ok=True)
        create_reference_w_paralogs(
            args.data_folder,
            all_hits_for_reference,
            args.minimum_divergence,
            args.maximum_divergence,
            blocklist,
            args.num_cores,
            log_queue,
        )

    logger.info("ParalogWizard cast_detect completed.")


@log_command
def run_cast_separate(args, log_queue: multiprocessing.Queue):
    """
    Execute the cast_separate step.
    This step aligns sequences and generates PSLX files.
    """
    from ParalogWizard.cast_separate import align, generate_pslx

    logger = setup_logging()
    probes_customized = os.path.abspath(args.probes_customized)
    if not os.path.isfile(probes_customized):
        logger.error("Probes customized file '%s' not found.", probes_customized)
        sys.exit(1)

    logger.info(
        "ParalogWizard cast_separate running with data_folder: %s, probes_customized: %s, min_identity: %s, "
        "num_cores: %d",
        args.data_folder,
        probes_customized,
        args.min_identity,
        args.num_cores,
    )
    redlist = set(args.redlist) if args.redlist else set()
    if redlist:
        logger.info("Taxa excluded from paralog separation: %s", ", ".join(redlist))
    else:
        logger.info("All taxa included in paralog separation.")

    generate_pslx(
        args.data_folder,
        probes_customized,
        args.min_identity,
        redlist,
        args.num_cores,
        log_queue,
    )
    align(args.data_folder, probes_customized, args.num_cores, log_queue)
    logger.info("ParalogWizard cast_separate completed.")


@log_command
def run_cast_extend(args, log_queue: multiprocessing.Queue = None):
    """
    Execute the cast_extend step.
    This step extends the bait sequences.
    """
    from ParalogWizard.cast_extend import extend

    logger = setup_logging()
    baitfile = os.path.abspath(args.baitfile)
    if not os.path.isfile(baitfile):
        logger.error("Bait file '%s' not found.", baitfile)
        sys.exit(1)

    logger.info("Running cast_extend with baitfile: %s", baitfile)
    extend(args.data_folder, baitfile, args.num_cores)
    logger.info("ParalogWizard cast_extend completed.")


@log_command
def run_cast_remap(args, log_queue: multiprocessing.Queue):
    """
    Execute the cast_remap step.
    This step remaps sequences using the provided customized probes.
    """
    from ParalogWizard.cast_remap import remap

    logger = setup_logging()
    probes_customized = os.path.abspath(args.probes_customized)
    if not os.path.isfile(probes_customized):
        logger.error("Probes customized file '%s' not found.", probes_customized)
        sys.exit(1)

    logger.info(
        "Running cast_remap with probes_customized: %s and exon_length: %d",
        probes_customized,
        args.exon_length,
    )
    remap(
        probes_customized, args.data_folder, args.num_cores, args.exon_length, log_queue
    )
    logger.info("ParalogWizard cast_remap completed.")


@log_command
def run_cast_call(args, log_queue: multiprocessing.Queue):
    """
    Execute the cast_call step.
    This step calls variants from the assembled data.
    """
    from ParalogWizard.cast_call import call_variants

    logger = setup_logging()
    logger.info("Running cast_call")
    call_variants(args.data_folder, args.num_cores, log_queue)
    logger.info("ParalogWizard cast_call completed.")


@log_command
def run_cast_polyploid(args, log_queue: multiprocessing.Queue = None):
    """
    Execute the cast_polyploid step.
    This step detects and processes polyploid samples.
    """
    from ParalogWizard.cast_polyploid import polyploid

    logger = setup_logging()
    probes_customized = os.path.abspath(args.probes_customized)
    if not os.path.isfile(probes_customized):
        logger.error("Probes customized file '%s' not found.", probes_customized)
        sys.exit(1)

    logger.info("Running cast_polyploid with probes_customized: %s", probes_customized)
    polyploids_dir = os.path.join(args.data_folder, "polyploids")
    os.makedirs(polyploids_dir, exist_ok=True)
    source_polyploids = os.path.join(
        args.data_folder, "10deduplicated_reads", "polyploids.txt"
    )
    if not os.path.isfile(source_polyploids):
        logger.error("Polyploids file '%s' not found.", source_polyploids)
        sys.exit(1)
    dest_polyploids = os.path.join(polyploids_dir, "polyploids.txt")
    try:
        shutil.copyfile(source_polyploids, dest_polyploids)
    except Exception as e:
        logger.error("Failed to copy polyploids file: %s", e)
        sys.exit(1)
    polyploid(args.data_folder, probes_customized, args.num_cores)
    logger.info("ParalogWizard cast_polyploid completed.")


@log_command
def run_cast_phase(args, log_queue: multiprocessing.Queue = None):
    """
    Execute the cast_phase step.
    This step phases sequences.
    """
    from ParalogWizard.cast_phase import phase

    logger = setup_logging()
    logger.info("Running cast_phase")
    phased_dir = os.path.join(args.data_folder, "101phased")
    os.makedirs(phased_dir, exist_ok=True)
    phase(args.data_folder, args.num_cores, logger)
    logger.info("ParalogWizard cast_phase completed.")


@log_command
def run_cast_ploidy(args, log_queue: multiprocessing.Queue = None):
    """
    Execute the cast_ploidy step.
    This step estimates ploidy.
    """
    from ParalogWizard.cast_ploidy import ploidy

    logger = setup_logging()
    logger.info("Running cast_ploidy")
    ploidy_dir = os.path.join(args.data_folder, "102ploidy")
    os.makedirs(ploidy_dir, exist_ok=True)
    ploidy(args.data_folder, args.num_cores)
    logger.info("ParalogWizard cast_ploidy completed.")


# -----------------------------------------------------------------------------
# Main Dispatcher
# -----------------------------------------------------------------------------
def main():
    parser = setup_argparser()
    args = parser.parse_args()

    # Create a log file name and set it in an environment variable so children can use it.
    log_file = (
        f"ParalogWizard_{args.command}_{datetime.now().strftime('%d.%b.%y_%H:%M')}.log"
    )
    os.environ["PARALOGWIZARD_LOGFILE"] = log_file

    # Set up the logging queue and start the listener process.
    log_queue = multiprocessing.Queue(-1)
    listener = multiprocessing.Process(
        target=listener_process, args=(log_queue, log_file)
    )
    listener.start()

    # Create a logger for the main process that sends logs to the queue.
    main_logger = logging.getLogger("ParalogWizard")
    main_logger.setLevel(logging.INFO)
    qh = logging.handlers.QueueHandler(log_queue)
    main_logger.handlers = []  # Clear any default handlers.
    main_logger.addHandler(qh)

    # --- Add a console (stream) handler here ---
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    formatter = logging.Formatter(
        "[%(asctime)s] [%(processName)s:%(process)d] [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    console_handler.setFormatter(formatter)
    main_logger.addHandler(console_handler)
    # ------------------------------------------------

    main_logger.info("Starting ParalogWizard command: %s", args.command)

    # Map subcommands to their respective functions.
    commands = {
        "cast_assemble": run_cast_assemble,
        "cast_retrieve": run_cast_retrieve,
        "cast_analyze": run_cast_analyze,
        "cast_detect": run_cast_detect,
        "cast_separate": run_cast_separate,
        "cast_extend": run_cast_extend,
        "cast_remap": run_cast_remap,
        "cast_call": run_cast_call,
        "cast_polyploid": run_cast_polyploid,
        "cast_phase": run_cast_phase,
        "cast_ploidy": run_cast_ploidy,
    }

    if args.command in commands:
        try:
            commands[args.command](args, log_queue)
        except Exception as e:
            main_logger.exception("Error executing command %s: %s", args.command, e)
            sys.exit(1)
    else:
        main_logger.error("Unrecognized command: %s", args.command)
        sys.exit(1)

    main_logger.info("ParalogWizard command '%s' completed successfully.", args.command)

    # Signal the listener to exit.
    time.sleep(1)
    log_queue.put(None)

    # Close the queue and wait for its internal thread to finish.
    log_queue.close()
    log_queue.join_thread()

    # Now join the listener process.
    listener.join()


if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()
