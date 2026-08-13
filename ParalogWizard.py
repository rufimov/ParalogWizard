#!/usr/bin/env python3
"""
ParalogWizard – HybSeq paralog detection, remapping, and polyploid phasing.

Commands (typical order):
  cast_assemble → cast_retrieve → cast_analyze → cast_detect → cast_separate
  → cast_remap → cast_call → cast_ploidy
  → cast_polyploid   (hybrid A/B phasing vs a backbone)
  → cast_phase       (VCF consensus / IterPol; optional)
  → cast_extend      (local exonerate supercontigs; optional)
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
        description=(
            "ParalogWizard: HybSeq paralog detection and polyploid phasing. "
            "Commands: assemble, retrieve, analyze, detect, separate, remap, "
            "call, ploidy, polyploid, phase."
        )
    )

    # Shared logging flags for every subcommand.
    log_parent = argparse.ArgumentParser(add_help=False)
    log_parent.add_argument(
        "--debug",
        action="store_true",
        help=(
            "Enable debug logging: DEBUG lines are written to the single *.log "
            "file and also shown on the console (INFO/ERROR always go to *.log)"
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

    # ---- cast_call -----------------------------------------------------------
    pcl = subparsers.add_parser(
        "cast_call", parents=[log_parent], help="Run the cast_call step"
    )
    pcl.add_argument("-d", "--data_folder", required=True, help="Main data folder")
    pcl.add_argument("-nc", "--num_cores", default=1, type=int, help="Cores (default 1)")

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
        "cast_polyploid",
        parents=[log_parent],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help=(
            "Phase polyploids vs a backbone → 101polyploid/ "
            "(lists, nearest-N ranks, IQ-TREE, AMAS concat)"
        ),
        description=POLYPLOID_DESCRIPTION,
        epilog=POLYPLOID_EPILOG,
    )
    ppp.add_argument(
        "-d",
        "--data_folder",
        required=True,
        help=(
            "Main data folder. Needs 41detected_par/all_paralogs_for_reference.tsv "
            "(from cast_detect -p) and 10deduplicated_reads/samples_list.txt"
        ),
    )
    ppp.add_argument(
        "-bb",
        "--backbone_list",
        default="backbone_list.txt",
        metavar="FILE",
        help=(
            "Diploid / main-phylogeny samples used as nearest-species references "
            "(one name per line). Default: backbone_list.txt in cwd or data folder"
        ),
    )
    ppp.add_argument(
        "-pl",
        "--polyploid_list",
        default="polyploid_list.txt",
        metavar="FILE",
        help=(
            "Polyploid/hybrid samples to phase (one name per line; several allowed). "
            "Must not overlap backbone_list. Default: polyploid_list.txt in cwd or "
            "data folder"
        ),
    )
    ppp.add_argument(
        "-ns",
        "--n_nearest",
        default=5,
        type=int,
        metavar="N",
        help=(
            "How many closest backbone samples to keep for each polyploid tip when "
            "building the rank matrix used to cluster phases A/B (default: 5). "
            "Clipped to the number of backbone samples if N is larger"
        ),
    )
    ppp.add_argument(
        "-nc",
        "--num_cores",
        default=1,
        type=int,
        help="CPU cores for MAFFT / IQ-TREE (default: 1)",
    )

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
    from ParalogWizard.cast_assemble import assemble

    assemble(
        data_folder=args.data_folder,
        baitfile=args.baitfile,
        num_cores=args.num_cores,
        log_queue=log_queue,
    )


@log_command
def run_cast_retrieve(args, log_queue):
    from ParalogWizard.cast_retrieve import retrieve

    retrieve(
        data_folder=os.path.abspath(args.data_folder),
        collect=bool(args.collect_contigs),
        probe_exons=os.path.abspath(args.probes_exons),
        num_cores=args.num_cores,
        length_cover=args.length_cover,
        spades_cover=args.spades_cover,
        log_queue=log_queue,
    )


@log_command
def run_cast_analyze(args, log_queue):
    from ParalogWizard.cast_analyze import analyze

    blocklist = set(args.blocklist) if args.blocklist else set()
    analyze(
        data_folder=args.data_folder,
        num_cores=args.num_cores,
        log_queue=log_queue,
        blocklist=blocklist,
    )


@log_command
def run_cast_detect(args, log_queue):
    from ParalogWizard.cast_detect import create_reference_w_paralogs, create_reference_wo_paralogs

    logger = logging.getLogger("ParalogWizard")

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

    required_cols = {
        "saccver",
        "pident",
        "qcovhsp",
        "evalue",
        "bitscore",
        "k-mer_cover",
        "exon",
        "locus",
        "sample",
        "sequence",
    }
    missing = required_cols - set(all_hits_for_reference.columns)
    if missing:
        fatal(f"all_hits.tsv missing required column(s): {', '.join(sorted(missing))}")

    blocklist = set(args.blocklist) if args.blocklist else set()
    if blocklist:
        logger.info(
            "Blocklist contains %d species: %s",
            len(blocklist),
            ", ".join(sorted(blocklist)),
        )

    if not args.paralogs:
        logger.info("Running WITHOUT paralog detection")
        create_reference_wo_paralogs(
            data_folder, all_hits_for_reference, blocklist, args.num_cores, log_queue
        )
        logger.info("Reference without paralogs created successfully")
    else:
        if args.minimum_divergence is None or args.maximum_divergence is None:
            fatal(
                "Minimum and maximum divergence must be specified when enabling "
                "paralog detection."
            )
        dist_path = os.path.join(data_folder, "40aln_orth_par", "pairwise_distances.tsv")
        if not os.path.isfile(dist_path):
            fatal(
                f"Required file '{dist_path}' not found. Run cast_analyze first."
            )
        logger.info(
            "Running WITH paralog detection: min_div=%.3f, max_div=%.3f",
            args.minimum_divergence,
            args.maximum_divergence,
        )
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
    
    remap(probes_customized, data_folder, args.num_cores, args.exon_length, log_queue)
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
    
    call_variants(data_folder, args.num_cores, log_queue)
    logger.info("ParalogWizard cast_call completed successfully.")


@log_command
def run_cast_polyploid(args, log_queue):
    from ParalogWizard.cast_polyploid import polyploid

    logger = logging.getLogger("ParalogWizard")

    data_folder = os.path.abspath(args.data_folder)
    if not os.path.isdir(data_folder):
        fatal(f"Data folder '{data_folder}' does not exist or is not a directory.")

    paralogs_tsv = os.path.join(
        data_folder, "41detected_par", "all_paralogs_for_reference.tsv"
    )
    if not os.path.isfile(paralogs_tsv):
        fatal(
            f"Required file '{paralogs_tsv}' not found. "
            "Run cast_detect with -p first."
        )

    if args.n_nearest < 1:
        fatal("-ns / --n_nearest must be an integer >= 1")

    logger.info(
        "Running cast_polyploid: data_folder=%s, backbone_list=%s, "
        "polyploid_list=%s, n_nearest=%d, cores=%d",
        data_folder,
        args.backbone_list,
        args.polyploid_list,
        args.n_nearest,
        args.num_cores,
    )
    polyploid(
        data_folder,
        args.num_cores,
        log_queue,
        backbone_list=args.backbone_list,
        polyploid_list=args.polyploid_list,
        n_nearest=args.n_nearest,
    )
    logger.info("ParalogWizard cast_polyploid completed successfully.")


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

    logger = logging.getLogger("ParalogWizard")

    data_folder = os.path.abspath(args.data_folder)
    if not os.path.isdir(data_folder):
        fatal(f"Data folder '{data_folder}' does not exist or is not a directory.")

    remapped = os.path.join(data_folder, "100remapped")
    if not os.path.isdir(remapped):
        fatal(f"Required directory '{remapped}' not found. Run cast_remap first.")

    logger.info(
        "Running cast_ploidy: data_folder=%s, cores=%d", data_folder, args.num_cores
    )
    ploidy(data_folder, args.num_cores, log_queue)
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

    # Queue + listener (single combined *.log; DEBUG included only with --debug)
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
    main_logger.info("Log file: %s", paths["info"])
    if debug:
        main_logger.debug("Debug mode enabled (DEBUG lines included in log file)")
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
        "cast_polyploid": run_cast_polyploid,
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
