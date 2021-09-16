#!/usr/bin/env python3

import argparse
import logging
import os
import sys
from datetime import datetime
from glob import glob
from typing import Set, Dict, List, Union


class ParsedArgs:
    def __init__(self):
        parser = argparse.ArgumentParser(
            usage="""ParalogWizard <command> [<args>]
The ParalogWizard commands are:
    cast_assemble
    cast_retrieve
    cast_analyze
    cast_detect
    cast_separate    
Use ParalogWizard <command> -h for help with arguments of the command of interest
"""
        )
        parser.add_argument(
            "command",
            help="Subcommand to run",
        )
        self.args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, self.args.command):
            print("Unrecognized command")
            parser.print_help()
            exit(1)

    def common_args(self, parser):
        parser.add_argument("-d", "--data_folder", required=True)
        parser.add_argument(
            "-nc",
            "--num_cores",
            default=1,
            type=int,
            help="Number of cores used. Default: 1",
        )

    def cast_assemble(self):
        parser = argparse.ArgumentParser()
        parser.add_argument(
            "-pr",
            "--baitfile",
            help="FASTA file containing bait sequences for each gene. Each probe sequence is concatenated exons. If "
            "there are multiple baits for a gene, "
            "the id must be of the form: >Taxon-geneName",
            required=True,
        )
        parser.add_argument(
            "-np",
            "--no-parallel",
            dest="parallel",
            action="store_true",
            default=False,
        )
        self.common_args(parser)
        args = parser.parse_args(sys.argv[2:])
        return args

    def cast_retrieve(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("-pe", "--probes_exons", required=True)
        parser.add_argument("-l", "--length_cover", default=None, type=float)
        parser.add_argument("-s", "--spades_cover", default=None, type=float)
        parser.add_argument(
            "-c",
            "--collect_contigs",
            action="store_true",
            default=False,
        )
        self.common_args(parser)
        args = parser.parse_args(sys.argv[2:])
        return args

    def cast_analyze(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("-b", "--blocklist", nargs="+", required=False)
        self.common_args(parser)
        args = parser.parse_args(sys.argv[2:])
        if args.blocklist is None:
            args.blocklist = set()
        else:
            args.blocklist = set(args.blocklist)
        return args

    def cast_detect(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("-b", "--blocklist", nargs="+", required=False)
        parser.add_argument("-pe", "--probes_exons", required=True)
        parser.add_argument(
            "-p",
            "--paralogs",
            default=False,
            action="store_true",
            help="Paralog detection. Default: off",
        )
        parser.add_argument(
            "-mi",
            "--minimum_divergence",
            required="--paralogs" in sys.argv,
            type=float,
        )
        parser.add_argument(
            "-ma",
            "--maximum_divergence",
            required="--paralogs" in sys.argv,
            type=float,
        )
        self.common_args(parser)
        args = parser.parse_args(sys.argv[2:])
        if args.blocklist is None:
            args.blocklist = set()
        else:
            args.blocklist = set(args.blocklist)
        if args.paralogs is True and (
            args.minimum_divergence is None or args.maximum_divergence is None
        ):
            print("Minimum and maximum divergence is not specified")
            parser.print_help()
            exit(1)
        return args

    def cast_separate(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("-r", "--redlist", nargs="+", required=False)
        parser.add_argument("-i", "--min_identity", required=True, type=float)
        parser.add_argument("-pc", "--probes_customized", required=True)
        self.common_args(parser)
        args = parser.parse_args(sys.argv[2:])
        if args.redlist is None:
            args.redlist = set()
        else:
            args.redlist = set(args.redlist)
        return args

    def cast_align(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("-r", "--redlist", nargs="+", required=False)
        parser.add_argument("-i", "--min_identity", required=True, type=float)
        parser.add_argument("-pc", "--probes_customized", required=True)
        self.common_args(parser)
        args = parser.parse_args(sys.argv[2:])

    def get_args_dict(self):
        command = self.args.__dict__
        arguments = getattr(self, self.args.command)().__dict__
        argument_dictionary = command.copy()
        argument_dictionary.update(arguments)
        return argument_dictionary


def main():
    arguments = ParsedArgs().get_args_dict()
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    log_handler = logging.FileHandler(
        f'ParalogWizard_{arguments["command"]}_{datetime.now().strftime("%d.%b.%y_%H.%M")}.log',
        "w",
    )
    log_handler.setLevel(logging.INFO)
    log_formatter = logging.Formatter(
        fmt="%(asctime)s - %(message)s", datefmt="%d-%b-%y %H:%M:%S"
    )
    log_handler.setFormatter(log_formatter)
    logger.addHandler(log_handler)
    if arguments["command"] == "cast_assemble":
        from ParalogWizard.HybPiper_spades import spades, clean_up
        from ParalogWizard.HybPiper_bwa_distribute import bwa
        from ParalogWizard.HybPiper_bwa_distribute import distribute_bwa

        os.makedirs(arguments["data_folder"], exist_ok=True)
        os.chdir(arguments["data_folder"])

        with open(
            os.path.join("10deduplicated_reads", "samples_list.txt")
        ) as samples_list:
            samples = [x.strip() for x in samples_list.readlines()]
        os.makedirs("20assemblies", exist_ok=True)
        os.chdir("20assemblies")
        for sample in samples:
            readfiles = glob(
                f"{os.path.join('..', '10deduplicated_reads', sample)}.*.fastq"
            )
            # Generate directory
            os.makedirs(sample, exist_ok=True)
            os.chdir(sample)
            readfiles = [
                os.path.join(
                    "..", "..", "10deduplicated_reads", os.path.basename(readfiles[0])
                ),
                os.path.join(
                    "..", "..", "10deduplicated_reads", os.path.basename(readfiles[1])
                ),
            ]
            levels_to_workdir = "../" * (arguments["data_folder"].count("/") + 1)
            baitfile = os.path.join(
                levels_to_workdir, "..", "..", arguments["baitfile"]
            )
            bamfile = bwa(
                readfiles,
                baitfile,
                sample,
                cpu=arguments["num_cores"],
            )
            if not bamfile:
                print("ERROR: Something went wrong with the BWA step, exiting!")
                return

            pre_existing_fastas = glob("./*/*_interleaved.fasta")
            for fn in pre_existing_fastas:
                os.remove(fn)

            exitcode = distribute_bwa(bamfile, readfiles)
            if exitcode:
                sys.exit(1)
            genes = [
                x
                for x in os.listdir(".")
                if os.path.isfile(os.path.join(x, x + "_interleaved.fasta"))
            ]
            if not arguments["parallel"]:
                spades(
                    readfiles=readfiles,
                    genes=genes,
                    cpu=arguments["num_cores"],
                )
            else:
                spades(
                    readfiles=readfiles,
                    genes=genes,
                    cpu=arguments["num_cores"],
                    parallel=False,
                )

            os.chdir("..")
            clean_up(sample)
            os.chdir("..")
    elif arguments["command"] == "cast_retrieve":
        from ParalogWizard.cast_retrieve import (
            collect_contigs,
            prepare_contigs,
            create_hit_tables,
            correct_contgis,
            write_stats,
            rename_contigs,
            clean,
        )

        logger.info(
            f"""ParalogWizard cast_retrieve running with the following settings
            main data folder - {arguments["data_folder"]}
            collect contigs - {arguments['collect_contigs']}
            probe file - {arguments["probes_exons"]}
            filter for blast length cover - {arguments["length_cover"]}
            k-mer cover threshold for spades contigs - {arguments["spades_cover"]}
            number of used cores - {arguments["num_cores"]}"""
        )
        logger.info("Retrieving data...\n")
        if not arguments["collect_contigs"] and not os.path.isdir(
            os.path.join(arguments["data_folder"], "30raw_contigs")
        ):
            print(
                "ERROR: No raw contigs found. Run ParalogWizard cast_retrieve with -c specified."
            )
            exit(1)
        elif arguments["collect_contigs"]:
            collect_contigs(arguments["data_folder"], logger)
        statistics: Dict[str, Dict[str, Union[Dict[str, List[str]], int]]] = dict()
        all_hits_for_reference: List[str] = list()
        main_path = os.path.join(arguments["data_folder"], "31exonic_contigs")
        os.makedirs(main_path, exist_ok=True)
        prepare_contigs(
            main_path,
            logger,
        )
        if arguments["length_cover"]:
            create_hit_tables(
                main_path,
                arguments["probes_exons"],
                arguments["num_cores"],
                logger,
                arguments["length_cover"],
            )
        else:
            create_hit_tables(
                main_path,
                arguments["probes_exons"],
                arguments["num_cores"],
                logger,
            )
        if arguments["spades_cover"]:
            correct_contgis(
                main_path,
                statistics,
                all_hits_for_reference,
                logger,
                arguments["spades_cover"],
            )
        else:
            correct_contgis(
                main_path,
                statistics,
                all_hits_for_reference,
                logger,
            )
        write_stats(main_path, arguments["probes_exons"], statistics, logger)
        rename_contigs(main_path, logger)
        clean(main_path, logger)
        logger.info("Data was successfully retrieved!")
    elif arguments["command"] == "cast_analyze":
        from ParalogWizard.cast_analyze import estimate_divergence, build_alignments

        if len(arguments["blocklist"]) > 0:
            blocklist_string = ", ".join(sp for sp in list(arguments["blocklist"]))
            logger.info(
                f"""ParalogWizard cast_collect running with the following settings
            main data folder - {arguments["data_folder"]}
            species not taken to paralogs divergency estimation - {blocklist_string}
            number of used cores - {arguments["num_cores"]}"""
            )
        else:
            logger.info(
                f"""ParalogWizard cast_collect running with the following settings
            main data folder - {arguments["data_folder"]}
            all species taken to paralogs divergency estimation
            number of used cores - {arguments["num_cores"]}"""
            )
        build_alignments(arguments["data_folder"], arguments["num_cores"], logger)
        estimate_divergence(arguments["data_folder"], arguments["blocklist"], logger)
    elif arguments["command"] == "cast_detect":
        from ParalogWizard.cast_detect import (
            score_samples,
            create_reference_wo_paralogs,
            create_reference_w_paralogs,
            refine_phasing,
            write_paralog_stats,
        )

        if not arguments["paralogs"]:
            os.makedirs(os.path.join(arguments["data_folder"], "41without_par"))
            if len(arguments["blocklist"]) > 0:
                blocklist_string = ", ".join(sp for sp in list(arguments["blocklist"]))
                logger.info(
                    f"""ParalogWizard cast_collect running with the following settings
                main data folder - {arguments["data_folder"]}
                paralogs are not being searched
                species not taken to paralogs divergency estimation - {blocklist_string}"""
                )
            else:
                logger.info(
                    f"""ParalogWizard cast_collect running with the following settings
                main data folder - {arguments["data_folder"]}
                paralogs are not being searched
                all species taken to paralogs divergency estimation"""
                )
            with open(
                os.path.join(
                    arguments["data_folder"], "31exonic_contigs", "all_hits.txt"
                )
            ) as all_hits:
                all_hits_for_reference: List[str] = [
                    x[:-1] for x in all_hits.readlines()
                ]
            all_hits_for_reference_scored = score_samples(all_hits_for_reference)
            create_reference_wo_paralogs(
                arguments["data_folder"],
                all_hits_for_reference_scored,
                arguments["blocklist"],
                logger,
            )
        else:
            if len(arguments["blocklist"]) > 0:
                blocklist_string = ", ".join(sp for sp in list(arguments["blocklist"]))
                logger.info(
                    f"""ParalogWizard cast_collect running with the following settings
                main data folder - {arguments["data_folder"]}
                paralogs are searched with minium {arguments["minimum_divergence"]} \
                and maximum {arguments["maximum_divergence"]} divergence 
                species not taken to paralogs divergency estimation - {blocklist_string}"""
                )
            else:
                logger.info(
                    f"""ParalogWizard cast_collect running with the following settings
                main data folder - {arguments["data_folder"]}
                paralogs are searched with minium {arguments["minimum_divergence"]} \
                and maximum {arguments["maximum_divergence"]} divergence 
                all species taken to paralogs divergency estimation"""
                )
            with open(
                os.path.join(
                    arguments["data_folder"], "31exonic_contigs", "all_hits.txt"
                )
            ) as all_hits:
                all_hits_for_reference: List[str] = [
                    x[:-1] for x in all_hits.readlines()
                ]
            all_hits_for_reference_scored = score_samples(all_hits_for_reference)
            paralog_statistic: Dict[str, Set[str]] = dict()
            os.makedirs(os.path.join(arguments["data_folder"], "41detected_par"))
            create_reference_w_paralogs(
                arguments["data_folder"],
                all_hits_for_reference_scored,
                paralog_statistic,
                arguments["minimum_divergence"],
                arguments["maximum_divergence"],
                arguments["blocklist"],
                logger,
            )
            refine_phasing(
                arguments["data_folder"],
                arguments["minimum_divergence"],
                arguments["maximum_divergence"],
                logger,
            )
            write_paralog_stats(
                arguments["data_folder"],
                arguments["minimum_divergence"],
                arguments["maximum_divergence"],
                paralog_statistic,
                arguments["probes_exons"],
                logger,
            )
    elif arguments["command"] == "cast_separate":

        from ParalogWizard.cast_separate import run_blat, correct, align

        if len(arguments["redlist"]) > 0:
            redlist_string = ", ".join(sp for sp in list(arguments["redlist"]))
            logger.info(
                f"""ParalogWizard cast_collect running with the following settings
                main data folder - {arguments["data_folder"]}
                probe file with separated paralogs - {arguments["probes_paralogs"]}
                minimum identity for BLAT - {arguments["min_identity"]}
                list of taxa excluded from paralogs separation - {redlist_string}"""
            )
        else:
            logger.info(
                f"""ParalogWizard cast_collect running with the following settings
                main data folder - {arguments["data_folder"]}
                probe file with separated paralogs - {arguments["probes_paralogs"]}
                minimum identity for BLAT - {arguments["min_identity"]}
                all taxa included to paralogs separation"""
            )
        run_blat(
            arguments["data_folder"],
            arguments["probes_customized"],
            arguments["min_identity"],
            logger,
        )
        if os.path.exists(
            os.path.join(
                arguments["data_folder"], "50pslx", "corrected", "list_pslx.txt"
            )
        ):
            os.remove(
                os.path.join(
                    arguments["data_folder"], "50pslx", "corrected", "list_pslx.txt"
                )
            )
        correct(
            arguments["data_folder"],
            arguments["redlist"],
            logger,
        )
        align(
            arguments["data_folder"],
            arguments["probes_customized"],
            arguments["num_cores"],
        )


if __name__ == "__main__":
    main()
