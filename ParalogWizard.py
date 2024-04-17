#!/usr/bin/ python3

"""
ParalogWizard is a pipeline for paralog detection in HybPiper output. It is based on the HybPiper pipeline and uses
the same input files. The pipeline consists of 5 steps: cast_assemble, cast_retrieve, cast_analyze, cast_detect,
and cast_separate. The first step is the assembly using SPAdes. The second step is the retrieval of
contigs from the assembly. The third step is the estimation of paralog divergence. The fourth step is the detection of
paralogs. The fifth step is the separation of paralogs to alignemnts. The pipeline is run using the command line
interface. The pipeline is described in the README file.
"""

import argparse
import logging
import multiprocessing
import os
import shutil
import sys
from datetime import datetime
from glob import glob
import pandas


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
        parser.add_argument("-l", "--length_cover", default=75, type=float)
        parser.add_argument("-s", "--spades_cover", default=5, type=float)
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

    def cast_extend(self):
        parser = argparse.ArgumentParser()
        parser.add_argument(
            "-pr",
            "--baitfile",
            required=True,
        )
        self.common_args(parser)
        args = parser.parse_args(sys.argv[2:])
        return args

    def cast_flow(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("-pc", "--probes_customized", required=True)
        parser.add_argument("-o", "--outgroup", required=True)
        parser.add_argument("-e", "--exon_length", type=int, default=150)
        self.common_args(parser)
        args = parser.parse_args(sys.argv[2:])
        return args

    def cast_ploidy(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("-pc", "--probes_customized", required=True)
        parser.add_argument("-e", "--exon_length", type=int, default=150)
        self.common_args(parser)
        args = parser.parse_args(sys.argv[2:])
        return args

    def cast_hybrid(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("-pc", "--probes_customized", required=True)
        self.common_args(parser)
        args = parser.parse_args(sys.argv[2:])
        return args

    def get_args_dict(self):
        command = self.args.__dict__
        arguments = getattr(self, self.args.command)().__dict__
        argument_dictionary = command.copy()
        argument_dictionary.update(arguments)
        return argument_dictionary


def create_logger(log_file):
    logger = multiprocessing.get_logger()
    logger.setLevel(logging.INFO)
    log_handler_info = logging.FileHandler(log_file)
    log_handler_info.setLevel(logging.INFO)
    log_formatter_info = logging.Formatter(
        fmt="%(asctime)s - %(message)s", datefmt="%d-%b-%y %H:%M:%S"
    )
    log_handler_info.setFormatter(log_formatter_info)
    logger.addHandler(log_handler_info)
    return logger


def main():
    arguments = ParsedArgs().get_args_dict()
    log_file = f'ParalogWizard_{arguments["command"]}_{datetime.now().strftime("%d.%b.%y_%H:%M")}.log'
    logger = create_logger(log_file)
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
        from ParalogWizard.cast_retrieve import retrieve

        logger.info(
            f"""ParalogWizard cast_retrieve running with the following settings
            main data folder - {arguments["data_folder"]}
            collect contigs - {arguments['collect_contigs']}
            probe file - {arguments["probes_exons"]}
            filter for blast length cover - {arguments["length_cover"]}
            k-mer cover threshold for spades contigs - {arguments["spades_cover"]}
            number of used cores - {arguments["num_cores"]}"""
        )
        retrieve(
            arguments["data_folder"],
            arguments["collect_contigs"],
            arguments["probes_exons"],
            arguments["num_cores"],
            arguments["length_cover"],
            arguments["spades_cover"],
            log_file,
        )
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
        estimate_divergence(
            arguments["data_folder"],
            arguments["blocklist"],
            arguments["num_cores"],
            logger,
        )
    elif arguments["command"] == "cast_detect":
        from ParalogWizard.cast_detect import (
            create_reference_w_paralogs,
            create_reference_wo_paralogs,
        )

        folder_31 = os.path.join(arguments["data_folder"], "31exonic_contigs")
        all_hits_for_reference = pandas.read_csv(
            os.path.join(folder_31, "all_hits.tsv"), sep="\t"
        )
        if not arguments["paralogs"] and len(arguments["blocklist"]) > 0:
            blocklist_string = ", ".join(sp for sp in list(arguments["blocklist"]))
            logger.info(
                f"""ParalogWizard cast_collect running with the following settings
                    main data folder - {arguments["data_folder"]}
                    paralogs are not being searched
                    species not taken to paralogs divergency estimation - {blocklist_string}
                    number of used cores - {arguments["num_cores"]}"""
            )
        elif not arguments["paralogs"] and len(arguments["blocklist"]) == 0:
            logger.info(
                f"""ParalogWizard cast_collect running with the following settings
                    main data folder - {arguments["data_folder"]}
                    paralogs are not being searched
                    all species taken to paralogs divergency estimation
                    number of used cores - {arguments["num_cores"]}"""
            )
        elif len(arguments["blocklist"]) > 0:
            blocklist_string = ", ".join(sp for sp in list(arguments["blocklist"]))
            logger.info(
                f"""ParalogWizard cast_collect running with the following settings
                    main data folder - {arguments["data_folder"]}
                    paralogs are searched with minium {arguments["minimum_divergence"]} \
                    and maximum {arguments["maximum_divergence"]} divergence 
                    species not taken to paralogs divergency estimation - {blocklist_string}
                    number of used cores - {arguments["num_cores"]}"""
            )
        elif len(arguments["blocklist"]) == 0:
            logger.info(
                f"""ParalogWizard cast_collect running with the following settings
                    main data folder - {arguments["data_folder"]}
                    paralogs are searched with minium {arguments["minimum_divergence"]} \
                    and maximum {arguments["maximum_divergence"]} divergence 
                    all species taken to paralogs divergency estimation
                    number of used cores - {arguments["num_cores"]}"""
            )
        if not arguments["paralogs"]:
            folder_41 = os.path.join(arguments["data_folder"], "41without_par")
            os.makedirs(folder_41, exist_ok=True)
            create_reference_wo_paralogs(
                arguments["data_folder"],
                all_hits_for_reference,
                arguments["blocklist"],
                arguments["num_cores"],
                log_file,
            )
        else:
            folder_41 = os.path.join(arguments["data_folder"], "41detected_par")
            os.makedirs(folder_41, exist_ok=True)
            create_reference_w_paralogs(
                arguments["data_folder"],
                all_hits_for_reference,
                arguments["minimum_divergence"],
                arguments["maximum_divergence"],
                arguments["blocklist"],
                arguments["num_cores"],
                log_file,
            )
    elif arguments["command"] == "cast_separate":
        from ParalogWizard.cast_separate import align, generate_pslx

        if len(arguments["redlist"]) > 0:
            redlist_string = ", ".join(sp for sp in list(arguments["redlist"]))
            logger.info(
                f"""ParalogWizard cast_collect running with the following settings
                main data folder - {arguments["data_folder"]}
                probe file with separated paralogs - {arguments["probes_customized"]}
                minimum identity for BLAT - {arguments["min_identity"]}
                list of taxa excluded from paralogs separation - {redlist_string}"""
            )
        else:
            logger.info(
                f"""ParalogWizard cast_collect running with the following settings
                main data folder - {arguments["data_folder"]}
                probe file with separated paralogs - {arguments["probes_customized"]}
                minimum identity for BLAT - {arguments["min_identity"]}
                all taxa included to paralogs separation"""
            )
        generate_pslx(
            arguments["data_folder"],
            arguments["probes_customized"],
            arguments["min_identity"],
            arguments["redlist"],
            arguments["num_cores"],
            log_file,
        )
        align(
            arguments["data_folder"],
            arguments["probes_customized"],
            arguments["num_cores"],
            log_file,
        )
    elif arguments["command"] == "cast_extend":
        from ParalogWizard.cast_extend import extend

        extend(arguments["data_folder"], arguments["baitfile"], arguments["num_cores"])

    elif arguments["command"] == "cast_ploidy":
        from ParalogWizard.cast_ploidy import ploidy

        ploidy(
            arguments["probes_customized"],
            arguments["data_folder"],
            arguments["num_cores"],
            arguments["exon_length"],
        )
    elif arguments["command"] == "cast_flow":
        from ParalogWizard.cast_flow import flow

        flow(
            arguments["probes_customized"],
            arguments["data_folder"],
            arguments["num_cores"],
            arguments["outgroup"],
            arguments["exon_length"],
        )
    elif arguments["command"] == "cast_hybrid":
        logger.info(
            f"""ParalogWizard cast_hybrid running with the following settings
                        main data folder - {arguments["data_folder"]}
                        probe file with separated paralogs - {arguments["probes_customized"]}"""
        )
        from ParalogWizard.cast_hybrid import hybrid

        os.makedirs(os.path.join(arguments["data_folder"], "hybrids"), exist_ok=True)
        shutil.copyfile(
            os.path.join(
                arguments["data_folder"], "10deduplicated_reads", "hybrids.txt"
            ),
            os.path.join(arguments["data_folder"], "hybrids", "hybrids.txt"),
        )
        hybrid(
            arguments["data_folder"],
            arguments["probes_customized"],
            arguments["num_cores"],
            log_file,
        )


if __name__ == "__main__":
    main()
