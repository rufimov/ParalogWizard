import glob
import os
import re
import shutil
from typing import Set, Dict, List

import Bio.SeqRecord
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastnCommandline


def sort_hit_table_cover(hittable_as_list: List[str], primary_field: str):
    """"""
    hittable_as_list.sort(
        key=lambda x: float(x.split()[1].split("_")[-1]), reverse=True
    )
    hittable_as_list.sort(key=lambda x: float(x.split()[5]), reverse=True)
    hittable_as_list.sort(key=lambda x: float(x.split()[4]))
    hittable_as_list.sort(key=lambda x: float(x.split()[2]), reverse=True)
    hittable_as_list.sort(key=lambda x: float(x.split()[3]), reverse=True)
    hittable_as_list.sort(key=globals()[primary_field])


def exon(string: str) -> str:
    """"""
    return string.split()[0].split("-")[1]


def locus(string: str) -> str:
    """"""
    return string.split()[0].split("-")[1].split("_")[0]


def contig(string: str) -> str:
    """"""
    return string.split()[1]


def contig_locus(string: str) -> str:
    """"""
    return string.split()[1].split("_N_")[0]


def slicing(
    dictionary: Dict[str, Bio.SeqRecord.SeqRecord],
    current_string: str,
    rev: bool,
) -> str:
    """"""
    sequence = dictionary[current_string.split()[1]].seq
    if rev:
        start = int(current_string.split()[7])
        end = int(current_string.split()[6])
        result = str(sequence.reverse_complement())[-start : -end - 1 : -1][::-1]
    else:
        start = int(current_string.split()[6])
        end = int(current_string.split()[7])
        result = str(sequence)[start - 1 : end]
    return result


def collect_contigs(path_to_data, logger):
    """"""

    import os

    os.makedirs(os.path.join(path_to_data, "30raw_contigs"), exist_ok=True)
    path_to_assemblies = os.path.join(path_to_data, "20assemblies")
    for folder in os.listdir(path_to_assemblies):
        if os.path.isdir(os.path.join(path_to_assemblies, folder)):
            logger.info(f"Processing {folder}")
            with open(
                os.path.join(path_to_data, "30raw_contigs", f"{folder}_contigs.fasta"),
                "w",
            ) as contigs:
                for locus in os.listdir(os.path.join(path_to_assemblies, folder)):
                    if os.path.exists(
                        os.path.join(
                            path_to_assemblies,
                            folder,
                            locus,
                            f"{locus}_contigs.fasta",
                        )
                    ):
                        logger.info(f"\tProcessing {locus}")
                        with open(
                            os.path.join(
                                path_to_assemblies,
                                folder,
                                locus,
                                f"{locus}_contigs.fasta",
                            )
                        ) as locus_contigs:
                            file_content = locus_contigs.read().replace(
                                ">", f">{locus}_"
                            )
                            contigs.write(file_content)
                        logger.info("\tOK")
            if (
                os.stat(
                    os.path.join(
                        path_to_data, "30raw_contigs", f"{folder}_contigs.fasta"
                    )
                ).st_size
                == 0
            ):
                os.remove(
                    os.path.join(
                        path_to_data, "30raw_contigs", f"{folder}_contigs.fasta"
                    )
                )
            logger.info("Ok")


def prepare_contigs(main_path, logger):
    """"""
    logger.info("Preparing congits...")
    for file in glob.glob(
        os.path.join(os.path.dirname(main_path), "30raw_contigs", "*contigs.fasta")
    ):
        shutil.copy(file, os.path.join(main_path, os.path.basename(file)))
    for file in glob.glob(os.path.join(main_path, "*.fasta")):
        with open(file, "r") as fasta:
            lines: List[str] = fasta.readlines()
        with open(file, "w") as fasta:
            for line in lines:
                line = line.replace("NODE", "N")
                line = re.sub(
                    r"length_([0-9]+)_cov_([0-9]+\.[0-9][0-9]).*", r"\1_c_\2", line
                )
                fasta.write(line)
        name_of_file: str = "_".join(
            os.path.basename(file).split(".")[0].split("_")[0:2]
        )
        os.rename(file, os.path.join(os.path.dirname(file), f"{name_of_file}.fasta"))
    logger.info("Done\n")


def create_hit_tables(main_path, probe_exons, n_cpu, logger, length_cover=75):
    """"""
    logger.info("Creating hit tables...")
    for file in glob.glob(os.path.join(main_path, "*.fasta")):
        file: str = os.path.basename(file)
        sample: str = file[:-6]
        logger.info(f"\tProcessing {sample}")
        NcbimakeblastdbCommandline(
            dbtype="nucl",
            input_file=os.path.join(main_path, file),
            out=os.path.join(main_path, sample),
            parse_seqids=True,
        )()
        logger.info("\tRunning BLAST...")
        NcbiblastnCommandline(
            task="blastn",
            query=probe_exons,
            db=os.path.join(main_path, sample),
            out=os.path.join(main_path, f"reference_in_{sample}_contigs.txt"),
            qcov_hsp_perc=length_cover,
            num_threads=n_cpu,
            outfmt="6 qaccver saccver pident qcovhsp evalue bitscore sstart send",
        )()
        logger.info("\tOK")
    logger.info("Done\n")


def correct_contgis(
    main_path, statistics, all_hits_for_reference, logger, spades_cover=5
):
    """"""
    logger.info("Correcting contigs...")
    for file in glob.glob(os.path.join(main_path, "*.fasta")):
        sample: str = os.path.basename(os.path.splitext(file)[0])
        logger.info(f" Processing {sample}")
        statistics[sample]: Dict[str, Dict[str, List[str]]] = dict()
        hits: List[str] = list()
        with open(
            os.path.join(main_path, f"reference_in_{sample}_contigs.txt")
        ) as blast_results:
            blast_results_as_list = [x[:-1] for x in blast_results.readlines()]
            for line in blast_results_as_list:
                if (
                    contig_locus(line) == locus(line)
                    and float(line.split()[1].split("_c_")[1]) >= spades_cover
                ):
                    hits.append(line)
        sort_hit_table_cover(hits, "exon")
        hits_exons_contigs: Set[str] = set()
        hits_dedup: List[str] = list()
        for hit in hits:
            if str(exon(hit) + contig(hit)) not in hits_exons_contigs:
                hits_dedup.append(hit)
            else:
                pass
            hits_exons_contigs.add(exon(hit) + contig(hit))
        with open(os.path.join(main_path, f"{sample}.fas"), "w") as result_fasta, open(
            file
        ) as contigs:
            contigs_fasta_parsed = SeqIO.to_dict(
                SeqIO.parse(contigs, "fasta", generic_dna)
            )
            for hit_dedup in hits_dedup:
                if int(hit_dedup.split()[6]) > int(hit_dedup.split()[7]):
                    sequence: str = slicing(contigs_fasta_parsed, hit_dedup, True)
                else:
                    sequence: str = slicing(contigs_fasta_parsed, hit_dedup, False)
                result_fasta.write(
                    f">{exon(hit_dedup)}_{'_'.join(contig(hit_dedup).split('_')[1:])}\n{sequence}\n"
                )
                all_hits_for_reference.append(f"{hit_dedup}\t{sample}\t{sequence}")
        sort_hit_table_cover(hits_dedup, "locus")
        hits_loci: Set[str] = set()
        hits_exons: Set[str] = set()
        for hit_dedup in hits_dedup:
            if locus(hit_dedup) not in hits_loci:
                statistics[sample][locus(hit_dedup)] = dict()
                if exon(hit_dedup) not in hits_exons:
                    statistics[sample][locus(hit_dedup)][exon(hit_dedup)] = [
                        contig(hit_dedup)
                    ]
                    hits_exons.add(exon(hit_dedup))
                else:
                    statistics[sample][locus(hit_dedup)][exon(hit_dedup)].append(
                        contig(hit_dedup)
                    )
                hits_loci.add(locus(hit_dedup))
            else:
                if exon(hit_dedup) not in hits_exons:
                    statistics[sample][locus(hit_dedup)][exon(hit_dedup)] = [
                        contig(hit_dedup)
                    ]
                    hits_exons.add(exon(hit_dedup))
                else:
                    statistics[sample][locus(hit_dedup)][exon(hit_dedup)].append(
                        contig(hit_dedup)
                    )
        for i in statistics[sample].keys():
            stat = []
            for j in statistics[sample][i].keys():
                stat.append(len(statistics[sample][i][j]))
            statistics[sample][i]: int = max(stat)
        with open(
            os.path.join(main_path, f"reference_against_{sample}_contigs.txt"), "w"
        ) as hittable:
            sort_hit_table_cover(hits, "locus")
            for hit in hits:
                hittable.write(f"{hit}\n")
    with open(os.path.join(main_path, "all_hits.txt"), "w") as all_hits_to_write:
        for hit in all_hits_for_reference:
            all_hits_to_write.write(f"{hit}\n")
    logger.info(" OK")
    logger.info("All contigs were successfully corrected!\n")


def write_stats(main_path, probe_exons, statistics, logger):
    """"""

    logger.info("Writing statistics...")
    with open(os.path.join(main_path, "statistics.tsv"), "w") as stats, open(
        probe_exons
    ) as reference:
        stats_dict: Dict[str, str] = dict([("gene\t", "")])
        loci: Set[str] = set()
        samples: List[str] = list()
        reference_as_list: List[str] = [x[:-1] for x in reference.readlines()]
        for line in reference_as_list:
            if line.startswith(">"):
                loci.add("_".join(line.split("-")[1].split("_")[:-2]))
        for key in statistics.keys():
            samples.append(key)
        samples.sort()
        for sample in samples:
            stats_dict["gene\t"]: str = stats_dict["gene\t"] + sample + "\t"
        for loc in loci:
            stats_dict[loc + "\t"]: str = ""
            for sample in samples:
                loci_in_sample: Set[str] = set(statistics[sample].keys())
                if loc in loci_in_sample:
                    stats_dict[loc + "\t"]: str = (
                        stats_dict[loc + "\t"] + str(statistics[sample][loc]) + "\t"
                    )
                else:
                    stats_dict[loc + "\t"]: str = stats_dict[loc + "\t"] + "NA" + "\t"
        stats.write("gene\t" + stats_dict["gene\t"] + "\n")
        del stats_dict["gene\t"]
        for key in sorted(list(stats_dict.keys())):
            stats.write(f"{key}{stats_dict[key]}\n")
    logger.info("Statistics file created!\n")


def rename_contigs(main_path, logger):
    """"""

    logger.info("Renaming contigs...")
    for file in glob.glob(os.path.join(main_path, "*.fas")):
        sample = os.path.basename(os.path.splitext(file)[0])
        logger.info(f" Processing {sample}")
        with open(file) as result_fasta:
            fasta_parsed = SeqIO.to_dict(
                SeqIO.parse(result_fasta, "fasta", generic_dna)
            )
            counter: int = 1
            fasta_to_write: List[str] = list()
            for line in sorted(fasta_parsed.keys()):
                fasta_to_write.append(
                    f">Contig{str(counter)}_{sample}-{line.replace('_', '-')}\n"
                )
                fasta_to_write.append(str(fasta_parsed[line].seq) + "\n")
                counter += 1
        with open(file, "w") as result_fasta:
            result_fasta.writelines(fasta_to_write)
        logger.info(" OK")
    logger.info("All contigs were successfully renamed!\n")


def clean(main_path, logger):
    """"""

    logger.info("Removing temporary files...")
    for file in glob.glob(os.path.join(main_path, "*.fasta")):
        os.remove(file)
    for file in glob.glob(os.path.join(main_path, "*.n*")):
        os.remove(file)
    for file in glob.glob(os.path.join(main_path, "reference_in*")):
        os.remove(file)
    logger.info("Done\n")
