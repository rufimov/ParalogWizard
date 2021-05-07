import copy
import os
from typing import Dict
from typing import List
from typing import Set

from Bio import SeqIO

from ParalogWizard.cast_analyze import sort_hit_table_ident
from ParalogWizard.cast_retrieve import exon, contig


def score_samples(list_with_hits: List[str]) -> List[str]:
    """"""
    sort_hit_table_ident(list_with_hits, "exon")
    # Add counters to hits for each exon
    counter: int = 0
    current_exon: str = ""
    for hit in list_with_hits:
        if exon(hit) != current_exon:
            counter = 1
            current_exon = exon(hit)
        else:
            counter += 1
        list_with_hits[list_with_hits.index(hit)] = f"{hit}\t{counter}"
    # Prepare dictionary for counting scores for contigs in each locus
    for_samples_scores = dict()
    samples: Set[str] = set()
    loci_samples: Set[str] = set()
    for hit in list_with_hits:
        sample: str = hit.split()[8]
        locus: str = exon(hit).split("_")[0]
        locus_sample: str = locus + sample
        spades_contig: str = contig(hit)
        if sample not in samples:
            for_samples_scores[sample] = dict()
            for_samples_scores[sample][locus] = [dict(), set()]
            for_samples_scores[sample][locus][0][spades_contig] = []
            samples.add(sample)
            loci_samples.add(locus_sample)
        else:
            if locus_sample not in loci_samples:
                for_samples_scores[sample][locus] = [dict(), set()]
                for_samples_scores[sample][locus][0][spades_contig] = []
                loci_samples.add(locus_sample)
            else:
                for_samples_scores[sample][locus][0][spades_contig] = []
    # Count scores for contigs in each locus
    for hit in list_with_hits:
        sample: str = hit.split()[8]
        locus: str = exon(hit).split("_")[0]
        spades_contig: str = contig(hit)
        for_samples_scores[sample][locus][0][spades_contig].append(int(hit.split()[-1]))
        for_samples_scores[sample][locus][1].add(exon(hit))
    samples_scores = copy.deepcopy(for_samples_scores)
    for sample in for_samples_scores.keys():
        for locus in for_samples_scores[sample].keys():
            samples_scores[sample][locus][1]: int = len(
                for_samples_scores[sample][locus][1]
            )
            sample_locus_scores: List[float] = []
            for spades_contig in for_samples_scores[sample][locus][0].keys():
                samples_scores[sample][locus][0][spades_contig] = sum(
                    for_samples_scores[sample][locus][0][spades_contig]
                ) / len(for_samples_scores[sample][locus][0][spades_contig])
                sample_locus_scores.append(
                    samples_scores[sample][locus][0][spades_contig]
                )
            sample_locus_score: float = sum(sample_locus_scores) / len(
                sample_locus_scores
            )
            samples_scores[sample][locus].append(sample_locus_score)
    # Adding scores to each hit
    for i in range(0, len(list_with_hits)):
        hit: str = list_with_hits[i]
        locus: str = exon(hit).split("_")[0]
        spades_contig: str = contig(hit)
        sample: str = hit.split()[8]
        list_with_hits[i] = (
            f"{hit}\t{str(samples_scores[sample][locus][1])}\t{str(samples_scores[sample][locus][2])}"
            f"\t{str(samples_scores[sample][locus][0][spades_contig])}"
        )
    list_with_hits.sort(key=lambda x: float(x.split()[-1]))
    list_with_hits.sort(key=lambda x: x.split()[8])
    list_with_hits.sort(key=lambda x: float(x.split()[-2]))
    list_with_hits.sort(key=lambda x: float(x.split()[-3]), reverse=True)
    list_with_hits.sort(key=lambda x: exon(x))
    return ["\t".join(x.split()[:10]) for x in list_with_hits]


def create_reference_wo_paralogs(
    path_to_data, all_hits_for_reference_scored, blocklist, logger
):
    """"""

    logger.info("Creating new reference...")
    exons: Set[str] = set()
    with open(
        os.path.join(
            path_to_data, "41without_par", "new_reference_for_HybPhyloMaker.fas"
        ),
        "w",
    ) as new_reference_HPM, open(
        os.path.join(
            path_to_data,
            "41without_par",
            "new_reference_for_HybPiper_separate_exons.fas",
        ),
        "w",
    ) as fasta_to_concatenate:
        for hit in all_hits_for_reference_scored:
            sample: str = hit.split()[8]
            if sample not in blocklist:
                if exon(hit) not in exons:
                    name_of_locus: str = (
                        exon(hit)
                        .replace("exon", "Contig")
                        .replace("Exon", "Contig")
                        .replace("contig", "Contig")
                        .replace("_", "")
                        .replace("Contig", "_Contig_")
                    )
                    new_reference_HPM.write(
                        f">Assembly_{name_of_locus}_{sample}_{contig(hit)}\n{hit.split()[9]}\n"
                    )
                    fasta_to_concatenate.write(
                        f">{sample.replace('-', '_')}_{contig(hit)}-{name_of_locus.replace('Contig', 'exon')}\n"
                        f"{hit.split()[9]}\n"
                    )
                    exons.add(exon(hit))
    with open(
        os.path.join(
            path_to_data,
            "41without_par",
            "new_reference_for_HybPiper_separate_exons.fas",
        )
    ) as fasta_to_concatenate, open(
        os.path.join(
            path_to_data,
            "41without_par",
            "new_reference_for_HybPiper_concatenated_exons.fas",
        ),
        "w",
    ) as concatenated_fasta:
        fasta_parsed = SeqIO.to_dict(
            SeqIO.parse(fasta_to_concatenate, "fasta")
        )
        current_locus = ""
        list_of_keys = list(fasta_parsed.keys())
        list_of_keys.sort(key=lambda x: int(x.split("-")[1].split("_")[2]))
        list_of_keys.sort(key=lambda x: x.split("-")[1].split("_")[0])
        count = 1
        for key in list_of_keys:
            locus = key.split("-")[1].split("_")[0]
            if count == 1:
                if locus != current_locus:
                    concatenated_fasta.write(
                        ">"
                        + "_".join(key.split("-")[0].split("_")[0:3])
                        + "-"
                        + locus
                        + "\n"
                        + str(fasta_parsed[key].seq)
                    )
                else:
                    concatenated_fasta.write(str(fasta_parsed[key].seq))
            else:
                if locus != current_locus:
                    concatenated_fasta.write(
                        "\n>"
                        + "_".join(key.split("-")[0].split("_")[0:3])
                        + "-"
                        + locus
                        + "\n"
                        + str(fasta_parsed[key].seq)
                    )
                else:
                    concatenated_fasta.write(str(fasta_parsed[key].seq))
            current_locus = locus
            count += 1

    logger.info("New reference created!\n")


def create_reference_w_paralogs(
    path_to_data,
    all_hits_for_reference_scored,
    paralog_statistic,
    paralog_min_divergence,
    paralog_max_divergence,
    blocklist,
    logger,
):
    """"""

    with open(
        os.path.join(path_to_data, "40aln_orth_par", "pairwise_distances.txt")
    ) as distances:
        pairwise_distances: Dict[str, float] = dict()
        for line in distances.read().splitlines():
            pairwise_distances[f"{line.split()[0]}_{line.split()[2]}"] = float(
                line.split()[1]
            )
            pairwise_distances[f"{line.split()[2]}_{line.split()[0]}"] = float(
                line.split()[1]
            )
    logger.info("Creating new reference...")
    all_paralogs_for_reference = []
    exons: Set[str] = set()
    count: int = 0
    paralog_found: bool = bool()
    current_best: str = ""
    current_locus: Dict[str, str] = dict()
    samples_current_locus: Set[str] = set()
    samples_with_paralogs: Set[str] = set()
    set_of_samples = set()
    for hit in all_hits_for_reference_scored:
        sample: str = hit.split()[8]
        if sample not in set_of_samples:
            set_of_samples.add(sample)
            paralog_statistic[sample] = set()
        if exon(hit) not in exons:
            if not paralog_found and count != 0:
                logger.info(
                    f"No paralog found for {current_best.split()[0].split('-')[1]}"
                )
                for sample_wo_paralog in current_locus.keys():
                    all_paralogs_for_reference.append(current_locus[sample_wo_paralog])
            current_locus: Dict[str, str] = dict()
            samples_current_locus: Set[str] = set()
            samples_with_paralogs: Set[str] = set()
            current_best: str = hit
            current_locus[sample] = hit
            samples_current_locus.add(sample)
            paralog_found: bool = False
            exons.add(exon(hit))
        else:
            if sample in samples_current_locus and sample not in samples_with_paralogs:
                current_exonic_contig = (
                    f"{hit.split()[0].split('-')[1]}_"
                    f"{'_'.join(hit.split()[1].split('_')[1:])}_{sample}"
                )
                exonic_contigs_to_compare = (
                    f"{current_locus[sample].split()[0].split('-')[1]}_"
                    f"{'_'.join(current_locus[sample].split()[1].split('_')[1:])}"
                    f"_{sample}"
                )
                if (
                    paralog_max_divergence
                    > pairwise_distances[
                        f"{current_exonic_contig}_{exonic_contigs_to_compare}"
                    ]
                    > paralog_min_divergence
                ):
                    logger.info(
                        f"Paralog detected for {hit.split()[0].split('-')[1]} in {sample}"
                    )
                    paralog_statistic[sample].add(hit.split()[1].split("_")[0])
                    all_paralogs_for_reference.append(current_locus[sample])
                    name_of_locus_para: str = "_".join(
                        [exon(hit).split("_")[0] + "_para"] + exon(hit).split("_")[1:]
                    )
                    name_of_locus_para = (
                        hit.split()[0].split("-")[0] + "-" + name_of_locus_para
                    )
                    hit = "\t".join([name_of_locus_para] + hit.split()[1:])
                    all_paralogs_for_reference.append(hit)
                    samples_with_paralogs.add(sample)
                    paralog_found: bool = True
            else:
                current_locus[sample] = hit
                samples_current_locus.add(sample)
        count += 1
    all_paralogs_for_reference = score_samples(all_paralogs_for_reference)
    exons: Set[str] = set()
    with open(
        os.path.join(
            path_to_data,
            "41detected_par",
            f"new_reference_for_HybPhyloMaker_div_{paralog_min_divergence}_{paralog_max_divergence}.fas",
        ),
        "w",
    ) as new_reference:
        for hit in all_paralogs_for_reference:
            sample = hit.split()[8]
            if sample not in blocklist:
                if exon(hit) not in exons:
                    name_of_locus = (
                        exon(hit)
                        .replace("exon", "Contig")
                        .replace("Exon", "Contig")
                        .replace("contig", "Contig")
                        .replace("_", "")
                        .replace("Contig", "_Contig_")
                    )
                    new_reference.write(
                        f">Assembly_{name_of_locus}_{sample}_{contig(hit)}\n{hit.split()[9]}\n"
                    )
                    exons.add(exon(hit))
        logger.info("New reference created!\n")


def refine_phasing(
    path_to_data, paralog_min_divergence, paralog_max_divergence, logger
):
    """"""

    with open(
        os.path.join(
            path_to_data,
            "41detected_par",
            f"new_reference_for_HybPhyloMaker_div_{paralog_min_divergence}_{paralog_max_divergence}.fas",
        ),
    ) as reference_to_check:
        new_ref_parsed = SeqIO.to_dict(
            SeqIO.parse(reference_to_check, "fasta")
        )
        all_seqs = set(new_ref_parsed.keys())
        initial_dict_to_check = dict()
        for seq in all_seqs:
            locus = seq.split("_")[1]
            exon = "_".join(seq.split("_")[2:4])
            name_contig = "_".join(seq.split("_")[4:])
            contig = (
                name_contig,
                new_ref_parsed[f"Assembly_{locus}_{exon}_{name_contig}"],
            )
            if locus not in initial_dict_to_check.keys():
                initial_dict_to_check[locus] = dict()
                initial_dict_to_check[locus][exon] = contig
            else:
                initial_dict_to_check[locus][exon] = contig
    final_dict_to_check = dict()
    for key1 in initial_dict_to_check.keys():
        if key1[-4:] != "para":
            final_dict_to_check[key1] = dict()
            final_dict_to_check[key1]["non-para"] = dict()
            final_dict_to_check[key1]["non-para"]["all"] = list()
            for key2 in initial_dict_to_check[key1]:
                final_dict_to_check[key1]["non-para"][key2] = initial_dict_to_check[
                    key1
                ][key2]
                final_dict_to_check[key1]["non-para"]["all"].append(
                    initial_dict_to_check[key1][key2][0]
                )
            if f"{key1}para" in initial_dict_to_check.keys():
                final_dict_to_check[key1]["para"] = dict()
                final_dict_to_check[key1]["para"]["all"] = list()

                for key2 in initial_dict_to_check[f"{key1}para"]:
                    final_dict_to_check[key1]["para"][key2] = initial_dict_to_check[
                        f"{key1}para"
                    ][key2]
                    final_dict_to_check[key1]["para"]["all"].append(
                        initial_dict_to_check[f"{key1}para"][key2][0]
                    )
    checked_dict = copy.deepcopy(final_dict_to_check)
    warn = list()
    for key in final_dict_to_check.keys():
        for item in [
            x for x in final_dict_to_check[key]["non-para"].keys() if x != "all"
        ]:
            if "para" in final_dict_to_check[key].keys():
                if (
                    final_dict_to_check[key]["non-para"][item][0]
                    in final_dict_to_check[key]["para"]["all"]
                ):
                    swap1 = final_dict_to_check[key]["non-para"][item]
                    checked_dict[key]["para"][item] = swap1
                    if item not in final_dict_to_check[key]["para"].keys():
                        del checked_dict[key]["non-para"][item]
                        warn.append(
                            f"Moving {swap1[0]} to a paralog for {key} {item}\n"
                        )
                    else:
                        swap2 = final_dict_to_check[key]["para"][item]
                        checked_dict[key]["non-para"][item] = swap2
                        warn.append(
                            f"Swaping {swap1[0]} and {swap2[0]} for {key} {item}\n"
                        )

    if len(warn) > 0:
        with open(
            os.path.join(path_to_data, "41detected_par", "warnings.txt"), "w"
        ) as warnings:
            warnings.write(
                "Following genes seems to be phased improperly based on similarity to reference. "
                "Refining attempted.\n"
            )
            for x in warn:
                warnings.write(x)
        with open(
            os.path.join(
                path_to_data,
                "41detected_par",
                f"refined_new_reference_for_HybPhyloMaker_div_{paralog_min_divergence}_{paralog_max_divergence}.fas",
            ),
            "w",
        ) as reference_to_write:
            for key in sorted(list(checked_dict.keys())):
                for item in [
                    x
                    for x in sorted(list(checked_dict[key]["non-para"].keys()))
                    if x != "all"
                ]:
                    name = f"Assembly_{key}_{item}_{checked_dict[key]['non-para'][item][0]}"
                    sequence = str(checked_dict[key]["non-para"][item][1].seq)
                    reference_to_write.write(f">{name}\n{sequence}\n")
                if "para" in checked_dict[key].keys():
                    for item in [
                        x
                        for x in sorted(list(checked_dict[key]["para"].keys()))
                        if x != "all"
                    ]:
                        name = f"Assembly_{key}para_{item}_{checked_dict[key]['para'][item][0]}"
                        sequence = str(checked_dict[key]["para"][item][1].seq)
                        reference_to_write.write(f">{name}\n{sequence}\n")
        with open(
            os.path.join(
                path_to_data,
                "41detected_par",
                f"refined_new_reference_for_HybPhyloMaker_div_{paralog_min_divergence}_{paralog_max_divergence}.fas",
            ),
        ) as reference_to_check:
            new_ref_parsed = SeqIO.to_dict(
                SeqIO.parse(reference_to_check, "fasta")
            )
            all_seqs = set(new_ref_parsed.keys())
            initial_dict_to_check = dict()
            for seq in all_seqs:
                locus = seq.split("_")[1]
                exon = "_".join(seq.split("_")[2:4])
                name_contig = "_".join(seq.split("_")[4:])
                contig = (
                    name_contig,
                    new_ref_parsed[f"Assembly_{locus}_{exon}_{name_contig}"],
                )
                if locus not in initial_dict_to_check.keys():
                    initial_dict_to_check[locus] = dict()
                    initial_dict_to_check[locus][exon] = contig
                else:
                    initial_dict_to_check[locus][exon] = contig
        final_dict_to_check = dict()
        for key1 in initial_dict_to_check.keys():
            if key1[-4:] != "para":
                final_dict_to_check[key1] = dict()
                final_dict_to_check[key1]["non-para"] = dict()
                final_dict_to_check[key1]["non-para"]["all"] = list()
                for key2 in initial_dict_to_check[key1]:
                    final_dict_to_check[key1]["non-para"][key2] = initial_dict_to_check[
                        key1
                    ][key2]
                    final_dict_to_check[key1]["non-para"]["all"].append(
                        initial_dict_to_check[key1][key2][0]
                    )
                if f"{key1}para" in initial_dict_to_check.keys():
                    final_dict_to_check[key1]["para"] = dict()
                    final_dict_to_check[key1]["para"]["all"] = list()

                    for key2 in initial_dict_to_check[f"{key1}para"]:
                        final_dict_to_check[key1]["para"][key2] = initial_dict_to_check[
                            f"{key1}para"
                        ][key2]
                        final_dict_to_check[key1]["para"]["all"].append(
                            initial_dict_to_check[f"{key1}para"][key2][0]
                        )
        checked_dict = copy.deepcopy(final_dict_to_check)
        warn = list()
        for key in final_dict_to_check.keys():
            for item in [
                x for x in final_dict_to_check[key]["non-para"].keys() if x != "all"
            ]:
                if "para" in final_dict_to_check[key].keys():
                    if (
                        final_dict_to_check[key]["non-para"][item][0]
                        in final_dict_to_check[key]["para"]["all"]
                    ):
                        swap1 = final_dict_to_check[key]["non-para"][item]
                        checked_dict[key]["para"][item] = swap1
                        if item not in final_dict_to_check[key]["para"].keys():
                            del checked_dict[key]["non-para"][item]
                            warn.append(f"{key} {item}\n")
                        else:
                            swap2 = final_dict_to_check[key]["para"][item]
                            checked_dict[key]["non-para"][item] = swap2
                            warn.append(f"{key} {item}\n")
        if len(warn) > 0:
            with open(
                os.path.join(path_to_data, "41detected_par", "warnings.txt"), "a"
            ) as warnings:
                warnings.write(
                    "\nFollowing genes and exons seems to be impossible to phase properly with the current algorithms. "
                    "Consider removing.\n"
                )
                for x in warn:
                    warnings.write(x)
                warnings.write(
                    "Pay closer attention to the loci mentioned above. It is worth having a look at the "
                    "alignments for particular exons in aln_orth_par/"
                )
        with open(
            os.path.join(
                path_to_data,
                "41detected_par",
                f"refined_new_reference_for_HybPhyloMaker_div_{paralog_min_divergence}_{paralog_max_divergence}.fas",
            ),
            "w",
        ) as reference_to_write:
            for key in sorted(list(checked_dict.keys())):
                for item in [
                    x
                    for x in sorted(list(checked_dict[key]["non-para"].keys()))
                    if x != "all"
                ]:
                    name = f"Assembly_{key}_{item}_{checked_dict[key]['non-para'][item][0]}"
                    sequence = str(checked_dict[key]["non-para"][item][1].seq)
                    reference_to_write.write(f">{name}\n{sequence}\n")
                if "para" in checked_dict[key].keys():
                    for item in [
                        x
                        for x in sorted(list(checked_dict[key]["para"].keys()))
                        if x != "all"
                    ]:
                        name = f"Assembly_{key}para_{item}_{checked_dict[key]['para'][item][0]}"
                        sequence = str(checked_dict[key]["para"][item][1].seq)
                        reference_to_write.write(f">{name}\n{sequence}\n")
    else:
        with open(
            os.path.join(path_to_data, "41detected_par", "warnings.txt"), "w"
        ) as warnings:
            warnings.write("No warnings generated.\n")


def write_paralog_stats(
    path_to_data,
    paralog_min_divergence,
    paralog_max_divergence,
    paralog_statistic,
    probes,
    logger,
):
    """"""

    with open(
        os.path.join(
            path_to_data,
            "41detected_par",
            f"paralog_statistics_div_{paralog_min_divergence}_{paralog_max_divergence}.tsv",
        ),
        "w",
    ) as par_stat:
        for sample in sorted(list(paralog_statistic.keys())):
            par_stat.write(f"{sample}\t{str(len(paralog_statistic[sample]))}\n")
    with open(probes) as probe_loci:
        all_probe_exons = SeqIO.to_dict(
            SeqIO.parse(probe_loci, "fasta")
        ).keys()
        all_loci: Set[str] = set()
        for key in all_probe_exons:
            all_loci.add(key.split("-")[1].split("_")[0])
    with open(
        os.path.join(
            path_to_data,
            "41detected_par",
            f"locus_statistics_div_{paralog_min_divergence}_{paralog_max_divergence}.tsv",
        ),
        "w",
    ) as loci_par_stat:
        loci_par_stat.write(r"sample\locus")
        for locus in sorted(list(all_loci)):
            loci_par_stat.write(f"\t{locus}")
        loci_par_stat.write("\n")
        for sample in sorted(list(paralog_statistic.keys())):
            loci_par_stat.write(sample)
            for locus in sorted(list(all_loci)):
                if locus in paralog_statistic[sample]:
                    loci_par_stat.write("\tYes")
                else:
                    loci_par_stat.write("\tN/A")
            loci_par_stat.write("\n")
