import sys
import copy
from typing import Set, Dict, List

from ParalogWizard_Functions import sort_hit_table_ident, exon, contig


def score_samples(list_with_hits: List[str]) -> List[str]:
    sort_hit_table_ident(list_with_hits, 'exon')
    # Add counters to hits for each exon
    counter: int = 0
    current_exon: str = ''
    for hit in list_with_hits:
        if exon(hit) != current_exon:
            counter = 1
            current_exon = exon(hit)
        else:
            counter += 1
        list_with_hits[list_with_hits.index(hit)] = f'{hit}\t{counter}'
    # Prepare dictionary for counting scores for contigs in each locus
    for_samples_scores = dict()
    samples: Set[str] = set()
    loci_samples: Set[str] = set()
    for hit in list_with_hits:
        sample: str = hit.split()[8]
        locus: str = exon(hit).split('_')[0]
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
        locus: str = exon(hit).split('_')[0]
        spades_contig: str = contig(hit)
        for_samples_scores[sample][locus][0][spades_contig].append(int(hit.split()[-1]))
        for_samples_scores[sample][locus][1].add(exon(hit))
    samples_scores = copy.deepcopy(for_samples_scores)
    for sample in for_samples_scores.keys():
        for locus in for_samples_scores[sample].keys():
            samples_scores[sample][locus][1]: int = len(for_samples_scores[sample][locus][1])
            sample_locus_scores: List[float] = []
            for spades_contig in for_samples_scores[sample][locus][0].keys():
                samples_scores[sample][locus][0][spades_contig] = \
                    sum(for_samples_scores[sample][locus][0][spades_contig]) / \
                    len(for_samples_scores[sample][locus][0][spades_contig])
                sample_locus_scores.append(samples_scores[sample][locus][0][spades_contig])
            sample_locus_score: float = sum(sample_locus_scores) / len(sample_locus_scores)
            samples_scores[sample][locus].append(sample_locus_score)
    # Adding scores to each hit
    for i in range(0, len(list_with_hits)):
        hit: str = list_with_hits[i]
        locus: str = exon(hit).split('_')[0]
        spades_contig: str = contig(hit)
        sample: str = hit.split()[8]
        list_with_hits[i] = f'{hit}\t{str(samples_scores[sample][locus][1])}\t{str(samples_scores[sample][locus][2])}' \
                            f'\t{str(samples_scores[sample][locus][0][spades_contig])}'
    list_with_hits.sort(key=lambda x: float(x.split()[-1]))
    list_with_hits.sort(key=lambda x: x.split()[8])
    list_with_hits.sort(key=lambda x: float(x.split()[-2]))
    list_with_hits.sort(key=lambda x: float(x.split()[-3]), reverse=True)
    list_with_hits.sort(key=lambda x: exon(x))
    return ['\t'.join(x.split()[:10]) for x in list_with_hits]


def create_reference_wo_paralogs():
    print('Creating new reference...')
    exons: Set[str] = set()
    with open(f'{path_to_data_HPM}/exons/new_reference_for_HybPhyloMaker.fas', 'w') as new_reference:
        for hit in all_hits_for_reference_scored:
            sample: str = hit.split()[8]
            if sample not in blacklist and float(hit.split()[3]) >= 75:
                if exon(hit) not in exons:
                    name_of_locus: str = exon(hit).replace('exon', 'Contig').replace('Exon', 'Contig') \
                        .replace('contig', 'Contig').replace('_', '').replace('Contig', '_Contig_')
                    new_reference.write(f'>Assembly_{name_of_locus}_{sample}_{contig(hit)}\n{hit.split()[9]}\n')
                    exons.add(exon(hit))
    print('New reference created!\n')


def create_reference_w_paralogs():
    with open(f'{path_to_data_HPM}/exons/aln_orth_par/pairwise_distances.txt') as distances:
        pairwise_distances: Dict[str, float] = dict()
        for line in distances.read().splitlines():
            pairwise_distances[f'{line.split()[0]}_{line.split()[2]}'] = float(line.split()[1])
            pairwise_distances[f'{line.split()[2]}_{line.split()[0]}'] = float(line.split()[1])
    print('Creating new reference...')
    all_paralogs_for_reference = []
    exons: Set[str] = set()
    count: int = 0
    paralog_found: bool = bool()
    current_best: str = ''
    current_locus: Dict[str, str] = dict()
    samples_current_locus: Set[str] = set()
    samples_with_paralogs: Set[str] = set()
    paralog_statistic = dict()
    set_of_samples = set()
    for hit in all_hits_for_reference_scored:
        sample: str = hit.split()[8]
        if sample not in set_of_samples:
            set_of_samples.add(sample)
            paralog_statistic[sample] = set()
        if sample not in blacklist and float(hit.split()[3]) >= 75:
            if exon(hit) not in exons:
                if not paralog_found and count != 0:
                    print(f"No paralog found for {current_best.split()[0].split('-')[1]}")
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
                    current_exonic_contig = f"{hit.split()[0].split('-')[1]}_" \
                                            f"{'_'.join(hit.split()[1].split('_')[1:])}_{sample}"
                    exonic_contigs_to_compare = f"{current_locus[sample].split()[0].split('-')[1]}_" \
                                                f"{'_'.join(current_locus[sample].split()[1].split('_')[1:])}" \
                                                f"_{sample}"
                    if pairwise_distances[f'{current_exonic_contig}_{exonic_contigs_to_compare}'] > \
                            paralog_min_divergence:
                        print(f"Paralog detected for {hit.split()[0].split('-')[1]} in {sample}")
                        paralog_statistic[sample].add(hit.split()[1].split('_')[0])
                        all_paralogs_for_reference.append(current_locus[sample])
                        name_of_locus_para: str = '_'.join([exon(hit).split('_')[0] + '_para'] +
                                                           exon(hit).split('_')[1:])
                        name_of_locus_para = hit.split()[0].split('-')[0] + '-' + name_of_locus_para
                        hit = '\t'.join([name_of_locus_para] + hit.split()[1:])
                        all_paralogs_for_reference.append(hit)
                        samples_with_paralogs.add(sample)
                        paralog_found: bool = True
                else:
                    current_locus[sample] = hit
                    samples_current_locus.add(sample)
            count += 1
    all_paralogs_for_reference = score_samples(all_paralogs_for_reference)
    exons: Set[str] = set()
    with open(f'{path_to_data_HPM}/exons/new_reference_for_HybPhyloMaker_div_{paralog_min_divergence}.fas', 'w') \
            as new_reference:
        for hit in all_paralogs_for_reference:
            sample = hit.split()[8]
            if sample not in blacklist and float(hit.split()[3]) >= 75:
                if exon(hit) not in exons:
                    name_of_locus = exon(hit).replace('exon', 'Contig').replace('Exon', 'Contig') \
                        .replace('contig', 'Contig').replace('_', '').replace('Contig', '_Contig_')
                    new_reference.write(
                        f'>Assembly_{name_of_locus}_{sample}_{contig(hit)}\n{hit.split()[9]}\n')
                    exons.add(exon(hit))
    with open(f'{path_to_data_HPM}/exons/paralog_statistics_div_{paralog_min_divergence}.tsv', 'w') as par_stat:
        for sample in sorted(list(paralog_statistic.keys())):
            par_stat.write(f'{sample}\t{str(len(paralog_statistic[sample]))}\n')
    print('New reference created!\n')


if __name__ == '__main__':
    path_to_data_HPM = sys.argv[1]
    blacklist: set = set([x.strip() for x in sys.argv[2].split(',')])
    paralogs_bool = sys.argv[3]
    paralog_min_divergence: float = float(sys.argv[4].strip())
    with open(f'{path_to_data_HPM}/exons/all_hits.txt') as all_hits:
        all_hits_for_reference: List[str] = [x[:-1] for x in all_hits.readlines()]
    all_hits_for_reference_scored = score_samples(all_hits_for_reference)
    if paralogs_bool != 'yes':
        create_reference_wo_paralogs()
    else:
        create_reference_w_paralogs()
