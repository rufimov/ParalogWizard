import glob
import os
import sys
from typing import List, Dict, Set

import Bio
import Bio.Application
import matplotlib.pyplot
import numpy
from Bio import SeqRecord, SeqIO
from Bio.Align.Applications import MafftCommandline
from Bio.Alphabet import generic_dna
from Bio.Phylo.Applications import _Fasttree
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ParalogWizard_Functions import (
    sort_hit_table_ident,
    exon,
    contig,
    get_distance_matrix,
    get_plot,
)


def build_alignments():
    print("Building individual exon alignments...")
    with open(path_to_data_HPM + "/exons/all_hits.txt") as all_hits:
        all_hits_for_reference: List[str] = [x[:-1] for x in all_hits.readlines()]
    sort_hit_table_ident(all_hits_for_reference, "exon")
    exons: Dict[str, List[SeqRecord]] = dict()
    for hit in all_hits_for_reference:
        if exon(hit) not in set(exons.keys()):
            exons[exon(hit)]: List[SeqRecord] = [
                SeqRecord(
                    Seq(hit.split()[-1], generic_dna),
                    id=contig(hit) + "_" + hit.split()[-2],
                    description="",
                )
            ]
        else:
            new_list: List[SeqRecord] = exons[exon(hit)]
            new_list.append(
                SeqRecord(
                    Seq(hit.split()[-1], generic_dna),
                    id=contig(hit) + "_" + hit.split()[-2],
                    description="",
                )
            )
            exons[exon(hit)] = new_list
    os.makedirs(path_to_data_HPM + "/exons/aln_orth_par", exist_ok=True)
    for key in exons.keys():
        for record in exons[key]:
            new_id = f"{key}_N_{record.id.split('_N_')[1]}"
            record.id = new_id
    for key in exons.keys():
        SeqIO.write(
            exons[key], f"{path_to_data_HPM}/exons/aln_orth_par/{key}.fasta", "fasta"
        )
        stdout, stderr = MafftCommandline(
            input=f"{path_to_data_HPM}/exons/aln_orth_par/{key}.fasta",
            adjustdirectionaccurately=True,
            auto=True,
            thread=4,
        )()
        with open(
            f"{path_to_data_HPM}/exons/aln_orth_par/{key}.mafft.fasta", "w"
        ) as aligned:
            aligned.write(stdout.replace(">_R_", ">"))
        try:
            _Fasttree.FastTreeCommandline(
                "fasttreemp",
                nt=True,
                gtr=True,
                input=f"{path_to_data_HPM}/exons/aln_orth_par/{key}.mafft.fasta",
                out=f"{path_to_data_HPM}/exons/aln_orth_par/{key}.mafft.fasta.tre",
            )()
        except Bio.Application.ApplicationError:
            _Fasttree.FastTreeCommandline(
                "fasttree",
                nt=True,
                gtr=True,
                input=f"{path_to_data_HPM}/exons/aln_orth_par/{key}.mafft.fasta",
                out=f"{path_to_data_HPM}/exons/aln_orth_par/{key}.mafft.fasta.tre",
            )()

    print("Done\n")


def estimate_divergence():
    print("Estimating divergence of paralogs...")
    divergency_distribution: List[float] = []
    divergencies_to_write: List[str] = []
    matplotlib.use("Agg")
    fig, axis = matplotlib.pyplot.subplots(figsize=(15, 10))
    for file in sorted(
        glob.glob(path_to_data_HPM + "/exons/aln_orth_par/*.mafft.fasta")
    ):
        get_distance_matrix(
            file,
            divergency_distribution,
            divergencies_to_write,
            blacklist,
            axis,
        )
    end = axis.get_xlim()[1]
    end = numpy.round(end, 0)
    axis.xaxis.set_ticks(numpy.arange(0, end, 1))
    axis.set_xlabel("Divergence (%)")
    axis.set_ylabel("Frequency")
    fig.savefig(
        f"{path_to_data_HPM}/exons/aln_orth_par/individual_distributions.png",
        dpi=300,
        format="png",
    )
    fig.savefig(
        f"{path_to_data_HPM}/exons/aln_orth_par/individual_distributions.svg",
        dpi=300,
        format="svg",
    )
    matplotlib.pyplot.close(fig)
    with open(
        path_to_data_HPM + "/exons/aln_orth_par/pairwise_distances.txt", "w"
    ) as divergency_distribution_to_write:
        for i in divergencies_to_write:
            divergency_distribution_to_write.write(str(i) + "\n")
    divergency_distribution_array: numpy.ndarray = numpy.array(
        [[x] for x in divergency_distribution]
    )
    get_plot(
        path_to_data_HPM + "/exons/aln_orth_par/",
        "pairwise_distances_distribution_1_comp",
        divergency_distribution_array,
        1,
    )
    get_plot(
        path_to_data_HPM + "/exons/aln_orth_par/",
        "pairwise_distances_distribution_2_comp",
        divergency_distribution_array,
        2,
    )
    get_plot(
        path_to_data_HPM + "/exons/aln_orth_par/",
        "pairwise_distances_distribution_3_comp",
        divergency_distribution_array,
        3,
    )

    print("Done\n")


if __name__ == "__main__":
    path_to_data_HPM: str = sys.argv[1].strip()
    blacklist: Set[str] = set([x.strip() for x in sys.argv[2].split(",")])
    build_alignments()
    estimate_divergence()
