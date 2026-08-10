# -------------------------------------------------------------------------------
#
# Filename: IterPol.py
#
# Version: 0.2
#
# Copyright 2024 Étienne Léveillé-Bourret under the terms of the GNU
# General Public License as published by the Free Software Foundation,
# either version 3 of the License, or (at your option) any later version.
#
#
# This script is used to phase the alleles of hybrids and allopolyploids
# represented by a single consensus sequence (with IUPAC ambiguities at
# sites where the alleles/subgenomes differ) in an alignment of one or
# many concatenated loci. At least one of the two parental lineages of
# each hybrid or allopolyploid need to be present in the alignment for
# maximal accuracy. The script can use external programs for phylogenetic
# inference (PAUP*, RAxML and ASTRAL).
#
# Usage: python3 IterPol.py --help (will show all available options)
#
# -------------------------------------------------------------------------------


import argparse
import subprocess
import sys
import os
import glob
from Bio import AlignIO, Phylo
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.BaseTree import Clade
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
import numpy as np

# Define IUPAC ambiguity codes
IUPAC_AMBIGUITIES = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"G", "C"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},
    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},
    "N": set(),  # Handle Ns as empty sets since they usually mean "no data"
    "-": set(),  # Handle gaps as empty sets
}

# Create a reverse lookup dictionary for IUPAC ambiguity codes
IUPAC_REVERSE_LOOKUP = {"".join(sorted(v)): k for k, v in IUPAC_AMBIGUITIES.items()}


def proportional_match(set_a, set_b):
    intersection = set_a.intersection(set_b)
    if not intersection:
        return 0
    return len(intersection) / (len(set_a) * len(set_b))


def proportional_distance_function(seq1, seq2):
    total_distance = 0
    valid_pairs = 0
    for a, b in zip(seq1, seq2):
        set_a = IUPAC_AMBIGUITIES.get(a, {a})
        set_b = IUPAC_AMBIGUITIES.get(b, {b})
        if set_a and set_b:  # Exclude gaps or undefined characters
            match_score = proportional_match(set_a, set_b)
            total_distance += 1 - match_score
            valid_pairs += 1
    return total_distance / valid_pairs if valid_pairs > 0 else 0


class ProportionalDistanceCalculator:
    def get_distance(self, msa):
        n = len(msa)
        matrix = []
        for i in range(n):
            row = []
            for j in range(i + 1):  # Only fill lower-triangular matrix
                if i == j:
                    row.append(0.0)
                else:
                    dist = proportional_distance_function(msa[i], msa[j])
                    row.append(dist)
            matrix.append(row)
        return DistanceMatrix(names=[record.id for record in msa], matrix=matrix)


def read_alignment(phylip_file):
    # Read the alignment from the Phylip file and convert to uppercase
    alignment = AlignIO.read(phylip_file, "phylip-relaxed")
    for record in alignment:
        record.seq = record.seq.upper()

    # Extract list of sequence names
    sequence_names = [seq.id for seq in alignment]

    return alignment, sequence_names


def estimate_nj_tree(phylip_file, tree_file):
    # Read the alignment from the phylip file
    alignment, _ = read_alignment(phylip_file)

    # Calculate the distance matrix using the custom distance calculator
    calculator = ProportionalDistanceCalculator()
    distance_matrix = calculator.get_distance(alignment)

    # Construct the tree
    constructor = DistanceTreeConstructor()
    nj_tree = constructor.nj(distance_matrix)

    # Write the tree to a file
    Phylo.write(nj_tree, tree_file, "newick")

    print(f"Neighbor-Joining tree written to {tree_file}")


def run_paup(phylip_file, tree_file, paup_path=None):
    # Create a PAUP* command script
    paup_script = f"""#NEXUS
    begin PAUP;
    toNexus format=relPhylip fromFile={phylip_file} toFile=tmp.IterPol.alignment.nex;
    execute tmp.IterPol.alignment.nex;
    hs mulTrees=NO;
    saveTrees trees=firstOnly brLens=YES file={tree_file} format=newick replace=YES;
    quit;
    end;
    """

    paup_script_file = "tmp.IterPol.paup_script.nex"
    with open(paup_script_file, "w") as f:
        f.write(paup_script)

    # Run PAUP*
    process = subprocess.Popen(
        [paup_path, paup_script_file],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    # Print stdout and stderr in real-time
    for stdout_line in iter(process.stdout.readline, ""):
        print(stdout_line.strip())
    for stderr_line in iter(process.stderr.readline, ""):
        print(stderr_line.strip(), file=sys.stderr)

    sys.stdout.flush()

    # Wait for the process to complete
    process.communicate()

    print(f"PAUP* parsimony tree written to {tree_file}")


def run_raxml(phylip_file, tree_file, threads, raxml_path=None):
    raxml_command = [
        raxml_path,
        "-f",
        "d",
        "-T",
        str(threads),
        "-s",
        phylip_file,
        "-n",
        "tmp.IterPol",
        "-m",
        "GTRCAT",  # Change this according to your model of choice
        "-p",
        "666",  # Random seed for parsimony inferences
    ]

    # run RAxML
    process = subprocess.Popen(
        raxml_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )

    # Print stdout and stderr in real-time
    for stdout_line in iter(process.stdout.readline, ""):
        print(stdout_line.strip())
    for stderr_line in iter(process.stderr.readline, ""):
        print(stderr_line.strip(), file=sys.stderr)

    sys.stdout.flush()

    # Wait for the process to complete
    process.communicate()

    # Rename the output tree file
    os.rename("RAxML_bestTree.tmp.IterPol", tree_file)
    print(f"RAxML tree generation complete. Output saved to {tree_file}")

    # Remove temporary files
    files_to_remove = glob.glob("*.tmp.IterPol")
    for file in files_to_remove:
        os.remove(file)


def run_raxml_ng(phylip_file, tree_file, threads, raxml_path=None):
    raxml_command = [
        raxml_path,
        "--threads",
        str(threads),
        "--msa",
        phylip_file,
        "--prefix",
        "tmp.IterPol",
        "--model",
        "GTR+G",  # Change this according to your model of choice
        "--seed",
        "666",  # Random seed for parsimony inferences
    ]

    # run RAxML
    process = subprocess.Popen(
        raxml_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )

    # Print stdout and stderr in real-time
    for stdout_line in iter(process.stdout.readline, ""):
        print(stdout_line.strip())
    for stderr_line in iter(process.stderr.readline, ""):
        print(stderr_line.strip(), file=sys.stderr)

    sys.stdout.flush()

    # Wait for the process to complete
    process.communicate()

    # Rename the output tree file
    os.rename("tmp.IterPol.raxml.bestTree", tree_file)
    print(f"RAxML tree generation complete. Output saved to {tree_file}")

    # Remove temporary files
    files_to_remove = glob.glob("tmp.IterPol*")
    for file in files_to_remove:
        os.remove(file)


def run_caster(phylip_file, tree_file, threads, caster_path=None):
    caster_command = [
        caster_path,
        "-t",
        str(threads),
        "-i",
        phylip_file,
        "-o",
        "tmp.IterPol.tre",
    ]

    # run CASTER
    process = subprocess.Popen(
        caster_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )

    # Print stdout and stderr in real-time
    for stdout_line in iter(process.stdout.readline, ""):
        print(stdout_line.strip())
    for stderr_line in iter(process.stderr.readline, ""):
        print(stderr_line.strip(), file=sys.stderr)

    sys.stdout.flush()

    # Wait for the process to complete
    process.communicate()

    # Rename the output tree file
    os.rename("tmp.IterPol.tre", tree_file)
    print(f"RAxML tree generation complete. Output saved to {tree_file}")

    # Remove temporary files
    files_to_remove = glob.glob("tmp.IterPol*")
    for file in files_to_remove:
        os.remove(file)


def estimate_tree(
    phylip_file,
    tree_file,
    method,
    threads=1,
    paup_path=None,
    raxml_path=None,
    caster_path=None,
):
    if method == "nj":
        estimate_nj_tree(phylip_file, tree_file)
    elif method == "paup":
        run_paup(phylip_file, tree_file, paup_path)
    elif method == "raxml":
        run_raxml(phylip_file, tree_file, threads, raxml_path)
    elif method == "raxml-ng":
        run_raxml_ng(phylip_file, tree_file, threads, raxml_path)
    elif method == "caster":
        run_caster(phylip_file, tree_file, threads, caster_path)
    else:
        raise SystemExit(
            f"ERROR: Unsupported method '{method}'. Choose 'nj', 'paup' or 'raxml'."
        )


def transform_to_cladogram(tree):
    for clade in tree.find_clades():
        if clade.is_terminal():
            clade.branch_length = 0  # Set external branches to zero length
        else:
            clade.branch_length = (
                1  # Internal nodes can have branch length 1 or as needed
            )
    return tree


def generate_patristic_matrix(tree, sequence_names, matrix_output_file):
    # Generate patristic matrix in the order of sequence_names
    patristic_matrix = []
    for seq_id1 in sequence_names:
        row = []
        for seq_id2 in sequence_names:
            # Calculate distance between seq_id1 and seq_id2 along the tree branches
            distance = tree.distance(seq_id1, seq_id2)
            row.append(distance)
        patristic_matrix.append(row)

    # Write the new distance matrix to a file
    with open(matrix_output_file, "w") as f:
        f.write("\t" + "\t".join(sequence_names) + "\n")
        for i, row in enumerate(patristic_matrix):
            f.write(sequence_names[i] + "\t" + "\t".join(map(str, row)) + "\n")
    print(f"Patristic matrix written to {matrix_output_file}")

    return patristic_matrix


def polarize_sequence(reference, sequence):
    polarized_sequence = []
    for ref_base, seq_base in zip(reference, sequence):
        if ref_base in {"N", "-"} or seq_base in {"N", "-"}:
            polarized_sequence.append(seq_base)
        else:
            ref_set = IUPAC_AMBIGUITIES.get(ref_base, {ref_base})
            seq_set = IUPAC_AMBIGUITIES.get(seq_base, {seq_base})
            # Get the bases in the sequence set that are not in the reference set
            non_ref_bases = seq_set - ref_set
            if non_ref_bases:
                non_ref_bases_str = "".join(sorted(non_ref_bases))
                polarized_sequence.append(
                    IUPAC_REVERSE_LOOKUP.get(non_ref_bases_str, seq_base)
                )
            else:
                polarized_sequence.append(seq_base)
    return "".join(polarized_sequence)


def polarize_alignment(alignment, closest_indices):
    polarized_sequences = []
    reverse_polarized_sequences = []

    for i, record in enumerate(alignment):
        original_sequence = record.seq
        polarized_seq = original_sequence

        # Successively polarize using closest sequences until fully polarized
        for index in closest_indices[i]:
            reference = alignment[index].seq
            polarized_seq = polarize_sequence(reference, polarized_seq)

        polarized_sequences.append(polarized_seq)

        # Using the fully polarized sequences (=allele1), polarize original sequences to allele2
        reverse_polarized_seq = polarize_sequence(polarized_seq, original_sequence)
        reverse_polarized_sequences.append(reverse_polarized_seq)

    return polarized_sequences, reverse_polarized_sequences


def write_polarized_alignment(
    alignment,
    polarized_sequences,
    reverse_polarized_sequences,
    output_file,
    reverse=False,
):
    with open(output_file, "w") as f:
        # Write the number of sequences and their length
        if reverse:
            f.write(f" {len(alignment) * 2} {len(alignment[0].seq)}")
        else:
            f.write(f" {len(alignment)} {len(alignment[0].seq)}")

        # Write the sequence data
        for record, polarized_seq, reverse_polarized_seq in zip(
            alignment, polarized_sequences, reverse_polarized_sequences
        ):
            f.write(f"\n{record.id}_allele1".ljust(50) + f" {polarized_seq}")
            if reverse:
                f.write(
                    f"\n{record.id}_allele2".ljust(50) + f" {reverse_polarized_seq}"
                )

        # Add a final linebreak at the end
        f.write(f"\n")

    print(f"Polarized alignment written to {output_file}")


def import_and_prune_allele_tree(tree_file, keep_allele1=True):
    # Read the tree from the file
    tree = Phylo.read(tree_file, "newick")

    # List to store terminals to prune
    terminals_to_prune = []

    # Determine which terminals to prune based on keep_allele1
    for terminal in tree.get_terminals():
        if keep_allele1:
            if terminal.name.endswith("_allele1"):
                terminal.name = terminal.name[:-8]
            elif terminal.name.endswith("_allele2"):
                terminals_to_prune.append(terminal)
        else:
            if terminal.name.endswith("_allele1"):
                terminals_to_prune.append(terminal)
            elif terminal.name.endswith("_allele2"):
                terminal.name = terminal.name[:-8]

    # Prune the terminals marked for removal
    for terminal in terminals_to_prune:
        if terminal in tree.get_terminals():
            tree.prune(terminal)

    Phylo.write(tree, f"{tree_file}.pruned", "newick")

    return tree


def dichotomize_tree(tree):
    def dichotomize_clade(clade):
        if len(clade.clades) > 2:
            # Create new internal nodes to convert polychotomies to dichotomies
            while len(clade.clades) > 2:
                # Create a new internal node
                new_internal = Clade()
                # Attach two clades to this new internal node
                new_internal.clades.append(clade.clades.pop(0))
                new_internal.clades.append(clade.clades.pop(0))
                # Insert this new internal node back into the clade
                clade.clades.insert(0, new_internal)

        # Recursively process all child clades
        for child_clade in clade.clades:
            if not child_clade.is_terminal():
                dichotomize_clade(child_clade)

    # Start the dichotomization process from the root
    dichotomize_clade(tree.root)
    return tree


def BACKUP_merge_sister_alleles(
    tree, polarized_alignment_file, raw_alignment_file, new_alignment_file
):
    terminals_to_prune = set()
    new_polarized_seqs = []

    # Read the alignments
    polarized_seqs = AlignIO.read(polarized_alignment_file, "phylip-relaxed")
    original_seqs = AlignIO.read(raw_alignment_file, "phylip-relaxed")

    # Create a dictionary for quick lookup of original sequences
    original_seq_dict = {record.id: record.seq for record in original_seqs}

    # Make sure the tree is fully dichotomous
    dichotomous_tree = dichotomize_tree(tree)

    # Make sure each branch has a length
    for clade in tree.find_clades():
        if clade.branch_length is None:
            clade.branch_length = 0

    # Root at midpoint
    dichotomous_tree.root_at_midpoint()

    # Iterate through all clades in the tree
    for clade in dichotomous_tree.find_clades():
        if clade.is_terminal():
            # Get the sister clade of the terminal clade
            path = dichotomous_tree.get_path(clade)

            if len(path) > 1:
                parent_clade = path[-2]

                if len(parent_clade.clades) == 2:
                    sister_clade = (
                        parent_clade.clades[1]
                        if parent_clade.clades[0] == clade
                        else parent_clade.clades[0]
                    )

                    # Extract names of the terminals
                    clade_name, sister_name = clade.name, sister_clade.name

                    if clade_name and sister_name:
                        # Check if names differ only by suffixes '_allele1' and '_allele2'
                        print(f"Clade name is {clade_name}")
                        print(f"Sister name is {sister_name}")
                        if (
                            clade_name.endswith("_allele1")
                            and sister_name.endswith("_allele2")
                            and clade_name[:-8] == sister_name[:-8]
                        ):
                            # Use the raw sequence for the "_allele1" and set the name to "_unphased"
                            unphased_name = clade_name[:-8] + "_unphased"
                            raw_name = clade_name[:-8]
                            if raw_name in original_seq_dict:
                                new_polarized_seqs.append(
                                    SeqRecord(
                                        Seq(str(original_seq_dict[raw_name])),
                                        id=unphased_name,
                                        description="",
                                    )
                                )
                            terminals_to_prune.add(sister_clade)

                        elif (
                            clade_name.endswith("_allele2")
                            and sister_name.endswith("_allele1")
                            and clade_name[:-8] == sister_name[:-8]
                        ):
                            pass

                        else:
                            # Else if the allele1 and allele2 of a sample are not sisters, check how distant they are
                            distance_between_alleles = tree.distance(
                                clade_name[:-8] + "_allele1",
                                clade_name[:-8] + "_allele2",
                            )

                            if distance_between_alleles > 0:
                                for record in polarized_seqs:
                                    if record.id == clade_name:
                                        new_polarized_seqs.append(record)

                            elif clade_name.endswith("_allele1"):
                                # Use the raw sequence for the "_allele1" and set the name to "_unphased"
                                unphased_name = clade_name[:-8] + "_unphased"
                                raw_name = clade_name[:-8]
                                if raw_name in original_seq_dict:
                                    new_polarized_seqs.append(
                                        SeqRecord(
                                            Seq(str(original_seq_dict[raw_name])),
                                            id=unphased_name,
                                            description="",
                                        )
                                    )
                                terminals_to_prune.add(sister_clade)

    # Prune the terminals marked for removal
    for terminal in terminals_to_prune:
        if terminal in dichotomous_tree.get_terminals():
            dichotomous_tree.prune(terminal)

    # Convert new_polarized_seqs to a MultipleSeqAlignment object
    new_polarized_alignment = MultipleSeqAlignment(new_polarized_seqs)

    # Write the updated polarized alignment to a new file
    with open(new_alignment_file, "w") as output_handle:
        output_handle.write(
            f"{len(new_polarized_alignment)} {len(new_polarized_alignment[0].seq)}\n"
        )
        for record in new_polarized_alignment:
            output_handle.write(
                f"{record.id}".ljust(50) + f" {str(record.seq).replace(' ', '')}\n"
            )

    return tree


def merge_sister_alleles(
    tree, polarized_alignment_file, raw_alignment_file, new_alignment_file
):
    new_polarized_seqs = []

    # Read the alignments
    polarized_seqs = AlignIO.read(polarized_alignment_file, "phylip-relaxed")
    original_seqs = AlignIO.read(raw_alignment_file, "phylip-relaxed")

    # Create a dictionary for quick lookup of original sequences
    original_seq_dict = {record.id: record.seq for record in original_seqs}

    # Make sure the tree is fully dichotomous
    dichotomous_tree = dichotomize_tree(tree)

    # Make sure each branch has a length
    for clade in tree.find_clades():
        if clade.branch_length is None:
            clade.branch_length = 0

    # Iterate through all clades in the tree
    for clade in dichotomous_tree.find_clades():
        if clade.is_terminal():
            # !!!!!! CONTINUE FROM HERE !!!!!
            # Get the closest terminal to the target terminal
            path = dichotomous_tree.get_path(clade)

            if len(path) > 1:
                parent_clade = path[-2]

                if len(parent_clade.clades) == 2:
                    sister_clade = (
                        parent_clade.clades[1]
                        if parent_clade.clades[0] == clade
                        else parent_clade.clades[0]
                    )

                    # Extract names of the terminals
                    clade_name, sister_name = clade.name, sister_clade.name

                    if clade_name.endswith("_allele1"):
                        # Check if alleles are sisters
                        if clade_name and sister_name:
                            print(f"Clade name is {clade_name}.")
                            print(f"Sister name is {sister_name}.")
                            # Check if names differ only by suffixes '_allele1' and '_allele2'
                            if (
                                clade_name.endswith("_allele1")
                                and sister_name.endswith("_allele2")
                                and clade_name[:-8] == sister_name[:-8]
                            ):
                                print(
                                    f"The alleles of {clade_name[:-8]} are sisters, so keeping the unphased sequence."
                                )
                                # Use the raw sequence for the "_allele1" and set the name to "_unphased"
                                unphased_name = clade_name[:-8] + "_unphased"
                                raw_name = clade_name[:-8]
                                if raw_name in original_seq_dict:
                                    new_polarized_seqs.append(
                                        SeqRecord(
                                            Seq(str(original_seq_dict[raw_name])),
                                            id=unphased_name,
                                            description="",
                                        )
                                    )

                            else:
                                print(
                                    f"The alleles of {clade_name[:-8]} are NOT sisters, so calculating distance between them."
                                )
                                # Else if the allele1 and allele2 of a sample are not sisters, check how distant they are
                                distance_between_alleles = tree.distance(
                                    clade_name[:-8] + "_allele1",
                                    clade_name[:-8] + "_allele2",
                                )
                                print(
                                    f"Distance between the alleles of {clade_name[:-8]} is {distance_between_alleles}."
                                )

                                if distance_between_alleles > 0:
                                    print(
                                        f"The distance is greater than zero, so keeping the phased alleles of {clade_name[:-8]}."
                                    )
                                    # Write the polarized sequences (alleles 1 and 2) to the new_polarized_seqs object
                                    for record in polarized_seqs:
                                        if (
                                            record.id == clade_name[:-8] + "_allele1"
                                            or record.id == clade_name[:-8] + "_allele2"
                                        ):
                                            new_polarized_seqs.append(record)

                                else:
                                    print(
                                        f"The distance is zero, so keeping the unphased sequence for {clade_name[:-8]}."
                                    )
                                    # Use the raw sequence for the "_allele1" and set the name to "_unphased"
                                    unphased_name = clade_name[:-8] + "_unphased"
                                    raw_name = clade_name[:-8]
                                    if raw_name in original_seq_dict:
                                        new_polarized_seqs.append(
                                            SeqRecord(
                                                Seq(str(original_seq_dict[raw_name])),
                                                id=unphased_name,
                                                description="",
                                            )
                                        )

                        else:
                            # Alleles cannot be sister since the sister clade does not have a name
                            print(
                                f"The alleles of {clade_name[:-8]} are NOT sisters, so calculating distance between them."
                            )

                            # Else if the allele1 and allele2 of a sample are not sisters, check how distant they are
                            distance_between_alleles = tree.distance(
                                clade_name[:-8] + "_allele1",
                                clade_name[:-8] + "_allele2",
                            )
                            print(
                                f"Distance between the alleles of {clade_name[:-8]} is {distance_between_alleles}."
                            )

                            if distance_between_alleles > 0:
                                print(
                                    f"The distance is greater than zero, so keeping the phased alleles of {clade_name[:-8]}."
                                )
                                # Write the polarized sequences (alleles 1 and 2) to the new_polarized_seqs object
                                for record in polarized_seqs:
                                    if (
                                        record.id == clade_name[:-8] + "_allele1"
                                        or record.id == clade_name[:-8] + "_allele2"
                                    ):
                                        new_polarized_seqs.append(record)

                            else:
                                print(
                                    f"The distance is zero, so keeping the unphased sequence for {clade_name[:-8]}."
                                )
                                # Use the raw sequence for the "_allele1" and set the name to "_unphased"
                                unphased_name = clade_name[:-8] + "_unphased"
                                raw_name = clade_name[:-8]
                                if raw_name in original_seq_dict:
                                    new_polarized_seqs.append(
                                        SeqRecord(
                                            Seq(str(original_seq_dict[raw_name])),
                                            id=unphased_name,
                                            description="",
                                        )
                                    )

    # Convert new_polarized_seqs to a MultipleSeqAlignment object
    new_polarized_alignment = MultipleSeqAlignment(new_polarized_seqs)

    # Write the updated polarized alignment to a new file
    with open(new_alignment_file, "w") as output_handle:
        output_handle.write(
            f"{len(new_polarized_alignment)} {len(new_polarized_alignment[0].seq)}\n"
        )
        for record in new_polarized_alignment:
            output_handle.write(
                f"{record.id}".ljust(50) + f" {str(record.seq).replace(' ', '')}\n"
            )

    return tree


def iterative_polarization(
    raw_alignment_file,
    out_prefix,
    iterations,
    method,
    no_branch_length,
    threads=1,
    paup_path=None,
    raxml_path=None,
    caster_path=None,
):
    # Generate the starting tree and read the alignment
    print("Creating a tree based on the raw alignment.")
    raw_tree_file = f"{out_prefix}_raw.tre"
    alignment, sequence_names = read_alignment(raw_alignment_file)
    estimate_tree(
        raw_alignment_file,
        raw_tree_file,
        method,
        threads,
        paup_path,
        raxml_path,
        caster_path,
    )
    tree = Phylo.read(raw_tree_file, "newick")

    print(f"\nStarting iteration 0 of polarization and tree inference.")

    # Transform the tree to a cladogram if no_branch_length is True
    if no_branch_length:
        tree = transform_to_cladogram(tree)

    # Generate the patristic matrix from the cladogram tree
    patristic_matrix_file = f"{out_prefix}_iter0.distmat"
    patristic_matrix = generate_patristic_matrix(
        tree, sequence_names, patristic_matrix_file
    )

    # Find the closest sequences based on the patristic matrix
    closest_indices = []
    for row in patristic_matrix:
        indices = sorted(range(len(row)), key=lambda i: row[i])
        closest_indices.append(indices)

    # Polarize the alignment using the closest sequences
    polarized_sequences, reverse_polarized_sequences = polarize_alignment(
        alignment, closest_indices
    )

    # Write the polarized alignment to a file in relaxed Phylip format
    polarized_alignment_file = f"{out_prefix}_iter0.polarized"
    write_polarized_alignment(
        alignment,
        polarized_sequences,
        reverse_polarized_sequences,
        polarized_alignment_file,
        reverse=False,
    )

    # Generate a tree on polarized alignment of iteration 0
    tree_file = f"{out_prefix}_iter0.tre"
    estimate_tree(
        polarized_alignment_file,
        tree_file,
        method,
        threads,
        paup_path,
        raxml_path,
        caster_path,
    )
    tree = Phylo.read(tree_file, "newick")

    for i in range(1, iterations):
        print(f"\nStarting iteration {i} of polarization and tree inference.")

        # Prune the tree from previous iteration
        pruned_tree = import_and_prune_allele_tree(tree_file, keep_allele1=True)

        # Transform the tree to a cladogram if no_branch_length is True
        if no_branch_length:
            pruned_tree = transform_to_cladogram(pruned_tree)

        # Generate the patristic matrix from the pruned tree
        patristic_matrix_file = f"{out_prefix}_iter{i}.distmat"
        patristic_matrix = generate_patristic_matrix(
            pruned_tree, sequence_names, patristic_matrix_file
        )

        # Find the closest sequences based on the patristic matrix
        closest_indices = []
        for row in patristic_matrix:
            indices = sorted(range(len(row)), key=lambda i: row[i])
            closest_indices.append(indices)

        # Polarize the alignment using the closest sequences
        if i == iterations - 1:
            polarized_sequences, reverse_polarized_sequences = polarize_alignment(
                alignment, closest_indices
            )
        else:
            polarized_sequences, reverse_polarized_sequences = polarize_alignment(
                alignment, closest_indices
            )

        # Write the new polarized alignment to a file in relaxed Phylip format
        polarized_alignment_file = f"{out_prefix}_iter{i}.polarized"
        if i == iterations - 1:
            write_polarized_alignment(
                alignment,
                polarized_sequences,
                reverse_polarized_sequences,
                polarized_alignment_file,
                reverse=True,
            )
        else:
            write_polarized_alignment(
                alignment,
                polarized_sequences,
                reverse_polarized_sequences,
                polarized_alignment_file,
                reverse=False,
            )

        # Generate the new tree on polarized alignment of the current iteration
        tree_file = f"{out_prefix}_iter{i}.tre"
        estimate_tree(
            polarized_alignment_file,
            tree_file,
            method,
            threads,
            paup_path,
            raxml_path,
            caster_path,
        )
        tree = Phylo.read(tree_file, "newick")

    print(f"\nIterative polarization finished.")

    # Merge sister alleles
    print(f"Now merging sister alleles and alleles with 0 distance.")
    final_polarized_alignment_file = f"{out_prefix}_final.polarized"
    merged_alleles_tree = merge_sister_alleles(
        tree,
        polarized_alignment_file,
        raw_alignment_file,
        final_polarized_alignment_file,
    )

    # Create a final tree
    print(f"Estimating the final tree with phased alleles.")
    final_tree_file = f"{out_prefix}_final.tre"
    estimate_tree(
        final_polarized_alignment_file,
        final_tree_file,
        method,
        threads,
        paup_path,
        raxml_path,
        caster_path,
    )

    print(f"\nIterative polarization completed.")
    print(f"Final polarized alignment: {final_polarized_alignment_file}")
    print(f"Final tree: {final_tree_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Phase sequences by iterative polarization."
    )
    parser.add_argument(
        "--infile", type=str, required=True, help="Input Phylip alignment file."
    )
    parser.add_argument(
        "--out_prefix",
        type=str,
        required=True,
        help="Prefix for the output files (.tre for tree, .txt for matrix, and .polarized for polarized alignment).",
    )
    parser.add_argument(
        "--iterations",
        type=int,
        default=3,
        help="Number of iterations of polarization and tree inference. Default 3.",
    )
    parser.add_argument(
        "--method",
        type=str,
        default="nj",
        help="Method for tree inference. Either 'nj' for Neighbor-Joining or 'paup' for parsimony analysis using PAUP*.",
    )
    parser.add_argument(
        "--no_branch_length",
        action="store_true",
        default=False,
        help="Discard branch length information from guide tree, and instead assume all branches are of equal length.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of threads to use for RAxML analyses only. Default 1.",
    )
    parser.add_argument(
        "--paup_path",
        type=str,
        default=None,
        help="Path to PAUP* command line executable. Use forward slashes.",
    )
    parser.add_argument(
        "--raxml_path",
        type=str,
        default=None,
        help="Path to RAxML executable. Use forward slashes.",
    )
    parser.add_argument(
        "--caster_path",
        type=str,
        default=None,
        help="Path to CASTER executable. Use forward slashes.",
    )

    args = parser.parse_args()

    iterative_polarization(
        args.infile,
        args.out_prefix,
        args.iterations,
        args.method,
        args.no_branch_length,
        args.threads,
        args.paup_path,
        args.raxml_path,
        args.caster_path,
    )
