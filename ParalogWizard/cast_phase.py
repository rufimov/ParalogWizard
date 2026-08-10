# -------------------------------------------------------------------------------
# Filename: cast_phase.py
#
# Version: 0.1
#
# Copyright 2024 Roman Ufimov
#
# This module is part of the ParalogWizard pipeline for paralog detection. It processes
# VCF variant files to generate phased consensus sequences for each sample and contig.
# The module performs the following steps:
#
#   - Copies and concatenates reference exon FASTA files.
#   - Compresses and indexes the VCF file.
#   - Extracts sample names and contig identifiers from the VCF.
#   - Splits the VCF into individual files per sample and contig.
#   - Generates consensus FASTA sequences using bcftools consensus.
#   - Concatenates, trims, and converts the alignment using AMAS.py.
#   - Runs the IterPol phasing pipeline.
#
# Detailed logging is provided for each step to facilitate debugging and traceability.
#
# Usage: python3 ParalogWizard.py cast_phase <data_folder> <num_cores>
#
# -------------------------------------------------------------------------------

import os
import subprocess
import sys


def phase(data_folder, num_cores, logger):
    """
    Process VCF variant files to generate phased consensus sequences for each sample
    and contig, then concatenate, trim, convert the alignment, and run an iterative
    phasing pipeline.

    Parameters:
        data_folder (str): Path to the folder containing input and output data.
        num_cores (int): Number of CPU cores to use.
        logger (logging.Logger): Logger object for recording log messages.
    """
    logger.info("Starting phasing pipeline.")

    try:
        # ----------------------------------------------------------------------
        # 1. Copy and combine reference exon files.
        # ----------------------------------------------------------------------
        logger.info("Copying reference exon FASTA files...")
        ref_source = os.path.join(
            data_folder, "100remapped", "exons*", "reference_exons*.fas"
        )
        ref_dest = os.path.join(data_folder, "101phased")
        cp_command = f"cp {ref_source} {ref_dest}"
        logger.debug(f"Running command: {cp_command}")
        subprocess.run(cp_command, shell=True, check=True)
        logger.info("Reference exon FASTA files copied successfully.")

        logger.info("Concatenating reference FASTA files into 'reference.fas'...")
        concat_ref_command = (
            f"cat {os.path.join(data_folder, '101phased', 'reference_exons*.fas')} > "
            f"{os.path.join(data_folder, '101phased', 'reference.fas')}"
        )
        logger.debug(f"Running command: {concat_ref_command}")
        subprocess.run(concat_ref_command, shell=True, check=True)
        logger.info("Reference FASTA files concatenated successfully.")

        # ----------------------------------------------------------------------
        # 2. Stage compressed/indexed VCF under 101phased.
        # ----------------------------------------------------------------------
        vcf_source_gz = os.path.join(
            data_folder, "100remapped", "all_variants.vcf.gz"
        )
        vcf_source = os.path.join(data_folder, "100remapped", "all_variants.vcf")
        vcf_compressed = os.path.join(data_folder, "101phased", "all_variants.vcf.gz")

        if os.path.isfile(vcf_source_gz):
            logger.info("Copying compressed VCF from cast_call output...")
            subprocess.run(
                ["cp", vcf_source_gz, vcf_compressed], check=True
            )
            # Refresh index next to the copy (tbi may be absent or stale).
            logger.info("Indexing compressed VCF file...")
            subprocess.run(
                ["bcftools", "index", "-f", "-t", vcf_compressed], check=True
            )
        elif os.path.isfile(vcf_source):
            logger.info("Compressing legacy uncompressed VCF file...")
            compress_command = f"bgzip -c {vcf_source} > {vcf_compressed}"
            logger.debug(f"Running command: {compress_command}")
            subprocess.run(compress_command, shell=True, check=True)
            logger.info("VCF file compressed successfully.")
            logger.info("Indexing compressed VCF file...")
            index_command = f"tabix -p vcf {vcf_compressed}"
            logger.debug(f"Running command: {index_command}")
            subprocess.run(index_command, shell=True, check=True)
            logger.info("VCF file indexed successfully.")
        else:
            raise FileNotFoundError(
                "Neither all_variants.vcf.gz nor all_variants.vcf found in "
                f"{os.path.join(data_folder, '100remapped')}"
            )

        base_dir = os.path.join(data_folder, "101phased")
        vcf_file = os.path.join(base_dir, "all_variants.vcf.gz")

        # ----------------------------------------------------------------------
        # 3. Extract sample names and contig identifiers from the VCF.
        # ----------------------------------------------------------------------
        logger.info("Extracting sample names and contig identifiers from VCF.")

        def get_samples(vcf):
            logger.info("Querying sample names using bcftools.")
            result_get_samples = subprocess.run(
                ["bcftools", "query", "-l", vcf],
                stdout=subprocess.PIPE,
                text=True,
                check=True,
            )
            samples_list = result_get_samples.stdout.strip().splitlines()
            logger.debug(f"Samples found: {samples_list}")
            return samples_list

        def get_contigs(vcf):
            logger.info("Querying contig identifiers using bcftools.")
            result_get_contigs = subprocess.run(
                ["bcftools", "query", "-f", "%CHROM\n", vcf],
                stdout=subprocess.PIPE,
                text=True,
                check=True,
            )
            contigs_list = sorted(set(result_get_contigs.stdout.strip().splitlines()))
            logger.debug(f"Contigs found: {contigs_list}")
            return contigs_list

        samples = get_samples(vcf_file)
        contigs = get_contigs(vcf_file)
        logger.info(f"Extracted {len(samples)} samples and {len(contigs)} contigs.")

        # ----------------------------------------------------------------------
        # 4. Create individual VCF files for each sample and contig.
        # ----------------------------------------------------------------------
        logger.info("Creating individual VCF files for each sample and contig.")
        for sample in samples:
            for contig in contigs:
                output_vcf = os.path.join(base_dir, f"{sample}.{contig}.vcf.gz")
                logger.info(f"Creating VCF for sample '{sample}', contig '{contig}'.")
                view_command = [
                    "bcftools",
                    "view",
                    "-s",
                    sample,
                    "-r",
                    contig,
                    "-Oz",
                    "-o",
                    output_vcf,
                    vcf_file,
                ]
                logger.debug(f"Running command: {' '.join(view_command)}")
                subprocess.run(view_command, check=True)
                logger.info(f"VCF file {output_vcf} created successfully.")

                logger.info(f"Indexing VCF file {output_vcf}.")
                subprocess.run(["tabix", "-p", "vcf", output_vcf], check=True)
                logger.info(f"Indexed VCF file {output_vcf} successfully.")

        # ----------------------------------------------------------------------
        # 5. Generate consensus FASTA files and combine them per contig.
        # ----------------------------------------------------------------------
        logger.info("Generating consensus FASTA files for each contig.")
        for contig in contigs:
            combined_fasta = os.path.join(base_dir, f"combined_{contig}.fasta")
            fasta_files = []
            logger.info(f"Processing contig: {contig}")

            for sample in samples:
                in_vcf = os.path.join(base_dir, f"{sample}.{contig}.vcf.gz")
                reference_fasta = os.path.join(base_dir, f"reference_{contig}.fas")
                sample_fasta = os.path.join(base_dir, f"{sample}.{contig}.fa")
                logger.info(
                    f"Generating consensus sequence for sample '{sample}', contig '{contig}'."
                )

                # Capture stdout and stderr from bcftools consensus.
                consensus_cmd = [
                    "bcftools",
                    "consensus",
                    "-I",  # Include indels.
                    "-s",
                    sample,  # Use the specified sample.
                    "-f",
                    reference_fasta,
                    in_vcf,
                ]
                logger.debug("Running command: " + " ".join(consensus_cmd))
                result = subprocess.run(
                    consensus_cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=True,
                )
                # Write the consensus sequence to the designated FASTA file.
                with open(sample_fasta, "w") as out_f:
                    out_f.write(result.stdout)
                logger.info(f"Consensus FASTA generated at {sample_fasta}.")

                # Log any stderr output from bcftools consensus.
                if result.stderr:
                    logger.info(
                        f"bcftools consensus stderr for sample '{sample}', contig '{contig}': {result.stderr}"
                    )

                # Update FASTA header to include only the sample name.
                with open(sample_fasta, "r") as f:
                    lines = f.readlines()
                if lines:
                    lines[0] = f">{sample}\n"
                with open(sample_fasta, "w") as f:
                    f.writelines(lines)
                logger.info(f"Updated FASTA header for sample '{sample}'.")
                fasta_files.append(sample_fasta)

            # Concatenate individual FASTA files into one combined file.
            logger.info(
                f"Combining FASTA files for contig '{contig}' into {combined_fasta}."
            )
            with open(combined_fasta, "w") as outfile:
                for fasta in fasta_files:
                    with open(fasta, "r") as infile:
                        outfile.write(infile.read())
            logger.info(
                f"Combined FASTA file for contig '{contig}' created successfully."
            )

        # ----------------------------------------------------------------------
        # 6. Concatenate, trim, and convert the alignment using AMAS.py.
        # ----------------------------------------------------------------------
        python_executable = sys.executable
        combined_fasta_path = os.path.join(data_folder, "101phased", "combined.fasta")
        partitions_path = os.path.join(
            data_folder, "101phased", "combined_partitions.txt"
        )

        # --- Concatenation using AMAS.py ---
        logger.info("Concatenating exon FASTA files using AMAS.py.")
        concat_cmd = (
            f"{python_executable} {os.path.join('ParalogWizard', 'AMAS.py')} concat "
            f"-i {os.path.join(data_folder, '101phased', 'combined_exons*.fasta')} "
            f"-f fasta -d dna -t {combined_fasta_path} "
            f"-p {partitions_path}"
        )
        logger.debug(f"Running command: {concat_cmd}")
        result_concat = subprocess.run(
            concat_cmd, shell=True, check=True, capture_output=True, text=True
        )
        logger.info("Concatenation with AMAS.py completed successfully.")
        logger.info(f"AMAS concat output: {result_concat.stdout}")
        if result_concat.stderr:
            logger.error(f"AMAS concat error: {result_concat.stderr}")

        # # Replace ambiguous bases 'N' with gaps '-'
        # logger.info(
        #     "Replacing ambiguous bases ('N') with gaps ('-') in the combined FASTA."
        # )
        # sed_cmd = f"sed -i 's/N/-/g' {combined_fasta_path}"
        # logger.debug(f"Running command: {sed_cmd}")
        # subprocess.run(sed_cmd, shell=True, check=True)
        # logger.info("Ambiguous bases replaced successfully.")

        # --- Trimming using AMAS.py ---
        logger.info("Trimming the combined alignment using AMAS.py.")
        trim_output = os.path.join(data_folder, "101phased", "combined_trimmed.fasta")
        trim_cmd = (
            f"{python_executable} {os.path.join('ParalogWizard', 'AMAS.py')} trim -t 1 "
            f"-i {combined_fasta_path} -f fasta -d dna -o {trim_output}"
        )
        logger.debug(f"Running command: {trim_cmd}")
        result_trim = subprocess.run(
            trim_cmd, shell=True, check=True, capture_output=True, text=True
        )
        logger.info("Trimming completed successfully.")
        logger.info(f"AMAS trim output: {result_trim.stdout}")
        if result_trim.stderr:
            logger.error(f"AMAS trim error: {result_trim.stderr}")

        # --- Conversion using AMAS.py ---
        logger.info("Converting trimmed alignment to PHYLIP format using AMAS.py.")
        convert_cmd = (
            f"{python_executable} {os.path.join('ParalogWizard', 'AMAS.py')} convert -u phylip "
            f"-i {trim_output} -f fasta -d dna"
        )
        logger.debug(f"Running command: {convert_cmd}")
        result_convert = subprocess.run(
            convert_cmd, shell=True, check=True, capture_output=True, text=True
        )
        logger.info("Conversion to PHYLIP format completed successfully.")
        logger.info(f"AMAS conversion output: {result_convert.stdout}")
        if result_convert.stderr:
            logger.error(f"AMAS conversion error output: {result_convert.stderr}")

        phy_file = os.path.join(data_folder, "101phased", "combined_trimmed.phy")
        mv_cmd = f"mv combined_trimmed.fasta-out.phy {phy_file}"
        logger.debug(f"Running command: {mv_cmd}")
        subprocess.run(mv_cmd, shell=True, check=True)
        logger.info(f"Renamed PHYLIP file to {phy_file}.")

        # ----------------------------------------------------------------------
        # 7. Run the IterPol phasing pipeline.
        # ----------------------------------------------------------------------
        logger.info("Running IterPol phasing pipeline.")
        iterpol_cmd = (
            f"{python_executable} {os.path.join('ParalogWizard', 'IterPol_v0.4.py')} "
            f"--infile {phy_file} --out_prefix phased "
            f"--threads {num_cores} --raxml_path raxmlHPC-PTHREADS-SSE3 --method raxml"
        )
        logger.debug(f"Running command: {iterpol_cmd}")
        output_iterpol = subprocess.run(
            iterpol_cmd, shell=True, check=True, capture_output=True, text=True
        )
        logger.info("IterPol pipeline completed successfully.")
        logger.info(f"IterPol output: {output_iterpol.stdout}")
        if output_iterpol.stderr:
            logger.error(f"IterPol error output: {output_iterpol.stderr}")

        logger.info("Moving IterPol output files into the phased folder.")
        mv_iterpol_cmd = f"mv phased* {os.path.join(data_folder, '101phased')}"
        logger.debug(f"Running command: {mv_iterpol_cmd}")
        subprocess.run(mv_iterpol_cmd, shell=True, check=True)
        logger.info("IterPol output files moved successfully.")

    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed: {e.cmd}")
        logger.error(f"Return code: {e.returncode}")
        logger.error(f"Output: {e.output}")
        raise
    except Exception as e:
        logger.exception(
            f"An unexpected error occurred during the phasing pipeline: {e}"
        )
        raise

    logger.info("Phasing pipeline completed successfully.")
