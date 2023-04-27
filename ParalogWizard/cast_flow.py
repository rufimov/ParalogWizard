import fileinput
import itertools
import multiprocessing
import os.path
import re
from glob import glob
from math import ceil

import pandas
from Bio import SeqRecord
from Bio import SeqIO
from Bio.Align.Applications import MafftCommandline
from Bio.Sequencing.Applications import BwaIndexCommandline
from Bio.Sequencing.Applications import BwaMemCommandline
from Bio.Sequencing.Applications import SamtoolsIndexCommandline
from Bio.Sequencing.Applications import SamtoolsVersion1xSortCommandline
from Bio.Sequencing.Applications import SamtoolsViewCommandline


def trim_ends(seq1, seq2):
    while seq1[0] == "-" or seq2[0] == "-":
        seq1 = seq1[1:]
        seq2 = seq2[1:]
    while seq1[-1] == "-" or seq2[-1] == "-":
        seq1 = seq1[:-1]
        seq2 = seq2[:-1]
    return seq1, seq2


def mafft_align(file):
    stdout, stderr = MafftCommandline(
        input=file,
        auto=True,
    )()
    with open(f"{os.path.splitext(file)[0]}.mafft.fasta", "w") as aligned:
        aligned.write(stdout)


def bwa_index(reference_to_map_to):
    BwaIndexCommandline(infile=reference_to_map_to, a="is")()


def bwa_map(exon, sample, main_data_folder):
    output_path = os.path.join(main_data_folder, "ABBA_BABA", exon)
    reference_to_map_to = os.path.join(
        main_data_folder, "ABBA_BABA", exon, f"reference_{exon}.fas"
    )
    stdout, stderr = BwaMemCommandline(
        reference=reference_to_map_to,
        read_file1=os.path.join(
            main_data_folder, "10deduplicated_reads", sample + ".R1.fastq.gz"
        ),
        read_file2=os.path.join(
            main_data_folder, "10deduplicated_reads", sample + ".R2.fastq.gz"
        ),
    )()
    with open(os.path.join(output_path, f"{sample}.sam"), "w") as mapped_sample:
        mapped_sample.write(stdout)
    SamtoolsViewCommandline(
        input_file=os.path.join(output_path, f"{sample}.sam"),
        h=True,
        b=True,
        o=os.path.join(output_path, f"{sample}.bam"),
    )()
    SamtoolsViewCommandline(
        input_file=os.path.join(output_path, f"{sample}.bam"),
        h=True,
        b=True,
        q=3,
        o=os.path.join(output_path, f"{sample}_filtered.bam"),
    )()
    SamtoolsViewCommandline(
        input_file=os.path.join(output_path, f"{sample}_filtered.bam"),
        h=True,
        b=True,
        # f=3,
        F=0x90C,
        o=os.path.join(output_path, f"{sample}_filtered_uniq.bam"),
    )()
    SamtoolsVersion1xSortCommandline(
        input=os.path.join(output_path, f"{sample}_filtered_uniq.bam"),
        o=os.path.join(output_path, f"{sample}_filtered_uniq_sorted.bam"),
    )()
    SamtoolsIndexCommandline(
        input_bam=os.path.join(output_path, f"{sample}_filtered_uniq_sorted.bam")
    )()


def variant_call(exon, main_data_folder, outgroup):
    bam_files = " ".join(
        glob(
            os.path.join(
                main_data_folder, "ABBA_BABA", exon, "*_filtered_uniq_sorted.bam"
            )
        )
    )
    variants_vcf_file = os.path.join(
        main_data_folder, "ABBA_BABA", exon, f"{exon}_variants.vcf"
    )
    reference_file = os.path.join(
        main_data_folder, "ABBA_BABA", exon, f"reference_{exon}.fas"
    )
    os.popen(
        f"bcftools mpileup -Ou -I -d 2000 -A -f {reference_file} {bam_files} | "
        f"bcftools call -Ov -mv > {variants_vcf_file}"
    ).read()
    variants_vcf_file_filtered = os.path.join(
        main_data_folder, "ABBA_BABA", exon, exon + "_variants_filtered.vcf"
    )
    os.popen(
        f" bcftools filter -e 'QUAL<20' {variants_vcf_file} > {variants_vcf_file_filtered}"
    ).read()
    with open(variants_vcf_file) as vcf_file:
        vcf_file_lines = vcf_file.readlines()
        num_lines = len(vcf_file_lines)
        for i in range(num_lines):
            line = vcf_file_lines[i]
            if line.startswith("#CHROM"):
                corrected_line = re.sub(
                    r"%s.+?/%s/" % (main_data_folder, exon), r"", line
                )
                corrected_line = re.sub(
                    r"_filtered_uniq_sorted.bam", r"", corrected_line
                )
                corrected_line = re.sub(f"{outgroup}", r"Outgroup", corrected_line)
                vcf_file_lines[i] = corrected_line
                break
        with open(
            f"{os.path.join(main_data_folder, 'ABBA_BABA', exon, exon + '_variants_filtered_corrected.vcf')}",
            "w",
        ) as vcf_file_corrected:
            for line in vcf_file_lines:
                vcf_file_corrected.write(line)


def split_to_chunks(list_a, chunk_size):
    for i in range(0, len(list_a), chunk_size):
        yield list_a[i : i + chunk_size]


def split_reference(reference_dict, data_folder, exons_allowed, chunk_size):
    exons_chunks = split_to_chunks(exons_allowed, chunk_size)
    count = 1
    return_chunks = []
    for chunk in exons_chunks:
        to_write = ""
        chunk_count = 1
        for item in chunk:
            exon = item[0]
            if item[1]:
                exon_para = (
                    f"{exon.split('_')[0] + 'para'}_{'_'.join(exon.split('_')[1:])}"
                )
                if (exon in reference_dict.keys()) and (
                    exon_para in reference_dict.keys()
                ):
                    seq = reference_dict[exon].seq
                    seq_para = reference_dict[exon_para].seq
                    with open(f"{exon}.fasta", "w") as file:
                        SeqIO.write(
                            SeqRecord.SeqRecord(seq, exon, exon, exon),
                            file,
                            "fasta-2line",
                        )
                        SeqIO.write(
                            SeqRecord.SeqRecord(
                                seq_para, exon_para, exon_para, exon_para
                            ),
                            file,
                            "fasta-2line",
                        )
                    mafft_align(f"{exon}.fasta")
                    fasta = list(SeqIO.parse(f"{exon}.mafft.fasta", "fasta"))
                    seq = fasta[0].seq
                    seq_para = fasta[1].seq
                    seq_trimmed, seq_para_trimmed = trim_ends(seq, seq_para)
                    if chunk_count == 1:
                        to_write = (
                            to_write
                            + seq_trimmed.replace("-", "")
                            + "N" * 200
                            + seq_para_trimmed.replace("-", "")
                        )
                    else:
                        to_write = (
                            to_write
                            + "N" * 200
                            + seq_trimmed.replace("-", "")
                            + "N" * 200
                            + seq_para_trimmed.replace("-", "")
                        )
            else:
                seq = reference_dict[exon].seq
                if chunk_count == 1:
                    to_write = to_write + seq
                else:
                    to_write = to_write + "N" * 200 + seq
            chunk_count += 1
        seq_record_to_write = SeqRecord.SeqRecord(
            to_write, f"exons{str(count)}", f"exons{str(count)}", f"exons{str(count)}"
        )
        os.makedirs(
            os.path.join(data_folder, "ABBA_BABA", f"exons{str(count)}"), exist_ok=True
        )
        with open(
            os.path.join(
                data_folder,
                "ABBA_BABA",
                f"exons{str(count)}",
                f"reference_exons{str(count)}.fas",
            ),
            "w",
        ) as exon_ref_file:
            SeqIO.write(seq_record_to_write, exon_ref_file, "fasta")
        return_chunks.append(f"exons{str(count)}")
        count += 1
    return return_chunks


def flow(reference, data_folder, num_cores, outgroup, read_length):
    with open(
        os.path.join(data_folder, "10deduplicated_reads", "samples_list.txt")
    ) as samples:
        sample_list = samples.read().splitlines()
    reference_dict = SeqIO.to_dict(SeqIO.parse(reference, "fasta"))
    reference_dict_corrected = dict()
    for exon in reference_dict.keys():
        exon_short = "_".join(exon.split("_")[1:4]).replace("Contig", "exon")
        sequence = reference_dict[exon].seq.ungap("-")
        seq_record = reference_dict[exon]
        seq_record.seq = sequence
        exon_length = len(sequence)
        if exon_length >= read_length:
            reference_dict_corrected[exon_short] = seq_record

    exons_in_ref = reference_dict_corrected.keys()
    all_paralogs_for_reference = pandas.read_csv(
        os.path.join(data_folder, "41detected_par", "all_paralogs_for_reference.tsv"),
        sep="\t",
    )
    from ParalogWizard.cast_detect import exon_stats

    all_paralogs_for_reference_exon_stat = exon_stats(all_paralogs_for_reference)
    all_paralogs_for_reference_exon_stat.to_csv(
        os.path.join(
            data_folder,
            "41detected_par",
            "exon_statistics.tsv",
        ),
        na_rep="NaN",
        sep="\t",
    )
    allowed_exons = []
    for column in all_paralogs_for_reference_exon_stat:
        num_yes = 0
        num_no = 0
        try:
            num_yes = all_paralogs_for_reference_exon_stat[column].value_counts()["Yes"]
        except KeyError:
            try:
                num_no = all_paralogs_for_reference_exon_stat[column].value_counts()[
                    "No"
                ]
            except KeyError:
                pass
        perc_yes = num_yes / len(sample_list)
        perc_no = num_no / len(sample_list)
        if perc_no >= 0.9 and (column in exons_in_ref):
            allowed_exons.append((column, False))
        elif (
            perc_yes >= 0.5
            and (column in exons_in_ref)
            and (
                f"{column.split('_')[0] + 'para'}_{'_'.join(column.split('_')[1:])}"
                in exons_in_ref
            )
        ):
            allowed_exons.append((column, True))
    chunk_size = ceil(len(allowed_exons) / 10)
    exons = split_reference(
        reference_dict_corrected, data_folder, allowed_exons, chunk_size
    )
    # with multiprocessing.Pool(processes=num_cores) as pool_index:
    #     pool_index.map(
    #         bwa_index,
    #         [
    #             os.path.join(data_folder, "ABBA_BABA", x, f"reference_{x}.fas")
    #             for x in exons
    #         ],
    #     )
    #
    # args_map = []
    # for pair in list(itertools.product(exons, sample_list)):
    #     exon_map = pair[0]
    #     sample_map = pair[1]
    #     args_map.append((exon_map, sample_map, data_folder))
    # with multiprocessing.Pool(processes=num_cores) as pool_bwa:
    #     pool_bwa.starmap(bwa_map, args_map, chunksize=100)
    # args_call = list(zip(exons, [data_folder] * len(exons), [outgroup] * len(exons)))
    # with multiprocessing.Pool(processes=num_cores) as pool_call:
    #     pool_call.starmap(variant_call, args_call, chunksize=100)
    vcf_to_concat = [
        os.path.join(
            data_folder, "ABBA_BABA", x, x + "_variants_filtered_corrected.vcf"
        )
        for x in exons
    ]

    with open(
        os.path.join(data_folder, "10deduplicated_reads", "samples_list.txt")
    ) as sample_file, open(
        os.path.join(data_folder, "ABBA_BABA", "samples_list.txt"), "w"
    ) as corrected_sample_file:
        for line in sample_file.readlines():
            corrected_line = re.sub(outgroup, "Outgroup", line)
            corrected_sample_file.write(corrected_line)

    for file in vcf_to_concat:
        folder = os.path.dirname(file)
        os.popen(
            f"bcftools view -S {os.path.join(data_folder, 'ABBA_BABA', 'samples_list.txt')} "
            f"{file} > {os.path.join(folder + '_variants_filtered_corrected.vcf')}"
        ).read()
    vcf_to_concat = " ".join(
        [
            os.path.join(
                data_folder, "ABBA_BABA", x + "_variants_filtered_corrected.vcf"
            )
            for x in exons
        ]
    )
    os.popen(
        f"bcftools concat -o {os.path.join(data_folder, 'ABBA_BABA', 'all_variants.vcf')} {vcf_to_concat}"
    ).read()
