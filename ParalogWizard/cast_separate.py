import glob
import multiprocessing
import os
import re
import shutil
import subprocess
import fileinput

from ParalogWizard.cast_analyze import mafft_align


def aln_similarity(blat_hit_string: str) -> float:
    """"""
    return (
        float(blat_hit_string.split()[0])
        / (float(blat_hit_string.split()[0]) + float(blat_hit_string.split()[1]))
    ) * 100


def number_locus_for_sort(blat_hit_string: str) -> str:
    """"""
    return blat_hit_string.split()[13].split("_")[1]


def number_exon_for_sort(blat_hit_string: str) -> int:
    """"""
    return int(blat_hit_string.split()[13].split("_")[3])


def number_contig_for_sort(blat_hit_string: str) -> int:
    """"""
    return int(blat_hit_string.split()[9].split("_")[0][6:])


def run_blat(path_to_data, probes, minident, logger):
    """"""

    logger.info("Generating pslx files using BLAT...\n")
    os.makedirs(os.path.join(path_to_data, "50pslx"), exist_ok=True)
    for contigfile in glob.glob(
        os.path.join(path_to_data, "31exonic_contigs", "*.fas")
    ):
        file = os.path.basename(contigfile)
        if file != probes:
            logger.info(f"Processing {file}...")
            blat_cmd_output = os.popen(
                f"blat -t=DNA -q=DNA -out=pslx -minIdentity={minident} {probes} {contigfile} \
                {os.path.join(path_to_data, '50pslx', f'{os.path.basename(contigfile)}.pslx')}"
            ).read()
            logger.info(blat_cmd_output)
            logger.info("Done")


def correct(path_to_data, redlist, logger):
    """"""

    os.makedirs(os.path.join(path_to_data, "50pslx", "corrected"), exist_ok=True)
    for file in glob.glob(os.path.join(path_to_data, "50pslx", "*.pslx")):
        with open(
            os.path.join(path_to_data, "50pslx", "corrected", "list_pslx.txt"), "a"
        ) as list:
            list.write(
                f"{os.path.join(os.path.dirname(file),'corrected',os.path.basename(file))}\n"
            )
        with open(file) as pslx_file, open(
            os.path.join(path_to_data, "50pslx", "corrected", os.path.basename(file)),
            "w",
        ) as corrected_pslx_file:
            file = pslx_file.readlines()
            head = file[0:5]
            list_to_work = file[5:]
            list_to_work.sort(key=aln_similarity, reverse=True)
            list_to_work.sort(key=number_exon_for_sort)
            list_to_work.sort(key=number_locus_for_sort)
            list_to_work_cleaned1 = []
            hits1 = set()
            for line in list_to_work:
                if line.split()[13] not in hits1:
                    list_to_work_cleaned1.append(line)
                    hits1.add(line.split()[13])
            list_to_work_cleaned1.sort(key=aln_similarity, reverse=True)
            list_to_work_cleaned1.sort(key=number_contig_for_sort)
            list_to_work_cleaned2 = []
            hits2 = set()
            for line in list_to_work_cleaned1:
                sample = f"{line.split()[9].split('-')[0].split('_')[1]}-{line.split()[9].split('-')[1]}"
                if sample not in redlist:
                    if line.split()[9] not in hits2:
                        list_to_work_cleaned2.append(line)
                        hits2.add(line.split()[9])
                else:
                    list_to_work_cleaned2.append(line)
            for line in head:
                corrected_pslx_file.write(line)
            for line in list_to_work_cleaned2:
                corrected_pslx_file.write(line)


def align(path_to_data, probes, n_cpu):
    shutil.rmtree(os.path.join(path_to_data, "60mafft"), ignore_errors=True)
    subprocess.call(
        f"python3 ParalogWizard/assembled_exons_to_fastas.py \
-l {os.path.join(path_to_data, '50pslx', 'corrected', 'list_pslx.txt')} -f {probes} \
-d {os.path.join(path_to_data, '60mafft')}",
        shell=True,
    )
    all_loci = set()
    pool_aln = multiprocessing.Pool(processes=n_cpu)
    for file in glob.glob(os.path.join(path_to_data, "60mafft", "*.fasta")):
        with fileinput.FileInput(file, inplace=True) as file_to_correct:
            for line in file_to_correct:
                print(re.sub(r">.+/", ">", line), end="")
        locus = os.path.basename(file).split("_")[3]
        all_loci.add(locus)
        pool_aln.apply_async(mafft_align, (file,))
    pool_aln.close()
    pool_aln.join()
    os.makedirs(
        os.path.join(path_to_data, "70concatenated_exon_alignments"), exist_ok=True
    )
    amas_ex = os.path.abspath("ParalogWizard/AMAS.py")
    os.chdir(os.path.join(path_to_data, "70concatenated_exon_alignments"))
    for locus in all_loci:
        loci_to_concat = glob.glob(
            os.path.join("..", "60mafft", f"*{locus}_*.mafft.fasta")
        )
        loci_to_concat.sort(key=lambda x: int(x.split("_")[5]))
        line_loci_to_concat = " ".join(loci_to_concat)
        subprocess.call(
            f"python3 {amas_ex} concat -i {line_loci_to_concat} -f fasta -d dna -t Assembly_{locus}.fas \
-p Assembly_{locus}.part",
            shell=True,
        )
        subprocess.call(
            f"python3 {amas_ex} convert -i Assembly_{locus}.fas -f fasta -d dna -u phylip",
            shell=True,
        )
        os.rename(f"Assembly_{locus}.fas-out.phy", f"Assembly_{locus}.phy")
