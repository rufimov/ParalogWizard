def aln_similarity(blat_hit_string: str) -> float:
    """

    :param blat_hit_string:
    :type blat_hit_string:
    :return:
    :rtype:
    """
    return (
        float(blat_hit_string.split()[0])
        / (float(blat_hit_string.split()[0]) + float(blat_hit_string.split()[1]))
    ) * 100


def number_locus_for_sort(blat_hit_string: str) -> str:
    """

    :rtype: object
    """
    return blat_hit_string.split()[13].split("_")[1]


def number_exon_for_sort(blat_hit_string: str) -> int:
    """

    :param blat_hit_string:
    :type blat_hit_string:
    :return:
    :rtype:
    """
    return int(blat_hit_string.split()[13].split("_")[3])


def number_contig_for_sort(blat_hit_string: str) -> int:
    """

    :param blat_hit_string:
    :type blat_hit_string:
    :return:
    :rtype:
    """
    return int(blat_hit_string.split()[9].split("_")[0][6:])


def run_blat(path_to_data, probes, minident, logger):
    """"""

    import glob
    import os
    import shutil

    logger.info("Generating pslx files using BLAT...\n")
    os.makedirs(f"{path_to_data}/exons/50pslx", exist_ok=True)
    for contigfile in glob.glob(f"{path_to_data}/exons/40exonic_contigs/*.fas"):
        file = os.path.basename(contigfile)
        if file != probes:
            logger.info(f"Processing {file}...")
            blat_cmd_output = os.popen(
                f"blat -t=DNA -q=DNA -out=pslx -minIdentity={minident} {probes} {contigfile} {contigfile}.pslx"
            ).read()
            logger.info(blat_cmd_output)
            shutil.move(
                f"{path_to_data}/exons/40exonic_contigs/{file}.pslx",
                f"{path_to_data}/exons/50pslx/{file}.pslx",
            )
            logger.info("Done")


def correct(path_to_data, probes, redlist, logger):
    """"""

    import glob
    import os

    os.makedirs(f"{path_to_data}/exons/50pslx/corrected", exist_ok=True)
    for file in glob.glob(f"{path_to_data}/exons/50pslx/*.pslx"):
        if os.path.basename(file) != f"{probes}.pslx":
            with open(file) as pslx_file, open(
                f"{path_to_data}/exons/50pslx/corrected/{os.path.basename(file)}",
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
