import sys
import os
import glob
import shutil


def aln_similarity(y: str) -> float:
    return (
        float(y.split("\t")[0]) / (float(y.split("\t")[0]) + float(y.split("\t")[1]))
    ) * 100


def number_locus_for_sort(z: str) -> str:
    return z.split()[13].split("_")[1]


def number_exon_for_sort(a: str) -> int:
    return int(a.split()[13].split("_")[3])


def number_contig_for_sort(b: str) -> int:
    return int(b.split()[9].split("_")[0][6:])


def run_blat():
    print("Generating pslx files using BLAT...\n")
    os.makedirs(f"{path_to_data_HPM}/exons/50pslx", exist_ok=True)
    for contigfile in glob.glob(f"{path_to_data_HPM}/exons/40contigs/*.fas"):
        file = contigfile.split("/")[-1]
        if file != probes:
            print(f"Processing {file}...")
            os.system(
                f"blat -t=DNA -q=DNA -out=pslx -minIdentity={minident} {probes} {contigfile} {contigfile}.pslx"
            )
            shutil.move(
                f"{path_to_data_HPM}/exons/40contigs/{file}.pslx",
                f"{path_to_data_HPM}/exons/50pslx/{file}.pslx",
            )
            print("Done")


def correct():
    os.makedirs(f"{path_to_data_HPM}/exons/50pslx/corrected", exist_ok=True)
    for file in glob.glob(f"{path_to_data_HPM}/exons/50pslx/*.pslx"):
        if file.split("/")[-1] != f"{probes}.pslx":
            with open(file) as pslx_file, open(
                f"{path_to_data_HPM}/exons/50pslx/corrected/{file.split('/')[-1]}", "w"
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
                    if sample not in whitelist:
                        if line.split()[9] not in hits2:
                            list_to_work_cleaned2.append(line)
                            hits2.add(line.split()[9])
                    else:
                        list_to_work_cleaned2.append(line)
                for line in head:
                    corrected_pslx_file.write(line)
                for line in list_to_work_cleaned2:
                    corrected_pslx_file.write(line)


if __name__ == "__main__":
    path_to_data_HPM = sys.argv[1]
    probes = sys.argv[2]
    minident = float(sys.argv[3].strip())
    whitelist = set([x.strip() for x in sys.argv[4].split(",")])
    run_blat()
    correct()
