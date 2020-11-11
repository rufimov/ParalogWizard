import sys
import os

path_to_data_HP: str = sys.argv[1].strip()
path_to_data_HPM: str = sys.argv[2].strip()
os.makedirs(f"{path_to_data_HPM}/HybPiper_contigs", exist_ok=True)
for folder in os.listdir(path_to_data_HP):
    if os.path.isdir(f"{path_to_data_HP}/{folder}"):
        print(f"Processing {folder}")
        with open(
            f"{path_to_data_HPM}/HybPiper_contigs/{folder}_contigs.fasta", "w"
        ) as contigs:
            for locus in os.listdir(f"{path_to_data_HP}/{folder}"):
                if os.path.isdir(
                    f"{path_to_data_HP}/{folder}/{locus}"
                ) and os.path.exists(
                    f"{path_to_data_HP}/{folder}/{locus}/{locus}_contigs.fasta"
                ):
                    print(f"\tProcessing {locus}")
                    with open(
                        f"{path_to_data_HP}/{folder}/{locus}/{locus}_contigs.fasta"
                    ) as locus_contigs:
                        file_content = locus_contigs.read().replace(">", f">{locus}_")
                        contigs.write(file_content)
                    print("\tOK")
        if (
            os.stat(
                f"{path_to_data_HPM}/HybPiper_contigs/{folder}_contigs.fasta"
            ).st_size
            == 0
        ):
            os.remove(f"{path_to_data_HPM}/HybPiper_contigs/{folder}_contigs.fasta")
        print("Ok")
