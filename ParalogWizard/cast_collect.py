def collect_contigs(path_to_data, logger):
    """"""

    import os

    os.makedirs(os.path.join(path_to_data, "30raw_contigs"), exist_ok=True)
    path_to_assemblies = os.path.join(path_to_data, "20assemblies")
    for folder in os.listdir(path_to_assemblies):
        if os.path.isdir(os.path.join(path_to_assemblies, folder)):
            logger.info(f"Processing {folder}")
            with open(
                os.path.join(path_to_data, "30raw_contigs", f"{folder}_contigs.fasta"),
                "w",
            ) as contigs:
                for locus in os.listdir(os.path.join(path_to_assemblies, folder)):
                    if os.path.isdir(
                        os.path.join(path_to_assemblies, folder, "locus")
                    ) and os.path.exists(
                        os.path.join(
                            path_to_assemblies,
                            folder,
                            "locus",
                            f"{locus}_contigs.fasta",
                        )
                    ):
                        logger.info(f"\tProcessing {locus}")
                        with open(
                            os.path.join(
                                path_to_assemblies,
                                folder,
                                locus,
                                f"{locus}_contigs.fasta",
                            )
                        ) as locus_contigs:
                            file_content = locus_contigs.read().replace(
                                ">", f">{locus}_"
                            )
                            contigs.write(file_content)
                        logger.info("\tOK")
            if (
                os.stat(
                    os.path.join(
                        path_to_data, "30raw_contigs", f"{folder}_contigs.fasta"
                    )
                ).st_size
                == 0
            ):
                os.remove(
                    os.path.join(
                        path_to_data, "30raw_contigs", f"{folder}_contigs.fasta"
                    )
                )
            logger.info("Ok")
