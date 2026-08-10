#!/usr/bin/env python
import multiprocessing
import os
import re
import shutil
import subprocess
from glob import glob

from ParalogWizard import worker_initializer
from ParalogWizard.cast_analyze import setup_logging

logger = setup_logging()


# -----------------------------------------------------------------------------
# Logging Decorator (with detailed context)
# -----------------------------------------------------------------------------
def log_exceptions(func):
    """
    Decorator that logs function entry, exit, and any exceptions.
    Logs the function's name along with its positional and keyword arguments.
    Instead of calling sys.exit(), it re-raises the exception so that the main process
    can catch the error and shut down the pool gracefully.
    """
    from functools import wraps

    @wraps(func)
    def wrapper(*args, **kwargs):
        try:
            logger.debug(
                "Entering %s with args: %s, kwargs: %s", func.__name__, args, kwargs
            )
            result = func(*args, **kwargs)
            logger.debug("Exiting %s", func.__name__)
            return result
        except Exception as e:
            arg_str = ", ".join(str(arg) for arg in args)
            kwarg_str = ", ".join(f"{k}={v}" for k, v in kwargs.items())
            logger.exception(
                f"Exception in {func.__name__} (args: {arg_str}; kwargs: {kwarg_str}): {e}"
            )
            raise

    return wrapper


# -----------------------------------------------------------------------------
# Variant Calling Functions
# -----------------------------------------------------------------------------
@log_exceptions
def variant_call(exon, main_data_folder):
    """
    Calls variants using bcftools, filters and adjusts them, outputting a VCF file.
    :param exon: the exon to call variants for.
    :param main_data_folder: the main data folder.
    :return: None
    """
    # Build file paths.
    bam_files_pattern = os.path.join(
        main_data_folder, "100remapped", exon, "*_filtered_uniq_sorted.bam"
    )
    bam_files = " ".join(glob(bam_files_pattern))
    variants_vcf_file = os.path.join(
        main_data_folder, "100remapped", exon, f"{exon}_variants.vcf"
    )
    reference_file = os.path.join(
        main_data_folder, "100remapped", exon, f"reference_{exon}.fas"
    )

    # Run bcftools mpileup and call.
    mpileup_call_cmd = (
        f"bcftools mpileup -Ou -I -d 2000 -A -f {reference_file} {bam_files} | "
        f"bcftools call -Ov -mv > {variants_vcf_file}"
    )
    logger.info("Running bcftools mpileup and call for exon %s", exon)
    try:
        subprocess.run(
            mpileup_call_cmd, shell=True, check=True, capture_output=True, text=True
        )
    except subprocess.CalledProcessError as e:
        logger.error("Variant calling failed for exon %s: %s", exon, e.stderr)
        raise

    # Filter the raw VCF.
    variants_vcf_file_filtered = os.path.join(
        main_data_folder, "100remapped", exon, f"{exon}_variants_filtered.vcf"
    )
    filter_cmd = f"bcftools filter -e 'QUAL<20' {variants_vcf_file} > {variants_vcf_file_filtered}"
    logger.info("Filtering variants for exon %s", exon)
    try:
        subprocess.run(
            filter_cmd, shell=True, check=True, capture_output=True, text=True
        )
    except subprocess.CalledProcessError as e:
        logger.error("Filtering failed for exon %s: %s", exon, e.stderr)
        raise

    # Adjust the VCF header.
    corrected_vcf_file = os.path.join(
        main_data_folder, "100remapped", exon, f"{exon}_variants_filtered_corrected.vcf"
    )
    try:
        with open(variants_vcf_file_filtered, "r") as vcf_file_filtered:
            vcf_file_filtered_lines = vcf_file_filtered.readlines()
    except Exception as e:
        logger.error("Failed to read VCF file for exon %s: %s", exon, e)
        raise

    for i, line in enumerate(vcf_file_filtered_lines):
        if line.startswith("#CHROM"):
            # Remove file path details from the header.
            pattern = r"%s.+?/%s/" % (re.escape(main_data_folder), re.escape(exon))
            corrected_line = re.sub(pattern, "", line)
            corrected_line = re.sub(r"_filtered_uniq_sorted.bam", "", corrected_line)
            vcf_file_filtered_lines[i] = corrected_line
            break

    try:
        with open(corrected_vcf_file, "w") as vcf_file_corrected:
            vcf_file_corrected.writelines(vcf_file_filtered_lines)
    except Exception as e:
        logger.error("Failed to write corrected VCF for exon %s: %s", exon, e)
        raise
    logger.info("Ordering samples in VCF for exon %s", exon)
    samples_list_dst = os.path.join(main_data_folder, "100remapped", "samples_list.txt")
    folder = os.path.dirname(corrected_vcf_file)
    temp_vcf = os.path.join(folder, "temp_variants.vcf")
    view_order_sample_cmd = (
        f"bcftools view -S {samples_list_dst} {corrected_vcf_file} > "
        f"{os.path.join(folder, 'temp_variants.vcf')}"
    )
    try:
        subprocess.run(
            view_order_sample_cmd,
            shell=True,
            check=True,
            capture_output=True,
            text=True,
        )
        # Replace the original file with the temporary output.
        os.replace(temp_vcf, corrected_vcf_file)
    except subprocess.CalledProcessError as e:
        logger.error(
            "bcftools sample order failed for %s: %s", corrected_vcf_file, e.stderr
        )
        raise

    logger.info("Variant calling and filtering completed for exon %s", exon)


@log_exceptions
def call_variants(data_folder, num_cores, log_queue):
    """
    Calls variants using bcftools for all exons.
    :param data_folder: main data folder.
    :param num_cores: number of cores to use.
    :param log_queue: logging queue.
    :return: None
    """
    logger.info("Calling variants")

    # Copy samples list.
    samples_list_src = os.path.join(
        data_folder, "10deduplicated_reads", "samples_list.txt"
    )
    shutil.copy(samples_list_src, os.path.join(data_folder, "100remapped"))
    exons_dir_pattern = os.path.join(data_folder, "100remapped", "exons*")
    exons = [os.path.basename(x) for x in glob(exons_dir_pattern) if os.path.isdir(x)]
    args_call = list(zip(exons, [data_folder] * len(exons)))

    # Use asynchronous multiprocessing
    pool_call = multiprocessing.Pool(
        processes=num_cores, initializer=worker_initializer, initargs=(log_queue,)
    )
    async_results = []
    try:
        for args in args_call:
            async_results.append(pool_call.apply_async(variant_call, args))
        for args, async_result in zip(args_call, async_results):
            exon, _ = args
            try:
                async_result.get()
            except Exception as e:
                logger.error("Error during variant calling for exon %s: %s", exon, e)
                pool_call.terminate()
                raise
    finally:
        pool_call.close()
        pool_call.join()

    # Concatenate all corrected VCF files into one.
    vcf_files_concat = " ".join(
        [
            os.path.join(
                data_folder,
                "100remapped",
                exon,
                f"{exon}_variants_filtered_corrected.vcf",
            )
            for exon in exons
        ]
    )
    concat_cmd = f"bcftools concat -o {os.path.join(data_folder, '100remapped', 'all_variants.vcf')} {vcf_files_concat}"
    logger.info("Concatenating VCF files")
    try:
        subprocess.run(
            concat_cmd, shell=True, check=True, capture_output=True, text=True
        )
    except subprocess.CalledProcessError as e:
        logger.error("bcftools concat failed: %s", e.stderr)
        raise

    logger.info("Variant calling completed for all exons.")
