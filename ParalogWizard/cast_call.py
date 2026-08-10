#!/usr/bin/env python
"""
------------------------------------------------------------------------------------------------------------------------
Copyright 2024 Roman Ufimov under the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

cast_call — variant calling on remapped BAMs under 100remapped/.

Requires cast_remap output: for every sample in samples_list.txt and every
chunk directory exons*, a {sample}_filtered_uniq_sorted.bam must exist.

  1. Copy samples_list.txt into 100remapped/ and abort if any BAM is missing.
  2. Per chunk (in parallel): bcftools mpileup | call | filter → compressed VCF,
     then rename/reorder sample columns to match samples_list.txt and index.
  3. bcftools concat all chunk VCFs → 100remapped/all_variants.vcf.gz (+ index).

num_cores (-nc) is split into chunk workers × bcftools threads for local or
cluster use. Complete per-chunk outputs are skipped unless --force.
------------------------------------------------------------------------------------------------------------------------
"""

import gzip
import logging
import os
import shutil
import subprocess
from glob import glob
from shutil import which

from ParalogWizard import log_exceptions, managed_pool

logger = logging.getLogger("ParalogWizard")

BAM_SUFFIX = "_filtered_uniq_sorted.bam"
CORRECTED_VCF_NAME = "{exon}_variants_filtered_corrected.vcf.gz"
ALL_VARIANTS_VCF = "all_variants.vcf.gz"
REQUIRED_TOOLS = ("bcftools", "bgzip", "tabix")


def allocate_workers_and_threads(num_cores, n_tasks):
    """
    Split num_cores into (n_workers, threads_per_worker) for a pool of n_tasks.

    n_workers = min(tasks, cores); threads_per_worker = cores // n_workers.
    Used by cast_call and cast_remap so -nc behaves the same locally and on a cluster.
    """
    if n_tasks < 1:
        raise ValueError("n_tasks must be >= 1")
    n_workers = max(1, min(int(n_tasks), int(num_cores)))
    threads_per_worker = max(1, int(num_cores) // n_workers)
    return n_workers, threads_per_worker


def require_file(path, description="file"):
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Missing {description}: {path}")
    size = os.path.getsize(path)
    if size == 0:
        raise ValueError(f"Empty {description}: {path}")
    logger.debug("OK %s: %s (%d bytes)", description, path, size)
    return path


def require_dir(path, description="directory"):
    if not os.path.isdir(path):
        raise FileNotFoundError(f"Missing {description}: {path}")
    logger.debug("OK %s: %s", description, path)
    return path


def require_tools(tools):
    missing = [name for name in tools if which(name) is None]
    if missing:
        raise EnvironmentError(
            "Required executable(s) not found on PATH: " + ", ".join(missing)
        )
    for name in tools:
        logger.debug("Found tool on PATH: %s -> %s", name, which(name))


def _clean_sample_name(raw_name):
    """
    Turn a VCF sample column (path, BAM basename, or SM tag) into a samples_list ID.
    Strips directory, _filtered_uniq_sorted.bam, and a trailing .bam if present.
    """
    name = os.path.basename(raw_name.strip())
    if name.endswith(BAM_SUFFIX):
        name = name[: -len(BAM_SUFFIX)]
    if name.endswith(".bam"):
        name = name[: -len(".bam")]
    if name != raw_name.strip():
        logger.debug("Cleaned sample name: %r -> %r", raw_name.strip(), name)
    return name


def _read_samples_list(samples_list_path):
    require_file(samples_list_path, "samples list")
    with open(samples_list_path) as fh:
        samples = [line.strip() for line in fh if line.strip()]
    if not samples:
        raise ValueError(f"No samples found in '{samples_list_path}'.")
    duplicates = sorted({s for s in samples if samples.count(s) > 1})
    if duplicates:
        raise ValueError(
            f"Duplicate sample name(s) in {samples_list_path}: {', '.join(duplicates)}"
        )
    logger.debug("samples_list (%d): %s", len(samples), ", ".join(samples))
    return samples


def _bam_path(main_data_folder, exon, sample):
    return os.path.join(
        main_data_folder, "100remapped", exon, f"{sample}{BAM_SUFFIX}"
    )


def _corrected_vcf_path(main_data_folder, exon):
    return os.path.join(
        main_data_folder,
        "100remapped",
        exon,
        CORRECTED_VCF_NAME.format(exon=exon),
    )


def _missing_samples(desired_samples, available_samples):
    return [s for s in desired_samples if s not in available_samples]


def _run_shell(cmd, log_file):
    """Run a shell command; write combined stdout/stderr to log_file."""
    logger.debug("Shell command -> %s\n%s", log_file, cmd)
    with open(log_file, "w") as log_fh:
        completed = subprocess.run(
            cmd,
            shell=True,
            check=False,
            stdout=log_fh,
            stderr=subprocess.STDOUT,
        )
    logger.debug(
        "Shell finished rc=%d log=%s", completed.returncode, log_file
    )
    if completed.returncode != 0:
        tail = ""
        try:
            with open(log_file) as fh:
                lines = fh.readlines()
                tail = "".join(lines[-40:])
        except OSError:
            pass
        raise subprocess.CalledProcessError(
            completed.returncode, cmd, output=tail
        )


def check_all_samples_have_bams(data_folder, exons, desired_samples):
    """
    Require every samples_list entry to have a non-empty remapped BAM (and .bai)
    in every chunk dir. Raises RuntimeError with a per-chunk summary.
    """
    problems = []
    for exon in exons:
        logger.debug("Checking BAMs for chunk %s", exon)
        missing = []
        empty = []
        no_index = []
        for sample in desired_samples:
            path = _bam_path(data_folder, exon, sample)
            if not os.path.isfile(path):
                missing.append(sample)
                logger.debug("  missing BAM: %s", path)
            elif os.path.getsize(path) == 0:
                empty.append(sample)
                logger.debug("  empty BAM: %s", path)
            elif not os.path.isfile(path + ".bai"):
                no_index.append(sample)
                logger.debug("  missing index: %s.bai", path)
            else:
                logger.debug(
                    "  OK %s (%d bytes)", path, os.path.getsize(path)
                )
        parts = []
        if missing:
            parts.append(f"missing BAM for {len(missing)}: {', '.join(missing)}")
        if empty:
            parts.append(f"empty BAM for {len(empty)}: {', '.join(empty)}")
        if no_index:
            parts.append(f"missing .bai for {len(no_index)}: {', '.join(no_index)}")
        if parts:
            problems.append(f"{exon}: " + "; ".join(parts))
        else:
            logger.debug("Chunk %s: all %d BAM(s) OK", exon, len(desired_samples))

    if problems:
        details = "\n  - ".join(problems)
        logger.error(
            "Remapped BAM check failed. Fix cast_remap output or edit "
            "samples_list.txt:\n  - %s",
            details,
        )
        raise RuntimeError(
            "cast_call aborted: remapped BAM check failed. See log for details."
        )


def _vcf_samples(vcf_path):
    result = subprocess.run(
        ["bcftools", "query", "-l", vcf_path],
        capture_output=True,
        text=True,
        check=True,
    )
    samples = [line.strip() for line in result.stdout.splitlines() if line.strip()]
    logger.debug("bcftools query -l %s -> %d sample(s)", vcf_path, len(samples))
    return samples


def _output_is_complete(vcf_path, desired_samples):
    """True if vcf_path exists and its sample order exactly matches desired_samples."""
    if not os.path.isfile(vcf_path) or os.path.getsize(vcf_path) == 0:
        logger.debug("Output incomplete (missing/empty): %s", vcf_path)
        return False
    try:
        samples = _vcf_samples(vcf_path)
    except (subprocess.CalledProcessError, OSError) as e:
        logger.warning("Could not read samples from %s (%s); will redo", vcf_path, e)
        return False
    ok = samples == desired_samples
    logger.debug(
        "Output complete check %s: %s (have %d, want %d)",
        vcf_path,
        ok,
        len(samples),
        len(desired_samples),
    )
    if not ok:
        logger.debug("  have: %s", ", ".join(samples[:20]))
        logger.debug("  want: %s", ", ".join(desired_samples[:20]))
    return ok


def _open_vcf(path, mode):
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t")
    return open(path, mode)


def _reorder_vcf_samples(vcf_path, samples_list_path, exon, out_path=None):
    """
    Normalize sample column names and write a BGZF VCF ordered as samples_list.txt.

    Raises RuntimeError if any listed sample is absent from the input VCF.
    out_path defaults to overwriting vcf_path; cast_call writes the corrected
    *_variants_filtered_corrected.vcf.gz here.
    """
    require_file(vcf_path, f"filtered VCF for {exon}")
    desired_samples = _read_samples_list(samples_list_path)
    out_path = out_path or vcf_path
    logger.debug(
        "Reordering VCF samples: in=%s out=%s chunk=%s",
        vcf_path,
        out_path,
        exon,
    )

    with _open_vcf(vcf_path, "r") as fh:
        lines = fh.readlines()
    logger.debug("Read %d line(s) from %s", len(lines), vcf_path)

    header_idx = None
    for i, line in enumerate(lines):
        if line.startswith("#CHROM"):
            header_idx = i
            break
    if header_idx is None:
        raise ValueError(f"No #CHROM header found in {vcf_path}")

    header_fields = lines[header_idx].rstrip("\n").split("\t")
    fixed_cols = header_fields[:9]
    raw_sample_cols = header_fields[9:]
    if not raw_sample_cols:
        raise ValueError(f"No sample columns in VCF header: {vcf_path}")
    logger.debug(
        "Chunk %s: %d raw VCF sample column(s); first=%r",
        exon,
        len(raw_sample_cols),
        raw_sample_cols[0] if raw_sample_cols else None,
    )

    clean_to_raw = {}
    for raw in raw_sample_cols:
        clean = _clean_sample_name(raw)
        if clean in clean_to_raw:
            logger.warning(
                "Chunk %s: duplicate sample name after cleanup: %s", exon, clean
            )
        clean_to_raw[clean] = raw

    missing = _missing_samples(desired_samples, clean_to_raw)
    unexpected = sorted(set(clean_to_raw) - set(desired_samples))
    logger.debug(
        "Chunk %s sample map: cleaned=%d missing=%d unexpected=%d",
        exon,
        len(clean_to_raw),
        len(missing),
        len(unexpected),
    )

    if missing:
        logger.error(
            "Chunk %s: %d sample(s) from samples_list.txt missing in VCF: %s",
            exon,
            len(missing),
            ", ".join(missing),
        )
        raise RuntimeError(
            f"Chunk {exon}: missing samples in VCF: {', '.join(missing)}"
        )
    if unexpected:
        logger.warning(
            "Chunk %s: %d VCF sample(s) not in samples_list.txt (ignored): %s%s",
            exon,
            len(unexpected),
            ", ".join(unexpected[:10]),
            " ..." if len(unexpected) > 10 else "",
        )

    raw_index = {raw: idx for idx, raw in enumerate(raw_sample_cols)}
    new_header = fixed_cols + desired_samples
    out_lines = lines[:header_idx] + ["\t".join(new_header) + "\n"]

    n_sites = 0
    for line in lines[header_idx + 1 :]:
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        site_cols = fields[:9]
        sample_cols = fields[9:]
        if len(sample_cols) != len(raw_sample_cols):
            raise ValueError(
                f"Malformed VCF record in {vcf_path}: expected "
                f"{len(raw_sample_cols)} sample fields, got {len(sample_cols)}"
            )
        ordered = [
            sample_cols[raw_index[clean_to_raw[sample]]] for sample in desired_samples
        ]
        out_lines.append("\t".join(site_cols + ordered) + "\n")
        n_sites += 1

    logger.info(
        "Chunk %s: writing corrected VCF with %d sample(s), %d site(s)",
        exon,
        len(desired_samples),
        n_sites,
    )

    folder = os.path.dirname(out_path)
    temp_plain = os.path.join(folder, f"temp_variants_{exon}.vcf")
    with open(temp_plain, "w") as fh:
        fh.writelines(out_lines)
    logger.debug("Wrote plain VCF temp %s (%d bytes)", temp_plain, os.path.getsize(temp_plain))
    # bgzip → BGZF (required by bcftools index / tabix; Python gzip is not enough).
    try:
        subprocess.run(
            ["bgzip", "-f", temp_plain],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True,
        )
    except subprocess.CalledProcessError as e:
        logger.error("bgzip failed for %s: %s", temp_plain, e.stderr)
        raise
    temp_gz = temp_plain + ".gz"
    os.replace(temp_gz, out_path)
    require_file(out_path, f"corrected VCF for {exon}")
    logger.debug("Corrected VCF ready: %s", out_path)


def _index_vcf(vcf_path):
    logger.debug("Indexing VCF: %s", vcf_path)
    try:
        subprocess.run(
            ["bcftools", "index", "-f", "-t", vcf_path],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True,
        )
    except subprocess.CalledProcessError as e:
        logger.error("bcftools index failed for %s: %s", vcf_path, e.stderr)
        raise
    logger.debug("Indexed VCF: %s", vcf_path)


# -----------------------------------------------------------------------------
# Per-chunk and all-chunk calling
# -----------------------------------------------------------------------------
@log_exceptions
def variant_call(exon, main_data_folder, n_threads=1, force=False):
    """
    Call and filter variants for one remapped chunk (e.g. exons2).

    Uses bam_list.txt in samples_list order, runs
    bcftools mpileup | call | filter → {exon}_variants_filtered.vcf.gz, then
    writes {exon}_variants_filtered_corrected.vcf.gz with samples reordered and
    indexed. Skips if that corrected VCF already has the expected sample set
    unless force=True.
    """
    samples_list_dst = os.path.join(main_data_folder, "100remapped", "samples_list.txt")
    desired_samples = _read_samples_list(samples_list_dst)
    corrected_vcf_file = _corrected_vcf_path(main_data_folder, exon)
    exon_dir = require_dir(
        os.path.join(main_data_folder, "100remapped", exon),
        f"chunk directory for {exon}",
    )
    logger.debug(
        "variant_call start chunk=%s force=%s threads=%d out=%s",
        exon,
        force,
        n_threads,
        corrected_vcf_file,
    )

    if not force and _output_is_complete(corrected_vcf_file, desired_samples):
        logger.info(
            "Skipping chunk %s: complete output already exists (%s)",
            exon,
            corrected_vcf_file,
        )
        return

    reference_file = require_file(
        os.path.join(exon_dir, f"reference_{exon}.fas"),
        f"reference FASTA for {exon}",
    )

    bam_files = []
    missing_bams = []
    for sample in desired_samples:
        path = _bam_path(main_data_folder, exon, sample)
        if not os.path.isfile(path) or os.path.getsize(path) == 0:
            missing_bams.append(sample)
        else:
            bam_files.append(path)
            logger.debug("Using BAM %s (%d bytes)", path, os.path.getsize(path))
    if missing_bams:
        logger.error(
            "Chunk %s: %d sample(s) lack usable remapped BAMs: %s",
            exon,
            len(missing_bams),
            ", ".join(missing_bams),
        )
        raise RuntimeError(
            f"Chunk {exon}: missing remapped BAMs for: {', '.join(missing_bams)}"
        )

    bam_list_file = os.path.join(exon_dir, "bam_list.txt")
    with open(bam_list_file, "w") as fh:
        fh.write("\n".join(bam_files) + "\n")
    logger.debug("Wrote bam_list.txt with %d path(s): %s", len(bam_files), bam_list_file)
    logger.info(
        "Chunk %s: %d BAM(s), reference=%s, threads=%d",
        exon,
        len(bam_files),
        reference_file,
        n_threads,
    )

    filtered_vcf = os.path.join(exon_dir, f"{exon}_variants_filtered.vcf.gz")
    call_log = os.path.join(exon_dir, f"{exon}_bcftools_call.log")

    mpileup_call_filter_cmd = (
        f"bcftools mpileup -Ou -I -d 2000 -A --threads {n_threads} "
        f"-f {reference_file} -b {bam_list_file} | "
        f"bcftools call -Ou -mv --threads {n_threads} | "
        f"bcftools filter -e 'QUAL<20' -Oz -o {filtered_vcf}"
    )
    logger.info("Running bcftools mpileup|call|filter for chunk %s", exon)
    try:
        _run_shell(mpileup_call_filter_cmd, call_log)
    except subprocess.CalledProcessError as e:
        logger.error(
            "Variant calling failed for chunk %s. See %s\n%s",
            exon,
            call_log,
            e.output or "",
        )
        raise

    require_file(filtered_vcf, f"filtered VCF for {exon}")
    logger.debug(
        "Filtered VCF size for %s: %d bytes",
        exon,
        os.path.getsize(filtered_vcf),
    )
    logger.info("Ordering samples in VCF for chunk %s", exon)
    try:
        _reorder_vcf_samples(
            filtered_vcf, samples_list_dst, exon, out_path=corrected_vcf_file
        )
        _index_vcf(corrected_vcf_file)
    except Exception as e:
        logger.error("Sample ordering/index failed for %s: %s", corrected_vcf_file, e)
        raise

    logger.info("Variant calling completed for chunk %s -> %s", exon, corrected_vcf_file)


@log_exceptions
def call_variants(data_folder, num_cores, log_queue, force=False):
    """
    Run cast_call for all 100remapped/exons* chunks, then concat to all_variants.vcf.gz.

    :param data_folder: pipeline data root
    :param num_cores: CPUs from -nc (split across chunk workers and bcftools threads)
    :param log_queue: multiprocessing logging queue
    :param force: redo chunks even when corrected VCFs are already complete
    """
    logger.info("Starting cast_call")
    logger.debug(
        "call_variants args: data_folder=%s num_cores=%d force=%s",
        data_folder,
        num_cores,
        force,
    )
    require_tools(REQUIRED_TOOLS)
    require_dir(data_folder, "data folder")

    remapped_dir = require_dir(
        os.path.join(data_folder, "100remapped"), "100remapped directory"
    )
    samples_list_src = require_file(
        os.path.join(data_folder, "10deduplicated_reads", "samples_list.txt"),
        "samples_list.txt",
    )
    samples_list_dst = os.path.join(remapped_dir, "samples_list.txt")
    shutil.copy(samples_list_src, samples_list_dst)
    logger.debug("Copied samples list -> %s", samples_list_dst)
    desired_samples = _read_samples_list(samples_list_dst)
    logger.info("Loaded %d sample(s) from %s", len(desired_samples), samples_list_src)

    exons_dir_pattern = os.path.join(remapped_dir, "exons*")
    exons = sorted(
        os.path.basename(x)
        for x in glob(exons_dir_pattern)
        if os.path.isdir(x)
    )
    if not exons:
        raise FileNotFoundError(
            f"No exon chunk directories matching {exons_dir_pattern}. "
            "Run cast_remap first."
        )
    logger.info("Found %d remapped chunk(s): %s", len(exons), ", ".join(exons))
    logger.debug("Chunk dirs: %s", exons)

    logger.info(
        "Checking remapped BAMs for %d sample(s) across %d chunk(s)",
        len(desired_samples),
        len(exons),
    )
    check_all_samples_have_bams(data_folder, exons, desired_samples)
    logger.info("All samples have remapped BAMs (+ indexes) for all chunks")

    n_workers, threads_per_worker = allocate_workers_and_threads(num_cores, len(exons))
    logger.info(
        "Parallelism: %d chunk worker(s) x %d bcftools thread(s) "
        "(requested cores=%d, chunks=%d, force=%s)",
        n_workers,
        threads_per_worker,
        num_cores,
        len(exons),
        force,
    )
    logger.debug(
        "Pool plan: workers=%d threads_per_worker=%d tasks=%d",
        n_workers,
        threads_per_worker,
        len(exons),
    )

    with managed_pool(n_workers, log_queue) as pool_call:
        async_results = [
            (
                exon,
                pool_call.apply_async(
                    variant_call,
                    (exon, data_folder),
                    {"n_threads": threads_per_worker, "force": force},
                ),
            )
            for exon in exons
        ]
        logger.debug("Submitted %d chunk job(s) to pool", len(async_results))
        for i, (exon, async_result) in enumerate(async_results, start=1):
            try:
                async_result.get()
                logger.info("Finished chunk %s (%d/%d)", exon, i, len(exons))
            except Exception as e:
                logger.error(
                    "Error during variant calling for chunk %s (%d/%d): %s",
                    exon,
                    i,
                    len(exons),
                    e,
                )
                raise

    vcf_files = [_corrected_vcf_path(data_folder, exon) for exon in exons]
    for path in vcf_files:
        require_file(path, "corrected chunk VCF")
        logger.debug("Concat input: %s", path)

    all_variants = os.path.join(remapped_dir, ALL_VARIANTS_VCF)
    concat_log = os.path.join(remapped_dir, "bcftools_concat.log")
    concat_cmd = (
        f"bcftools concat -Oz -o {all_variants} {' '.join(vcf_files)} && "
        f"bcftools index -f -t {all_variants}"
    )
    logger.info("Concatenating %d VCF(s) -> %s", len(vcf_files), all_variants)
    try:
        _run_shell(concat_cmd, concat_log)
    except subprocess.CalledProcessError as e:
        logger.error(
            "bcftools concat failed. See %s\n%s", concat_log, e.output or ""
        )
        raise
    require_file(all_variants, "concatenated VCF")
    logger.debug(
        "Final VCF %s (%d bytes)", all_variants, os.path.getsize(all_variants)
    )

    logger.info("cast_call completed: %s", all_variants)
