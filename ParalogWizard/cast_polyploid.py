#!/usr/bin/env python
"""
------------------------------------------------------------------------------------------------------------------------
Copyright 2024 Roman Ufimov under the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later version.

cast_polyploid — export, filter, align, tree, and unsupervised hybrid phasing.

  1. Read 41detected_par/all_paralogs_for_reference.tsv.
  2. Keep only samples in backbone_list + polyploid_list; write
     101polyploid/{exon}_{main|para}.fasta.
  3. Delete FASTAs that lack any polyploid sample (or have < MIN_SEQS_FOR_TREE).
  4. MAFFT-align remaining FASTAs → *.fasta.mafft.
  5. Build a simple IQ-TREE tree per alignment → *.treefile
     (auto-picks highest available of iqtree3 / iqtree3-mpi / iqtree2 /
     iqtree2-mpi / iqtree, plus any other iqtree* on PATH).
  6. Unsupervised phasing of each polyploid (several allowed in one run):
       - Distances / top-N nearest ranks use backbone tips only
         (other polyploids are ignored for ranking; N default 5, -ns).
       - Cluster each polyploid's tips into A / B.
       - Enforce same-contig consistency within each polyploid gene
         (main/para = separate genes).
       - Write tables + phased FASTAs (phase letter on species name;
         one sequence per sample; longest contig if several).
  7. Concatenate phased exons by gene with AMAS into
     101polyploid/concatenated/Assembly_{gene}.fasta|.part|.phylip
     (does not touch 70concatenated_exon_alignments/).

Lists default to backbone_list.txt / polyploid_list.txt in the working
directory (or under -d). Each run clears and rebuilds 101polyploid/.
------------------------------------------------------------------------------------------------------------------------
"""

from __future__ import annotations

import logging
import os
import re
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from glob import glob
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from Bio import Phylo
from sklearn.cluster import KMeans
from sklearn.metrics.pairwise import cosine_distances

from ParalogWizard import log_exceptions, managed_pool
from ParalogWizard.cast_analyze import mafft_align
from ParalogWizard.cast_call import (
    allocate_workers_and_threads,
    require_dir,
    require_file,
    require_tools,
)

logger = logging.getLogger("ParalogWizard")

POLYPLOID_DIRNAME = "101polyploid"
DETECTED_DIRNAME = "41detected_par"
PARALOGS_TSV = "all_paralogs_for_reference.tsv"
REQUIRED_COLUMNS = ("qaccver", "sample", "sequence", "contig", "copy")
COPY_VALUES = ("main", "para")
REQUIRED_TOOLS = ("mafft",)
# Serial names first so they win a version tie over MPI/OpenMP builds.
IQTREE_CANDIDATES = (
    "iqtree3",
    "iqtree3-mpi",
    "iqtree3-omp",
    "iqtree2",
    "iqtree2-mpi",
    "iqtree2-omp",
    "iqtree",
    "iqtree-mpi",
    "iqtree-omp",
)
DEFAULT_REQUIRED_SAMPLE = "necopinnata"  # legacy; prefer polyploid_list.txt
MIN_SEQS_FOR_TREE = 4
TOP_N_NEAREST = 5
# Cosine distance on rank-score vectors: tip pairs above this within a locus
# are forced onto opposite phases when global clustering left them together.
WITHIN_LOCUS_SPLIT_COSINE = 0.35
BACKBONE_LIST_DEFAULT = "backbone_list.txt"
POLYPLOID_LIST_DEFAULT = "polyploid_list.txt"


def _exon_from_qaccver(qaccver: str) -> str:
    """Return the exon id: everything after the first '-' in qaccver."""
    text = str(qaccver).strip()
    if "-" not in text:
        raise ValueError(f"qaccver has no '-': {qaccver!r}")
    exon = text.split("-", 1)[1].strip()
    if not exon:
        raise ValueError(f"Empty exon id derived from qaccver: {qaccver!r}")
    return exon


def _seq_header(sample: str, contig: str) -> str:
    return f"{sample}_{contig}"


def _load_paralogs_tsv(path: str) -> pd.DataFrame:
    require_file(path, PARALOGS_TSV)
    df = pd.read_csv(path, sep="\t")
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(
            f"{path} missing column(s): {', '.join(missing)}. "
            f"Found: {list(df.columns)}"
        )
    if df.empty:
        raise ValueError(f"{path} is empty")

    before = len(df)
    df = df.dropna(subset=list(REQUIRED_COLUMNS)).copy()
    if df.empty:
        raise ValueError(f"{path}: all rows missing required fields")
    if len(df) < before:
        logger.warning(
            "Dropped %d row(s) with missing required fields from %s",
            before - len(df),
            path,
        )

    df["copy"] = df["copy"].astype(str).str.strip().str.lower()
    bad_copy = sorted(set(df["copy"]) - set(COPY_VALUES))
    if bad_copy:
        raise ValueError(
            f"{path}: unexpected copy value(s) {bad_copy}; expected {list(COPY_VALUES)}"
        )

    df["exon_id"] = df["qaccver"].map(_exon_from_qaccver)
    df["seq_id"] = [
        _seq_header(str(s), str(c)) for s, c in zip(df["sample"], df["contig"])
    ]
    df["sequence"] = df["sequence"].astype(str).str.strip()
    empty_seq = df["sequence"].eq("") | df["sequence"].str.lower().eq("nan")
    if empty_seq.any():
        n = int(empty_seq.sum())
        raise ValueError(f"{path}: {n} row(s) have empty sequence")

    dups = df.duplicated(subset=["exon_id", "copy", "seq_id"], keep=False)
    if dups.any():
        examples = (
            df.loc[dups, ["exon_id", "copy", "seq_id"]]
            .drop_duplicates()
            .head(10)
            .to_dict("records")
        )
        raise ValueError(
            f"{path}: duplicate sequence names within the same exon+copy "
            f"(examples: {examples})"
        )

    logger.info(
        "Loaded %d contig hit(s) from %s (%d exon id(s), copy=%s)",
        len(df),
        path,
        df["exon_id"].nunique(),
        df["copy"].value_counts().to_dict(),
    )
    return df


def _write_fasta(path: str, records: Iterable[Tuple[str, str]]) -> int:
    n = 0
    with open(path, "w") as fh:
        for seq_id, sequence in records:
            fh.write(f">{seq_id}\n{sequence}\n")
            n += 1
    return n


def _fasta_headers(path: str) -> List[str]:
    headers = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                headers.append(line[1:].strip().split()[0])
    return headers


_IQTREE_CACHE: Optional[Tuple[str, Tuple[int, ...]]] = None


def _parse_iqtree_version(text: str) -> Optional[Tuple[int, ...]]:
    """Extract a version tuple from IQ-TREE --version / -version output."""
    if not text:
        return None
    m = re.search(
        r"(?:version\s+|IQ-TREE\s*)(\d+)(?:\.(\d+))?(?:\.(\d+))?",
        text,
        flags=re.IGNORECASE,
    )
    if not m:
        return None
    parts = [int(p) for p in m.groups() if p is not None]
    return tuple(parts) if parts else None


def _probe_iqtree_binary(name: str) -> Optional[Tuple[str, Tuple[int, ...]]]:
    """Return (absolute_path, version_tuple) if `name` is a usable IQ-TREE."""
    path = shutil.which(name)
    if not path:
        return None
    version: Optional[Tuple[int, ...]] = None
    for flag in ("--version", "-version", "-v"):
        try:
            result = subprocess.run(
                [path, flag],
                capture_output=True,
                text=True,
                check=False,
            )
        except OSError:
            return None
        blob = (result.stdout or "") + "\n" + (result.stderr or "")
        version = _parse_iqtree_version(blob)
        if version is not None:
            break
    if version is None:
        version = _version_from_binary_name(name)
        logger.warning(
            "Could not parse version for %s (%s); inferring %s from the name",
            name,
            path,
            ".".join(str(p) for p in version),
        )
    return path, version


def _version_from_binary_name(name: str) -> Tuple[int, ...]:
    """iqtree3-mpi → (3,), iqtree2 → (2,), unknown → (0,)."""
    base = os.path.basename(name)
    m = re.search(r"iqtree(\d+)", base, flags=re.IGNORECASE)
    if m:
        return (int(m.group(1)),)
    return (0,)


def _iqtree_names_on_path() -> List[str]:
    """Named candidates first, then any other executable whose name starts with iqtree."""
    names: List[str] = []
    seen = set()
    for name in IQTREE_CANDIDATES:
        if name not in seen:
            names.append(name)
            seen.add(name)
    for directory in os.environ.get("PATH", "").split(os.pathsep):
        if not directory or not os.path.isdir(directory):
            continue
        try:
            listing = os.listdir(directory)
        except OSError:
            continue
        for fname in sorted(listing):
            if not fname.startswith("iqtree") or fname in seen:
                continue
            full = os.path.join(directory, fname)
            if os.path.isfile(full) and os.access(full, os.X_OK):
                names.append(fname)
                seen.add(fname)
    return names


def resolve_iqtree(force: bool = False) -> str:
    """
    Pick the highest-version IQ-TREE on PATH.

    Tries iqtree3, iqtree3-mpi, iqtree2, iqtree2-mpi, iqtree, then any other
    iqtree* executable. Preference is by parsed version; at equal version a
    serial binary wins over MPI/OpenMP. Raises if none are found.
    """
    global _IQTREE_CACHE
    if _IQTREE_CACHE is not None and not force:
        return _IQTREE_CACHE[0]

    candidates = _iqtree_names_on_path()
    found: List[Tuple[Tuple[int, ...], str, str]] = []
    for name in candidates:
        probed = _probe_iqtree_binary(name)
        if probed is None:
            continue
        path, version = probed
        found.append((version, path, name))

    if not found:
        raise RuntimeError(
            "No IQ-TREE binary found on PATH. Tried: " + ", ".join(candidates)
        )

    # Highest version wins; tie-break by candidate order (serial iqtree3 first).
    order = {name: i for i, name in enumerate(candidates)}
    found.sort(key=lambda x: (x[0], -order[x[2]]), reverse=True)
    version, path, name = found[0]
    ver_str = ".".join(str(p) for p in version)
    logger.info("Using IQ-TREE binary %s (%s, version %s)", name, path, ver_str)
    for version_i, path_i, name_i in found[1:]:
        logger.debug(
            "Also found %s (%s, version %s)",
            name_i,
            path_i,
            ".".join(str(p) for p in version_i),
        )
    _IQTREE_CACHE = (path, version)
    return path


def _iqtree_major_version() -> int:
    resolve_iqtree()
    assert _IQTREE_CACHE is not None
    return int(_IQTREE_CACHE[1][0]) if _IQTREE_CACHE[1] else 0


@log_exceptions
def _iqtree_simple(alignment_path: str, n_threads: int = 1) -> str:
    """
    Simple IQ-TREE search (--fast / -fast, GTR+F+G, no bootstrap).

    Uses the highest available IQ-TREE on PATH (iqtree3, iqtree3-mpi, …).
    Writes <alignment_path>.treefile.
    """
    require_file(alignment_path, "alignment for IQ-TREE")
    n_seq = sum(1 for line in open(alignment_path) if line.startswith(">"))
    if n_seq < MIN_SEQS_FOR_TREE:
        raise ValueError(
            f"Need >= {MIN_SEQS_FOR_TREE} sequences for IQ-TREE, "
            f"got {n_seq} in {alignment_path}"
        )

    binary = resolve_iqtree()
    major = _iqtree_major_version()
    prefix = alignment_path
    treefile = f"{prefix}.treefile"
    threads = max(1, int(n_threads))

    if major >= 2:
        cmd = [
            binary,
            "-s",
            alignment_path,
            "--prefix",
            prefix,
            "-m",
            "GTR+F+G",
            "--fast",
            "-T",
            str(threads),
            "--quiet",
            "--redo",
        ]
    else:
        # IQ-TREE 1.x flag set
        cmd = [
            binary,
            "-s",
            alignment_path,
            "-pre",
            prefix,
            "-m",
            "GTR+F+G",
            "-fast",
            "-nt",
            str(threads),
            "-quiet",
            "-redo",
        ]

    logger.debug("IQ-TREE: %s", " ".join(cmd))
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        logger.error(
            "IQ-TREE failed for %s:\n%s\n%s",
            alignment_path,
            e.stdout or "",
            e.stderr or "",
        )
        raise
    if result.stderr:
        logger.debug("IQ-TREE stderr tail: %s", result.stderr[-500:])
    require_file(treefile, f"IQ-TREE treefile for {alignment_path}")
    logger.debug("IQ-TREE done: %s", treefile)
    return treefile


@log_exceptions
def _align_and_tree(fasta_path: str, n_threads: int = 1) -> Tuple[str, str]:
    """MAFFT + IQ-TREE for one FASTA. Returns (alignment, treefile)."""
    aln = mafft_align(fasta_path, n_threads=n_threads)
    tree = _iqtree_simple(aln, n_threads=n_threads)
    return aln, tree


def _export_fastas(df: pd.DataFrame, out_dir: str) -> int:
    groups: Dict[Tuple[str, str], List[Tuple[str, str]]] = defaultdict(list)
    for row in df.itertuples(index=False):
        groups[(row.exon_id, row.copy)].append((row.seq_id, row.sequence))

    written = 0
    for (exon_id, copy_label), records in sorted(groups.items()):
        out_path = os.path.join(out_dir, f"{exon_id}_{copy_label}.fasta")
        _write_fasta(out_path, records)
        written += 1
    n_main = sum(1 for (_, c) in groups if c == "main")
    n_para = sum(1 for (_, c) in groups if c == "para")
    logger.info(
        "Wrote %d FASTA file(s) (%d main, %d para)", written, n_main, n_para
    )
    return written


def _read_name_list(path: str, label: str) -> List[str]:
    require_file(path, label)
    with open(path) as fh:
        names = [ln.strip() for ln in fh if ln.strip() and not ln.strip().startswith("#")]
    if not names:
        raise ValueError(f"No names in {path}")
    dups = sorted({n for n in names if names.count(n) > 1})
    if dups:
        raise ValueError(f"Duplicate name(s) in {path}: {', '.join(dups)}")
    return names


def _resolve_list_path(arg: Optional[str], default_name: str, data_folder: str) -> str:
    """
    Resolve backbone/polyploid list path:
      1) explicit CLI path if it exists
      2) ./default_name (cwd)
      3) {data_folder}/default_name
    """
    candidates = []
    if arg:
        candidates.append(arg)
        if not os.path.isabs(arg):
            candidates.append(os.path.join(os.getcwd(), arg))
            candidates.append(os.path.join(data_folder, arg))
    candidates.extend(
        [
            os.path.join(os.getcwd(), default_name),
            os.path.join(data_folder, default_name),
        ]
    )
    seen = set()
    for path in candidates:
        path = os.path.abspath(path)
        if path in seen:
            continue
        seen.add(path)
        if os.path.isfile(path):
            return path
    raise FileNotFoundError(
        f"Could not find {default_name}. Tried: {', '.join(sorted(seen))}"
    )


def _read_samples_list(data_folder: str) -> List[str]:
    path = os.path.join(data_folder, "10deduplicated_reads", "samples_list.txt")
    require_file(path, "samples_list.txt")
    with open(path) as fh:
        samples = [ln.strip() for ln in fh if ln.strip()]
    if not samples:
        raise ValueError(f"No samples in {path}")
    return samples


def _tip_to_sample(tip: str, samples: Sequence[str]) -> Optional[str]:
    """Map a tip name '{sample}_{contig…}' back to a samples_list entry."""
    for sample in sorted(samples, key=len, reverse=True):
        if tip == sample or tip.startswith(sample + "_"):
            return sample
    return None


def _filter_df_to_allowed_samples(
    df: pd.DataFrame, allowed: Sequence[str]
) -> pd.DataFrame:
    allowed_set = set(allowed)
    before = len(df)
    out = df[df["sample"].astype(str).isin(allowed_set)].copy()
    logger.info(
        "Sample filter (backbone+polyploid): kept %d / %d contig hit(s); "
        "%d sample(s) allowed",
        len(out),
        before,
        len(allowed_set),
    )
    if out.empty:
        raise RuntimeError(
            "No contig hits left after filtering to backbone_list + polyploid_list"
        )
    return out


def _fasta_has_any_sample(path: str, samples: Sequence[str]) -> bool:
    headers = _fasta_headers(path)
    for tip in headers:
        if _tip_to_sample(tip, samples) is not None:
            return True
    # fallback: substring match for short needles
    lower_headers = [h.lower() for h in headers]
    for sample in samples:
        needle = sample.lower()
        if any(needle in h for h in lower_headers):
            return True
    return False


def _filter_fastas_missing_polyploid(
    out_dir: str, polyploid_samples: Sequence[str]
) -> List[str]:
    """
    Delete *.fasta lacking any polyploid sample, or with too few sequences
    for IQ-TREE. Return kept paths.
    """
    kept: List[str] = []
    removed_no_poly = 0
    removed_too_few = 0
    for path in sorted(glob(os.path.join(out_dir, "*.fasta"))):
        if path.endswith(".fasta.mafft"):
            continue
        if not _fasta_has_any_sample(path, polyploid_samples):
            os.remove(path)
            removed_no_poly += 1
            logger.debug(
                "Removed (no polyploid sample): %s", os.path.basename(path)
            )
            continue
        n_seq = len(_fasta_headers(path))
        if n_seq < MIN_SEQS_FOR_TREE:
            os.remove(path)
            removed_too_few += 1
            logger.info(
                "Removed (<%d seqs, n=%d): %s",
                MIN_SEQS_FOR_TREE,
                n_seq,
                os.path.basename(path),
            )
            continue
        kept.append(path)
    logger.info(
        "Polyploid filter (%d sample(s)): kept %d FASTA(s), "
        "deleted %d (no polyploid), %d (<%d seqs)",
        len(polyploid_samples),
        len(kept),
        removed_no_poly,
        removed_too_few,
        MIN_SEQS_FOR_TREE,
    )
    if not kept:
        raise RuntimeError(
            f"No FASTA files contain any polyploid sample under {out_dir}: "
            f"{list(polyploid_samples)}"
        )
    return kept

def _hybrid_distance_profiles(
    treefiles: Sequence[str],
    backbone_samples: Sequence[str],
    polyploid_samples: Sequence[str],
) -> pd.DataFrame:
    """
    One row per polyploid tip; columns = backbone samples only.
    Other polyploids are ignored when ranking nearest tips.
    Value = phylogenetic distance to the nearest tip of that backbone sample.
    """
    if not backbone_samples:
        raise RuntimeError("backbone_list is empty")
    if not polyploid_samples:
        raise RuntimeError("polyploid_list is empty")

    all_known = list(backbone_samples) + list(polyploid_samples)
    poly_set = set(polyploid_samples)
    back_set = set(backbone_samples)

    rows = []
    for treefile in treefiles:
        tree = Phylo.read(treefile, "newick")
        tips = [t.name for t in tree.get_terminals() if t.name]
        tip_sample = {tip: _tip_to_sample(tip, all_known) for tip in tips}

        by_backbone: Dict[str, List[str]] = defaultdict(list)
        poly_tips: List[Tuple[str, str]] = []  # (tip, polyploid_sample)
        for tip, sample in tip_sample.items():
            if sample is None:
                logger.debug(
                    "Ignoring unlisted tip %s in %s",
                    tip,
                    os.path.basename(treefile),
                )
                continue
            if sample in poly_set:
                poly_tips.append((tip, sample))
            elif sample in back_set:
                by_backbone[sample].append(tip)

        locus = os.path.basename(treefile).replace(".fasta.mafft.treefile", "")
        # n hybrid tips counted per polyploid sample at this locus
        n_by_poly = Counter(s for _t, s in poly_tips)
        for hybrid_tip, poly_sample in poly_tips:
            row = {
                "locus": locus,
                "treefile": os.path.basename(treefile),
                "hybrid_tip": hybrid_tip,
                "polyploid_sample": poly_sample,
                "n_hybrid_tips": int(n_by_poly[poly_sample]),
            }
            for sample in backbone_samples:
                refs = by_backbone.get(sample, [])
                if not refs:
                    row[sample] = np.nan
                else:
                    row[sample] = min(tree.distance(hybrid_tip, r) for r in refs)
            rows.append(row)

    if not rows:
        raise RuntimeError(
            f"No polyploid tips matching {list(polyploid_samples)!r} found in trees"
        )
    return pd.DataFrame(rows)


def _feature_cols(profiles: pd.DataFrame) -> List[str]:
    meta = {
        "locus",
        "treefile",
        "hybrid_tip",
        "n_hybrid_tips",
        "polyploid_sample",
        "phase",
    }
    return [c for c in profiles.columns if c not in meta]

def _nearest_species_table(
    profiles: pd.DataFrame,
    feature_cols: Sequence[str],
    top_n: int = TOP_N_NEAREST,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    From tip×backbone distances, build:
      - long table of top-N nearest backbone species (rank 1 = closest)
      - wide rank-score matrix (score = top_n - rank + 1; 0 if not in top-N)
    """
    long_rows: List[dict] = []
    rank_rows: List[dict] = []

    for row in profiles.to_dict(orient="records"):
        poly = row.get("polyploid_sample", "")
        dists = []
        for sample in feature_cols:
            d = row.get(sample)
            if d is None or (isinstance(d, float) and np.isnan(d)):
                continue
            dists.append((float(d), sample))
        dists.sort(key=lambda x: (x[0], x[1]))
        top = dists[:top_n]

        rank_scores = {sample: 0 for sample in feature_cols}
        nearest_names = []
        nearest_dists = []
        for rank, (dist, sample) in enumerate(top, start=1):
            score = top_n - rank + 1
            rank_scores[sample] = score
            nearest_names.append(sample)
            nearest_dists.append(dist)
            long_rows.append(
                {
                    "polyploid_sample": poly,
                    "locus": row["locus"],
                    "treefile": row["treefile"],
                    "hybrid_tip": row["hybrid_tip"],
                    "n_hybrid_tips": row["n_hybrid_tips"],
                    "rank": rank,
                    "species": sample,
                    "distance": dist,
                    "rank_score": score,
                }
            )

        rank_row = {
            "polyploid_sample": poly,
            "locus": row["locus"],
            "treefile": row["treefile"],
            "hybrid_tip": row["hybrid_tip"],
            "n_hybrid_tips": row["n_hybrid_tips"],
        }
        for i in range(top_n):
            rank_row["nearest_%d" % (i + 1)] = (
                nearest_names[i] if i < len(nearest_names) else ""
            )
            rank_row["dist_%d" % (i + 1)] = (
                nearest_dists[i] if i < len(nearest_dists) else np.nan
            )
        rank_row.update(rank_scores)
        rank_rows.append(rank_row)

    nearest_long = pd.DataFrame(long_rows)
    rank_matrix = pd.DataFrame(rank_rows)
    logger.info(
        "Nearest-species table: %d tip×rank row(s); rank matrix %d tip(s) × %d backbone",
        len(nearest_long),
        len(rank_matrix),
        len(feature_cols),
    )
    return nearest_long, rank_matrix


def _cluster_one_polyploid_rank_matrix(
    rank_matrix: pd.DataFrame,
    feature_cols: Sequence[str],
) -> Tuple[pd.DataFrame, np.ndarray]:
    """Cluster tips of one polyploid into A/B. Returns (assignment frame, phase array)."""
    meta_cols = [
        "polyploid_sample",
        "locus",
        "treefile",
        "hybrid_tip",
        "n_hybrid_tips",
    ]
    meta_cols.extend(
        c
        for c in rank_matrix.columns
        if c.startswith("nearest_") or c.startswith("dist_")
    )
    # tolerate missing polyploid_sample for legacy frames
    use_cols = [c for c in meta_cols if c in rank_matrix.columns]
    out = rank_matrix[use_cols].copy()

    X = rank_matrix[list(feature_cols)].to_numpy(dtype=float)
    if X.shape[0] == 1:
        out["phase"] = ["A"]
        out["cluster_id"] = [0]
        out["locus_pattern"] = ["single_or_same"]
        return out, np.array(["A"], dtype=object)

    norms = np.linalg.norm(X, axis=1, keepdims=True)
    norms[norms < 1e-12] = 1.0
    X_n = X / norms

    kmeans = KMeans(n_clusters=2, n_init=50, random_state=0)
    labels = kmeans.fit_predict(X_n)
    counts = pd.Series(labels).value_counts()
    major = int(counts.idxmax())
    label_map = {major: "A", 1 - major: "B"}
    phase = np.array([label_map[int(x)] for x in labels], dtype=object)

    loci = rank_matrix["locus"].to_numpy()
    for locus in sorted(set(loci)):
        idx = np.where(loci == locus)[0]
        if len(idx) < 2:
            continue
        if set(phase[idx]) == {"A", "B"}:
            continue
        D = cosine_distances(X_n[idx])
        i_local, j_local = np.unravel_index(np.argmax(D), D.shape)
        if D[i_local, j_local] < WITHIN_LOCUS_SPLIT_COSINE:
            continue
        axis = kmeans.cluster_centers_[1 - major] - kmeans.cluster_centers_[major]
        scores = X_n[idx] @ axis
        order = np.argsort(scores)
        phase[idx[order[0]]] = "A"
        phase[idx[order[-1]]] = "B"
        for k in order[1:-1]:
            phase[idx[k]] = "A" if scores[k] < 0 else "B"
        logger.info(
            "Locus %s (%s): split hybrid tips into A/B "
            "(max cosine distance on rank vectors=%.3f)",
            locus,
            out["polyploid_sample"].iloc[0] if "polyploid_sample" in out.columns else "?",
            float(D[i_local, j_local]),
        )

    out["phase"] = phase
    out["cluster_id"] = [0 if p == "A" else 1 for p in phase]
    both = {
        loc
        for loc, g in out.groupby("locus")
        if set(g["phase"]) == {"A", "B"}
    }
    out["locus_pattern"] = out["locus"].map(
        lambda loc: "both_parents" if loc in both else "single_or_same"
    )
    return out, phase


def _cluster_rank_matrix_into_ab(
    rank_matrix: pd.DataFrame,
    feature_cols: Sequence[str],
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Cluster each polyploid’s tips into A/B using backbone nearest-species ranks.
    """
    if rank_matrix.empty:
        raise RuntimeError("Empty rank matrix")

    if "polyploid_sample" not in rank_matrix.columns:
        rank_matrix = rank_matrix.copy()
        rank_matrix["polyploid_sample"] = "unknown"

    parts = []
    phase_all = []
    for poly, sub in rank_matrix.groupby("polyploid_sample", sort=True):
        logger.info(
            "Clustering polyploid %s (%d tip profile(s))", poly, len(sub)
        )
        part, phase = _cluster_one_polyploid_rank_matrix(sub, feature_cols)
        parts.append(part)
        phase_all.append(sub.assign(phase=phase))

    out = pd.concat(parts, ignore_index=True)
    ranked = pd.concat(phase_all, ignore_index=True)
    both_n = out.loc[out["locus_pattern"] == "both_parents", "locus"].nunique()
    logger.info(
        "Phased %d hybrid tip(s) across %d polyploid(s): A=%d, B=%d; "
        "loci with both A+B=%d / %d",
        len(out),
        int(out["polyploid_sample"].nunique()),
        int((out["phase"] == "A").sum()),
        int((out["phase"] == "B").sum()),
        both_n,
        int(out["locus"].nunique()),
    )
    return out, ranked

def _gene_from_locus(locus: str) -> str:
    """
    Locus gene id with main/para kept separate (different genes).
    4848_exon_3_main → 4848_main; 4848_exon_5_para → 4848_para.
    """
    m = re.match(r"^(.+)_exon_\d+_(main|para)$", str(locus))
    if m:
        return f"{m.group(1)}_{m.group(2)}"
    return str(locus)


def _assembled_contig_from_tip(tip: str) -> str:
    """
    Contig id after sample name (before 2nd '_' is sample).
    Full SPAdes id `{node}_{length}_c_{cov}` — length/cov are of the assembled
    contig, so the same id across exons means the exons were cut from one contig.
    """
    _sample, contig = _sample_and_contig(tip)
    return contig


def _contig_length_from_id(contig: str) -> int:
    m = re.match(r"^(\d+)_(\d+)_c_([0-9.]+)$", contig)
    return int(m.group(2)) if m else 0


def _enforce_shared_contig_phases(
    assignments: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    If several exons of one gene (main/para separate) for the same polyploid
    were cut from the same assembled contig, force those tips to the same phase.

    Conflict resolution: majority vote; on tie, phase from earliest locus name.
    """
    out = assignments.copy()
    if "polyploid_sample" not in out.columns:
        out["polyploid_sample"] = "unknown"
    out["gene"] = out["locus"].map(_gene_from_locus)
    out["assembled_contig"] = out["hybrid_tip"].map(_assembled_contig_from_tip)
    out["contig_length"] = out["assembled_contig"].map(_contig_length_from_id)
    out["phase_before_contig_link"] = out["phase"]

    link_rows: List[dict] = []
    n_corrected = 0
    group_cols = ["polyploid_sample", "gene", "assembled_contig"]

    for (poly, gene, contig), g in out.groupby(group_cols, sort=False):
        if not contig or g["locus"].nunique() < 2:
            continue

        phases = list(g["phase"])
        counts = Counter(phases)
        if len(counts) == 1:
            consensus = phases[0]
            conflict = False
        else:
            conflict = True
            top = counts.most_common()
            if len(top) >= 2 and top[0][1] == top[1][1]:
                first_locus = sorted(g["locus"].unique())[0]
                consensus = g.loc[g["locus"] == first_locus, "phase"].iloc[0]
            else:
                consensus = top[0][0]

        loci = sorted(g["locus"].unique())
        link_rows.append(
            {
                "polyploid_sample": poly,
                "gene": gene,
                "assembled_contig": contig,
                "contig_length": int(g["contig_length"].iloc[0]),
                "n_exons": len(loci),
                "exons": ",".join(loci),
                "n_tips": len(g),
                "phases_before": ",".join(sorted(set(phases))),
                "consensus_phase": consensus,
                "conflict": conflict,
            }
        )

        for idx in g.index:
            if out.at[idx, "phase"] != consensus:
                logger.info(
                    "Contig link %s / %s / %s: %s @ %s  %s → %s",
                    poly,
                    gene,
                    contig,
                    out.at[idx, "hybrid_tip"],
                    out.at[idx, "locus"],
                    out.at[idx, "phase"],
                    consensus,
                )
                out.at[idx, "phase"] = consensus
                n_corrected += 1

    out["cluster_id"] = [0 if p == "A" else 1 for p in out["phase"]]
    out["contig_phase_corrected"] = out["phase"] != out["phase_before_contig_link"]
    both = {
        (poly, loc)
        for (poly, loc), g in out.groupby(["polyploid_sample", "locus"])
        if set(g["phase"]) == {"A", "B"}
    }
    out["locus_pattern"] = [
        "both_parents" if (r.polyploid_sample, r.locus) in both else "single_or_same"
        for r in out.itertuples(index=False)
    ]

    links = pd.DataFrame(link_rows)
    if not links.empty:
        links = links.sort_values(
            ["conflict", "contig_length", "polyploid_sample", "gene"],
            ascending=[False, False, True, True],
        ).reset_index(drop=True)

    logger.info(
        "Shared-contig phase check: %d multi-exon contig link(s), "
        "%d tip phase correction(s), %d conflict(s) resolved "
        "(per polyploid; main/para separate genes)",
        len(links),
        n_corrected,
        int(links["conflict"].sum()) if len(links) else 0,
    )
    return out, links

def _phase_nearest_summary(
    nearest_long: pd.DataFrame, assignments: pd.DataFrame, top_n: int = TOP_N_NEAREST
) -> pd.DataFrame:
    """Per polyploid × phase: how often each backbone species appears at each rank."""
    cols = ["hybrid_tip", "treefile", "phase"]
    if "polyploid_sample" in assignments.columns:
        cols = ["polyploid_sample"] + cols
    merged = nearest_long.merge(
        assignments[cols],
        on=["hybrid_tip", "treefile"]
        + (["polyploid_sample"] if "polyploid_sample" in nearest_long.columns and "polyploid_sample" in assignments.columns else []),
        how="left",
        suffixes=("", "_asn"),
    )
    if "polyploid_sample" not in merged.columns:
        merged["polyploid_sample"] = "unknown"
    top_col = "times_in_top%d" % int(top_n)
    rows = []
    for (poly, phase), g in merged.groupby(["polyploid_sample", "phase"]):
        n_tips = g[["hybrid_tip", "treefile"]].drop_duplicates().shape[0]
        for species, sg in g.groupby("species"):
            rank_counts = Counter(sg["rank"].astype(int))
            rec = {
                "polyploid_sample": poly,
                "phase": phase,
                "species": species,
                "n_tips_in_phase": n_tips,
                "top_n": int(top_n),
                top_col: int(len(sg)),
            }
            for rank in range(1, int(top_n) + 1):
                rec["times_rank%d" % rank] = int(rank_counts.get(rank, 0))
            rec["mean_rank"] = float(sg["rank"].mean())
            rec["mean_distance"] = float(sg["distance"].mean())
            rows.append(rec)
    summary = pd.DataFrame(rows)
    if not summary.empty:
        summary = summary.sort_values(
            ["polyploid_sample", "phase", top_col, "times_rank1"],
            ascending=[True, True, False, False],
        ).reset_index(drop=True)
    return summary

def _tip_with_phase_letter(tip: str, phase: str) -> str:
    """
    Insert phase letter immediately after the species token.

    Tip format: '{Genus}-{species}_{rest…}'.
    Species = text after the first '-' and before the first '_'.
    Example: Crataegus-necopinnata_Cr710… + A → Crataegus-necopinnataA_Cr710…
    """
    letter = str(phase).strip()
    if not letter:
        return tip
    if "-" not in tip:
        return f"{tip}{letter}"
    prefix, after_dash = tip.split("-", 1)
    if "_" not in after_dash:
        return f"{prefix}-{after_dash}{letter}"
    species, rest = after_dash.split("_", 1)
    return f"{prefix}-{species}{letter}_{rest}"


def _sample_and_contig(tip: str) -> Tuple[str, str]:
    """
    Sample name = everything before the second '_'; contig = the rest.
    Example: Crataegus-necopinnataA_Cr710GrAl_1_1046_c_104.18
      → sample Crataegus-necopinnataA_Cr710GrAl, contig 1_1046_c_104.18
    """
    parts = tip.split("_", 2)
    if len(parts) < 3:
        return tip, ""
    return f"{parts[0]}_{parts[1]}", parts[2]


def _contig_rank_key(contig: str, aligned_seq: str) -> Tuple[int, float, int]:
    """
    Sort key for choosing one contig per sample (higher is better):
      1. contig length from name `{n}_{length}_c_{cov}` (fallback: ungapped bases)
      2. coverage from name
      3. lower SPAdes node number (prefer primary contig)
    """
    m = re.match(r"^(\d+)_(\d+)_c_([0-9.]+)$", contig)
    if m:
        node = int(m.group(1))
        length = int(m.group(2))
        cov = float(m.group(3))
    else:
        node = 10**9
        length = sum(1 for c in aligned_seq if c not in "- \t\n")
        cov = 0.0
    return (length, cov, -node)


def _read_fasta_pairs(path: str) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    name: Optional[str] = None
    chunks: List[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks)))
                name = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line)
        if name is not None:
            records.append((name, "".join(chunks)))
    return records


def _write_phased_alignments(
    out_dir: str, assignments: pd.DataFrame
) -> None:
    """
    Write phased/ alignments with:
      - hybrid tips renamed by appending the phase letter to the species name
      - one sequence per sample (sample = before second '_'; rest = contig)
      - if multiple contigs for a sample, keep the longest (then highest coverage)
    """
    phased_dir = os.path.join(out_dir, "phased")
    if os.path.isdir(phased_dir):
        shutil.rmtree(phased_dir)
    os.makedirs(phased_dir, exist_ok=True)

    tip_to_phase = {
        (row.treefile, row.hybrid_tip): str(row.phase)
        for row in assignments.itertuples(index=False)
    }

    choice_rows: List[dict] = []
    n_written = 0
    n_dropped = 0

    for treefile_name, _group in assignments.groupby("treefile"):
        aln_name = treefile_name.replace(".treefile", "")
        aln_path = os.path.join(out_dir, aln_name)
        if not os.path.isfile(aln_path):
            logger.warning("Alignment missing for %s", treefile_name)
            continue

        # Collect sequences keyed by sample; keep best contig per sample
        best: Dict[str, Tuple[Tuple[int, float, int], str, str, str]] = {}
        # sample -> (rank_key, contig, full_tip_used, sequence)
        order: List[str] = []

        for tip, seq in _read_fasta_pairs(aln_path):
            phase = tip_to_phase.get((treefile_name, tip))
            renamed = _tip_with_phase_letter(tip, phase) if phase else tip
            sample, contig = _sample_and_contig(renamed)
            key = _contig_rank_key(contig, seq)
            if sample not in best:
                order.append(sample)
                best[sample] = (key, contig, renamed, seq)
            else:
                prev_key, prev_contig, prev_tip, _prev_seq = best[sample]
                if key > prev_key:
                    choice_rows.append(
                        {
                            "locus": aln_name.replace(".fasta.mafft", ""),
                            "sample": sample,
                            "kept_contig": contig,
                            "dropped_contig": prev_contig,
                            "kept_tip": renamed,
                            "dropped_tip": prev_tip,
                            "reason": "longer_or_higher_cov",
                        }
                    )
                    best[sample] = (key, contig, renamed, seq)
                    n_dropped += 1
                else:
                    choice_rows.append(
                        {
                            "locus": aln_name.replace(".fasta.mafft", ""),
                            "sample": sample,
                            "kept_contig": prev_contig,
                            "dropped_contig": contig,
                            "kept_tip": prev_tip,
                            "dropped_tip": renamed,
                            "reason": "longer_or_higher_cov",
                        }
                    )
                    n_dropped += 1

        out_path = os.path.join(phased_dir, aln_name)
        with open(out_path, "w") as fout:
            for sample in order:
                _key, contig, _tip, seq = best[sample]
                fout.write(f">{sample}\n{seq}\n")
                n_written += 1

    choices_path = os.path.join(out_dir, "hybrid_contig_choices.tsv")
    pd.DataFrame(choice_rows).to_csv(choices_path, sep="\t", index=False)
    logger.info(
        "Wrote phased alignments under %s (%d sequence(s), dropped %d extra contig(s); "
        "rule=longest then highest coverage). Choices: %s",
        phased_dir,
        n_written,
        n_dropped,
        choices_path,
    )


def _amas_script_path() -> str:
    return os.path.abspath(os.path.join(os.path.dirname(__file__), "AMAS.py"))


def _parse_phased_exon_name(basename: str) -> Optional[Tuple[str, int, str]]:
    """
    4848_exon_3_main.fasta.mafft → (gene_id_for_Assembly, exon_num, copy).
    Main → Assembly_4848; para → Assembly_4848para.
    """
    m = re.match(
        r"^(.+)_exon_(\d+)_(main|para)\.fasta\.mafft$",
        basename,
    )
    if not m:
        return None
    gene, exon_num, copy = m.group(1), int(m.group(2)), m.group(3)
    assembly_id = gene if copy == "main" else f"{gene}para"
    return assembly_id, exon_num, copy


@log_exceptions
def _concat_phased_genes(out_dir: str) -> None:
    """
    AMAS-concatenate phased exon alignments into genes under
    101polyploid/concatenated/ (Assembly_{id}.fasta|.part|.phylip).
    Does not modify 70concatenated_exon_alignments/.
    """
    phased_dir = os.path.join(out_dir, "phased")
    require_dir(phased_dir, "phased alignments")
    concat_dir = os.path.join(out_dir, "concatenated")
    if os.path.isdir(concat_dir):
        shutil.rmtree(concat_dir)
    os.makedirs(concat_dir, exist_ok=True)

    amas = _amas_script_path()
    require_file(amas, "AMAS.py")
    py = sys.executable

    by_gene: Dict[str, List[Tuple[int, str]]] = defaultdict(list)
    for path in sorted(glob(os.path.join(phased_dir, "*.fasta.mafft"))):
        parsed = _parse_phased_exon_name(os.path.basename(path))
        if not parsed:
            logger.warning("Skip unrecognized phased file: %s", path)
            continue
        assembly_id, exon_num, _copy = parsed
        by_gene[assembly_id].append((exon_num, path))

    if not by_gene:
        raise RuntimeError(f"No phased exon alignments to concatenate under {phased_dir}")

    written = 0
    for assembly_id in sorted(by_gene):
        exons = sorted(by_gene[assembly_id], key=lambda x: x[0])
        exon_paths = [os.path.abspath(p) for _n, p in exons]
        fasta_name = f"Assembly_{assembly_id}.fasta"
        part_name = f"Assembly_{assembly_id}.part"
        phy_name = f"Assembly_{assembly_id}.phylip"
        out_fasta = os.path.join(concat_dir, fasta_name)
        out_part = os.path.join(concat_dir, part_name)
        out_phy = os.path.join(concat_dir, phy_name)

        concat_cmd = [
            py,
            amas,
            "concat",
            "-i",
            *exon_paths,
            "-f",
            "fasta",
            "-d",
            "dna",
            "-t",
            fasta_name,
            "-p",
            part_name,
        ]
        logger.info(
            "AMAS concat Assembly_%s (%d exon(s): %s)",
            assembly_id,
            len(exon_paths),
            ",".join(str(n) for n, _ in exons),
        )
        result = subprocess.run(
            concat_cmd,
            capture_output=True,
            text=True,
            check=False,
            cwd=concat_dir,
        )
        if result.returncode != 0 or not os.path.isfile(out_fasta):
            logger.error(
                "AMAS concat failed for Assembly_%s:\n%s\n%s",
                assembly_id,
                result.stdout,
                result.stderr,
            )
            raise RuntimeError(f"AMAS concat failed for Assembly_{assembly_id}")
        if result.stdout.strip():
            logger.debug("AMAS concat: %s", result.stdout.strip()[-500:])

        # Partitions: RAxML-style "DNA, …" and readable exon labels
        with open(out_part) as fh:
            part_lines = fh.readlines()
        cleaned = []
        for i, line in enumerate(part_lines, start=1):
            line = line.rstrip("\n")
            if not line.strip():
                continue
            if " = " in line:
                _name, coords = line.split(" = ", 1)
                exon_num = exons[i - 1][0] if i - 1 < len(exons) else i
                label = f"p{i}_Assembly_{assembly_id}_exon_{exon_num}"
                cleaned.append(f"DNA, {label} = {coords}\n")
            else:
                cleaned.append(
                    line if line.startswith("DNA, ") else f"DNA, {line}\n"
                )
        with open(out_part, "w") as fh:
            fh.writelines(cleaned)

        convert_cmd = [
            py,
            amas,
            "convert",
            "-i",
            fasta_name,
            "-f",
            "fasta",
            "-d",
            "dna",
            "-u",
            "phylip",
        ]
        result = subprocess.run(
            convert_cmd,
            capture_output=True,
            text=True,
            check=False,
            cwd=concat_dir,
        )
        phy_tmp = os.path.join(concat_dir, f"{fasta_name}-out.phy")
        if result.returncode != 0 or not os.path.isfile(phy_tmp):
            logger.error(
                "AMAS convert failed for Assembly_%s:\n%s\n%s",
                assembly_id,
                result.stdout,
                result.stderr,
            )
            raise RuntimeError(f"AMAS convert failed for Assembly_{assembly_id}")
        os.rename(phy_tmp, out_phy)
        written += 1
        logger.info(
            "Wrote Assembly_%s.fasta/.part/.phylip (%d exons)",
            assembly_id,
            len(exons),
        )

    logger.info(
        "AMAS gene concatenation done: %d gene(s) under %s", written, concat_dir
    )


@log_exceptions
def _phase_hybrid_from_trees(
    out_dir: str,
    backbone_samples: Sequence[str],
    polyploid_samples: Sequence[str],
    n_nearest: int = TOP_N_NEAREST,
) -> None:
    """
    Backbone-only nearest ranks → per-polyploid A/B → contig consistency →
    phased FASTAs → AMAS gene concat.
    """
    treefiles = sorted(glob(os.path.join(out_dir, "*.treefile")))
    if not treefiles:
        raise RuntimeError(f"No *.treefile under {out_dir}")

    top_n = int(n_nearest)
    logger.info(
        "Nearest-species phasing for %d polyploid(s) vs %d backbone sample(s) "
        "across %d tree(s) (top %d); other polyploids ignored in ranks",
        len(polyploid_samples),
        len(backbone_samples),
        len(treefiles),
        top_n,
    )
    profiles = _hybrid_distance_profiles(
        treefiles, backbone_samples, polyploid_samples
    )
    feature_cols = _feature_cols(profiles)
    if top_n > len(feature_cols):
        logger.warning(
            "n_nearest=%d is larger than backbone size %d; using %d",
            top_n,
            len(feature_cols),
            len(feature_cols),
        )
        top_n = len(feature_cols)
    nearest_long, rank_matrix = _nearest_species_table(
        profiles, feature_cols, top_n=top_n
    )
    assignments, rank_with_phase = _cluster_rank_matrix_into_ab(
        rank_matrix, feature_cols
    )
    assignments, contig_links = _enforce_shared_contig_phases(assignments)
    phase_lookup = {
        (r.treefile, r.hybrid_tip): r.phase
        for r in assignments.itertuples(index=False)
    }
    rank_with_phase = rank_with_phase.copy()
    rank_with_phase["phase"] = [
        phase_lookup.get((t, h), p)
        for t, h, p in zip(
            rank_with_phase["treefile"],
            rank_with_phase["hybrid_tip"],
            rank_with_phase["phase"],
        )
    ]
    summary = _phase_nearest_summary(nearest_long, assignments, top_n=top_n)

    dist_path = os.path.join(out_dir, "hybrid_distance_matrix.tsv")
    nearest_path = os.path.join(out_dir, "hybrid_nearest_species.tsv")
    rank_path = os.path.join(out_dir, "hybrid_neighbor_rank_matrix.tsv")
    assign_path = os.path.join(out_dir, "hybrid_phase_assignments.tsv")
    summary_path = os.path.join(out_dir, "hybrid_phase_nearest_summary.tsv")
    links_path = os.path.join(out_dir, "hybrid_contig_exon_links.tsv")

    profiles.to_csv(dist_path, sep="\t", index=False)
    nearest_long.to_csv(nearest_path, sep="\t", index=False)
    rank_with_phase.to_csv(rank_path, sep="\t", index=False)
    assignments.to_csv(assign_path, sep="\t", index=False)
    summary.to_csv(summary_path, sep="\t", index=False)
    contig_links.to_csv(links_path, sep="\t", index=False)

    for path in (
        dist_path,
        nearest_path,
        rank_path,
        assign_path,
        summary_path,
        links_path,
    ):
        logger.info("Wrote %s", path)

    for poly in sorted(summary["polyploid_sample"].unique()) if len(summary) else []:
        for phase in ("A", "B"):
            sub = summary[
                (summary["polyploid_sample"] == poly) & (summary["phase"] == phase)
            ].head(5)
            if sub.empty:
                continue
            top_col = "times_in_top%d" % top_n
            peek = ", ".join(
                "%s(r1=%s, top%d=%s)"
                % (r.species, r.times_rank1, top_n, getattr(r, top_col))
                for r in sub.itertuples(index=False)
            )
            logger.info("%s phase %s top backbone: %s", poly, phase, peek)

    if len(contig_links):
        for r in contig_links.itertuples(index=False):
            logger.info(
                "Contig %s (%s, len=%d) spans %s → phase %s%s",
                r.assembled_contig,
                r.polyploid_sample,
                r.contig_length,
                r.exons,
                r.consensus_phase,
                " [conflict resolved]" if r.conflict else "",
            )

    _write_phased_alignments(out_dir, assignments)
    _concat_phased_genes(out_dir)


@log_exceptions
def polyploid(
    data_folder: str,
    num_cores: int,
    log_queue,
    backbone_list: Optional[str] = None,
    polyploid_list: Optional[str] = None,
    n_nearest: int = TOP_N_NEAREST,
) -> None:
    """
    Rebuild 101polyploid/ using backbone_list + polyploid_list:
    export filtered FASTAs, align, tree, phase each polyploid vs backbone only,
    AMAS-concatenate genes.
    """
    require_tools(REQUIRED_TOOLS)
    resolve_iqtree()  # fail early if no IQ-TREE on PATH
    data_folder = os.path.abspath(data_folder)
    require_dir(data_folder, "data folder")

    bb_path = _resolve_list_path(backbone_list, BACKBONE_LIST_DEFAULT, data_folder)
    pl_path = _resolve_list_path(polyploid_list, POLYPLOID_LIST_DEFAULT, data_folder)
    backbone_samples = _read_name_list(bb_path, "backbone_list")
    polyploid_samples = _read_name_list(pl_path, "polyploid_list")

    overlap = sorted(set(backbone_samples) & set(polyploid_samples))
    if overlap:
        raise ValueError(
            "Samples listed in both backbone_list and polyploid_list: "
            + ", ".join(overlap)
        )

    if int(n_nearest) < 1:
        raise ValueError("n_nearest must be an integer >= 1")

    logger.info(
        "Starting cast_polyploid (backbone=%d from %s, polyploid=%d from %s, "
        "n_nearest=%d, cores=%d)",
        len(backbone_samples),
        bb_path,
        len(polyploid_samples),
        pl_path,
        int(n_nearest),
        num_cores,
    )

    tsv_path = os.path.join(data_folder, DETECTED_DIRNAME, PARALOGS_TSV)
    df = _load_paralogs_tsv(tsv_path)
    allowed = list(backbone_samples) + list(polyploid_samples)
    df = _filter_df_to_allowed_samples(df, allowed)

    out_dir = os.path.join(data_folder, POLYPLOID_DIRNAME)
    if os.path.isdir(out_dir):
        shutil.rmtree(out_dir)
        logger.debug("Removed previous %s", out_dir)
    os.makedirs(out_dir, exist_ok=True)

    # record which lists were used
    with open(os.path.join(out_dir, "backbone_list.used.txt"), "w") as fh:
        fh.write("\n".join(backbone_samples) + "\n")
    with open(os.path.join(out_dir, "polyploid_list.used.txt"), "w") as fh:
        fh.write("\n".join(polyploid_samples) + "\n")

    _export_fastas(df, out_dir)
    kept = _filter_fastas_missing_polyploid(out_dir, polyploid_samples)

    n_workers, threads = allocate_workers_and_threads(num_cores, len(kept))
    logger.info(
        "Align+tree: %d FASTA(s), %d worker(s) x %d thread(s)",
        len(kept),
        n_workers,
        threads,
    )

    failures: List[Tuple[str, Exception]] = []
    ok = 0
    with managed_pool(n_workers, log_queue) as pool:
        async_results = [
            (
                path,
                pool.apply_async(_align_and_tree, (path, threads)),
            )
            for path in kept
        ]
        total = len(async_results)
        for done, (path, async_result) in enumerate(async_results, start=1):
            try:
                aln, tree = async_result.get()
                ok += 1
                logger.debug(
                    "Done %d/%d: %s → %s",
                    done,
                    total,
                    os.path.basename(path),
                    os.path.basename(tree),
                )
                if done == total or done % 5 == 0:
                    logger.info("Align+tree progress: %d / %d", done, total)
            except Exception as e:
                logger.error("Align/tree failed for %s: %s", path, e)
                failures.append((path, e))

    if failures:
        preview = ", ".join(os.path.basename(p) for p, _ in failures[:10])
        raise RuntimeError(
            f"cast_polyploid aborted: {len(failures)} job(s) failed ({preview}). See log."
        )

    n_aln = len(glob(os.path.join(out_dir, "*.fasta.mafft")))
    n_tree = len(glob(os.path.join(out_dir, "*.treefile")))
    logger.info(
        "Align+tree done (%d FASTA(s), %d alignment(s), %d tree(s)); phasing…",
        ok,
        n_aln,
        n_tree,
    )
    _phase_hybrid_from_trees(
        out_dir, backbone_samples, polyploid_samples, n_nearest=int(n_nearest)
    )
    logger.info("cast_polyploid completed under %s", out_dir)
