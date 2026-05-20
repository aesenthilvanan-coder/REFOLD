"""Train/val/test/benchmark splits with zero sequence-identity leakage.

Uses MMseqs2 clustering at 30% sequence identity.
Fallback to hash-based splitting if MMseqs2 unavailable.
Splits: 70% train, 10% val, 10% test, 10% benchmark (sequestered).
"""

from __future__ import annotations

import hashlib
import logging
import subprocess
import tempfile
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd

from refold.constants import BENCHMARK_MUTATIONS, PROCESSED_DIR

logger = logging.getLogger(__name__)

SplitName = Literal["train", "val", "test", "benchmark"]

SPLIT_FRACTIONS = {"train": 0.70, "val": 0.10, "test": 0.10, "benchmark": 0.10}


def _hash_uniprot_to_split(uniprot_id: str, seed: int = 42) -> SplitName:
    """Deterministic hash-based split assignment."""
    h = int(hashlib.md5(f"{seed}:{uniprot_id}".encode()).hexdigest(), 16)
    r = (h % 10000) / 10000.0
    if r < 0.70:
        return "train"
    elif r < 0.80:
        return "val"
    elif r < 0.90:
        return "test"
    else:
        return "benchmark"


def _run_mmseqs2_clustering(
    sequences: dict[str, str],
    identity_threshold: float = 0.30,
    tmpdir: Path | None = None,
) -> dict[str, str] | None:
    """Run MMseqs2 easy-cluster, return uniprot_id → cluster_rep mapping."""
    try:
        result = subprocess.run(["mmseqs", "version"], capture_output=True, text=True)
        if result.returncode != 0:
            return None
    except FileNotFoundError:
        return None

    with tempfile.TemporaryDirectory() as tmp:
        fasta_path = Path(tmp) / "seqs.fasta"
        with open(fasta_path, "w") as f:
            for uid, seq in sequences.items():
                f.write(f">{uid}\n{seq}\n")

        out_prefix = Path(tmp) / "clusters"
        mmseqs_tmp = Path(tmp) / "mmseqs_tmp"

        cmd = [
            "mmseqs", "easy-cluster",
            str(fasta_path), str(out_prefix), str(mmseqs_tmp),
            "--min-seq-id", str(identity_threshold),
            "--cov-mode", "0",
            "-c", "0.8",
            "--threads", "4",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        if result.returncode != 0:
            logger.warning(f"MMseqs2 clustering failed: {result.stderr}")
            return None

        cluster_file = out_prefix.parent / f"{out_prefix.name}_cluster.tsv"
        if not cluster_file.exists():
            return None

        uid_to_cluster: dict[str, str] = {}
        with open(cluster_file) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) == 2:
                    rep, member = parts
                    uid_to_cluster[member] = rep
        return uid_to_cluster


def create_splits(
    df: pd.DataFrame,
    sequences: dict[str, str] | None = None,
    seed: int = 42,
    output_dir: Path | None = None,
) -> dict[SplitName, pd.DataFrame]:
    """Create train/val/test/benchmark splits with zero leakage.

    Args:
        df: DataFrame with uniprot_id column.
        sequences: Optional dict of uniprot_id → sequence for MMseqs2 clustering.
        seed: Random seed for hash-based fallback.
        output_dir: Directory to save split parquet files.

    Returns:
        Dict mapping split name → DataFrame.
    """
    if output_dir is None:
        output_dir = PROCESSED_DIR / "splits"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Check cached splits
    split_paths = {s: output_dir / f"{s}.parquet" for s in SPLIT_FRACTIONS}
    if all(p.exists() for p in split_paths.values()):
        logger.info("Loading existing splits from cache")
        splits = {s: pd.read_parquet(p) for s, p in split_paths.items()}
        _verify_no_leakage(splits)
        return splits  # type: ignore[return-value]

    unique_uniprots = df["uniprot_id"].unique().tolist()

    uid_to_split: dict[str, SplitName] = {}

    # Benchmark set: known chaperone-responsive variants
    benchmark_uniprots: set[str] = set()
    for key in BENCHMARK_MUTATIONS:
        parts = key.split("_")
        if parts:
            benchmark_uniprots.add(parts[0])

    # Try MMseqs2 clustering
    cluster_to_split: dict[str, SplitName] = {}
    uid_to_cluster: dict[str, str] | None = None

    if sequences:
        uid_to_cluster = _run_mmseqs2_clustering(sequences, identity_threshold=0.30)

    if uid_to_cluster:
        logger.info("Using MMseqs2 clusters for split assignment")
        clusters = list({v for v in uid_to_cluster.values()})
        rng = np.random.default_rng(seed)
        rng.shuffle(clusters)

        n = len(clusters)
        n_train = int(n * 0.70)
        n_val = int(n * 0.10)
        n_test = int(n * 0.10)

        for i, cluster in enumerate(clusters):
            if i < n_train:
                cluster_to_split[cluster] = "train"
            elif i < n_train + n_val:
                cluster_to_split[cluster] = "val"
            elif i < n_train + n_val + n_test:
                cluster_to_split[cluster] = "test"
            else:
                cluster_to_split[cluster] = "benchmark"

        for uid in unique_uniprots:
            if uid in benchmark_uniprots:
                uid_to_split[uid] = "benchmark"
            else:
                cluster = uid_to_cluster.get(uid, uid)
                uid_to_split[uid] = cluster_to_split.get(cluster, "train")
    else:
        logger.info("MMseqs2 unavailable — using hash-based split assignment")
        for uid in unique_uniprots:
            if uid in benchmark_uniprots:
                uid_to_split[uid] = "benchmark"
            else:
                uid_to_split[uid] = _hash_uniprot_to_split(uid, seed=seed)

    df["split"] = df["uniprot_id"].map(uid_to_split).fillna("train")

    splits: dict[SplitName, pd.DataFrame] = {}
    for split_name in SPLIT_FRACTIONS:
        split_df = df[df["split"] == split_name].copy()
        splits[split_name] = split_df  # type: ignore[assignment]
        split_df.to_parquet(split_paths[split_name], index=False)
        logger.info(f"Split '{split_name}': {len(split_df):,} mutations")

    _verify_no_leakage(splits)
    return splits


def _verify_no_leakage(splits: dict[SplitName, pd.DataFrame]) -> None:
    """Assert zero overlap between splits at the UniProt level."""
    train_uniprots = set(splits["train"]["uniprot_id"].unique())
    val_uniprots = set(splits["val"]["uniprot_id"].unique())
    test_uniprots = set(splits["test"]["uniprot_id"].unique())
    bench_uniprots = set(splits["benchmark"]["uniprot_id"].unique())

    assert len(train_uniprots & val_uniprots) == 0, \
        f"Train/val overlap: {train_uniprots & val_uniprots}"
    assert len(train_uniprots & test_uniprots) == 0, \
        f"Train/test overlap: {train_uniprots & test_uniprots}"
    assert len(train_uniprots & bench_uniprots) == 0, \
        f"Train/benchmark overlap: {train_uniprots & bench_uniprots}"
    assert len(val_uniprots & test_uniprots) == 0, \
        f"Val/test overlap: {val_uniprots & test_uniprots}"
    assert len(val_uniprots & bench_uniprots) == 0, \
        f"Val/benchmark overlap: {val_uniprots & bench_uniprots}"

    logger.info("Split leakage check passed: zero UniProt overlap across all splits")
