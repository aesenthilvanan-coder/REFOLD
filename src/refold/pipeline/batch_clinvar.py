"""Batch ClinVar scan pipeline — runs REFOLD over all pathogenic missense variants."""

from __future__ import annotations

import argparse
import gc
import logging
import time
from pathlib import Path
from typing import Any

import pandas as pd

from refold.constants import PROCESSED_DIR, RESULTS_DIR

logger = logging.getLogger(__name__)

PROGRESS_FILE = RESULTS_DIR / "clinvar_scan" / "progress.json"


def run_clinvar_scan(
    mutations_parquet: Path | None = None,
    stage: int = 1,
    resume: bool = True,
    n_molecules_per_pocket: int = 200,
    batch_size: int = 100,
    max_mutations: int | None = None,
    device: "torch.device | None" = None,
) -> pd.DataFrame:
    """Run multi-stage ClinVar scan.

    Stage 1: Classify only (ΔΔG + rescue probability) — fast.
    Stage 2: Stage 1 + pocket detection.
    Stage 3: Stage 2 + molecule generation.

    Args:
        mutations_parquet: Path to ClinVar missense mutations parquet.
        stage: Pipeline stage (1, 2, or 3).
        resume: Resume from previous checkpoint.
        n_molecules_per_pocket: Number of molecules per pocket (stage 3).
        batch_size: Mutations per batch for checkpointing.
        max_mutations: Maximum mutations to process (for testing).
        device: Compute device.

    Returns:
        DataFrame with all results.
    """
    from refold.utils.logging import setup_logging
    setup_logging(RESULTS_DIR / "clinvar_scan" / "logs")

    if mutations_parquet is None:
        mutations_parquet = PROCESSED_DIR / "mutations" / "clinvar_missense.parquet"

    if not mutations_parquet.exists():
        raise FileNotFoundError(
            f"ClinVar parquet not found: {mutations_parquet}\n"
            "Run: python scripts/preprocess_clinvar.py"
        )

    output_dir = RESULTS_DIR / "clinvar_scan"
    output_dir.mkdir(parents=True, exist_ok=True)
    results_path = output_dir / f"stage{stage}_results.parquet"

    df = pd.read_parquet(mutations_parquet)
    logger.info(f"Loaded {len(df):,} mutations from {mutations_parquet}")

    if max_mutations:
        df = df.head(max_mutations)

    if resume and PROGRESS_FILE.exists():
        import json
        with open(PROGRESS_FILE) as f:
            progress = json.load(f)
        completed_keys = set(progress.get("completed_keys", []))
        logger.info(f"Resuming from {len(completed_keys)} completed mutations")
    else:
        completed_keys = set()

    from refold.pipeline.single_mutation import REFOLDPipeline
    from refold.types import Mutation

    pipeline = REFOLDPipeline(
        device=device,
        n_conformations=50 if stage >= 2 else 1,
        n_molecules_per_pocket=n_molecules_per_pocket if stage >= 3 else 0,
    )

    all_results = []
    if results_path.exists():
        all_results = pd.read_parquet(results_path).to_dict("records")

    batch = []
    n_processed = 0

    for _, row in df.iterrows():
        key = f"{row['uniprot_id']}_{row['position']}_{row['wildtype_aa']}_{row['mutant_aa']}"
        if key in completed_keys:
            continue

        try:
            mutation = Mutation(
                uniprot_id=row["uniprot_id"],
                position=int(row["position"]),
                wildtype_aa=row["wildtype_aa"],
                mutant_aa=row["mutant_aa"],
                gene_name=row.get("gene_symbol", ""),
                disease=row.get("disease_name", ""),
                clinvar_id=str(row.get("clinvar_variation_id", "")),
            )

            result = pipeline.run(mutation, output_dir=output_dir / "individual")
            batch.append(result.to_dict())
            completed_keys.add(key)
            n_processed += 1

        except Exception as e:
            logger.warning(f"Failed for {key}: {e}")
            batch.append({
                "uniprot_id": row["uniprot_id"],
                "position": int(row["position"]),
                "wildtype_aa": row["wildtype_aa"],
                "mutant_aa": row["mutant_aa"],
                "error": str(e),
            })

        # MPS memory cleanup every 50 mutations
        if n_processed % 50 == 0:
            gc.collect()
            try:
                import torch
                if torch.backends.mps.is_available():
                    torch.mps.empty_cache()
            except Exception:
                pass

        if len(batch) >= batch_size:
            all_results.extend(batch)
            _flush_results(all_results, results_path)
            _save_progress(completed_keys)
            logger.info(f"Processed {n_processed} mutations — saved checkpoint")
            batch = []

    if batch:
        all_results.extend(batch)
        _flush_results(all_results, results_path)
        _save_progress(completed_keys)

    df_results = pd.DataFrame(all_results)
    _print_scan_summary(df_results)
    return df_results


def _flush_results(results: list[dict], path: Path) -> None:
    """Append results to parquet checkpoint."""
    path.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(results)
    df.to_parquet(path, index=False)


def _save_progress(completed_keys: set[str]) -> None:
    import json
    PROGRESS_FILE.parent.mkdir(parents=True, exist_ok=True)
    with open(PROGRESS_FILE, "w") as f:
        json.dump({"completed_keys": list(completed_keys)}, f)


def _print_scan_summary(df: pd.DataFrame) -> None:
    print(f"\nClinVar scan summary ({len(df):,} mutations processed)")
    if "mutation_class" in df.columns:
        print("\nMutation class distribution:")
        print(df["mutation_class"].value_counts().to_string())
    if "rescue_amenability" in df.columns:
        print("\nRescue amenability distribution:")
        print(df["rescue_amenability"].value_counts().to_string())
    if "ddg_predicted" in df.columns:
        ddg = df["ddg_predicted"].dropna()
        print(f"\nΔΔG: mean={ddg.mean():.2f} median={ddg.median():.2f} kcal/mol")
    if "n_passing_lipinski" in df.columns:
        n_with_cands = (df["n_passing_lipinski"] > 0).sum()
        print(f"\nMutations with viable candidates: {n_with_cands:,}")


def main() -> None:
    parser = argparse.ArgumentParser(description="REFOLD ClinVar scan")
    parser.add_argument("--mutations", type=Path, default=None)
    parser.add_argument("--stage", type=int, default=1, choices=[1, 2, 3])
    parser.add_argument("--no-resume", action="store_true")
    parser.add_argument("--n-molecules", type=int, default=200)
    parser.add_argument("--batch-size", type=int, default=100)
    parser.add_argument("--max-mutations", type=int, default=None)
    parser.add_argument("--device", default=None)
    args = parser.parse_args()

    import torch
    device = torch.device(args.device) if args.device else None

    run_clinvar_scan(
        mutations_parquet=args.mutations,
        stage=args.stage,
        resume=not args.no_resume,
        n_molecules_per_pocket=args.n_molecules,
        batch_size=args.batch_size,
        max_mutations=args.max_mutations,
        device=device,
    )


if __name__ == "__main__":
    main()
