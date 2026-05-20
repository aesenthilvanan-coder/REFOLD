#!/usr/bin/env python3
"""Run REFOLD pipeline across all pathogenic ClinVar variants for a full proteome scan."""

import argparse
import logging
from pathlib import Path

from refold.utils.logging import setup_logging
from refold.constants import PROCESSED_DIR, RESULTS_DIR


def main() -> None:
    parser = argparse.ArgumentParser(description="REFOLD proteome-wide scan")
    parser.add_argument("--stage", type=int, choices=[1, 2, 3], default=1)
    parser.add_argument("--n-conformations", type=int, default=10)
    parser.add_argument("--n-molecules", type=int, default=0)
    parser.add_argument("--device", default=None)
    parser.add_argument("--output-dir", type=Path, default=RESULTS_DIR / "proteome_scan")
    parser.add_argument(
        "--input", type=Path,
        default=PROCESSED_DIR / "clinvar_missense_filtered.parquet",
    )
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--n-workers", type=int, default=1)
    args = parser.parse_args()

    setup_logging()
    logger = logging.getLogger(__name__)

    if not args.input.exists():
        logger.error(f"Input file not found: {args.input}")
        logger.error("Run: make preprocess-clinvar preprocess-structures")
        return

    import pandas as pd
    from refold.types import Mutation

    df = pd.read_parquet(args.input)
    if args.limit:
        df = df.head(args.limit)

    logger.info(f"Running proteome scan on {len(df)} variants (stage {args.stage})")

    mutations = []
    for _, row in df.iterrows():
        try:
            mutations.append(Mutation(
                uniprot_id=str(row["uniprot_id"]),
                position=int(row["position"]),
                wildtype_aa=str(row["wildtype_aa"]),
                mutant_aa=str(row["mutant_aa"]),
                gene_name=str(row.get("gene_name", "")),
                disease=str(row.get("disease", "")),
                clinvar_id=str(row.get("clinvar_id", "")),
                source="clinvar",
            ))
        except Exception as e:
            logger.warning(f"Skipping row: {e}")

    if args.stage >= 2:
        n_conf = args.n_conformations
    else:
        n_conf = 0

    if args.stage >= 3:
        n_mol = args.n_molecules or 100
    else:
        n_mol = 0

    from refold.pipeline.inference.batch_runner import BatchRunner
    runner = BatchRunner(
        n_workers=args.n_workers,
        stage=args.stage,
        n_conformations=n_conf,
        n_molecules=n_mol,
        device_str=args.device,
        output_dir=args.output_dir / "per_mutation",
    )

    results = runner.run(mutations)

    # Write summary
    from refold.pipeline.inference.result_writer import write_results_to_parquet, write_csv_summary
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_results_to_parquet(results, args.output_dir / "results.parquet")
    write_csv_summary(results, args.output_dir / "summary.csv")

    n_amenable = sum(1 for r in results if r.is_rescue_amenable)
    logger.info(
        f"Proteome scan complete. {len(results)} variants processed. "
        f"{n_amenable} rescue-amenable."
    )


if __name__ == "__main__":
    main()
