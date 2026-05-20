#!/usr/bin/env python3
"""Download and parse AlphaFold structures for all UniProt IDs in the ClinVar dataset."""

import argparse
import logging
from pathlib import Path

from refold.utils.logging import setup_logging
from refold.constants import PROCESSED_DIR, RAW_DIR


def main() -> None:
    parser = argparse.ArgumentParser(description="Download and cache AlphaFold structures")
    parser.add_argument(
        "--input", type=Path,
        default=PROCESSED_DIR / "clinvar_mutations.parquet",
        help="Path to clinvar mutations parquet",
    )
    parser.add_argument("--max-workers", type=int, default=4)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--limit", type=int, default=None, help="Limit number of structures")
    args = parser.parse_args()

    setup_logging()
    logger = logging.getLogger(__name__)

    if not args.input.exists():
        logger.error(f"Input file not found: {args.input}. Run 'make preprocess-clinvar' first.")
        return

    import pandas as pd
    df = pd.read_parquet(args.input)
    uniprot_ids = df["uniprot_id"].unique().tolist()

    if args.limit:
        uniprot_ids = uniprot_ids[: args.limit]

    logger.info(f"Downloading AlphaFold structures for {len(uniprot_ids)} proteins...")

    from refold.data.alphafold_fetcher import batch_fetch_structures
    results = batch_fetch_structures(uniprot_ids, max_workers=args.max_workers, force=args.force)

    n_ok = sum(1 for v in results.values() if v is not None)
    logger.info(f"Downloaded {n_ok}/{len(uniprot_ids)} structures")

    # Parse and validate
    from refold.structure.parser import parse_pdb_to_structure
    n_parsed = 0
    n_failed = 0

    for uid, pdb_path in results.items():
        if pdb_path is None:
            n_failed += 1
            continue
        try:
            structure = parse_pdb_to_structure(pdb_path, uid)
            n_parsed += 1
            if n_parsed % 100 == 0:
                logger.info(f"Parsed {n_parsed} structures...")
        except Exception as e:
            logger.warning(f"Parse failed for {uid}: {e}")
            n_failed += 1

    logger.info(f"Done. Parsed: {n_parsed}, Failed: {n_failed}")


if __name__ == "__main__":
    main()
