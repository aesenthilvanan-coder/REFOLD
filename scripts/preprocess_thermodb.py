#!/usr/bin/env python3
"""Download and preprocess stability databases (ThermoMutDB + ProTherm)."""

import argparse
import logging
from pathlib import Path

from refold.utils.logging import setup_logging
from refold.data.thermomutdb_loader import build_stability_dataset


def main() -> None:
    parser = argparse.ArgumentParser(description="Preprocess stability databases")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    setup_logging()
    logger = logging.getLogger(__name__)

    logger.info("Building stability training dataset...")
    df = build_stability_dataset(force=args.force)
    logger.info(f"Done. {len(df):,} stability measurements.")
    print(df.describe())


if __name__ == "__main__":
    main()
