#!/usr/bin/env python3
"""
Extract viral candidates (tiers A, B, C) and save as summary table.

DESCRIPTION:
    Filters a length-filtered master contig table to viral confidence tiers
    (A_high-confidence_viral, B_medium-confidence_viral, C_probable_viral)
    and saves the result as a clean summary table.

USAGE:
    python build_viral_table.py \\
        --master_table master_contig_table_1000bp.tsv \\
        --output viral_candidates_ABC_1000bp.tsv

INPUT:
    master_contig_table_*bp.tsv    — Filtered master table from extract_viral_contigs.py
                                     (contains contig_id, sample_id, confidence_tier,
                                      and all tool output columns)

OUTPUT:
    viral_candidates_ABC_*bp.tsv   — Subset to tiers A, B, C only
                                     Same columns as input, just filtered

REQUIREMENTS:
    pandas

AUTHOR:
    Eddie Fuques

"""

import argparse
import sys
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument("--master_table", required=True,
                        help="Path to length-filtered master table "
                             "(from extract_viral_contigs.py)")
    parser.add_argument("--output", required=True,
                        help="Path for output viral candidates TSV")

    args = parser.parse_args()

    master_path = Path(args.master_table)
    output_path = Path(args.output)

    print("=" * 60)
    print("Build viral candidates table")
    print("=" * 60)
    print(f"Input:  {master_path}")
    print(f"Output: {output_path}")
    print("")

    # Validate input
    if not master_path.exists():
        sys.exit(f"ERROR: Master table not found: {master_path}")

    # Load table
    print("Loading master table...")
    master = pd.read_csv(master_path, sep="\t", low_memory=False)
    print(f"  Total contigs: {len(master):,}")

    # Filter to viral tiers
    print("\nFiltering to viral confidence tiers (A, B, C)...")
    viral_tiers = [
        "A_high-confidence_viral",
        "B_medium-confidence_viral",
        "C_probable_viral"
    ]
    viral = master[master["confidence_tier"].isin(viral_tiers)].copy()
    print(f"  Viral candidates: {len(viral):,}")

    # Summary by tier
    print("\n── Tier breakdown ──")
    for tier in viral_tiers:
        n = (viral["confidence_tier"] == tier).sum()
        pct = 100 * n / len(viral) if len(viral) > 0 else 0
        print(f"  {tier:<30} {n:>8,}  ({pct:.1f}%)")

    # Print column summary
    print(f"\nColumns in output: {len(viral.columns)}")
    print(f"  {', '.join(viral.columns[:5])}, ...")

    # Save output
    viral.to_csv(output_path, sep="\t", index=False)
    print(f"\nSaved: {output_path}")
    print("=" * 60)


if __name__ == "__main__":
    main()
