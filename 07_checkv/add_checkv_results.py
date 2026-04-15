#!/usr/bin/env python3
"""
Merge CheckV output into the viral candidate table.

DESCRIPTION:
    Reads viral_candidates_ABC_*bp.tsv and the two standard CheckV output
    files (quality_summary.tsv and contamination.tsv) and produces an
    annotated table with CheckV quality metrics.

    CheckV output files expected at:
      {CHECKV_DIR}/quality_summary.tsv
      {CHECKV_DIR}/contamination.tsv   (optional — only for proviruses)

KEY COLUMNS ADDED FROM CHECKV:
    checkv_quality              — Complete / High-quality / Medium-quality /
                                  Low-quality / Not-determined
    checkv_completeness         — Estimated genome completeness (%)
    checkv_completeness_method  — Method used (AAI-based / HMM-based / etc.)
    checkv_contamination        — Host contamination estimate (%)
    checkv_provirus             — True/False — is it an integrated provirus?
    checkv_warnings             — Any CheckV warning flags
    checkv_viral_genes          — Number of viral protein matches
    checkv_host_genes           — Number of host protein matches
    checkv_gene_count           — Total ORF count
    checkv_region_types         — Provirus region classification (if available)
    checkv_proviral_length      — Predicted viral sequence length (if provirus)
    checkv_host_length          — Predicted host contamination length (if provirus)

USAGE:
    python add_checkv_results.py \\
        --input viral_candidates_ABC_1kb.tsv \\
        --checkv checkv_output \\
        --output viral_candidates_ABC_1kb_checkv.tsv

    # With defaults (uses current directory):
    python add_checkv_results.py

INPUT:
    viral_candidates_ABC_*bp.tsv  — From build_viral_table.py
    {CHECKV_DIR}/quality_summary.tsv
    {CHECKV_DIR}/contamination.tsv (optional)

OUTPUT:
    viral_candidates_ABC_*bp_checkv.tsv  — Master table with CheckV results merged

REQUIREMENTS:
    pandas

AUTHOR:
    Viral Community Analysis Pipeline 2

"""

import argparse
import sys
from pathlib import Path

import pandas as pd


# CheckV quality tier order (for sorting / filtering later)
CHECKV_QUALITY_ORDER = [
    "Complete",
    "High-quality",
    "Medium-quality",
    "Low-quality",
    "Not-determined",
]


# ─────────────────────────────────────────────────────────────────────────────
# HELPERS
# ─────────────────────────────────────────────────────────────────────────────

def load_quality_summary(checkv_dir: Path) -> pd.DataFrame:
    """
    Parse CheckV quality_summary.tsv.

    Standard columns (CheckV >= 1.0):
      contig_id, contig_length, genome_copies, gene_count,
      viral_genes, host_genes, checkv_quality, miuvig_quality,
      completeness, completeness_method, contamination, warnings
    """
    f = checkv_dir / "quality_summary.tsv"
    if not f.exists():
        sys.exit(f"ERROR: CheckV quality_summary.tsv not found at {f}")

    df = pd.read_csv(f, sep="\t")

    # Rename to our namespace
    rename = {
        "contig_id":           "contig_id",          # key for merge
        "checkv_quality":      "checkv_quality",
        "completeness":        "checkv_completeness",
        "completeness_method": "checkv_completeness_method",
        "contamination":       "checkv_contamination",
        "warnings":            "checkv_warnings",
        "viral_genes":         "checkv_viral_genes",
        "host_genes":          "checkv_host_genes",
        "gene_count":          "checkv_gene_count",
        "provirus":            "checkv_provirus",
    }
    # Only keep columns that exist in this CheckV version
    present = {k: v for k, v in rename.items() if k in df.columns}
    df = df.rename(columns=present)
    keep = list(present.values())
    df = df[[c for c in keep if c in df.columns]].drop_duplicates("contig_id")

    print(f"  CheckV quality_summary: {len(df):,} entries")
    if "checkv_quality" in df.columns:
        print(f"\n  CheckV quality breakdown:")
        for q in CHECKV_QUALITY_ORDER:
            n = (df["checkv_quality"] == q).sum()
            if n:
                print(f"    {q:<20} {n:>6,}")

    return df


def load_contamination(checkv_dir: Path) -> pd.DataFrame:
    """
    Parse CheckV contamination.tsv (only present for provirus contigs).
    Adds provirus-specific trim coordinates if available.
    """
    f = checkv_dir / "contamination.tsv"
    if not f.exists():
        print("  contamination.tsv not found — no provirus trim info added")
        return pd.DataFrame(columns=["contig_id"])

    df = pd.read_csv(f, sep="\t")
    rename = {
        "contig_id":        "contig_id",
        "region_types":     "checkv_region_types",
        "proviral_length":  "checkv_proviral_length",
        "host_length":      "checkv_host_length",
    }
    present = {k: v for k, v in rename.items() if k in df.columns}
    df = df.rename(columns=present)
    keep = list(present.values())
    df = df[[c for c in keep if c in df.columns]].drop_duplicates("contig_id")
    print(f"  CheckV contamination.tsv: {len(df):,} provirus entries")
    return df


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input",   default="viral_candidates_ABC_1kb.tsv",
                        help="Path to viral_candidates_ABC_*.tsv")
    parser.add_argument("--checkv",  default="checkv_output",
                        help="Path to CheckV output directory")
    parser.add_argument("--output",  default="viral_candidates_ABC_1kb_checkv.tsv",
                        help="Path for output TSV")
    args = parser.parse_args()

    in_path  = Path(args.input)
    cv_dir   = Path(args.checkv)
    out_path = Path(args.output)

    print("=" * 60)
    print("Merging CheckV results into viral candidate table")
    print("=" * 60)

    # ── Load candidate table ──────────────────────────────────────────────────
    if not in_path.exists():
        sys.exit(f"ERROR: Input table not found: {in_path}\n"
                 f"       Run build_viral_table.py first.")
    print(f"\nReading candidate table: {in_path}")
    viral = pd.read_csv(in_path, sep="\t", low_memory=False)
    print(f"  {len(viral):,} viral candidate contigs")

    # ── Load CheckV outputs ───────────────────────────────────────────────────
    print(f"\nReading CheckV output from: {cv_dir}")
    quality = load_quality_summary(cv_dir)
    contam  = load_contamination(cv_dir)

    # ── Merge ─────────────────────────────────────────────────────────────────
    print("\nMerging...")
    merged = viral.merge(quality, on="contig_id", how="left")
    if len(contam.columns) > 1:
        merged = merged.merge(contam, on="contig_id", how="left")

    n_matched = merged["checkv_quality"].notna().sum()
    n_missing = merged["checkv_quality"].isna().sum()
    print(f"  Contigs with CheckV result: {n_matched:,}")
    if n_missing:
        print(f"  WARNING: {n_missing:,} contigs not found in CheckV output")
        print(f"           (may be below CheckV's minimum length or had no ORFs)")

    # ── Summary by tier + CheckV quality ─────────────────────────────────────
    print("\n── CheckV quality by confidence tier ──")
    tier_order = [
        "A_high-confidence_viral",
        "B_medium-confidence_viral",
        "C_probable_viral",
    ]
    for tier in tier_order:
        sub = merged[merged["confidence_tier"] == tier]
        if sub.empty:
            continue
        print(f"\n  {tier}  (n={len(sub):,}):")
        if "checkv_quality" in sub.columns:
            vc = sub["checkv_quality"].value_counts()
            for q in CHECKV_QUALITY_ORDER:
                n = vc.get(q, 0)
                pct = 100 * n / len(sub)
                if n:
                    print(f"    {q:<22} {n:>5,}  ({pct:.1f}%)")
            nd = vc.get("Not-determined", 0)
            if nd:
                print(f"    {'Not-determined':<22} {nd:>5,}  ({100*nd/len(sub):.1f}%)")

    # ── Completeness summary ──────────────────────────────────────────────────
    if "checkv_completeness" in merged.columns:
        comp = pd.to_numeric(merged["checkv_completeness"], errors="coerce")
        print(f"\n── Completeness (non-null, n={comp.notna().sum():,}) ──")
        print(f"  Mean:   {comp.mean():.1f}%")
        print(f"  Median: {comp.median():.1f}%")
        bins = [0, 50, 75, 90, 100]
        labels = ["<50%", "50–75%", "75–90%", ">90%"]
        cut = pd.cut(comp.dropna(), bins=bins, labels=labels, right=True)
        for lbl, n in cut.value_counts().sort_index().items():
            print(f"  {lbl:<10} {n:>6,}")

    # ── Save ─────────────────────────────────────────────────────────────────
    merged.to_csv(out_path, sep="\t", index=False)
    print(f"\nSaved: {out_path}")
    print(f"Columns in output: {len(merged.columns)}")
    print("=" * 60)


if __name__ == "__main__":
    main()
