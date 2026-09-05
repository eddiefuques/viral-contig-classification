#!/usr/bin/env python3
"""
Extract viral contigs by confidence tier and filter by length.

DESCRIPTION:
    Filters the master contig table to a minimum length threshold and extracts
    corresponding contig sequences from FASTA files into per-tier output FASTAs.

WORKFLOW:
    1. Load master_contig_table.tsv (from build_master_table.py)
    2. Filter contigs to >= MIN_LENGTH bp
    3. Save filtered master table as master_contig_table_{min_length}bp.tsv
    4. Build tier lookup: {contig_id -> tier}
    5. Stream through input FASTA files:
       - Classify each sequence by tier
       - Write to per-tier FASTA: viral_tier_A_{min_length}bp.fasta, etc.
       - Write to combined FASTA: viral_ABC_{min_length}bp.fasta
    6. Print summary statistics

USAGE:
    python extract_viral_contigs.py \\
        --master_table master_contig_table.tsv \\
        --contigs_dir /path/to/contigs \\
        --output_dir ./

    # With custom minimum length:
    python extract_viral_contigs.py \\
        --master_table master_contig_table.tsv \\
        --contigs_dir /path/to/contigs \\
        --output_dir ./ \\
        --min_length 2000

INPUT:
    master_contig_table.tsv    — From build_master_table.py (contains contig_id,
                                  sample_id, confidence_tier)
    Contigs FASTA files        — {contigs_dir}/*{fasta_suffix}

OUTPUT:
    master_contig_table_{X}bp.tsv      — Filtered master table (contigs >= X bp)
    viral_tier_A_{X}bp.fasta           — Tier A contigs (high-confidence viral)
    viral_tier_B_{X}bp.fasta           — Tier B contigs (medium-confidence viral)
    viral_tier_C_{X}bp.fasta           — Tier C contigs (probable viral)
    viral_ABC_{X}bp.fasta              — All viral tiers combined

REQUIREMENTS:
    pandas

AUTHOR:
    Eddie Fuques

"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict

import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument("--master_table", required=True,
                        help="Path to master_contig_table.tsv (from build_master_table.py)")
    parser.add_argument("--contigs_dir", required=True,
                        help="Directory with contig FASTA files")
    parser.add_argument("--output_dir", required=True,
                        help="Directory for output FASTA and TSV files")
    parser.add_argument("--min_length", type=int, default=1000,
                        help="Minimum contig length (bp) to extract (default: 1000)")
    parser.add_argument("--fasta_suffix",
                        default="_rnaspades_min500bp_transcripts.fasta",
                        help="FASTA file suffix (default: _rnaspades_min500bp_transcripts.fasta)")

    args = parser.parse_args()

    master_path = Path(args.master_table)
    contigs_dir = Path(args.contigs_dir)
    output_dir = Path(args.output_dir)

    print("=" * 60)
    print("Extract viral contigs by tier and length filter")
    print("=" * 60)
    print(f"Master table:    {master_path}")
    print(f"Contigs dir:     {contigs_dir}")
    print(f"Output dir:      {output_dir}")
    print(f"Min length:      {args.min_length} bp")
    print("")

    # Validate inputs
    if not master_path.exists():
        sys.exit(f"ERROR: Master table not found: {master_path}")
    if not contigs_dir.exists():
        sys.exit(f"ERROR: Contigs directory not found: {contigs_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)

    # Load master table
    print("Loading master table...")
    master = pd.read_csv(master_path, sep="\t", low_memory=False)
    print(f"  Total contigs in table: {len(master):,}")

    # Filter by length
    print(f"\nFiltering by length >= {args.min_length} bp...")
    filtered = master[master["length_bp"] >= args.min_length].copy()
    print(f"  Contigs >= {args.min_length} bp: {len(filtered):,}")

    # Save filtered table
    filtered_table_path = output_dir / f"master_contig_table_{args.min_length}bp.tsv"
    filtered.to_csv(filtered_table_path, sep="\t", index=False)
    print(f"  Saved: {filtered_table_path}")

    # Build tier lookup, keyed by (sample_id, contig_id).
    # Assemblers emit generic contig names (NODE_1_length_..._cov_...) that are
    # identical across samples, so contig_id alone is not a unique key.
    tier_lookup = {(row["sample_id"], row["contig_id"]): row["confidence_tier"]
                   for _, row in filtered.iterrows()}
    tier_groups = filtered["confidence_tier"].value_counts().to_dict()

    print(f"\n  Tier breakdown (>= {args.min_length} bp):")
    for tier in ["A_high-confidence_viral", "B_medium-confidence_viral", "C_probable_viral"]:
        n = tier_groups.get(tier, 0)
        print(f"    {tier:<30} {n:>6,}")

    # Open output FASTA files
    print(f"\nExtracting sequences by tier...")
    tier_files = {}
    for tier in ["A_high-confidence_viral", "B_medium-confidence_viral", "C_probable_viral"]:
        tier_fasta = output_dir / f"viral_tier_{tier.split('_')[0]}_{args.min_length}bp.fasta"
        tier_files[tier] = open(tier_fasta, "w")

    combined_fasta = output_dir / f"viral_ABC_{args.min_length}bp.fasta"
    combined_file = open(combined_fasta, "w")

    ab_fasta = output_dir / f"viral_AB_{args.min_length}bp.fasta"
    ab_file = open(ab_fasta, "w")

    # Track sequences written
    tier_counts = defaultdict(int)

    # Stream through FASTA files and extract sequences
    for fasta_path in sorted(contigs_dir.glob(f"*{args.fasta_suffix}")):
        print(f"  Processing: {fasta_path.name}")
        pool = fasta_path.name.replace(args.fasta_suffix, "")

        with open(fasta_path) as fh:
            current_id = None
            current_seq = []

            for line in fh:
                line = line.rstrip("\n")

                if line.startswith(">"):
                    # Process previous sequence
                    if current_id and (pool, current_id) in tier_lookup:
                        seq = "".join(current_seq)
                        tier = tier_lookup[(pool, current_id)]

                        # Write to tier file. The combined files must only get
                        # contigs that belong to a viral tier -- writing them
                        # outside this guard puts every contig >= min_length
                        # (tier D included) into viral_ABC_*.fasta.
                        tier_fh = tier_files.get(tier)
                        if tier_fh:
                            tier_fh.write(f">{current_id}\n{seq}\n")
                            tier_counts[tier] += 1
                            combined_file.write(f">{current_id}\n{seq}\n")
                            if tier.startswith(("A_", "B_")):
                                ab_file.write(f">{current_id}\n{seq}\n")

                    # Start new sequence
                    current_id = line.lstrip(">").split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)

            # Process final sequence
            if current_id and (pool, current_id) in tier_lookup:
                seq = "".join(current_seq)
                tier = tier_lookup[(pool, current_id)]
                tier_fh = tier_files.get(tier)
                if tier_fh:
                    tier_fh.write(f">{current_id}\n{seq}\n")
                    tier_counts[tier] += 1
                    combined_file.write(f">{current_id}\n{seq}\n")
                    if tier.startswith(("A_", "B_")):
                        ab_file.write(f">{current_id}\n{seq}\n")

    # Close files
    for fh in tier_files.values():
        fh.close()
    combined_file.close()
    ab_file.close()

    # Summary
    print(f"\n" + "=" * 60)
    print("OUTPUT FILES")
    print("=" * 60)
    print(f"Filtered table: {filtered_table_path}")

    for tier in ["A_high-confidence_viral", "B_medium-confidence_viral", "C_probable_viral"]:
        tier_letter = tier.split("_")[0]
        tier_fasta = output_dir / f"viral_tier_{tier_letter}_{args.min_length}bp.fasta"
        n = tier_counts[tier]
        print(f"Tier {tier_letter}:         {tier_fasta.name} ({n:,} contigs)")

    n_ab = tier_counts["A_high-confidence_viral"] + tier_counts["B_medium-confidence_viral"]
    print(f"Combined ABC:   {combined_fasta.name} ({sum(tier_counts.values()):,} contigs)")
    print(f"Combined AB:    {ab_fasta.name} ({n_ab:,} contigs)")
    print("=" * 60)


if __name__ == "__main__":
    main()
