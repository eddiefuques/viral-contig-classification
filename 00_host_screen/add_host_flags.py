#!/usr/bin/env python3
"""
Merge the host screen into the master contig table.

DESCRIPTION:
    Attaches the per-contig host alignment evidence produced by
    run_host_screen.sh to the master table from build_master_table.py, so that a
    viral call which also maps to the host genome is visible rather than
    silently trusted.

    Nothing is removed. Endogenous viral elements -- endogenous retroviruses in
    particular -- are viral by sequence, so no classifier can separate them from
    a circulating virus. Flagging them lets you make that judgement yourself,
    with the evidence in the same table as the tier.

INPUT:
    --master_table   master_contig_table.tsv from build_master_table.py
    --host_dir       directory of {sample}_host_screen.tsv files

ADDED COLUMNS:
    host_best_frac   fraction of the contig covered by its best host alignment
    host_identity    identity of that alignment (nmatch / alignment length).
                     A long alignment at low identity suggests an ancient
                     endogenous element rather than recent host transcript.
    host_aln_bp      aligned length in bp
    host_mapq        mapping quality of the best alignment
    host_n_hits      number of host alignments for this contig
    host_target      host sequence hit
    host_flag        host | host_partial | none

OUTPUT:
    <master_table>_hostflagged.tsv (or --output)
    host_flagged_viral_candidates.tsv — viral-tier calls that map to the host,
                                        written next to the output table

USAGE:
    python add_host_flags.py \\
        --master_table master_contig_table.tsv \\
        --host_dir     host_screen \\
        --output       master_contig_table_hostflagged.tsv

    Run this between build_master_table.py and extract_viral_contigs.py, and
    pass the flagged table to extract_viral_contigs.py so the host columns
    propagate into every candidate table.

REQUIREMENTS:
    pandas

AUTHOR:
    Eddie Fuques
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

# Tool results are keyed per sample: assemblers emit generic contig names that
# repeat across assemblies, so contig_id alone is not unique.
KEYS = ["sample_id", "contig_id"]

HOST_COLS = ["host_best_frac", "host_identity", "host_aln_bp",
             "host_mapq", "host_n_hits", "host_target"]


def load_host_screen(host_dir: Path, samples=None) -> pd.DataFrame:
    """Read every {sample}_host_screen.tsv and tag its rows with the sample."""
    files = sorted(host_dir.glob("*_host_screen.tsv"))
    if not files:
        sys.exit(f"ERROR: no *_host_screen.tsv files in {host_dir}\n"
                 f"       Run 00_host_screen/run_host_screen.sh first.")

    frames = []
    for f in files:
        sample = f.name.replace("_host_screen.tsv", "")
        if samples and sample not in samples:
            continue
        df = pd.read_csv(f, sep="\t")
        if df.empty:
            print(f"  {sample}: no host alignments")
            continue
        df["sample_id"] = sample
        frames.append(df)
        print(f"  {sample}: {len(df):,} contigs with at least one host alignment")

    if not frames:
        return pd.DataFrame(columns=KEYS + HOST_COLS)

    host = pd.concat(frames, ignore_index=True)
    # the screen carries its own length column; drop it to avoid a clash
    return host.drop(columns=[c for c in ["length_bp"] if c in host.columns])


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--master_table", required=True,
                   help="Master table from build_master_table.py")
    p.add_argument("--host_dir", required=True,
                   help="Directory with {sample}_host_screen.tsv files")
    p.add_argument("--output", default=None,
                   help="Output TSV (default: <master_table>_hostflagged.tsv)")
    p.add_argument("--frac_full", type=float, default=0.5,
                   help="host_best_frac >= this -> host (default: 0.5)")
    p.add_argument("--frac_partial", type=float, default=0.2,
                   help="host_best_frac >= this -> host_partial (default: 0.2)")
    args = p.parse_args()

    master_path = Path(args.master_table)
    if not master_path.exists():
        sys.exit(f"ERROR: master table not found: {master_path}")
    out_path = (Path(args.output) if args.output
                else master_path.with_name(master_path.stem + "_hostflagged.tsv"))

    print("=" * 60)
    print("Adding host flags to the master table")
    print("=" * 60)

    master = pd.read_csv(master_path, sep="\t", low_memory=False)
    for col in KEYS:
        if col not in master.columns:
            sys.exit(f"ERROR: master table has no '{col}' column")
    print(f"Master table: {len(master):,} contigs, "
          f"{master['sample_id'].nunique()} samples")

    print("\nReading host screen:")
    host = load_host_screen(Path(args.host_dir), set(master["sample_id"].unique()))

    dup = host.duplicated(KEYS).sum()
    if dup:
        sys.exit(f"ERROR: {dup:,} duplicate (sample_id, contig_id) rows in the host screen")

    merged = master.merge(host, on=KEYS, how="left")
    if len(merged) != len(master):
        sys.exit(f"ERROR: merge changed the row count "
                 f"({len(master):,} -> {len(merged):,})")

    merged["host_best_frac"] = pd.to_numeric(
        merged.get("host_best_frac"), errors="coerce").fillna(0.0)
    merged["host_n_hits"] = pd.to_numeric(
        merged.get("host_n_hits"), errors="coerce").fillna(0).astype(int)

    merged["host_flag"] = "none"
    merged.loc[merged["host_best_frac"] >= args.frac_partial, "host_flag"] = "host_partial"
    merged.loc[merged["host_best_frac"] >= args.frac_full, "host_flag"] = "host"

    merged.to_csv(out_path, sep="\t", index=False)

    # ── Report ──────────────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("HOST FLAGS")
    print("=" * 60)
    for flag in ["host", "host_partial", "none"]:
        n = int((merged["host_flag"] == flag).sum())
        print(f"  {flag:<14} {n:>8,}  ({100 * n / len(merged):.2f}%)")

    if "confidence_tier" in merged.columns:
        print("\n── Host flag x confidence tier ──")
        ct = pd.crosstab(merged["confidence_tier"], merged["host_flag"])
        for col in ["host", "host_partial", "none"]:
            if col not in ct.columns:
                ct[col] = 0
        print(ct[["host", "host_partial", "none"]].to_string())

        viral = merged[
            merged["confidence_tier"].astype(str).str.startswith(("A_", "B_", "C_")) &
            (merged["host_flag"] != "none")
        ].copy()

        sus_path = out_path.with_name("host_flagged_viral_candidates.tsv")
        cols = [c for c in ["contig_id", "sample_id", "length_bp", "confidence_tier",
                            "host_flag", "host_best_frac", "host_identity",
                            "host_target", "vs2_max_score", "vs2_group",
                            "genomad_virus_score", "genomad_taxonomy",
                            "deep6_top_class", "viralm_score"]
                if c in viral.columns]
        viral.sort_values(["confidence_tier", "host_best_frac"],
                          ascending=[True, False])[cols].to_csv(
            sus_path, sep="\t", index=False)

        print(f"\n  {len(viral):,} viral-tier contigs also map to the host genome")
        print(f"  -> {sus_path}")
        if len(viral):
            print("\n  Highest host coverage among them:")
            for _, r in viral.sort_values("host_best_frac", ascending=False).head(10).iterrows():
                tax = str(r.get("genomad_taxonomy", ""))[:40]
                print(f"    {str(r['contig_id'])[:46]:<46} "
                      f"{str(r['confidence_tier'])[:1]}  "
                      f"host={r['host_best_frac']:.2f}  {tax}")
            print("\n  Retroviral or retro-like calls in this list are very likely")
            print("  endogenous elements of the host genome, not circulating viruses.")

    print(f"\nSaved: {out_path}")
    print("=" * 60)


if __name__ == "__main__":
    main()
