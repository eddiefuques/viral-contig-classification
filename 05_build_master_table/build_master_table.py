#!/usr/bin/env python3
"""
Build a comprehensive master contig classification table from multi-tool outputs.

DESCRIPTION:
    Integrates results from four viral classification tools (VirSorter2, geNomad,
    Deep6, and ViraLM) into a single master table with consensus scoring,
    confidence tiers, and per-tool agreement metrics.

    Confidence tiers (assigned by tool agreement and score thresholds):
      A_high-confidence_viral       — Consensus viral calls (2+ tools agree)
      B_medium-confidence_viral     — One backbone + ≥ 1 support tool
      C_probable_viral              — Partial evidence (Deep6 + ViraLM only)
      D_unknown                     — Ambiguous or no viral signal

TOOLS & INPUTS:
    VirSorter2      — {vs2_dir}/{sample}_virsorter2_out/final-viral-score.tsv
    geNomad         — {genomad_dir}/{sample}_genomad_out/{sample}_summary/*_virus_summary.tsv
    Deep6           — {deep6_dir}/{sample}_predict_deep6.txt
    ViraLM          — {viralm_dir}/{sample}_viralm_out/result_{sample}.csv

USAGE:
    python build_master_table.py \\
        --contigs_dir /path/to/contigs \\
        --vs2_dir /path/to/virsorter2_outputs \\
        --genomad_dir /path/to/genomad_outputs \\
        --deep6_dir /path/to/deep6_outputs \\
        --viralm_dir /path/to/viralm_outputs \\
        --output master_contig_table.tsv

OPTIONAL ARGUMENTS:
    --fasta_suffix              File suffix (default: _rnaspades_min500bp_transcripts.fasta)
    --vs2_min                   VirSorter2 score threshold (default: 0.5)
    --vs2_high                  VirSorter2 high-confidence threshold (default: 0.9)
    --genomad_min               geNomad virus score threshold (default: 0.7)
    --genomad_plasmid           geNomad plasmid score threshold (default: 0.7)
    --deep6_factor              Deep6 top-score multiplier for viral call (default: 1.25)
    --deep6_min_score           Deep6 minimum absolute score (default: 0.70)
    --viralm_min                ViraLM score threshold (default: 0.70)
    --length_short              Minimum length for tier B support (default: 2500 bp)
    --require_all_tools         Filter to samples with all 4 tools completed (default: True)

OUTPUT:
    Master TSV with columns:
      contig_id, sample_id, length_bp,
      vs2_max_score, vs2_group, vs2_hallmark, vs2_viral_frac, vs2_cellular_frac, vs2_called,
      genomad_virus_score, genomad_fdr, genomad_n_hallmarks, genomad_marker_enrichment,
        genomad_taxonomy, genomad_n_genes, genomad_topology, genomad_called,
      genomad_plasmid_score, genomad_plasmid_flag, genomad_conjugation_genes, genomad_amr_genes,
      deep6_top_class, deep6_top_score, deep6_is_viral,
      viralm_score, viralm_called,
      confidence_tier

REQUIREMENTS:
    pandas, numpy

AUTHOR:
    Eddie Fuques
    
"""

import argparse
import sys
from pathlib import Path
from itertools import combinations

import numpy as np
import pandas as pd

# Merge key. Assemblers emit generic contig names (NODE_1_length_..._cov_...)
# that repeat across samples, so contig_id alone is NOT unique in a multi-sample
# run -- merging on it assigns one sample's scores to another sample's contigs.
KEYS = ["sample_id", "contig_id"]


# ─────────────────────────────────────────────────────────────────────────────
# PARSING FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────

def get_pools(contigs_dir, fasta_suffix):
    """
    Identify sample pool names from FASTA files in contigs_dir.
    """
    pattern = str(Path(contigs_dir) / f"*{fasta_suffix}")
    import glob
    fastas = sorted(glob.glob(pattern))
    if not fastas:
        sys.exit(f"ERROR: No FASTA files found matching *{fasta_suffix} in {contigs_dir}")
    pools = [Path(f).name.replace(fasta_suffix, "") for f in fastas]
    print(f"  Found {len(pools)} pools")
    return pools


def parse_contigs(pools, contigs_dir, fasta_suffix):
    """
    Extract contig IDs and lengths from FASTA files.
    """
    records = []
    for pool in pools:
        fasta = Path(contigs_dir) / f"{pool}{fasta_suffix}"
        if not fasta.exists():
            continue
        with open(fasta) as fh:
            for line in fh:
                if not line.startswith(">"):
                    continue
                contig_id = line.strip().lstrip(">").split()[0]
                try:
                    length = int(contig_id.split("_length_")[1].split("_")[0])
                except (IndexError, ValueError):
                    length = None
                records.append({
                    "contig_id": contig_id,
                    "sample_id": pool,
                    "length_bp": length
                })
    df = pd.DataFrame(records)
    print(f"  Total contigs: {len(df):,}")
    return df


def parse_virsorter2(pools, vs2_dir, vs2_min):
    """
    Parse VirSorter2 final-viral-score.tsv outputs.
    """
    dfs = []
    missing = []

    for pool in pools:
        f = Path(vs2_dir) / f"{pool}_virsorter2_out" / "final-viral-score.tsv"
        if not f.exists():
            missing.append(pool)
            continue
        df = pd.read_csv(f, sep="\t")
        df["contig_id"] = df["seqname"].str.split("||", regex=False).str[0]
        df["sample_id"] = pool
        dfs.append(df)

    if missing:
        print(f"  WARNING: VS2 missing for {len(missing)} pools")

    if not dfs:
        return pd.DataFrame(columns=["contig_id"])

    vs2 = pd.concat(dfs, ignore_index=True)
    vs2 = vs2.rename(columns={
        "max_score": "vs2_max_score",
        "max_score_group": "vs2_group",
        "hallmark": "vs2_hallmark",
        "viral": "vs2_viral_frac",
        "cellular": "vs2_cellular_frac"
    })

    # Call viral if score >= threshold OR has group assignment
    vs2["vs2_called"] = (
        (vs2["vs2_max_score"] >= vs2_min) |
        (vs2["vs2_max_score"].isna() & vs2["vs2_group"].notna())
    )

    keep = [
        "sample_id", "contig_id", "vs2_max_score", "vs2_group", "vs2_hallmark",
        "vs2_viral_frac", "vs2_cellular_frac", "vs2_called"
    ]
    vs2 = vs2[[c for c in keep if c in vs2.columns]].drop_duplicates(KEYS)

    print(f"  VS2: {len(vs2):,} entries | {vs2['vs2_called'].sum():,} viral")
    return vs2


def parse_genomad(pools, genomad_dir, fasta_suffix, genomad_virus_thresh, genomad_plasmid_thresh):
    """
    Parse geNomad virus_summary.tsv and plasmid_summary.tsv outputs.
    """
    v_dfs = []
    p_dfs = []
    missing = []

    for pool in pools:
        summary_dir = Path(genomad_dir) / f"{pool}_genomad_out" / \
                      f"{pool}{fasta_suffix.replace('.fasta', '')}_summary"
        prefix = f"{pool}{fasta_suffix.replace('.fasta', '')}"
        vf = summary_dir / f"{prefix}_virus_summary.tsv"
        pf = summary_dir / f"{prefix}_plasmid_summary.tsv"

        if not vf.exists() and not pf.exists():
            missing.append(pool)
            continue

        if vf.exists():
            _v = pd.read_csv(vf, sep="\t")
            _v["sample_id"] = pool
            v_dfs.append(_v)
        if pf.exists():
            _p = pd.read_csv(pf, sep="\t")
            _p["sample_id"] = pool
            p_dfs.append(_p)

    if missing:
        print(f"  WARNING: geNomad missing for {len(missing)} pools")

    result = pd.DataFrame()

    if v_dfs:
        virus = pd.concat(v_dfs, ignore_index=True)
        id_col = "seq_name" if "seq_name" in virus.columns else virus.columns[0]
        virus = virus.rename(columns={
            id_col: "contig_id",
            "virus_score": "genomad_virus_score",
            "fdr": "genomad_fdr",
            "n_hallmarks": "genomad_n_hallmarks",
            "marker_enrichment": "genomad_marker_enrichment",
            "taxonomy": "genomad_taxonomy",
            "n_genes": "genomad_n_genes",
            "topology": "genomad_topology"
        })

        virus["genomad_called"] = virus["genomad_virus_score"] >= genomad_virus_thresh

        keep_v = [
            "sample_id", "contig_id", "genomad_virus_score", "genomad_fdr", "genomad_n_hallmarks",
            "genomad_marker_enrichment", "genomad_taxonomy", "genomad_n_genes",
            "genomad_topology", "genomad_called"
        ]
        result = virus[[c for c in keep_v if c in virus.columns]].drop_duplicates(KEYS)
        print(f"  geNomad virus: {len(result):,} | {result['genomad_called'].sum():,} viral")

    if p_dfs:
        plasmid = pd.concat(p_dfs, ignore_index=True)
        id_col = "seq_name" if "seq_name" in plasmid.columns else plasmid.columns[0]
        plasmid = plasmid.rename(columns={
            id_col: "contig_id",
            "plasmid_score": "genomad_plasmid_score"
        })

        plasmid["genomad_plasmid_flag"] = plasmid["genomad_plasmid_score"] >= genomad_plasmid_thresh

        for col in ["conjugation_genes", "amr_genes"]:
            if col in plasmid.columns:
                plasmid = plasmid.rename(columns={col: f"genomad_{col}"})

        keep_p = [
            "sample_id", "contig_id", "genomad_plasmid_score", "genomad_plasmid_flag",
            "genomad_conjugation_genes", "genomad_amr_genes"
        ]
        plasmid = plasmid[[c for c in keep_p if c in plasmid.columns]].drop_duplicates(KEYS)
        print(f"  geNomad plasmid: {len(plasmid):,} | {plasmid['genomad_plasmid_flag'].sum():,} flagged")

        result = result.merge(plasmid, on=KEYS, how="outer") if not result.empty else plasmid

    return result


def parse_deep6(pools, deep6_dir, deep6_factor, deep6_min_score):
    """
    Parse Deep6 prediction files and compute viral calls.
    """
    SCORE_COLS = ["duplo", "euk", "mono", "pro", "ribo", "vari"]
    VIRAL_COLS = {"duplo", "mono", "ribo", "vari"}

    dfs = []
    missing = []

    for pool in pools:
        f = Path(deep6_dir) / f"{pool}_predict_deep6.txt"
        if not f.exists():
            missing.append(pool)
            continue
        _d = pd.read_csv(f, sep="\t")
        _d["sample_id"] = pool
        dfs.append(_d)

    if missing:
        print(f"  WARNING: Deep6 missing for {len(missing)} pools")

    if not dfs:
        return pd.DataFrame(columns=["contig_id"])

    deep6 = pd.concat(dfs, ignore_index=True)

    # Extract scores and compute top classification
    score_data = deep6[SCORE_COLS].astype(float)
    medians = score_data.median(axis=1)
    top_idx = score_data.values.argmax(axis=1)
    top_scores = score_data.values[range(len(score_data)), top_idx]
    top_groups = [SCORE_COLS[i] for i in top_idx]

    is_viral = (
        pd.Series([g in VIRAL_COLS for g in top_groups]) &
        (pd.Series(top_scores) >= medians.values * deep6_factor) &
        (pd.Series(top_scores) >= deep6_min_score)
    )

    result = pd.DataFrame({
        "sample_id": deep6["sample_id"].values,
        "contig_id": deep6["name"].values,
        "deep6_top_class": top_groups,
        "deep6_top_score": top_scores,
        "deep6_is_viral": is_viral.values
    }).drop_duplicates(KEYS)

    print(f"  Deep6: {len(result):,} entries | {result['deep6_is_viral'].sum():,} viral")
    return result


def parse_viralm(pools, viralm_dir, fasta_suffix, viralm_min):
    """
    Parse ViraLM result CSV files.
    """
    dfs = []
    missing = []

    for pool in pools:
        fname = f"result_{pool}{fasta_suffix.replace('.fasta', '.csv')}"
        f = Path(viralm_dir) / f"{pool}_viralm_out" / fname
        if not f.exists():
            missing.append(pool)
            continue
        _v = pd.read_csv(f)
        _v["sample_id"] = pool
        dfs.append(_v)

    if missing:
        print(f"  WARNING: ViraLM missing for {len(missing)} pools")

    if not dfs:
        return pd.DataFrame(columns=["contig_id"])

    viralm = pd.concat(dfs, ignore_index=True)
    viralm = viralm.rename(columns={
        "seq_name": "contig_id",
        "virus_score": "viralm_score"
    })

    viralm["viralm_called"] = viralm["viralm_score"] >= viralm_min

    keep = ["sample_id", "contig_id", "viralm_score", "viralm_called"]
    viralm = viralm[[c for c in keep if c in viralm.columns]].drop_duplicates(KEYS)

    print(f"  ViraLM: {len(viralm):,} entries | {viralm['viralm_called'].sum():,} viral")
    return viralm


# ─────────────────────────────────────────────────────────────────────────────
# TIER ASSIGNMENT
# ─────────────────────────────────────────────────────────────────────────────

def assign_tiers(master, vs2_high, length_short):
    """
    Assign confidence tiers based on tool agreement and score thresholds.
    """
    def tier(row):
        vs2 = bool(row.get("vs2_called", False))
        gmd = bool(row.get("genomad_called", False))
        d6 = bool(row.get("deep6_is_viral", False))
        vlm = bool(row.get("viralm_called", False))
        plasmid = bool(row.get("genomad_plasmid_flag", False))

        vs2_sc = float(row.get("vs2_max_score", 0) or 0)
        gmd_sc = float(row.get("genomad_virus_score", 0) or 0)
        vs2_h = float(row.get("vs2_hallmark", 0) or 0)
        gmd_h = float(row.get("genomad_n_hallmarks", 0) or 0)
        length = float(row.get("length_bp", 0) or 0)

        support_any = d6 or vlm
        support_both = d6 and vlm

        # No viral signal
        if not vs2 and not gmd and not d6 and not vlm:
            return "D_unknown"

        # Flagged as plasmid
        if plasmid and not (vs2 and gmd):
            return "D_plasmid_flagged"

        # High confidence: consensus between major tools
        if vs2 and gmd:
            return "A_high-confidence_viral"

        # High confidence: geNomad with strong hallmarks
        if gmd and gmd_sc >= 0.9 and gmd_h >= 1:
            return "A_high-confidence_viral"

        # High confidence: VirSorter2 with strong hallmarks
        if vs2 and vs2_sc >= vs2_high and vs2_h >= 1:
            return "A_high-confidence_viral"

        # Medium confidence: one major tool + secondary support
        if (vs2 or gmd) and support_any:
            return "B_medium-confidence_viral"

        # Medium confidence: one major tool + length support
        if (vs2 or gmd) and length >= length_short:
            return "B_medium-confidence_viral"

        # Medium confidence: one major tool alone
        if vs2 or gmd:
            return "B_medium-confidence_viral"

        # Probable: secondary tools agree
        if support_both:
            return "C_probable_viral"

        return "D_unknown"

    master["confidence_tier"] = master.apply(tier, axis=1)
    return master


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # Required arguments
    parser.add_argument("--contigs_dir", required=True,
                        help="Directory with *_rnaspades_min500bp_transcripts.fasta files")
    parser.add_argument("--vs2_dir", required=True,
                        help="Directory with VirSorter2 outputs ({sample}_virsorter2_out/)")
    parser.add_argument("--genomad_dir", required=True,
                        help="Directory with geNomad outputs ({sample}_genomad_out/)")
    parser.add_argument("--deep6_dir", required=True,
                        help="Directory with Deep6 outputs ({sample}_predict_deep6.txt)")
    parser.add_argument("--viralm_dir", required=True,
                        help="Directory with ViraLM outputs ({sample}_viralm_out/)")
    parser.add_argument("--output", required=True,
                        help="Output TSV file path")

    # Optional arguments
    parser.add_argument("--fasta_suffix",
                        default="_rnaspades_min500bp_transcripts.fasta",
                        help="FASTA file suffix (default: _rnaspades_min500bp_transcripts.fasta)")
    parser.add_argument("--vs2_min", type=float, default=0.5,
                        help="VirSorter2 score threshold (default: 0.5)")
    parser.add_argument("--vs2_high", type=float, default=0.9,
                        help="VirSorter2 high-confidence threshold (default: 0.9)")
    parser.add_argument("--genomad_min", type=float, default=0.7,
                        help="geNomad virus score threshold (default: 0.7)")
    parser.add_argument("--genomad_plasmid", type=float, default=0.7,
                        help="geNomad plasmid score threshold (default: 0.7)")
    parser.add_argument("--deep6_factor", type=float, default=1.25,
                        help="Deep6 top-score multiplier (default: 1.25)")
    parser.add_argument("--deep6_min_score", type=float, default=0.70,
                        help="Deep6 minimum absolute score (default: 0.70)")
    parser.add_argument("--viralm_min", type=float, default=0.70,
                        help="ViraLM score threshold (default: 0.70)")
    parser.add_argument("--length_short", type=int, default=2500,
                        help="Minimum length for tier B support (default: 2500 bp)")
    parser.add_argument("--require_all_tools", dest="require_all_tools",
                        action="store_true", default=True,
                        help="Only build the table from samples that have output from "
                             "all four tools (default)")
    parser.add_argument("--allow_missing_tools", dest="require_all_tools",
                        action="store_false",
                        help="Build the table even when a sample is missing one or more "
                             "tools. Missing tools score as not-called, which lowers the "
                             "tier a contig can reach -- interpret tiers accordingly.")

    args = parser.parse_args()

    print("=" * 60)
    print("Building master contig table")
    print("=" * 60)

    # Get pool list
    pools = get_pools(args.contigs_dir, args.fasta_suffix)

    # Which tools produced output for which pool?
    def tool_outputs(pool):
        stem = f"{pool}{args.fasta_suffix.replace('.fasta', '')}"
        return {
            "VirSorter2": Path(args.vs2_dir) / f"{pool}_virsorter2_out" / "final-viral-score.tsv",
            "geNomad":    Path(args.genomad_dir) / f"{pool}_genomad_out" / f"{stem}_summary" / f"{stem}_virus_summary.tsv",
            "Deep6":      Path(args.deep6_dir) / f"{pool}_predict_deep6.txt",
            "ViraLM":     Path(args.viralm_dir) / f"{pool}_viralm_out" / f"result_{stem}.csv",
        }

    complete, incomplete = [], {}
    for pool in pools:
        absent = [t for t, f in tool_outputs(pool).items() if not f.exists()]
        if absent:
            incomplete[pool] = absent
        else:
            complete.append(pool)

    if incomplete:
        print(f"\n  {len(incomplete)} sample(s) are missing tool output:")
        for pool, absent in list(incomplete.items())[:10]:
            print(f"    {pool}: missing {', '.join(absent)}")
        if len(incomplete) > 10:
            print(f"    ... and {len(incomplete) - 10} more")

    if args.require_all_tools:
        if not complete:
            sys.exit("ERROR: no sample has output from all four tools.\n"
                     "       Re-run the missing classifiers, or pass --allow_missing_tools "
                     "to build the table anyway (tiers will be conservative).")
        if incomplete:
            print(f"  Using the {len(complete)} sample(s) with all four tools. "
                  f"Pass --allow_missing_tools to include the rest.")
        pools = complete
    elif incomplete:
        print("  --allow_missing_tools: keeping every sample; absent tools count as not-called.")

    # Parse all inputs
    print("\nParsing contig lengths...")
    master = parse_contigs(pools, args.contigs_dir, args.fasta_suffix)

    print("\nParsing VirSorter2...")
    vs2 = parse_virsorter2(pools, args.vs2_dir, args.vs2_min)

    print("\nParsing geNomad...")
    genomad = parse_genomad(pools, args.genomad_dir, args.fasta_suffix,
                            args.genomad_min, args.genomad_plasmid)

    print("\nParsing Deep6...")
    deep6 = parse_deep6(pools, args.deep6_dir, args.deep6_factor, args.deep6_min_score)

    print("\nParsing ViraLM...")
    viralm = parse_viralm(pools, args.viralm_dir, args.fasta_suffix, args.viralm_min)

    # Merge all results
    for df in [vs2, genomad, deep6, viralm]:
        if not df.empty and "contig_id" in df.columns and len(df.columns) > 2:
            master = master.merge(df, on=KEYS, how="left")

    # Fill boolean columns
    bool_cols = ["vs2_called", "genomad_called", "deep6_is_viral",
                 "genomad_plasmid_flag", "viralm_called"]
    for col in bool_cols:
        if col in master.columns:
            master[col] = master[col].fillna(False).astype(bool)

    # Assign confidence tiers
    master = assign_tiers(master, args.vs2_high, args.length_short)

    # Print summary
    print("\n" + "=" * 60)
    print("RESULTS SUMMARY")
    print("=" * 60)
    print(f"Total contigs: {len(master):,}")
    print(f"Unique pools:  {master['sample_id'].nunique()}")

    print("\n── Confidence tiers ──")
    for tier, n in master["confidence_tier"].value_counts().items():
        pct = 100 * n / len(master)
        print(f"  {tier:<30} {n:>8,}  ({pct:.1f}%)")

    # Per-tool summary
    tool_cols = {
        "VirSorter2": "vs2_called",
        "geNomad": "genomad_called",
        "Deep6": "deep6_is_viral",
        "ViraLM": "viralm_called"
    }
    present = {name: col for name, col in tool_cols.items() if col in master.columns}

    print("\n── Per-tool viral calls ──")
    for name, col in present.items():
        n = master[col].sum()
        pct = 100 * n / len(master)
        print(f"  {name:<12} {n:>8,}  ({pct:.1f}%)")

    # Tool agreement
    if len(present) >= 2:
        cols = list(present.values())
        names = list(present.keys())
        print("\n── Tool agreement ──")
        for r in range(2, len(cols) + 1):
            for combo in combinations(range(len(cols)), r):
                mask = master[cols[combo[0]]].copy()
                label = names[combo[0]]
                for i in combo[1:]:
                    mask = mask & master[cols[i]]
                    label += " + " + names[i]
                print(f"  {label:<40} {mask.sum():>8,}")

    # Save output
    master.to_csv(args.output, sep="\t", index=False)
    print(f"\nSaved: {args.output}")
    print("=" * 60)


if __name__ == "__main__":
    main()
