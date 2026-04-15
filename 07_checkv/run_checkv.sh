#!/bin/bash
################################################################################
# CheckV: Assess viral genome completeness and contamination
################################################################################
#
# DESCRIPTION:
#   Runs CheckV (https://bitbucket.org/berkeleylab/checkv) to assess viral
#   contig quality by estimating genome completeness, detecting integrated
#   proviruses, and identifying potential host contamination.
#
# REQUIREMENTS & INSTALLATION:
#   1. Install CheckV via conda:
#      mamba create -n checkv -c bioconda checkv
#   2. Download the CheckV database:
#      checkv download-database /path/to/checkv-db-v1.5
#      (store this path and use it in CHECKV_DB below)
#
# USAGE:
#   # Set the database path and run:
#   export CHECKV_DB="/path/to/checkv-db-v1.5"
#   ./run_checkv.sh
#   # or with custom settings:
#   INPUT_FASTA="viral_ABC.fasta" OUTPUT_DIR="checkv_output" ./run_checkv.sh
#
# INPUT:
#   Input FASTA file (default: viral_ABC_1kb.fasta)
#   Should contain viral sequences from extract_viral_contigs.py
#
# OUTPUT:
#   OUTPUT_DIR/
#     quality_summary.tsv       — Completeness and quality estimates
#     contamination.tsv         — Provirus integration info (if detected)
#     protein_calls/            — ORF predictions
#     checkv_results.txt        — Summary report
#
# CONFIGURATION:
#   Edit the variables below or set via environment before calling this script.
#
################################################################################

# ─────────────────────────────────────────────────────────────────────────────
# CONFIGURATION — EDIT THESE OR SET VIA ENVIRONMENT
# ─────────────────────────────────────────────────────────────────────────────

# Input FASTA file containing viral contigs
INPUT_FASTA="${INPUT_FASTA:-viral_ABC_1kb.fasta}"

# Output directory for CheckV results
OUTPUT_DIR="${OUTPUT_DIR:-checkv_output}"

# Path to CheckV database (must be set by user)
CHECKV_DB="${CHECKV_DB:-}"

# Number of threads
THREADS="${THREADS:-12}"

# ─────────────────────────────────────────────────────────────────────────────
# MAIN WORKFLOW
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

echo "========================================================================"
echo "CheckV — Viral genome completeness and quality assessment"
echo "========================================================================"
echo "Input FASTA:    $INPUT_FASTA"
echo "Output dir:     $OUTPUT_DIR"
echo "CheckV DB:      $CHECKV_DB"
echo "Threads:        $THREADS"
echo ""

# Validate inputs
if [ ! -f "$INPUT_FASTA" ]; then
    echo "ERROR: Input FASTA not found: $INPUT_FASTA"
    exit 1
fi

if [ -z "$CHECKV_DB" ] || [ ! -d "$CHECKV_DB" ]; then
    echo "ERROR: CHECKV_DB not set or directory not found: $CHECKV_DB"
    echo "       1. Download: checkv download-database /path/to/checkv-db-v1.5"
    echo "       2. Set via environment: export CHECKV_DB='/path/to/checkv-db-v1.5'"
    exit 1
fi

echo "Running CheckV..."
checkv end_to_end \
    "$INPUT_FASTA" \
    "$OUTPUT_DIR" \
    -d "$CHECKV_DB" \
    --remove_tmp \
    -t "$THREADS"

echo ""
echo "========================================================================"
echo "CheckV analysis complete"
echo "Output directory: $OUTPUT_DIR"
echo ""
echo "Key output files:"
echo "  $OUTPUT_DIR/quality_summary.tsv (use with add_checkv_results.py)"
echo "  $OUTPUT_DIR/contamination.tsv (provirus info, if present)"
echo "========================================================================"
