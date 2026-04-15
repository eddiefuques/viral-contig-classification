#!/bin/bash
################################################################################
# Deep6: Realm classification via deep learning
################################################################################
#
# DESCRIPTION:
#   Runs Deep6 (https://github.com/janfelix/Deep6), a deep learning classifier that categorizes sequences into
#   six realms: duplo, euk, mono, pro, ribo, and vari.
#
# REQUIREMENTS & INSTALLATION:
#   Deep6 requires MANUAL installation from source:
#   1. Clone the repository:
#      git clone https://github.com/janteuni/Deep6.git
#   2. Install dependencies:
#      cd Deep6
#      pip install -r requirements.txt
#   3. Download or obtain the trained models directory
#
#   The DEEP6_SCRIPT variable below must point to:
#      {Deep6_repo}/Master/deep6.py
#
#   The MODELS_DIR variable below must point to:
#      {Deep6_repo}/Models
#
# USAGE:
#   # Set paths before running:
#   export DEEP6_SCRIPT="/path/to/Deep6/Master/deep6.py"
#   export MODELS_DIR="/path/to/Deep6/Models"
#   export CONTIGS_DIR="/my/contigs"
#   export OUTDIR="/my/output"
#   ./run_deep6.sh
#
# INPUT:
#   FASTA files matching: ${CONTIGS_DIR}/*${FASTA_SUFFIX}
#   Default: *_rnaspades_min500bp_transcripts.fasta
#
# OUTPUT:
#   For each sample: ${OUTDIR}/${sample}_predict_deep6.txt
#   Format: TSV with columns [name, duplo, euk, mono, pro, ribo, vari]
#
# CONFIGURATION:
#   Edit the variables below or set via environment before calling this script.
#
################################################################################

# ─────────────────────────────────────────────────────────────────────────────
# CONFIGURATION — EDIT THESE OR SET VIA ENVIRONMENT
# ─────────────────────────────────────────────────────────────────────────────

# Directory containing input contig FASTA files
CONTIGS_DIR="${CONTIGS_DIR:-.}"

# Output directory for Deep6 results
OUTDIR="${OUTDIR:-.}"

# Path to Deep6 main script (must be set by user)
DEEP6_SCRIPT="${DEEP6_SCRIPT:-}"

# Path to Deep6 models directory (must be set by user)
MODELS_DIR="${MODELS_DIR:-}"

# FASTA file suffix
FASTA_SUFFIX="${FASTA_SUFFIX:-_rnaspades_min500bp_transcripts.fasta}"

# Minimum contig length (bp) — must match VirSorter2 preprocessing
MIN_LENGTH="${MIN_LENGTH:-500}"

# ─────────────────────────────────────────────────────────────────────────────
# MAIN WORKFLOW
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

echo "========================================================================"
echo "Deep6 — Contigs classification via deep learning"
echo "========================================================================"
echo "Contigs dir:       $CONTIGS_DIR"
echo "Output dir:        $OUTDIR"
echo "FASTA suffix:      $FASTA_SUFFIX"
echo "Min length:        $MIN_LENGTH bp"
echo ""

# Validate required variables
if [ -z "$DEEP6_SCRIPT" ] || [ ! -f "$DEEP6_SCRIPT" ]; then
    echo "ERROR: DEEP6_SCRIPT not set or file not found: $DEEP6_SCRIPT"
    echo "       Set via environment: export DEEP6_SCRIPT='/path/to/Deep6/Master/deep6.py'"
    exit 1
fi

if [ -z "$MODELS_DIR" ] || [ ! -d "$MODELS_DIR" ]; then
    echo "ERROR: MODELS_DIR not set or directory not found: $MODELS_DIR"
    echo "       Set via environment: export MODELS_DIR='/path/to/Deep6/Models'"
    exit 1
fi

if [ ! -d "$CONTIGS_DIR" ]; then
    echo "ERROR: CONTIGS_DIR does not exist: $CONTIGS_DIR"
    exit 1
fi

mkdir -p "$OUTDIR"

# Process each FASTA file
cd "$CONTIGS_DIR"
count=0
for fasta in *"$FASTA_SUFFIX"; do
    [ -f "$fasta" ] || continue

    pool_id="${fasta%${FASTA_SUFFIX}}"
    outfile="${OUTDIR}/${pool_id}_predict_deep6.txt"

    # Skip if already completed
    if [ -f "$outfile" ]; then
        echo "[SKIP] $pool_id (already completed)"
        continue
    fi

    echo "[RUN] $pool_id ..."

    # Run Deep6
    python "$DEEP6_SCRIPT" \
        -i "$fasta" \
        -l "$MIN_LENGTH" \
        -m "$MODELS_DIR" \
        -o "$OUTDIR"

    # Deep6 outputs to OUTDIR with a specific naming convention;
    # rename to our standard format
    temp_out="${OUTDIR}/$(basename "$fasta")_predict_${MIN_LENGTH}bp_deep6.txt"
    if [ -f "$temp_out" ]; then
        mv "$temp_out" "$outfile"
    fi

    echo "[DONE] $pool_id"
    ((count++))
done

echo ""
echo "========================================================================"
echo "Deep6 analysis complete"
echo "Processed: $count samples"
echo "Output directory: $OUTDIR"
echo "========================================================================"
