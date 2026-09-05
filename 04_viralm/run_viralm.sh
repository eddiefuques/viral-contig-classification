#!/bin/bash
################################################################################
# ViraLM: Viral classification via language model
################################################################################
#
# DESCRIPTION:
#   Runs ViraLM (https://github.com/ChengPENG-wolf/ViraLM),
#   a large language model-based approach for viral sequence classification.
#   ViraLM leverages pre-trained transformer models to identify and classify
#   viral sequences with high accuracy.
#
# REQUIREMENTS & INSTALLATION:
#   ViraLM requires MANUAL installation from source. Follow the official
#   repository instructions at:
#   https://github.com/ChengPENG-wolf/ViraLM
#
#   After installation, locate:
#     VIRALM_SCRIPT: {ViraLM_repo}/viralm.py or equivalent entry point
#     VIRALM_MODEL:  {ViraLM_repo}/model or path to trained model weights
#
# USAGE:
#   # Set paths before running:
#   export VIRALM_SCRIPT="/path/to/ViraLM/viralm.py"
#   export VIRALM_MODEL="/path/to/ViraLM/model"
#   ./run_viralm.sh
#   # or with custom settings:
#   CONTIGS_DIR=/my/contigs OUTDIR_BASE=/my/output ./run_viralm.sh
#
# INPUT:
#   FASTA files matching: ${CONTIGS_DIR}/*${FASTA_SUFFIX}
#   Default: *_rnaspades_min500bp_transcripts.fasta
#
# OUTPUT:
#   For each sample: ${OUTDIR_BASE}/{sample}_viralm_out/
#   Key file: result_{sample}.csv (used by build_master_table.py)
#   Format: CSV with columns [seq_name, virus_score, ...]
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

# Base output directory (will create {sample}_viralm_out subdirs)
OUTDIR_BASE="${OUTDIR_BASE:-.}"

# Path to ViraLM main script (must be set by user)
VIRALM_SCRIPT="${VIRALM_SCRIPT:-}"

# Path to ViraLM trained model (must be set by user)
VIRALM_MODEL="${VIRALM_MODEL:-}"

# Number of threads per sample
THREADS_PER_SAMPLE="${THREADS_PER_SAMPLE:-12}"

# Maximum number of samples to run in parallel
MAX_PARALLEL="${MAX_PARALLEL:-3}"

# FASTA file suffix
FASTA_SUFFIX="${FASTA_SUFFIX:-_rnaspades_min500bp_transcripts.fasta}"

# Minimum contig length passed to ViraLM (must match the other tools)
MIN_LENGTH="${MIN_LENGTH:-500}"

# ─────────────────────────────────────────────────────────────────────────────
# MAIN WORKFLOW
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

echo "========================================================================"
echo "ViraLM — Viral classification (batch mode)"
echo "========================================================================"
echo "Contigs dir:        $CONTIGS_DIR"
echo "Output dir:         $OUTDIR_BASE"
echo "ViraLM script:      $VIRALM_SCRIPT"
echo "ViraLM model:       $VIRALM_MODEL"
echo "Threads/sample:     $THREADS_PER_SAMPLE"
echo "Max parallel:       $MAX_PARALLEL"
echo "FASTA suffix:       $FASTA_SUFFIX"
echo ""

# Validate required variables
if [ -z "$VIRALM_SCRIPT" ] || [ ! -f "$VIRALM_SCRIPT" ]; then
    echo "ERROR: VIRALM_SCRIPT not set or file not found: $VIRALM_SCRIPT"
    echo "       Set via environment: export VIRALM_SCRIPT='/path/to/ViraLM/viralm.py'"
    exit 1
fi

if [ -z "$VIRALM_MODEL" ] || [ ! -d "$VIRALM_MODEL" ]; then
    echo "ERROR: VIRALM_MODEL not set or directory not found: $VIRALM_MODEL"
    echo "       Set via environment: export VIRALM_MODEL='/path/to/ViraLM/model'"
    exit 1
fi

if [ ! -d "$CONTIGS_DIR" ]; then
    echo "ERROR: CONTIGS_DIR does not exist: $CONTIGS_DIR"
    exit 1
fi

mkdir -p "$OUTDIR_BASE"

# Build sample list. Written to the OUTPUT directory, never to CONTIGS_DIR --
# the input directory may be read-only or shared, and it is not ours to litter.
SAMPLE_LIST="${OUTDIR_BASE}/sample_list.txt"
echo "Building sample list..."
( cd "$CONTIGS_DIR" && ls *"$FASTA_SUFFIX" 2>/dev/null | sed "s/${FASTA_SUFFIX}//" ) > "$SAMPLE_LIST"
if [ ! -s "$SAMPLE_LIST" ]; then
    echo "ERROR: No FASTA files found matching *${FASTA_SUFFIX} in $CONTIGS_DIR"
    exit 1
fi
num_samples=$(wc -l < "$SAMPLE_LIST")
echo "Found $num_samples samples"
echo ""

# Export for parallel
export CONTIGS_DIR OUTDIR_BASE VIRALM_SCRIPT VIRALM_MODEL THREADS_PER_SAMPLE FASTA_SUFFIX

# Define worker function
run_viralm() {
    local sample="$1"
    local fasta="${CONTIGS_DIR}/${sample}${FASTA_SUFFIX}"
    local outdir="${OUTDIR_BASE}/${sample}_viralm_out"

    # Skip if already completed
    if [ -f "${outdir}/result_${sample}${FASTA_SUFFIX%.fasta}.csv" ]; then
        echo "[SKIP] $sample (already completed)"
        return 0
    fi

    echo "[RUN] $sample ..."

    # Do NOT create $outdir here. viralm.py exits 1 with "The output directory
    # already exists. Use -f or --force to overwrite." before it loads the model,
    # so pre-creating the directory guarantees failure. It creates the directory
    # itself; --force clears the remains of an interrupted run (we only reach
    # this point when the result CSV is absent).
    mkdir -p "$(dirname "$outdir")"

    python "$VIRALM_SCRIPT" \
        -i "$fasta" \
        -o "$outdir" \
        -d "$VIRALM_MODEL" \
        --len "$MIN_LENGTH" \
        --threads "$THREADS_PER_SAMPLE" \
        --force

    echo "[DONE] $sample"
}

export -f run_viralm

# Run in parallel
echo "Running ViraLM on $num_samples samples (max $MAX_PARALLEL in parallel)..."
parallel -j "$MAX_PARALLEL" run_viralm :::: "$SAMPLE_LIST"

echo ""
echo "========================================================================"
echo "ViraLM batch run complete"
echo "Output directory: $OUTDIR_BASE"
echo "========================================================================"
