#!/bin/bash
################################################################################
# VirSorter2: Viral classification via sequence scoring and hallmark genes
################################################################################
#
# DESCRIPTION:
#   Runs VirSorter2 on a batch of contig FASTA files using GNU parallel for
#   parallelization. VirSorter2 identifies viral sequences by examining
#   sequence composition, hallmark genes, and genome organization.
#
# REQUIREMENTS:
#   - VirSorter2 (installed via pixi: https://github.com/jiarong/VirSorter2)
#   - GNU parallel (https://www.gnu.org/software/parallel/)
#
# USAGE:
#   ./run_virsorter2.sh
#   # or with custom settings:
#   CONTIGS_DIR=/my/contigs OUTDIR_BASE=/my/output ./run_virsorter2.sh
#
# INPUT:
#   Expects FASTA files matching: ${CONTIGS_DIR}/*${FASTA_SUFFIX}
#   Default: *_rnaspades_min500bp_transcripts.fasta
#
# OUTPUT:
#   For each sample: ${OUTDIR_BASE}/{sample}_virsorter2_out/
#   Key file: final-viral-score.tsv (used by build_master_table.py)
#
# CONFIGURATION:
#   Edit the variables below to customize paths, thread allocation, and
#   parallelization settings.
#
################################################################################

# ─────────────────────────────────────────────────────────────────────────────
# CONFIGURATION — EDIT THESE
# ─────────────────────────────────────────────────────────────────────────────

# Directory containing input contig FASTA files
CONTIGS_DIR="${CONTIGS_DIR:-.}"

# Base output directory (will create {sample}_virsorter2_out subdirs)
OUTDIR_BASE="${OUTDIR_BASE:-.}"

# Number of threads per sample
THREADS_PER_SAMPLE="${THREADS_PER_SAMPLE:-8}"

# Maximum number of samples to run in parallel
MAX_PARALLEL="${MAX_PARALLEL:-5}"

# FASTA file suffix (files matching *${FASTA_SUFFIX} will be processed)
FASTA_SUFFIX="${FASTA_SUFFIX:-_rnaspades_min500bp_transcripts.fasta}"

# VirSorter2 minimum contig length (bp)
MIN_LENGTH="${MIN_LENGTH:-500}"

# Viral groups to include (comma-separated)
VIRAL_GROUPS="${VIRAL_GROUPS:-dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae}"

# ─────────────────────────────────────────────────────────────────────────────
# MAIN WORKFLOW
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

echo "========================================================================"
echo "VirSorter2 — Multi-sample viral classification (batch mode)"
echo "========================================================================"
echo "Contigs dir:        $CONTIGS_DIR"
echo "Output dir:         $OUTDIR_BASE"
echo "Threads/sample:     $THREADS_PER_SAMPLE"
echo "Max parallel:       $MAX_PARALLEL"
echo "FASTA suffix:       $FASTA_SUFFIX"
echo "Min length:         $MIN_LENGTH bp"
echo "Viral groups:       $VIRAL_GROUPS"
echo ""

# Verify directories
if [ ! -d "$CONTIGS_DIR" ]; then
    echo "ERROR: CONTIGS_DIR does not exist: $CONTIGS_DIR"
    exit 1
fi

mkdir -p "$OUTDIR_BASE"

# Build sample list
cd "$CONTIGS_DIR"
echo "Building sample list..."
ls *"$FASTA_SUFFIX" 2>/dev/null | sed "s/${FASTA_SUFFIX}//" > sample_list.txt || {
    echo "ERROR: No FASTA files found matching *${FASTA_SUFFIX} in $CONTIGS_DIR"
    exit 1
}
num_samples=$(wc -l < sample_list.txt)
echo "Found $num_samples samples"
echo ""

# Export for parallel
export CONTIGS_DIR OUTDIR_BASE THREADS_PER_SAMPLE FASTA_SUFFIX MIN_LENGTH VIRAL_GROUPS

# Define worker function
run_virsorter() {
    local sample="$1"
    local contigs_file="${CONTIGS_DIR}/${sample}${FASTA_SUFFIX}"
    local outdir="${OUTDIR_BASE}/${sample}_virsorter2_out"

    # Skip if already completed
    if [ -f "$outdir/final-viral-score.tsv" ]; then
        echo "[SKIP] $sample (already completed)"
        return 0
    fi

    echo "[RUN] $sample ..."
    mkdir -p "$outdir"

    pixi run virsorter run \
        -w "$outdir" \
        -i "$contigs_file" \
        --min-length "$MIN_LENGTH" \
        --include-groups "$VIRAL_GROUPS" \
        --prep-for-dramv \
        -j "$THREADS_PER_SAMPLE" \
        all

    echo "[DONE] $sample"
}

export -f run_virsorter

# Run in parallel
echo "Running VirSorter2 on $num_samples samples (max $MAX_PARALLEL in parallel)..."
parallel -j "$MAX_PARALLEL" run_virsorter :::: sample_list.txt

echo ""
echo "========================================================================"
echo "VirSorter2 batch run complete"
echo "Output directory: $OUTDIR_BASE"
echo "========================================================================"
