#!/bin/bash
################################################################################
# geNomad: Viral and plasmid identification via sequence composition and markers
################################################################################
#
# DESCRIPTION:
#   Runs geNomad (https://github.com/apcamargo/genomad) to identify viral and
#   plasmid sequences using a combination of sequence composition, open reading
#   frame analysis, and hidden Markov model-based marker detection.
#
# REQUIREMENTS & INSTALLATION:
#   1. Create conda environment:
#      mamba create -n genomad -c conda-forge -c bioconda genomad
#   2. Activate environment:
#      conda activate genomad
#   3. Download the default database:
#      genomad download-database /path/to/genomad_db
#      (store this path and use it in GENOMAD_DB below)
#
# USAGE:
#   # Set the database path and run:
#   export GENOMAD_DB="/path/to/genomad_db"
#   ./run_genomad.sh
#   # or with custom settings:
#   CONTIGS_DIR=/my/contigs OUTDIR_BASE=/my/output ./run_genomad.sh
#
# INPUT:
#   FASTA files matching: ${CONTIGS_DIR}/*${FASTA_SUFFIX}
#   Default: *_rnaspades_min500bp_transcripts.fasta
#
# OUTPUT:
#   For each sample: ${OUTDIR_BASE}/{sample}_genomad_out/
#   Key files:
#     {sample}_summary/{sample}..._virus_summary.tsv
#     {sample}_summary/{sample}..._plasmid_summary.tsv (if plasmids found)
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

# Base output directory (will create {sample}_genomad_out subdirs)
OUTDIR_BASE="${OUTDIR_BASE:-.}"

# Path to geNomad database (must be set — see REQUIREMENTS above)
GENOMAD_DB="${GENOMAD_DB:-}"

# Number of threads per sample
THREADS_PER_SAMPLE="${THREADS_PER_SAMPLE:-12}"

# Maximum number of samples to run in parallel
MAX_PARALLEL="${MAX_PARALLEL:-3}"

# FASTA file suffix
FASTA_SUFFIX="${FASTA_SUFFIX:-_rnaspades_min500bp_transcripts.fasta}"

# ─────────────────────────────────────────────────────────────────────────────
# MAIN WORKFLOW
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

echo "========================================================================"
echo "geNomad — Viral and plasmid identification (batch mode)"
echo "========================================================================"
echo "Contigs dir:        $CONTIGS_DIR"
echo "Output dir:         $OUTDIR_BASE"
echo "Database:           $GENOMAD_DB"
echo "Threads/sample:     $THREADS_PER_SAMPLE"
echo "Max parallel:       $MAX_PARALLEL"
echo "FASTA suffix:       $FASTA_SUFFIX"
echo ""

# Validate database
if [ -z "$GENOMAD_DB" ] || [ ! -d "$GENOMAD_DB" ]; then
    echo "ERROR: GENOMAD_DB not set or directory not found: $GENOMAD_DB"
    echo "       1. Download: genomad download-database /path/to/genomad_db"
    echo "       2. Set via environment: export GENOMAD_DB='/path/to/genomad_db'"
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
export CONTIGS_DIR OUTDIR_BASE GENOMAD_DB THREADS_PER_SAMPLE FASTA_SUFFIX

# Define worker function
run_genomad() {
    local sample="$1"
    local fasta="${CONTIGS_DIR}/${sample}${FASTA_SUFFIX}"
    local outdir="${OUTDIR_BASE}/${sample}_genomad_out"

    # Skip if already completed
    if [ -d "$outdir" ] && [ -f "${outdir}/${sample}${FASTA_SUFFIX%.fasta}_summary/${sample}${FASTA_SUFFIX%.fasta}_virus_summary.tsv" ]; then
        echo "[SKIP] $sample (already completed)"
        return 0
    fi

    echo "[RUN] $sample ..."
    mkdir -p "$outdir"

    genomad end-to-end \
        --threads "$THREADS_PER_SAMPLE" \
        --cleanup \
        "$fasta" \
        "$outdir" \
        "$GENOMAD_DB"

    echo "[DONE] $sample"
}

export -f run_genomad

# Run in parallel
echo "Running geNomad on $num_samples samples (max $MAX_PARALLEL in parallel)..."
parallel -j "$MAX_PARALLEL" run_genomad :::: "$SAMPLE_LIST"

echo ""
echo "========================================================================"
echo "geNomad batch run complete"
echo "Output directory: $OUTDIR_BASE"
echo "========================================================================"
