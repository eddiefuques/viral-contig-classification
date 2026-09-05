#!/bin/bash
################################################################################
# Host screen: flag contigs of host origin before trusting viral calls
################################################################################
#
# DESCRIPTION:
#   Aligns assembled contigs to a host reference genome with minimap2 and
#   records, for every contig, the fraction of its length covered by its best
#   host alignment. The result feeds add_host_flags.py, which merges the flags
#   into the master contig table.
#
# WHY THIS STEP EXISTS:
#   Read-level host removal with an unspliced aligner (bowtie2, bwa) misses host
#   mRNA reads that span exon junctions, so host transcripts survive filtering
#   and assemble. More importantly, vertebrate genomes carry endogenous viral
#   elements -- endogenous retroviruses above all -- that VirSorter2 and geNomad
#   score as genuinely viral. They are the standard false-positive class in
#   host-associated virome studies, and no viral classifier can distinguish them
#   from a circulating virus, because by sequence they are one.
#
#   Nothing is deleted here. A contig that is both host-mapping and called viral
#   is the one you want to look at by hand, not the one you want silently gone.
#
# REQUIREMENTS:
#   - minimap2 (https://github.com/lh3/minimap2)
#   - a host reference genome (FASTA, may be gzipped)
#
# USAGE:
#   export CONTIGS_DIR=/path/to/contigs
#   export HOST_GENOME=/path/to/host_genome.fna.gz
#   export OUTDIR=/path/to/host_screen
#   ./run_host_screen.sh
#
# INPUT:
#   FASTA files matching: ${CONTIGS_DIR}/*${FASTA_SUFFIX}
#
# OUTPUT (per sample, in ${OUTDIR}):
#   {sample}_host_screen.tsv    contig_id, length_bp, host_best_frac,
#                               host_identity, host_aln_bp, host_mapq,
#                               host_n_hits, host_target
#   {sample}_host_contigs.txt   contigs with host_best_frac >= FRAC_FULL
#   {sample}_host_partial.txt   contigs with FRAC_PARTIAL <= frac < FRAC_FULL
#   {sample}_host.paf           raw alignments, kept for auditing
#
# CONFIGURATION:
#   Edit the variables below or set them in the environment before calling.
#
################################################################################

# ─────────────────────────────────────────────────────────────────────────────
# CONFIGURATION — EDIT THESE OR SET VIA ENVIRONMENT
# ─────────────────────────────────────────────────────────────────────────────

# Directory containing input contig FASTA files
CONTIGS_DIR="${CONTIGS_DIR:-.}"

# Output directory for the screen
OUTDIR="${OUTDIR:-./host_screen}"

# Host reference genome (FASTA or FASTA.gz) — must be set
HOST_GENOME="${HOST_GENOME:-}"

# FASTA file suffix
FASTA_SUFFIX="${FASTA_SUFFIX:-_rnaspades_min500bp_transcripts.fasta}"

# Threads for minimap2
THREADS="${THREADS:-8}"

# minimap2 preset.
#   splice  — transcript assemblies vs. a genomic reference (RNA virome; this is
#             the point of the step, since it handles exon junctions)
#   asm10   — assembled contigs vs. a closely related genomic reference
#   asm20   — as above, more divergent
PRESET="${PRESET:-splice}"

# A contig is "host" when its best alignment covers >= FRAC_FULL of it, and
# "host_partial" from FRAC_PARTIAL up to that. The partial tier catches chimeric
# assemblies and divergent endogenous elements, which rarely align end to end.
FRAC_FULL="${FRAC_FULL:-0.5}"
FRAC_PARTIAL="${FRAC_PARTIAL:-0.2}"

# ─────────────────────────────────────────────────────────────────────────────
# MAIN WORKFLOW
# ─────────────────────────────────────────────────────────────────────────────

set -uo pipefail

echo "========================================================================"
echo "Host screen — minimap2 -x ${PRESET}"
echo "========================================================================"
echo "Contigs dir:   $CONTIGS_DIR"
echo "Host genome:   $HOST_GENOME"
echo "Output dir:    $OUTDIR"
echo "Thresholds:    host >= ${FRAC_FULL}, partial >= ${FRAC_PARTIAL}"
echo ""

command -v minimap2 >/dev/null 2>&1 || { echo "ERROR: minimap2 not on PATH"; exit 1; }

if [ -z "$HOST_GENOME" ] || [ ! -f "$HOST_GENOME" ]; then
    echo "ERROR: HOST_GENOME not set or not found: $HOST_GENOME"
    echo "       export HOST_GENOME=/path/to/host_genome.fna.gz"
    exit 1
fi

if [ ! -d "$CONTIGS_DIR" ]; then
    echo "ERROR: CONTIGS_DIR does not exist: $CONTIGS_DIR"
    exit 1
fi

mkdir -p "$OUTDIR"

n_done=0
for fasta in "${CONTIGS_DIR}"/*"${FASTA_SUFFIX}"; do
    [ -f "$fasta" ] || continue
    sample="$(basename "$fasta" "$FASTA_SUFFIX")"

    paf="${OUTDIR}/${sample}_host.paf"
    tsv="${OUTDIR}/${sample}_host_screen.tsv"
    full="${OUTDIR}/${sample}_host_contigs.txt"
    part="${OUTDIR}/${sample}_host_partial.txt"

    if [ -s "$paf" ] && [ "$paf" -nt "$fasta" ]; then
        echo "[SKIP] $sample (alignments up to date)"
    else
        echo "[RUN]  $sample"
        minimap2 -x "$PRESET" -t "$THREADS" "$HOST_GENOME" "$fasta" \
            2> "${OUTDIR}/${sample}_minimap2.log" > "$paf" || {
            echo "ERROR: minimap2 failed for $sample (see ${OUTDIR}/${sample}_minimap2.log)"
            exit 1
        }
    fi

    # PAF fields: 1 qname 2 qlen 3 qstart 4 qend 5 strand 6 tname 7 tlen
    #             8 tstart 9 tend 10 nmatch 11 alnlen 12 mapq
    # Keep the single best alignment per contig (largest query fraction).
    awk -v OFS='\t' '
        {
            frac = ($4 - $3) / $2
            if (frac > best[$1]) {
                best[$1] = frac; tgt[$1] = $6; nm[$1] = $10
                aln[$1] = $11; mapq[$1] = $12; qlen[$1] = $2
            }
            nhit[$1]++
        }
        END {
            for (c in best)
                print c, qlen[c], best[c], (aln[c] > 0 ? nm[c]/aln[c] : 0), \
                      aln[c], mapq[c], nhit[c], tgt[c]
        }
    ' "$paf" | sort -k3,3gr | cat <(printf 'contig_id\tlength_bp\thost_best_frac\thost_identity\thost_aln_bp\thost_mapq\thost_n_hits\thost_target\n') - > "$tsv"

    awk -v F="$FRAC_FULL" 'NR>1 && $3 >= F {print $1}' "$tsv" | sort -u > "$full"
    awk -v F="$FRAC_FULL" -v P="$FRAC_PARTIAL" \
        'NR>1 && $3 >= P && $3 < F {print $1}' "$tsv" | sort -u > "$part"

    total=$(grep -c "^>" "$fasta")
    nfull=$(wc -l < "$full")
    npart=$(wc -l < "$part")
    printf '       %-20s host=%-6s partial=%-6s of %s contigs\n' \
        "$sample" "$nfull" "$npart" "$total"
    n_done=$((n_done + 1))
done

if [ "$n_done" -eq 0 ]; then
    echo "ERROR: no FASTA files matching *${FASTA_SUFFIX} in $CONTIGS_DIR"
    exit 1
fi

echo ""
echo "========================================================================"
echo "Host screen complete — $n_done sample(s)"
echo "Output directory: $OUTDIR"
echo ""
echo "Next: merge the flags into the master table (after build_master_table.py):"
echo "  python 00_host_screen/add_host_flags.py \\"
echo "      --master_table master_contig_table.tsv \\"
echo "      --host_dir ${OUTDIR} \\"
echo "      --output master_contig_table_hostflagged.tsv"
echo "========================================================================"
