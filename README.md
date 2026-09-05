# Viral Contig Classification — Multi-tool Approach

A bioinformatics pipeline for identifying and quality-filtering viral contigs from assembled metagenomes, combining four independent viral classifiers with a confidence tier system. Works with both RNA and DNA viromes.

---

## Overview

**Starting point: assembled contigs (FASTA format).** This pipeline does not include a read processing or assembly step — it picks up from wherever you have assembled contigs. If you need to go from raw FASTQ reads to assembled contigs, see the companion repository: [rna-virome-blastx-pipeline](https://github.com/eddiefuques/rna-virome-blastx-pipeline), which includes QC, trimming, normalization, and rnaSPAdes assembly.

Viral identification is done by running four classifiers independently and then merging their results into a single master table. Every contig receives a confidence tier based on how many tools agree it is viral. The pipeline is particularly suited for discovering novel and divergent viruses that may have little or no similarity to sequences in reference databases — cases where a BLAST-based approach would miss them entirely.

---

## When to use this pipeline

- Your starting material is already assembled (any assembler, any library type)
- You want to identify viruses beyond those in reference databases
- You need orthogonal validation from multiple independent methods
- Your data includes DNA viruses, RNA viruses, or both

---

## Pipeline Overview

```
Assembled contigs (FASTA, ≥ 500 bp)
    │
    ├──── 00 Host screen    (optional) minimap2 vs. host genome → per-contig flags
    │
    ├──── 01 VirSorter2     Score-based | dsDNA, ssDNA, RNA, NCLDV, lavidaviridae
    ├──── 02 Deep6          Deep learning | 4 viral realms vs. bacteria + eukaryotes
    ├──── 03 geNomad        Marker genes + neural network | virus, plasmid, or chromosome
    └──── 04 ViraLM         Large language model | trained on viral sequences
    │
    ▼
05  build_master_table.py   Merge all 4 outputs; assign confidence tiers (A/B/C/D)
    │
05b add_host_flags.py       (optional) attach host-alignment evidence per contig
    │
    ▼
    Confidence tiers:
      A — High confidence   Both backbone tools (VS2 + geNomad)
                            OR solo call with score ≥ 0.9 + ≥ 1 hallmark gene
      B — Medium            One backbone + ≥ 1 support tool
                            OR one backbone + contig length ≥ 2500 bp
      C — Probable viral    Both support tools agree (Deep6 + ViraLM), no backbone
      D — Unknown/plasmid   ≤ 1 tool, or geNomad plasmid flag
    │
    ▼
06  extract_viral_contigs.py  Filter to ≥ 1 kb; extract per-tier FASTAs (A, B, C, AB, ABC)
    build_viral_table.py       Subset master table to viral candidates only
    │
    ▼
07  CheckV                  Estimate genome completeness and contamination
    add_checkv_results.py   Merge CheckV metrics into viral candidate table
    │
    ▼
    Final table: viral candidates with per-tool scores, confidence tiers,
                 and CheckV quality estimates
```

---

## Usage

### 1. Set up the environment

```bash
conda env create -f environment.yml
conda activate viral-metagenomics-p2
```

See `environment.yml` for notes on tools that require manual installation (VirSorter2, Deep6, ViraLM).

### 2. Run the four classifiers

Each tool runs independently. Edit the `CONFIGURATION` block at the top of each script to set your paths, then submit to your cluster or run locally.

```bash
# VirSorter2 (requires pixi installation — see script header)
bash 01_virsorter2/run_virsorter2.sh

# Deep6 (set DEEP6_SCRIPT and MODELS_DIR — see script header)
bash 02_deep6/run_deep6.sh

# geNomad (set GENOMAD_DB — see script header)
bash 03_genomad/run_genomad.sh

# ViraLM (set VIRALM_SCRIPT and VIRALM_MODEL — see script header)
bash 04_viralm/run_viralm.sh
```

### 2b. (Optional but recommended for host-associated samples) Host screen

Vertebrate and invertebrate genomes carry endogenous viral elements — endogenous
retroviruses above all — which VirSorter2 and geNomad score as viral because, by
sequence, they are. They are the standard false-positive class in host-associated
virome work, and read-level host removal does not catch them: unspliced aligners
(bowtie2, bwa) also miss host mRNA reads spanning exon junctions, so host
transcripts survive filtering and assemble.

This step aligns contigs to the host genome and records, per contig, the fraction
covered by its best alignment. It removes nothing — the flags travel with the
table so you can judge a host-mapping viral call yourself.

```bash
export CONTIGS_DIR=/path/to/contigs
export HOST_GENOME=/path/to/host_genome.fna.gz
export OUTDIR=/path/to/host_screen
bash 00_host_screen/run_host_screen.sh
```

Use `PRESET=splice` (the default) for transcript assemblies against a genomic
reference, `PRESET=asm10`/`asm20` for DNA assemblies. Then, after step 3, merge
the flags in and pass the flagged table to step 4:

```bash
python 00_host_screen/add_host_flags.py \
    --master_table master_contig_table.tsv \
    --host_dir     /path/to/host_screen \
    --output       master_contig_table_hostflagged.tsv
```

This adds `host_flag` (`host` / `host_partial` / `none`), `host_best_frac`,
`host_identity` and supporting columns, and writes
`host_flagged_viral_candidates.tsv` — the viral-tier calls that also map to the
host, which is the list worth reviewing by hand. A long alignment at low identity
points to an ancient endogenous element; a short high-identity one to recent host
transcript.

### 3. Build master contig table

```bash
python 05_build_master_table/build_master_table.py \
    --contigs_dir  /path/to/contigs \
    --vs2_dir      /path/to/virsorter2_outputs \
    --genomad_dir  /path/to/genomad_outputs \
    --deep6_dir    /path/to/deep6_outputs \
    --viralm_dir   /path/to/viralm_outputs \
    --output       master_contig_table.tsv
```

All classification thresholds are configurable (see `--help`). Defaults reflect validated settings from the original study:

| Parameter | Default | Description |
|---|---|---|
| `--vs2_min` | 0.5 | VirSorter2 min score |
| `--genomad_min` | 0.7 | geNomad min virus_score |
| `--deep6_factor` | 1.25 | Deep6 top score must be ≥ 1.25× median |
| `--deep6_min_score` | 0.70 | Deep6 absolute score floor |
| `--viralm_min` | 0.70 | ViraLM min virus_score (more stringent than tool default) |

By default the table is built only from samples that have output from all four
tools; samples missing a tool are listed and skipped. Pass `--allow_missing_tools`
to include them anyway — absent tools then count as not-called, so affected
contigs cannot reach the tiers that require those tools.

### 4. Extract viral candidates

```bash
# Filter to ≥ 1 kb and extract per-tier FASTAs
python 06_extract_viral_contigs/extract_viral_contigs.py \
    --master_table  master_contig_table.tsv \
    --contigs_dir   /path/to/contigs \
    --output_dir    ./viral_candidates

# Subset master table to tiers A + B + C
python 06_extract_viral_contigs/build_viral_table.py \
    --master_table  ./viral_candidates/master_contig_table_1kb.tsv \
    --output        viral_candidates_ABC_1kb.tsv
```

### 5. CheckV quality assessment

```bash
# Run CheckV (requires CheckV database — see environment.yml)
bash 07_checkv/run_checkv.sh \
    --input   ./viral_candidates/viral_ABC_1kb.fasta \
    --output  ./checkv_output \
    --db      /path/to/checkv_db

# Merge quality metrics into viral candidate table
python 07_checkv/add_checkv_results.py \
    --input   viral_candidates_ABC_1kb.tsv \
    --checkv  ./checkv_output \
    --output  viral_candidates_final.tsv
```

---

## Output

The final `viral_candidates_final.tsv` contains one row per viral candidate contig:

| Column | Source | Description |
|---|---|---|
| contig_id | assembly | Contig identifier (unique within a sample) |
| sample_id | assembly | Sample of origin |
| length_bp | assembly | Contig length (bp) |
| confidence_tier | pipeline | A / B / C / D |
| vs2_max_score | VirSorter2 | Max score across viral groups |
| vs2_group | VirSorter2 | Best-scoring viral group |
| vs2_hallmark | VirSorter2 | Hallmark gene count |
| genomad_virus_score | geNomad | Virus probability (0–1) |
| genomad_n_hallmarks | geNomad | Hallmark gene count |
| genomad_taxonomy | geNomad | Predicted taxonomic lineage |
| genomad_plasmid_flag | geNomad | Plasmid contamination flag |
| deep6_top_class | Deep6 | Top predicted class |
| deep6_top_score | Deep6 | Score for top class |
| viralm_score | ViraLM | Virus probability (0–1) |
| checkv_quality | CheckV | Complete / High / Medium / Low / Not-determined |
| checkv_completeness | CheckV | Estimated genome completeness (%) |
| checkv_contamination | CheckV | Host contamination estimate (%) |
| checkv_provirus | CheckV | Integrated provirus flag |
| host_flag | host screen | host / host_partial / none (optional step 00) |
| host_best_frac | host screen | Fraction of the contig covered by the host genome |
| host_identity | host screen | Identity of that alignment |

Per-tier FASTA files are also produced for downstream analyses (phylogenetics, genome annotation, etc.). The suffix is the `--min_length` you passed, so with the default 1000 you get:
- `viral_tier_A_1000bp.fasta` — High-confidence
- `viral_tier_B_1000bp.fasta` — Medium-confidence
- `viral_tier_C_1000bp.fasta` — Probable viral
- `viral_AB_1000bp.fasta`     — Tiers A + B (the usual working set)
- `viral_ABC_1000bp.fasta`    — Tiers A + B + C

---

## Assumptions and gotchas

**Each classifier needs its own environment.** VirSorter2, Deep6 and ViraLM have
conflicting dependencies, so the batch scripts deliberately do not activate
anything — activate the right environment, then run the script.

**Do not add `conda activate` inside these scripts.** They run under
`set -euo pipefail`, and MKL-linked environments source an `activate.d` hook that
reads an unset `$MKL_INTERFACE_LAYER`. Under `set -u` that is fatal, and the job
dies seconds in with a misleading "unbound variable" message. If you must
activate inside a script, wrap it: `set +u; conda activate <env>; set -u`.

**Contig names may repeat across samples.** Assemblers emit generic names
(`NODE_1_length_..._cov_...`) that are identical in every assembly, so tool
results are merged on `(sample_id, contig_id)`, not on `contig_id` alone. Contig
IDs only need to be unique *within* a sample. (Before this was fixed, one
sample's scores could be assigned to another sample's identically-named contigs,
and the per-tier FASTAs could come out empty.)

**GNU Parallel** is required by the VirSorter2, geNomad and ViraLM batch scripts.

**ViraLM** refuses to start if its output directory already exists; the batch
script passes `--force` and lets the tool create the directory itself.

**Re-running is safe.** Every batch script skips a sample whose expected output
file is already present, so an interrupted run resumes rather than restarts.

---

## Dependencies

| Tool | Install | Purpose |
|---|---|---|
| VirSorter2 | pixi (see script header) | Score-based viral detection |
| Deep6 | from source | Deep learning classification |
| geNomad | conda (bioconda) | Marker gene + neural network |
| ViraLM | from source | Language model classification |
| CheckV | conda (bioconda) | Genome completeness |
| minimap2 | conda (bioconda) | Host screen (optional step 00) |
| GNU Parallel | conda | Batch parallelization |
| Python ≥3.8 + numpy + pandas | conda | Data merging and analysis |

See `environment.yml` for the full install instructions and database download commands.

---

## Citations

If you use this pipeline, please also cite the individual tools:

- **VirSorter2:** Guo J. *et al.* (2021). VirSorter2: a multi-classifier, expert-guided approach to detect diverse DNA and RNA viruses. *Microbiome* 9:37. https://doi.org/10.1186/s40168-020-00990-y
- **Deep6:** Teunisse E. *et al.* (2023). Deep6: Classification of Metatranscriptomic Sequences into Cellular Empires and Viral Realms Using Deep Learning Models. *Microbiol Resour Announc* 12:e01079-22. https://doi.org/10.1128/mra.01079-22 | [GitHub](https://github.com/janfelix/Deep6)
- **geNomad:** Camargo A.P. *et al.* (2023). Identification of mobile genetic elements with geNomad. *Nature Biotechnology* 41:1303–1312. https://doi.org/10.1038/s41587-023-01953-y
- **ViraLM:** Peng C. *et al.* (2024). ViraLM: empowering virus discovery through the genome foundation model. *Bioinformatics* 40(12):btae704. https://doi.org/10.1093/bioinformatics/btae704 | [GitHub](https://github.com/ChengPENG-wolf/ViraLM)
- **CheckV:** Nayfach S. *et al.* (2021). CheckV assesses the quality and completeness of metagenome-assembled viral genomes. *Nature Biotechnology* 39:578–585. https://doi.org/10.1038/s41587-020-00774-7

---

## Contact

Eddie Fuques — eddiefuques@gmail.com
