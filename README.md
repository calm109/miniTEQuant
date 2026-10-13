# miniTEQuant — Long-Read RNA Isoform and TE Quantification

miniTEQuant is a long-read RNA-seq isoform quantification tool that supports single-cell (SC) and spatial transcriptomics (ST) analyses. It uses a community-based EM algorithm (expectation-maximization) to estimate isoform abundances from long-read alignments. It also provides transposable element (TE) quantification at the per-cell and per-spot level.

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Data Preparation](#data-preparation)
- [Quick Start](#quick-start)
- [Subcommand: `quantify`](#subcommand-quantify)
- [Subcommand: `cal_TE`](#subcommand-cal_te)
- [Module: `sp_tags` (high-resolution spatial barcodes and UMIs)](#module-sp_tags-high-resolution-spatial-barcodes-and-umis)
- [Running Modes](#running-modes)
- [Output Files](#output-files)

---

## Features

- **Multi-mode support**: single-cell (10x Chromium) and spatial transcriptomics (Visium / Visium HD / Stereo-seq)
- **High-resolution spatial preprocessing** (`sp_tags`): for PacBio HiFi and ONT reads, splits concatenated long reads into molecules, identifies spatial IDs (Stereo-seq CID) and UMIs (Stereo-seq MID), corrects them against a whitelist / by UMI clustering, optionally aggregates CIDs into cells (`--barcode2spot`), and writes reads that plug directly into `--st_mode`
- **EM-based isoform quantification**: community-based parallel EM on gene-level bipartite graphs
- **UMI deduplication**: exact-match UMI collapsing in SC/ST modes
- **TE quantification**: annotates transcripts by TE overlap and computes TE expression at locus / subfamily / family / class level for per-cell and per-spot data
![miniTEQuant Overview](miniTEQuant.png) 
---

## Installation

```bash
git clone https://github.com/calm109/miniTEQuant.git
cd miniTEQuant
conda env create -f environment.yml
conda activate miniTEQuant
conda install -c bioconda -c conda-forge flexiplex
conda install -c bioconda minimap2
conda install -c bioconda samtools
```

### Additional dependencies for SC / ST modes

| Tool | Purpose | Install |
|------|---------|---------|
| [minimap2](https://github.com/lh3/minimap2) | Long-read genome alignment | `conda install -c bioconda minimap2` |
| [flexiplex](https://davidsongroup.github.io/flexiplex/) | Barcode discovery and demultiplexing | `conda install -c bioconda -c conda-forge flexiplex` |
| [nailpolish](https://github.com/DavidsonGroup/nailpolish) | UMI-based read deduplication and consensus | See below |

> `minimap2` and `flexiplex` install cleanly via conda/bioconda, as shown above. `nailpolish` does not currently ship a conda package, so it is installed separately from a precompiled binary (below). All three are required for the SC/ST barcode → dedup → align pipeline described under [Data Preparation](#data-preparation).

```bash
# nailpolish (Linux precompiled binary)
wget https://github.com/DavidsonGroup/nailpolish/releases/download/nightly_develop/nailpolish
chmod +x nailpolish
mkdir -p ~/bin && mv nailpolish ~/bin/
echo 'export PATH=~/bin:$PATH' >> ~/.bashrc && source ~/.bashrc
```

---

## Data Preparation

### Long-read alignment (all modes)

Use `minimap2` to align long reads to the **reference genome** in spliced mode. Build a `.mmi` index first to speed up repeated runs.

```bash
# Build genome index (do once)
minimap2 -d genome.mmi genome.fa

# ONT cDNA / direct RNA
minimap2 -y -ax splice --MD --secondary=no -G 200k -t 16 genome.mmi reads.fastq.gz > LR.sam

# PacBio HiFi / Iso-Seq
minimap2 -y --secondary=no -ax splice:hq --MD -t 16 genome.mmi reads.fastq.gz \
    | samtools view -b | samtools sort -@ 8 -t CB -O SAM -o LR.sam
```

> **Note:** The `-y` flag propagates read tags (CB, UB) from the FASTQ header into the SAM output — required for SC/ST modes when barcodes are embedded in the read name.

### SC / ST preprocessing: barcode demultiplexing and UMI deduplication

Before alignment, extract and filter barcodes with flexiplex, then deduplicate UMIs with nailpolish.

```bash
# 1. Barcode discovery (generate whitelist)
gunzip -c reads.fastq.gz | flexiplex -d 10x3v3 -f 0 -p 16 > flexiplex.out
flexiplex-filter flexiplex_barcodes_counts.txt > flexiplex_barcodes_final.txt

# 2. Demultiplexing against the whitelist
# If you already have a barcode whitelist (e.g. a spatial barcode whitelist), skip step 1 (Barcode discovery),
# replace flexiplex_barcodes_final.txt with your own whitelist file,
# and run step 2 (Demultiplexing against the whitelist) directly.
gunzip -c reads.fastq.gz | flexiplex -d 10x3v3 -p 16 -k flexiplex_barcodes_final.txt > demux.fastq

# 3. UMI deduplication and consensus generation
nailpolish index --skip-unmatched demux.fastq
nailpolish summary demux.fastq
nailpolish consensus -t 16 demux.fastq > dedup.fastq
```

> For ST (Visium), a spatial barcode whitelist is typically available from Space Ranger output. Pass it directly to flexiplex with `-k` instead of the auto-discovered whitelist.

### High-resolution ST preprocessing (Stereo-seq): `sp_tags`

For Stereo-seq long-read libraries (PacBio HiFi or ONT), use the built-in `sp_tags` module. It splits concatenated reads into individual molecules, identifies the CID and MID of each molecule, corrects CIDs against the chip whitelist, corrects MID sequencing errors (pre-alignment, edit distance 1 within each CID), optionally aggregates CIDs into cells / spots (`--barcode2spot`), and writes a FASTQ ready for minimap2. After alignment, `run_umi_dedup.py` collapses MIDs per CID and gene before `quantify`.

```bash
# 1. Split concatenated reads, correct CID / MID
python isoform_quantification/sp_tags/run_sp_tags.py \
    -i hifi_reads.fastq.gz \
    -w cid_whitelist.txt \
    --data_type pacbio \
    -o sp_tags_out \
    -t 16
# --barcode2spot cell_map.tsv: ID_COL:SPOT_COL  -> quantify per cell / spot (omit to quantify per CID)

# 2. Align (-y copies the CB / UB / XC tags from the read name into the SAM)
minimap2 -y --secondary=no -ax splice:hq --MD -t 16 genome.mmi sp_tags_out/sp_tags.fastq \
    | samtools sort -@ 8 -o LR.sorted.bam -
samtools index LR.sorted.bam

# 3. Post-alignment MID deduplication per CID and gene (recommended; optional for HiFi)
python isoform_quantification/sp_tags/run_umi_dedup.py \
    -i LR.sorted.bam -g annotation.gtf -o LR.dedup.sam --umi_ed 4 -t 8

# 4. Quantify (barcodes are read from the CB / UB tags)
python isoform_quantification/main.py quantify \
    -gtf annotation.gtf -lrsam LR.dedup.sam \
    --st_mode --no_barcode_in_readname -t 16 -o quant_out
```

**ONT**: same steps, but run `sp_tags` with `--data_type nanopore` (see [ONT reads](#ont-reads-stereo-seq--nanopore)).

Notes:
- In cell / spot mode (`--barcode2spot`), CIDs that do not belong to any cell are dropped in step 1, so the aligned read count is lower than that of tools that keep every whitelisted CID.

See [Module: `sp_tags`](#module-sp_tags-high-resolution-spatial-barcodes-and-umis) for the processing steps, all arguments and output formats, [MID (UMI) handling](#mid-umi-handling) for the deduplication strategy, and [Full pipeline examples](#high-resolution-spatial-pipeline-stereo-seq-pacbio-hifi) for complete Stereo-seq scripts including `cal_TE`.

### TE annotation preparation (for `cal_TE`)

`cal_TE` requires a TE GTF (`-te_gtf`) built from a RepeatMasker `.out` file. Convert it once with `create_simple_TE_gtf.py`.

**`create_simple_TE_gtf.py`** — Converts RepeatMasker `.out` files to GTF format for TE annotations.

Purpose:
- Handles RepeatMasker `.out` column-width output format
- Generates a standardized GTF with TE hierarchy (`class_id` / `family_id` / `gene_id`=subfamily / `transcript_id`=locus name / `locus`)
- Calculates consensus TE length and truncation proportion from the `rep_start`/`rep_end`/`rep_left` columns

Command:

```bash
python isoform_quantification/TE_analyse/create_simple_TE_gtf.py \
    -i repeatmasker.fa.out \
    -o TE_annotation.gtf
```

Parameters:

| Argument | Description |
|---|---|
| `-i` / `--input` (required) | Input RepeatMasker `.out` file path |
| `-o` / `--output` (required) | Output GTF file path for TE annotations |

The resulting `TE_annotation.gtf` is the `-te_gtf` input for the [`cal_TE` subcommand](#subcommand-cal_te).

---

## Quick Start

### Single-cell mode

```bash
# Step 1–3: demux → dedup → align (see Data Preparation)

# Step 4: isoform quantification
python isoform_quantification/main.py quantify \
    -gtf annotation.gtf \
    -lrsam LR_dedup.sam \
    --sc_mode \
    --no_barcode_in_readname \
    --cb_tag CB --umi_tag UB \
    -t 16 \
    -o output/sc

# Step 5: TE quantification
python isoform_quantification/main.py cal_TE \
    -gtf annotation.gtf \
    -te_gtf TE_annotation.gtf \
    --sc_quant output/sc \
    -t 16 \
    -o output/sc/te
```

### Spatial transcriptomics mode

```bash
# Step 1–3: demux → dedup → align (see Data Preparation)

# Step 4: isoform quantification
python isoform_quantification/main.py quantify \
    -gtf annotation.gtf \
    -lrsam LR_dedup.sam \
    --st_mode \
    --no_barcode_in_readname \
    -t 16 \
    -o output/st

# Step 5: TE quantification
python isoform_quantification/main.py cal_TE \
    -gtf annotation.gtf \
    -te_gtf TE_annotation.gtf \
    --st_quant output/st \
    -t 16 \
    -o output/st/te
```

---

## Subcommand: `quantify`

```
python main.py quantify -gtf <GTF> -o <OUTPUT> [options]
```

### Required arguments

| Argument | Description |
|---|---|
| `-gtf` / `--gtf_annotation_path` | Path to transcript GTF annotation file |
| `-o` / `--output_path` | Output directory |

### Input data (required)

| Argument | Description |
|---|---|
| `-lrsam` / `--long_read_sam_path` | Path to long-read SAM file (genome-aligned) |

### General optional arguments

| Argument | Default | Description |
|---|---|---|
| `-t` / `--threads` | `1` | Number of parallel worker processes |
| `--EM_SR_num_iters` | `200` | Maximum EM iterations |
| `--isoform_start_end_site_tolerance` | `20` | Tolerance (bp) for matching LR read start/end to isoform boundaries |
| `--junction_site_tolerance` | `5` | Tolerance (bp) for matching splice junction sites |
| `--multi_mapping_filtering` | `best` | Multi-mapping filtering strategy: `best` (keep best-scoring alignment) or `unique_only` |

### Mode flags (mutually exclusive)

| Flag | Description |
|---|---|
| `--sc_mode` | Single-cell mode: extracts CB/UMI per read and outputs a cell × isoform count matrix |
| `--st_mode` | Spatial transcriptomics mode: extracts spot barcodes and outputs a spot × isoform count matrix |

### Single-cell mode arguments

| Argument | Default | Description |
|---|---|---|
| `--barcode_in_readname` | `True` | Barcode/UMI is embedded in the read name (flexiplex output) |
| `--no_barcode_in_readname` | — | Read barcode from SAM tags instead of the read name |
| `--cb_tag` | `CB` | SAM tag for cell barcode (used with `--no_barcode_in_readname`) |
| `--umi_tag` | `UB` | SAM tag for UMI (used with `--no_barcode_in_readname`) |

### Spatial transcriptomics mode arguments

| Argument | Default | Description |
|---|---|---|
| `--barcode_whitelist` | `None` | Path to spatial barcode whitelist file (one barcode per line); barcodes not in the list are discarded |

---

## Subcommand: `cal_TE`

Annotates each transcript by overlap with transposable elements and computes TE expression at locus / subfamily / family / class level for per-cell and/or per-spot data. The `-te_gtf` input is built from a RepeatMasker `.out` file — see [TE annotation preparation](#te-annotation-preparation-for-cal_te).

```
python main.py cal_TE -gtf <GTF> -te_gtf <TE_GTF> -o <OUTPUT> [options]
```

### Required arguments

| Argument | Description |
|---|---|
| `-gtf` / `--gtf_annotation_path` | Path to transcript GTF annotation file |
| `-te_gtf` / `--te_gtf_path` | Path to TE GTF annotation file |
| `-o` / `--output_path` | Output directory |

### General optional arguments

| Argument | Default | Description |
|---|---|---|
| `-t` / `--threads` | `1` | Number of parallel worker processes for the TE-overlap computation (transcripts are grouped by chromosome and distributed across workers; recommended value is the number of chromosomes or available CPU cores, whichever is smaller) |

### TE annotation thresholds

| Argument | Default | Description |
|---|---|---|
| `--first_exon_threshold` | `50.0` | First-exon TE proportion (%) threshold for TE-derived transcript detection |
| `--total_threshold` | `50.0` | Full-transcript TE proportion (%) threshold |
| `--te_overlap_threshold` | `10` | Minimum TE overlap length (bp) to distinguish gene-only from TE-containing transcripts |
| `--te_ratio_threshold` | `80.0` | TE proportion (%) in full transcript for TE-alone identification |
| `--te_feature_threshold` | `50` | TE overlap length (bp) within TSS/TES 200 bp window for sub-classification |

### Quantification inputs (at least one optional)

| Argument | Description |
|---|---|
| `--sc_quant` | Path to the quantification output directory from a prior miniTEQuant `quantify --sc_mode` run; triggers per-cell TE metrics, written to `<sc_quant>/SC_output/TE/` |
| `--st_quant` | Path to the quantification output directory from a prior miniTEQuant `quantify --st_mode` run; triggers per-spot TE metrics, written to `<st_quant>/ST_output/TE/` |
| `--custom_sc_quant` | Path to a per-cell isoform quantification directory produced by **any** tool (not necessarily miniTEQuant); triggers per-cell TE metrics, written to `<output_path>/Custom_SC_TE_output/` |
| `--custom_st_quant` | Path to a per-spot isoform quantification directory produced by **any** tool; triggers per-spot TE metrics, written to `<output_path>/Custom_ST_TE_output/` |

`--sc_quant` / `--st_quant` / `--custom_sc_quant` / `--custom_st_quant` all expect the same
directory layout: a standard 10x MEX triplet (`barcodes.tsv`, `features.tsv`, `matrix.mtx`),
found either directly in the given directory or under an `isoform/`, `SC_output/isoform/` or
`ST_output/isoform/` subdirectory. Each of the three MEX files may be plain text **or
gzip-compressed** (e.g. `barcodes.tsv.gz` — the default CellRanger layout); both are detected
automatically. Isoform IDs in `features.tsv` must match the `transcript_id` values produced
from `-gtf`. Matching tries an exact string match first, then falls back to comparing IDs with
any trailing Ensembl-style version suffix stripped (so `ENST00000448958` matches
`ENST00000448958.2`) — a version-only mismatch between `-gtf` and the quantification result is
therefore tolerated automatically. A low ID match rate is reported as a warning at runtime; if
it persists after the version-stripped fallback, it usually means `-gtf` does not match the
annotation used to produce the quantification result at all (e.g. novel transcripts from an
external tool's discovery mode, or a different transcript ID namespace such as RefSeq vs.
GENCODE).

Use `--custom_sc_quant` / `--custom_st_quant` instead of `--sc_quant` / `--st_quant` when the
isoform-level quantification was produced by a different tool — the TE aggregation logic is
identical, only the input provenance and output location differ (results are written under
`-o` instead of inside the input directory, so a read-only or third-party directory is never
written to).

**Example: TE quantification from another tool's per-cell results**

```bash
python main.py cal_TE \
    -gtf annotation.gtf \
    -te_gtf TE_annotation.gtf \
    --custom_sc_quant /path/to/other_tool/per_cell_isoform_mex_dir \
    -t 16 \
    -o output/te
```

Results are written to `output/te/Custom_SC_TE_output/`. `/path/to/other_tool/per_cell_isoform_mex_dir`
must contain `barcodes.tsv`/`features.tsv`/`matrix.mtx` (plain or `.gz`), with isoform IDs matching
`annotation.gtf`.

**Just want transcript structure information**

`--custom_sc_quant` (and the other `*_quant` options) is optional. If you only want to know how each
transcript in the annotation relates to TEs (Gene-alone / TE-Gene / TE-alone, TE overlap proportions,
etc.), drop `--custom_sc_quant` entirely:

```bash
python main.py cal_TE \
    -gtf annotation.gtf \
    -te_gtf TE_annotation.gtf \
    -t 16 \
    -o output/te
```

The per-transcript structure information is written to `output/te/transcript_quantification_with_TE.tsv`.

---

## Module: `sp_tags` (high-resolution spatial barcodes and UMIs)

`sp_tags` (in `isoform_quantification/sp_tags/`) turns Stereo-seq long reads (PacBio HiFi or ONT) into tagged reads for `quantify --st_mode`. The detection logic follows the Stereo-seq mode of [Spl-IsoQuant](https://github.com/ablab/IsoQuant). Run `python isoform_quantification/sp_tags/run_sp_tags.py -h` for all options.

### Stereo-seq read structure

```
PC1 primer | CID (25 nt) | linker | MID (10 nt) | polyT | cDNA (antisense) | TSO | PC1 primer | CID | ...
\__________________________ molecule 1 __________________________/  \___ molecule 2 ...
```

A read can contain several concatenated molecules, in either orientation; both strands are searched.

### Processing steps

1. **Detect each molecule**: polyT → linker (the anchor; no linker, no molecule) → PC1 primer → CID (between primer and linker) → MID (between linker and polyT) → TSO (end of the molecule). Adapters are found by k-mer seeding plus local alignment, so sequencing errors are tolerated.
2. **Correct the CID** against the whitelist (k-mer candidates, then alignment; accepted if score ≥ `--min_bc_score`, default CID length − 3).
3. **Split concatenated reads** (disable with `--no_split`); in multi-molecule reads only molecules closed by a TSO are written.
4. **Correct MIDs** within each CID (edit distance ≤ `--umi_ed`, default 1; see [MID (UMI) handling](#mid-umi-handling)).
5. **Write** the tagged reads, setting `CB` / `UB` according to the quantification unit below.

### Quantification unit: CID or spot

| Mode | `CB` | `UB` | Count matrix columns |
|---|---|---|---|
| default | CID | MID | one per CID (bin1) |
| `--barcode2spot FILE[:ID_COL:SPOT_COL]` | cell / spot of the CID | CID + MID | one per cell / spot |

- Molecules whose CID is not in the file (e.g. between cells) are dropped and counted as `Molecules dropped: spatial ID not in barcode2spot` in `sp_tags.stats.tsv`.

### Arguments

| Argument | Default | Description |
|---|---|---|
| `-i` / `--input` (required) | — | FASTQ / FASTA (optionally `.gz`) or unaligned BAM |
| `-o` / `--output_dir` (required) | — | Output directory |
| `-w` / `--whitelist` | — | CID whitelist: `CID [x y]` per line (`.gz` allowed). Strongly recommended; with coordinates, position files are also written |
| `--barcode2spot` | — | CID → cell / spot mapping (see above) |
| `--data_type` | `pacbio` | `pacbio` or `nanopore` (`ont`); `nanopore` uses k-mer 12 for the CID index |
| `--primer` / `--linker` / `--tso` | Stereo-seq defaults | Adapter sequences (5'→3', capture-oligo-first orientation); set only for a different kit |
| `--umi_ed` | `1` | Max edit distance for MID correction within a CID; `0` disables it |
| `--require_good_umi` | off | Drop molecules whose MID length is outside 10 ± 3 nt |
| `--min_bc_score` | CID length − 3 | Minimum alignment score to accept a whitelist CID |
| `--kmer_size` | `14` (`12` for nanopore) | K-mer size of the whitelist index |

### Output files

| File | Description |
|---|---|
| `sp_tags.fastq` | Tagged reads, one per molecule; input to minimap2 |
| `sp_tags.tags.tsv` | One row per molecule: read, CID, MID, CID score, MID validity, strand, adapter positions |
| `sp_tags.stats.tsv` | Detection, splitting, correction and filtering counts |
| `sp_tags.spatial_positions.csv` | CID mode with a coordinate whitelist: CID coordinates (`tissue_positions.csv` layout) |
| `sp_tags.spots.tsv` / `sp_tags.spot_positions.csv` | Spot mode: per-spot CIDs, reads and centroid |

Read header (tab-separated tags; minimap2 `-y` copies them into the SAM, then use `quantify --no_barcode_in_readname`):

```
@<read>_<start>_<end>_<strand>_<CID>_<MID>          CB:Z:<CID>   UB:Z:<MID>       CR:Z:<raw CID> UR:Z:<raw MID>
@<read>_<start>_<end>_<strand>_<spot>_<CID><MID>    CB:Z:<spot>  UB:Z:<CID><MID>  CR:Z:... UR:Z:... XC:Z:<CID>   # --barcode2spot
```

### MID (UMI) handling

| Stage | Tool | Scope | MID matching | Removes reads |
|---|---|---|---|---|
| 1. Correction | `sp_tags` (before alignment) | each CID | edit distance ≤ 1 | no |
| 2. Deduplication | `run_umi_dedup.py` (after alignment) | each CID **and gene** | edit distance ≤ 4 | yes, one read per molecule |
| 3. Counting | `quantify --st_mode` | each cell / CID, per gene | exact `CB` + `UB` | — |

- **Stage 1** is conservative: without gene information, a larger edit distance would merge distinct molecules of different genes in the same CID.
- **Stage 2** follows spl-IsoQuant's `UMIFilter` (code reimplemented). Each read is assigned to the gene it overlaps most (antisense strand preferred). MIDs within a (CID, gene) group are clustered greedily by read count. The kept read is chosen by unambiguous gene > more exons > more aligned bases > higher MAPQ. Grouping is by gene rather than isoform because truncated copies of one molecule may be assigned to different isoforms. Recommended for ONT, optional for PacBio HiFi.
- **Stage 3** is exact matching (`--umi_dedup_hamming` is not applied in ST / SC modes).

```bash
python isoform_quantification/sp_tags/run_umi_dedup.py \
    -i LR.sorted.bam -g annotation.gtf -o LR.dedup.sam [--umi_ed 4] [-t 8]
```

Input can be SAM, SAM.gz or BAM (no sorting needed). The statistics go to `<output>.umi_dedup.stats.tsv` (or `--stats`); pass the output SAM to `quantify -lrsam`.

### ONT reads (Stereo-seq + Nanopore)

Barcode detection is the same as for PacBio (as in spl-IsoQuant); three steps differ:

| Step | PacBio HiFi | ONT |
|---|---|---|
| `sp_tags` | `--data_type pacbio` | `--data_type nanopore` (CID k-mer 12) |
| minimap2 | `-ax splice:hq` | `-ax splice -k14` (build the index with `-x splice -k14`) |
| `run_umi_dedup.py` | optional | recommended |

### Notes

- **Speed / memory**: Large whitelists use about 2 GB plus ~56 bytes per CID, so crop the whitelist to the tissue region.
- **Kinnex / MAS-Seq**: run PacBio `skera split` first.
- **Coordinates**: `quantify` does not attach coordinates; join `sp_tags.spatial_positions.csv` / `sp_tags.spot_positions.csv` with `ST_output/isoform/barcodes.tsv` downstream.
- **Spot mode**: do not pass the CID whitelist to `quantify --barcode_whitelist`, because barcodes are now cell / spot names.

---

## Output Files

### `quantify --sc_mode` → `SC_output/`

```
SC_output/
├── isoform/
│   ├── barcodes.tsv         one cell barcode per line
│   ├── features.tsv         isoform_id  gene_id  Gene Expression
│   ├── matrix.mtx           sparse UMI count matrix  (isoforms × cells)
│   └── cpm_matrix.mtx       sparse CPM matrix        (isoforms × cells)
└── gene/
    ├── barcodes.tsv
    ├── features.tsv          gene_id  gene_id  Gene Expression
    ├── matrix.mtx            gene-level UMI count matrix (genes × cells)
    └── cpm_matrix.mtx        gene-level CPM matrix
```

All files use the standard MEX (Market Exchange Format), directly compatible with `Seurat::ReadMtx()` and `scanpy.read_mtx()`.

---

---

### `cal_TE` output

| Output | Condition | Description |
|---|---|---|
| `transcript_quantification_with_TE.tsv` | always | Per-transcript TE overlap annotation table (every transcript in `-gtf`) |
| `TE_derived_transcripts.tsv` | any transcripts meet the TE-derived thresholds | Subset of transcripts classified as TE-derived, with dominant TE info |
| `<sc_quant>/SC_output/TE/` | `--sc_quant` provided | Per-cell TE count tables, written **inside** the input quantification directory |
| `<st_quant>/ST_output/TE/` | `--st_quant` provided | Per-spot TE count tables, written **inside** the input quantification directory |
| `<output_path>/Custom_SC_TE_output/` | `--custom_sc_quant` provided | Per-cell TE count tables, for quantification results from any tool, written under `-o` |
| `<output_path>/Custom_ST_TE_output/` | `--custom_st_quant` provided | Per-spot TE count tables, for quantification results from any tool, written under `-o` |

**Per-cell / per-spot TE output files** (inside `SC_output/TE/`, `ST_output/TE/`,
`Custom_SC_TE_output/` or `Custom_ST_TE_output/` — `cell_*` for SC, `spot_*` for ST):

| File | Description |
|---|---|
| `{cell,spot}_transcript_classification_counts.tsv` | Per-cell/spot count breakdown by transcript classification (`Gene-alone` / `TE-alone` / `TE-Gene`) |
| `{cell,spot}_te_class_counts.tsv` | Per-cell/spot counts aggregated to TE class level (LINE, SINE, LTR, ...) |
| `{cell,spot}_te_family_counts.tsv` | Per-cell/spot counts aggregated to TE family level |
| `{cell,spot}_te_subfamily_counts.tsv` | Per-cell/spot counts aggregated to TE subfamily level |
| `{cell,spot}_te_loci_counts.tsv` (or `*_matrix.mtx` MEX triplet if very wide) | Per-cell/spot counts aggregated to individual TE locus — only with `--output_loci` |
| `README.txt` | Description of the files in this directory |

**TE annotation columns** (in `transcript_quantification_with_TE.tsv`):

| Column | Description |
|---|---|
| `transcript_id` | Transcript ID |
| `gene_id` | Gene ID |
| `transcript_TPM` | Transcript expression (TPM); not populated for SC/ST workflows (always `0`) |
| `transcript_length` | Transcript length (bp) |
| `num_TEs` | Number of distinct TE elements overlapping this transcript |
| `TE_class` | Deduplicated, semicolon-separated set of TE classes overlapping this transcript (unlike `class` below, each class appears at most once regardless of how many loci contribute it) |
| `transcript_classification` | TE classification: `Gene-alone`, `TE-alone`, `TE-Gene` |
| `transcript_sub_classification` | Sub-classification for `TE-Gene` transcripts, comma-separated if multiple: `TSS`, `TE-terminating`, `TE-exonization` |
| `first_exon_TE_proportion` | TE proportion in the first exon (%) |
| `total_transcript_TE_proportion` | TE proportion across the full transcript (%) |
| `tss_200_overlap` / `tes_200_overlap` | TE overlap length (bp) within the 200 bp TSS / TES window |
| `locus` / `subfamily` / `family` / `class` | Semicolon-separated TE locus / subfamily / family / class for each overlapping TE |
| `length` / `overlap_length` | Semicolon-separated TE length / overlap length (bp) for each overlapping TE |

---

## Full Pipeline Examples

### Single-cell pipeline (ONT, 10x 3' v3)

```bash
#!/bin/bash
SAMPLE="sample_name"
RAW_FQ="${SAMPLE}.fastq.gz"
REF="genome.fa"; INDEX="genome.mmi"; GTF="annotation.gtf"; TE_GTF="TE_annotation.gtf"
THREADS=16; OUT_DIR="output/sc"

# Step 1: barcode discovery & demultiplexing
gunzip -c "$RAW_FQ" | flexiplex -d 10x3v3 -f 0 -p $THREADS > flexiplex.out
flexiplex-filter flexiplex_barcodes_counts.txt > barcodes_whitelist.txt
gunzip -c "$RAW_FQ" | flexiplex -d 10x3v3 -p $THREADS -k barcodes_whitelist.txt > demux.fastq

# Step 2: UMI deduplication
nailpolish index --skip-unmatched demux.fastq
nailpolish summary demux.fastq
nailpolish consensus -t $THREADS demux.fastq > dedup.fastq && rm demux.fastq

# Step 3: alignment
[ ! -f "$INDEX" ] && minimap2 -d "$INDEX" "$REF"
minimap2 -y -ax splice --MD --secondary=no -G 200k -t $THREADS "$INDEX" dedup.fastq > LR.sam

# Step 4: isoform quantification
python isoform_quantification/main.py quantify \
    -gtf "$GTF" -lrsam LR.sam \
    --sc_mode --no_barcode_in_readname --cb_tag CB --umi_tag UB \
    -t $THREADS -o "$OUT_DIR"

# Step 5: TE quantification
python isoform_quantification/main.py cal_TE \
    -gtf "$GTF" -te_gtf "$TE_GTF" \
    --sc_quant "$OUT_DIR" \
    -t 16 \
    -o "$OUT_DIR/te"
```

### Spatial transcriptomics pipeline (Visium, ONT)

```bash
#!/bin/bash
SAMPLE="sample_name"
RAW_FQ="${SAMPLE}.fastq.gz"
REF="genome.fa"; INDEX="genome.mmi"; GTF="annotation.gtf"; TE_GTF="TE_annotation.gtf"
THREADS=16; OUT_DIR="output/st"

# Step 1: barcode discovery & demultiplexing
gunzip -c "$RAW_FQ" | flexiplex -d 10x3v3 -f 0 -p $THREADS > flexiplex.out
flexiplex-filter flexiplex_barcodes_counts.txt > barcodes_whitelist.txt
gunzip -c "$RAW_FQ" | flexiplex -d 10x3v3 -p $THREADS -k barcodes_whitelist.txt > demux.fastq

# Step 2: UMI deduplication
nailpolish index --skip-unmatched demux.fastq
nailpolish summary demux.fastq
nailpolish consensus -t $THREADS demux.fastq > dedup.fastq && rm demux.fastq

# Step 3: alignment
[ ! -f "$INDEX" ] && minimap2 -d "$INDEX" "$REF"
minimap2 -y -ax splice --MD --secondary=no -G 200k -t $THREADS "$INDEX" dedup.fastq > LR.sam

# Step 4: isoform quantification
python isoform_quantification/main.py quantify \
    -gtf "$GTF" -lrsam LR.sam \
    --st_mode --no_barcode_in_readname \
    -t $THREADS -o "$OUT_DIR"

# Step 5: TE quantification
python isoform_quantification/main.py cal_TE \
    -gtf "$GTF" -te_gtf "$TE_GTF" \
    --st_quant "$OUT_DIR" \
    -t 16 \
    -o "$OUT_DIR/te"
```

### High-resolution spatial pipeline (Stereo-seq, PacBio HiFi)

```bash
#!/bin/bash
SAMPLE="sample_name"
RAW_FQ="${SAMPLE}.hifi_reads.fastq.gz"     # or the unaligned HiFi BAM
CID_WHITELIST="chip_cid_whitelist.txt"     # CID [x y], exported from the chip mask
REF="genome.fa"; INDEX="genome.mmi"; GTF="annotation.gtf"; TE_GTF="TE_annotation.gtf"
THREADS=16; OUT_DIR="output/stereo"

# Step 1: split concatenated reads, identify and correct CID / MID (replaces flexiplex + nailpolish)
python isoform_quantification/sp_tags/run_sp_tags.py \
    -i "$RAW_FQ" -w "$CID_WHITELIST" \
    -o sp_tags_out -t $THREADS
# per cell instead of per CID: add  --barcode2spot cid2cell.tsv  (CID <tab> cell_id)

# Step 2: alignment (-y copies the CB/UB tags into the SAM)
[ ! -f "$INDEX" ] && minimap2 -d "$INDEX" "$REF"
minimap2 -y -ax splice:hq --MD --secondary=no -t $THREADS "$INDEX" sp_tags_out/sp_tags.fastq \
    | samtools sort -@ 8 -o LR.sorted.bam

# Step 3 (optional for HiFi): MID deduplication per CID and gene (edit distance 4, as spl-IsoQuant)
python isoform_quantification/sp_tags/run_umi_dedup.py \
    -i LR.sorted.bam -g "$GTF" -o LR.dedup.sam -t 8

# Step 4: isoform quantification (spot = CID and UMI = MID; with --barcode2spot: spot = cell, UMI = CID+MID)
# without Step 3, pass -lrsam LR.sorted.bam
python isoform_quantification/main.py quantify \
    -gtf "$GTF" -lrsam LR.dedup.sam \
    --st_mode --no_barcode_in_readname \
    -t $THREADS -o "$OUT_DIR"

# Step 5: TE quantification
python isoform_quantification/main.py cal_TE \
    -gtf "$GTF" -te_gtf "$TE_GTF" \
    --st_quant "$OUT_DIR" \
    -t 16 \
    -o "$OUT_DIR/te"

# Coordinates: sp_tags_out/sp_tags.spatial_positions.csv (per CID) or sp_tags.spot_positions.csv (with --barcode2spot);
# join on barcode with $OUT_DIR/ST_output/isoform/barcodes.tsv
```

### High-resolution spatial pipeline (Stereo-seq, ONT)

```bash
#!/bin/bash
SAMPLE="sample_name"
RAW_FQ="${SAMPLE}.ont_reads.fastq.gz"      # or the unaligned dorado BAM
CID_WHITELIST="chip_cid_whitelist.txt"     # CID [x y], exported from the chip mask
REF="genome.fa"; INDEX="genome.mmi"; GTF="annotation.gtf"; TE_GTF="TE_annotation.gtf"
THREADS=16; OUT_DIR="output/stereo_ont"

# Step 1: split concatenated reads, identify and correct CID / MID (ONT settings)
python isoform_quantification/sp_tags/run_sp_tags.py \
    -i "$RAW_FQ" -w "$CID_WHITELIST" --data_type nanopore \
    -o sp_tags_out -t $THREADS
# per cell instead of per CID: add  --barcode2spot cid2cell.tsv  (CID <tab> cell_id)

# Step 2: alignment with spl-IsoQuant ONT settings
[ ! -f "$INDEX" ] && minimap2 -d "$INDEX" "$REF"
minimap2 -y -ax splice --MD --secondary=no -G 200k -t $THREADS "$INDEX" sp_tags_out/sp_tags.fastq \
    | samtools sort -@ 8 -o LR.sorted.bam

# Step 3: UMI deduplication per CID and gene (edit distance 4, as spl-IsoQuant)
python isoform_quantification/sp_tags/run_umi_dedup.py \
    -i LR.sorted.bam -g "$GTF" -o LR.dedup.sam -t 8

# Step 4: isoform quantification
python isoform_quantification/main.py quantify \
    -gtf "$GTF" -lrsam LR.dedup.sam \
    --st_mode --no_barcode_in_readname \
    -t $THREADS -o "$OUT_DIR"

# Step 5: TE quantification
python isoform_quantification/main.py cal_TE \
    -gtf "$GTF" -te_gtf "$TE_GTF" \
    --st_quant "$OUT_DIR" \
    -t 16 \
    -o "$OUT_DIR/te"
```
