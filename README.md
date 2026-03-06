# WIPE

**Web of life Integrated Processing Ecosystem (WIPE)** is a lightweight CLI toolkit for processing large genome collections (e.g., WoL-style genome indices). It provides batch utilities for genome QC, annotation, summary compilation, and reference/functional database helpers.

> Designed for HPC-style workflows where per-genome outputs are written to organized directories and later “compiled” into summary tables.

---

## Installation

### Option A: Install from source (editable)
```bash
git clone https://github.com/sherlyn99/wipe.git
cd wipe
pip install -e .
```

### Option B: Conda environment (recommended on clusters)

#### Create from `wipe.yml` (recommended)
```bash
git clone https://github.com/sherlyn99/wipe.git
cd wipe
conda env create -f wipe.yml
conda activate wipe
pip install -e .
```

---

## Quickstart

Check that the CLI is available:
```bash
wipe --help
```

---

## Command overview

WIPE is a `click` CLI with top-level commands plus two subgroups: `gsearch` and `uniref`.

### Top-level commands

#### `wipe metadata`
Process/merge metadata tables (behavior depends on which optional file paths are provided).
```bash
wipe metadata -m metadata.tsv -o OUTDIR
```

#### `wipe qc`
Run CheckM2 QC over genomes described by a metadata TSV.
```bash
wipe qc \
  -m metadata.tsv \
  -o OUTDIR \
  -db /path/to/CheckM2_database \
  -t 16
```

#### `wipe compile`
Compile per-genome outputs into summary tables.
```bash
wipe compile \
  -i RESULTS_DIR \
  -o OUTDIR \
  --checkm2 --linearization --kofamscan --barrnap --proteins \
  -c
```

#### `wipe linearize`
Batch genome linearization / standardization utilities.
```bash
wipe linearize -o OUTDIR
```

Optional examples:
```bash
wipe linearize -o OUTDIR --metadata metadata.tsv
wipe linearize -o OUTDIR --gap "N*20" # stitch contigs together with 20 'N's
wipe linearize -o OUTDIR --filt "plasmid,phage" # filter out plasmids and phages
```

#### `wipe annotate`
Batch 16S and phylogenetic marker gene annotation using KO profiles + KO list (and other configured resources).
```bash
wipe annotate \
  --metadata metadata.tsv \
  -o OUTDIR \
  --rrna-cutoff 0.67 \
  --ko-profiles ko_profiles.tsv \
  --ko-list ko_list.tsv \
  --tmp-dir /scratch/tmp \
  --nthreads 16
```

#### `wipe annotate-ko`
Per-genome KO annotation (only) helper driven by `metadata.tsv`.
```bash
wipe annotate-ko \
  --metadata metadata.tsv \
  -o OUTDIR \
  --ko-profiles ko_profiles.tsv \
  --ko-list ko_list.tsv \
  --tmp-dir /scratch/tmp \
  --nthreads 16
```

#### `wipe gen-dm`
Generate a distance matrix with Bindash2.
```bash
wipe gen-dm -i GENOMES_DIR -o OUTDIR --nthreads 16
```

#### `wipe functiondb`
DIAMOND-based functional annotation (proteins → UniRef DB).
```bash
wipe functiondb \
  -i PROTEIN_DIR \
  -db /path/to/uniref.dmnd \
  -o OUTDIR \
  -t 16
```

---

## `gsearch` subgroup

[Gsearch](https://github.com/jean-pierreBoth/gsearch) is an ultra-fast and scalable genome search / classification tool for finding the closest matches of query genomes against very large reference collections (hundreds of thousands of genomes).

### `wipe gsearch create`
Create a genome-search database.
```bash
wipe gsearch create -i GENOMES_DIR -o DB_OUTDIR -t 8
```

### `wipe gsearch update`
Update an existing genome-search database.
```bash
wipe gsearch update \
  -i NEW_GENOMES_DIR \
  -db EXISTING_DB_DIR \
  -o DB_OUTDIR \
  -t 8 \
  -b BACKUP_DIR
```

---

## `uniref` subgroup

This subgropu of commands allow users to generate uniref annotations of genomes, which can then be used to generate genome-function stratified tables using [Woltka](https://github.com/qiyunzhu/woltka). 

### `wipe uniref download`
Download UniRef from UniProt.
```bash
wipe uniref download --level 90 -o uniref90_download/
```

### `wipe uniref build`
Build a DIAMOND database from UniRef FASTA.
```bash
wipe uniref build -i uniref90.fasta.gz -o db/uniref90.dmnd -t 16
```

### `wipe uniref process`
Process UniRef XML → TSV.
```bash
wipe uniref process -i uniref90.xml.gz -o uniref90.tsv
```

### `wipe uniref blastp`
Run DIAMOND blastp against the UniRef DB.
```bash
wipe uniref blastp -i PROTEIN_DIR -db db/uniref90.dmnd -o OUTDIR -t 16
```

### `wipe uniref merge-maps`
Merge UniRef90 and UniRef50 mapping files.
```bash
wipe uniref merge-maps \
  --uniref90-map uniref90.map.tsv \
  --uniref50-map uniref50.map.tsv \
  -o merged.map.tsv \
  --simplify
```

### `wipe uniref extract-names`
Extract UniRef names for IDs in a mapping file.
```bash
wipe uniref extract-names \
  --map-file merged.map.tsv \
  --uniref90-names uniref90_names.tsv \
  --uniref50-names uniref50_names.tsv \
  -o merged.names.tsv
```

---

## Input conventions

Many commands assume a metadata TSV with at least:
- `genome_id`
- `lgenome_path` (path to genome FASTA, often `.fna.gz`)

Additional required columns depend on the workflow stage and enabled modules.

---

## Output layout

Most commands write:
- per-genome outputs under `OUTDIR/<genome_id>/...`
- compiled summary tables in `OUTDIR/`
- logs and “failed genome IDs” lists when applicable

