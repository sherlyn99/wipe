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

WIPE is a `click` CLI with top-level commands plus three subgroups: `gsearch`, `uniref`, and `functional-db`.

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

## `functional-db` subgroup

End-to-end helpers for downloading and running functional annotation databases (UniRef and EggNOG). The typical workflow is:

1. **Download** the databases once with `functional-db download`
2. **Annotate** a protein FASTA with `functional-db build`

### `wipe functional-db download`

Download UniRef90 + UniRef50 FASTAs and build their DIAMOND databases, and/or download the EggNOG mapper database. Both are enabled by default.

```bash
# download and build everything into dbs/ (default location)
wipe functional-db download

# use more threads for the diamond makedb step
wipe functional-db download -t 32

# download into a custom directory
wipe functional-db download -o /path/to/dbs -t 32

# download and build only UniRef
wipe functional-db download --no-eggnog -t 32

# download only EggNOG
wipe functional-db download --no-uniref
```

Output files are placed under:
- `<outdir>/uniref/` — `uniref90.fasta.gz`, `uniref90.dmnd`, `uniref50.fasta.gz`, `uniref50.dmnd`, `uniref90.xml.gz`, `uniref50.xml.gz`, `uniref90_names.tsv`, `uniref50_names.tsv`
- `<outdir>/eggnog/` — EggNOG mapper database files

### `wipe functional-db build`

Annotate a `.faa` protein file against UniRef and/or EggNOG. At least one of `--uniref` or `--eggnog` must be specified. DB paths default to where `functional-db download` places them (`dbs/uniref` and `dbs/eggnog`).

```bash
# annotate against both databases (using default db paths)
wipe functional-db build -i all.faa -o annotation_out/ --uniref --eggnog

# UniRef only, with a custom db path
wipe functional-db build \
  -i all.faa \
  -o annotation_out/ \
  --uniref \
  --uniref-db /path/to/uniref \
  -t 32

# EggNOG only
wipe functional-db build \
  -i all.faa \
  -o annotation_out/ \
  --eggnog \
  --eggnog-db /path/to/eggnog \
  -t 32
```

**UniRef annotation** (`--uniref`) runs DIAMOND blastp against both uniref90 and uniref50, then merges the hits (preferring uniref90), producing:
- `annotation_out/uniref_map.txt.xz` — tab-separated `ORF_ID <tab> UniRef_ID`
- `annotation_out/uniref_names.txt` — `UniRef_ID <tab> cluster name` (only if `uniref90_names.tsv` / `uniref50_names.tsv` are present in the DB dir, as generated by `functional-db download`)

**EggNOG annotation** (`--eggnog`) runs `emapper.py` and extracts per-annotation mapping files, producing:
- `annotation_out/eggnog_map.tsv` — full emapper annotations table (ORFs as rows)
- `annotation_out/orf_to_go.tsv` — ORF → GO term
- `annotation_out/orf_to_ec.tsv` — ORF → EC number
- `annotation_out/orf_to_ko.tsv` — ORF → KEGG KO
- `annotation_out/orf_to_cazy.tsv` — ORF → CAZy family
- `annotation_out/orf_to_cog.tsv` — ORF → COG category
- `annotation_out/orf_to_pfam.tsv` — ORF → Pfam domain

All mapping files are tab-separated with one ORF–annotation pair per line (multi-valued fields are expanded into multiple rows).

### Running the tests

```bash
pytest tests/test_functional_db.py -v
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

