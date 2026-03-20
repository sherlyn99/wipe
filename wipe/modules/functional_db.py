import os
import shutil
import tempfile
import click
from wipe.modules.utils import run_command
from wipe.modules.functiondb import merge_uniref


def download_uniref(outdir, threads=4):
    """
    Download UniRef90 and UniRef50 FASTA files from the UniProt FTP server
    and build DIAMOND databases from each.

    Args:
        outdir (str): Destination directory; files go into <outdir>/uniref/.
        threads (int): Number of threads for diamond makedb.
    """
    uniref_dir = os.path.join(outdir, "uniref")
    os.makedirs(uniref_dir, exist_ok=True)

    for level in ("90", "50"):
        url = (
            f"ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/"
            f"uniref/uniref{level}/uniref{level}.fasta.gz"
        )
        fasta = os.path.join(uniref_dir, f"uniref{level}.fasta.gz")
        dmnd = os.path.join(uniref_dir, f"uniref{level}.dmnd")

        click.echo(f"Downloading UniRef{level} FASTA...")
        run_command(["wget", "-P", uniref_dir, url])

        click.echo(f"Building UniRef{level} DIAMOND database...")
        run_command([
            "diamond", "makedb",
            "--in", fasta,
            "--db", dmnd,
            "--threads", str(threads),
        ])
        click.echo(f"UniRef{level} database ready.")


def annotate_uniref(faa, uniref_db_dir, outdir, threads):
    """
    Annotate a .faa file against UniRef90 and UniRef50 using DIAMOND blastp,
    then merge the hits into a single ORF -> UniRef ID map.

    Expects uniref_db_dir to contain uniref90.dmnd and uniref50.dmnd.
    Produces outdir/uniref_map.txt.xz.

    Args:
        faa (str): Path to input protein FASTA file.
        uniref_db_dir (str): Directory containing uniref90.dmnd and uniref50.dmnd.
        outdir (str): Output directory.
        threads (int): Number of threads for DIAMOND.
    """
    os.makedirs(outdir, exist_ok=True)

    diamond_flags = [
        "--index-chunks", "1",
        "--id", "90",
        "--subject-cover", "80",
        "--query-cover", "80",
        "--max-target-seqs", "1",
    ]

    m8_paths = {}
    with tempfile.TemporaryDirectory() as tmpdir:
        for level in ("90", "50"):
            db = os.path.join(uniref_db_dir, f"uniref{level}.dmnd")
            m8 = os.path.join(tmpdir, f"uniref{level}.m8")
            click.echo(f"Running DIAMOND blastp against UniRef{level}...")
            run_command([
                "diamond", "blastp",
                "--threads", str(threads),
                "--db", db,
                "--query", faa,
                "--out", m8,
                "--tmpdir", tmpdir,
                *diamond_flags,
            ])
            m8_paths[level] = m8

        merged = os.path.join(outdir, "uniref_map.txt")
        click.echo("Merging UniRef90 and UniRef50 hits...")
        merge_uniref(m8_paths["90"], m8_paths["50"], merged, simplify=True)

    click.echo("Compressing uniref_map.txt...")
    run_command(["xz", merged])
    click.echo(f"Done. Output: {merged}.xz")


# (output filename, 0-based column index, split character-by-character)
_EGGNOG_MAPPINGS = [
    ("orf_to_go.tsv",   9,  False),  # GOs          — comma-separated
    ("orf_to_ec.tsv",   10, False),  # EC            — comma-separated
    ("orf_to_ko.tsv",   11, False),  # KEGG_ko       — comma-separated
    ("orf_to_cazy.tsv", 18, False),  # CAZy          — comma-separated
    ("orf_to_cog.tsv",  6,  True),   # COG_category  — one letter per category
    ("orf_to_pfam.tsv", 20, False),  # PFAMs         — comma-separated
]


def make_eggnog_mappings(annotations_path, outdir):
    """
    Parse an emapper annotations file and write individual ORF-to-annotation
    mapping files for GO, EC, KO, CAZy, COG, and Pfam.

    Each output file contains two tab-separated columns: ORF ID and annotation
    value, with one row per ORF-value pair (multi-valued fields are expanded).

    Args:
        annotations_path (str): Path to the .emapper.annotations file.
        outdir (str): Output directory for the mapping files.
    """
    handles = {
        fname: open(os.path.join(outdir, fname), "w")
        for fname, _, _ in _EGGNOG_MAPPINGS
    }
    try:
        with open(annotations_path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                orf = fields[0]
                for fname, col, by_char in _EGGNOG_MAPPINGS:
                    if col >= len(fields):
                        continue
                    val = fields[col]
                    if not val or val == "-":
                        continue
                    values = list(val) if by_char else val.split(",")
                    for v in values:
                        v = v.strip()
                        if v and v != "-":
                            handles[fname].write(f"{orf}\t{v}\n")
    finally:
        for h in handles.values():
            h.close()

    click.echo(f"  Wrote {len(_EGGNOG_MAPPINGS)} mapping files to {outdir}.")


def annotate_eggnog(faa, eggnog_db_dir, outdir, threads):
    """
    Annotate a .faa file against the EggNOG database using emapper.py, then
    extract individual ORF-to-annotation mapping files.

    Produces in outdir:
      eggnog_map.tsv      — full emapper annotations table
      orf_to_go.tsv       — ORF -> GO term
      orf_to_ec.tsv       — ORF -> EC number
      orf_to_ko.tsv       — ORF -> KEGG KO
      orf_to_cazy.tsv     — ORF -> CAZy family
      orf_to_cog.tsv      — ORF -> COG category
      orf_to_pfam.tsv     — ORF -> Pfam domain

    Args:
        faa (str): Path to input protein FASTA file.
        eggnog_db_dir (str): Directory containing the EggNOG mapper database.
        outdir (str): Output directory.
        threads (int): Number of threads for emapper.
    """
    os.makedirs(outdir, exist_ok=True)

    click.echo("Running EggNOG mapper...")
    run_command([
        "emapper.py",
        "-i", faa,
        "--output", "eggnog",
        "--output_dir", outdir,
        "--cpu", str(threads),
        "--data_dir", eggnog_db_dir,
    ])

    annotations = os.path.join(outdir, "eggnog.emapper.annotations")
    dest = os.path.join(outdir, "eggnog_map.tsv")
    if os.path.exists(annotations):
        click.echo("Extracting per-annotation mappings...")
        make_eggnog_mappings(annotations, outdir)
        shutil.move(annotations, dest)
    click.echo(f"Done. Output: {dest}")


_EGGNOG_BASE_URL = "http://eggnog5.embl.de/download/emapperdb-5.0.2/"

_EGGNOG_FILES = [
    "eggnog.db.gz",
    "eggnog.taxa.tar.gz",
    "eggnog_proteins.dmnd.gz",
    "mmseqs.tar.gz",
    "pfam.tar.gz",
]


def _file_exists_nonempty(path):
    return os.path.exists(path) and os.path.getsize(path) > 0


def _extract(filepath, dest_dir):
    if filepath.endswith(".tar.gz"):
        run_command(["tar", "-xzf", filepath, "-C", dest_dir])
    elif filepath.endswith(".gz"):
        run_command(["gunzip", filepath])


def download_eggnog(outdir):
    """
    Download individual EggNOG mapper database files and extract them.

    For each file in the required set, checks whether it already exists and
    is non-empty in the destination directory. If not, downloads and extracts it.

    Args:
        outdir (str): Destination directory; files go into <outdir>/eggnog/.
    """
    eggnog_dir = os.path.join(outdir, "eggnog")
    os.makedirs(eggnog_dir, exist_ok=True)

    for filename in _EGGNOG_FILES:
        filepath = os.path.join(eggnog_dir, filename)
        if _file_exists_nonempty(filepath):
            click.echo(f"  {filename} already exists, skipping.")
            continue
        click.echo(f"  Downloading {filename}...")
        run_command(["wget", "-P", eggnog_dir, _EGGNOG_BASE_URL + filename])
        click.echo(f"  Extracting {filename}...")
        _extract(filepath, eggnog_dir)

    click.echo("EggNOG download complete.")
