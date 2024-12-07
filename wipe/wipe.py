import os
import click
from os.path import join
from wipe.modules.constants import MSG_WELCOME
from wipe.modules.linearize import linearize_genomes
from wipe.modules.metadata import generate_metadata
from wipe.modules.prodigal import run_prodigal_batch
from wipe.modules.checkm2 import run_checkm2_batch
from wipe.modules.collection import compile_results
from wipe.modules.gsearch import create_db, update_db
from wipe.modules.annotate import annotate_multiple
from wipe.modules.barrnap import run_barrnap_single
from wipe.modules.recover_16s import extract_16s_from_tlp
from wipe.modules.utils import load_metadata


# takes in a directory of genomes
# linearize
# prodigal
# concatenate
# get coords.txt.gz
# get genome_length.gz
# get length.map.gz
# diamond
# get uniref.map.xz
# get uniref.names.xz


@click.group(help=MSG_WELCOME)
def wipe():
    pass


# fmt: off
@wipe.command()
@click.option("-m", "--metadata", 
              type=click.Path(exists=True), required=True,
              help="Path to metadata file containing genome information.")
@click.option("-o", "--outdir", required=True, help="Output directory for storing results.")
@click.option("-db", "--dbpath", 
              default="/home/y1weng/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd",
              help="Path to the CheckM2 database.")
@click.option("-t", "--threads", default=4, help="Number of threads to use.")
# fmt: on
def qc(metadata, outdir, dbpath, threads):
    """
    Run CheckM2 QC on a batch of genomes based on the metadata file.
    """
    run_checkm2_batch(metadata, outdir, dbpath, threads)


# fmt: off
@wipe.command()
@click.option("-i", "--indir", required=True, type=click.Path(exists=True),
              help="Input directory containing stats files.")
@click.option("-o", "--outdir", required=True,
              help="Directory storing compiled results.")
@click.option("-c", "--coords", is_flag=True, default=False,
              help="Genearte coords.txt.xz files.")
@click.option("--checkm2", is_flag=True, default=None, help="Compile checkm2 results")
# fmt: on
def compile(indir, outdir, coords, checkm2):
    compile_results(indir, outdir, checkm2=checkm2, coords=coords)


# look at distribution of contamination v.s. completion
# filter genomes

# annotate 16s

# annotate marker genes

# protocol selection
# - completion
# - contamination
# - filter and select based on marker gene content
# - bindash2 and selection

# tree building
# - marker gene extraction
# - multiple sequence alignment
# - build tree


### old code


@wipe.command()
@click.option("--metadata", type=click.Path(exists=True), required=True)
@click.option("--outdir", required=True)
def annotate_16s(metadata, outdir):
    df_md = load_metadata(metadata)
    infile = df_md["lin_path"]  # this needs to be verified
    name = "_".join(
        df_md["organism_name"].str.split(" ").str[:2]
    )  # this needs to be tested
    if metadata:
        out_file_tsv = run_barrnap_single(infile, outdir)

        if out_file_tsv and os.path.getsize(out_file_tsv) == 0:
            extract_16s_from_tlp(out_file_tsv, name, outdir)


# fmt: off
@wipe.command()
@click.option("-m", "--metadata", type=click.Path(exists=True),
              help="Tab-delimited mapping between genome filenames (without extension) and genome id (optional)")
@click.option("-e", "--ext",
              help="Filename extension following genome ID")
@click.option("-o", "--outdir", required=True,
              help="Output directory for multi-Fasta file")
@click.option("-g", "--gap", default=None,
              help='Fill sequence gaps with a string, use "*" to indicate repeats, e.g., "N*20"')
@click.option("-f", "--filt", default=None,
              help="Exclude sequences with any of the comma-delimited words in title, e.g., 'plasmid,phage'")
@click.option("-a", "--assembly", type=click.Path(exists=True))
# fmt: on
def linearize(metadata, ext, outdir, gap, filt, assembly):
    linearize_genomes(metadata, ext, outdir, gap, filt, assembly)


# fmt: off
@wipe.command()
@click.option('--metadata', required=True, type=click.Path(exists=True),
              help="Path to the metadata file (e.g., metadata_fqc.tsv).")
@click.option('--lin-dir', required=True, type=click.Path(exists=True),
              help="Directory containing linearized genomes.")
@click.option('--rrna-cutoff', type=float, default=0.67,
              help="rRNA cutoff value.")
@click.option('--ko-profiles', required=True, type=click.Path(exists=True),
              help="Path to KO profiles file.")
@click.option('--ko-list', required=True, type=click.Path(exists=True),
              help="Path to KO list file.")
@click.option('--tmp-dir', required=True, type=click.Path(),
              help="Temporary directory for intermediate files.")
@click.option('--nthreads', required=True, type=int,
              help="Number of threads to use.")
# fmt: on
def annotate(
    metadata, lin_dir, rrna_cutoff, ko_profiles, ko_list, tmp_dir, nthreads
):
    """
    Annotates genomes based on provided metadata and KO profiles.
    """
    annotate_multiple(
        metadata, lin_dir, rrna_cutoff, ko_profiles, ko_list, tmp_dir, nthreads
    )


# fmt: off
@wipe.command()
@click.option("-i", "--indir", required=True, type=click.Path(exists=True),
              help="Input directory containing all genome files.")
@click.option("-e", "--ext", required=True,
              help="Filename extension. e.g. 'fna' or '.fna'.")
@click.option("-o", "--outdir", required=True,
              help="Output directory for metadata.")
@click.option("-s", "--start-gid", required=False,
              help="Start genome gids, if they need to be specified.")
# fmt: on
def metadata(indir, ext, outdir, start_gid):
    md_df = generate_metadata(indir, ext, start_gid)
    md_df.to_csv(
        join(outdir, "metadata.tsv"),
        sep="\t",
        index=False,
        header=False,
    )


# # fmt: off
# @wipe.command()
# @click.option("-i", "--indir", required=True, type=click.Path(exists=True),
#               help="Input directory containing all genome files.")
# @click.option("-log", "--log_dir", required=True,
#               help="E.g. ./tests/data/out.")
# @click.option("-tmp", "--tmp_dir", required=True,
#               help="Tmp dir storing decompressed genome data.")
# # fmt: on
# def annotate(indir, log_dir, tmp_dir):
#     run_prodigal_batch(indir, log_dir, tmp_dir)


@wipe.group()
def gsearch():
    pass


# fmt: off
@gsearch.command()
@click.option("-i", "--indir", required=True, type=click.Path(exists=True),
              help="Input directory containing fasta or fna.gz files.")
@click.option("-o", "--outdir", required=True,
              help="Directory storing produced gsearch db files.")
@click.option("-t", "--nthreads", default=1, help="Number of threads to use (default is 1).")
# fmt: on
def create(indir, outdir, nthreads):
    create_db(indir, outdir, nthreads)


# fmt: off
@gsearch.command()
@click.option("-i", "--indir", required=True, type=click.Path(exists=True),
              help="Input directory containing new fasta or fna.gz files.")
@click.option("-db", "--existing-db", required=True, type=click.Path(exists=True),
              help="Path to the directory containing existing gsearch db.") 
@click.option("-o", "--outdir", required=True,
              help="Directory storing new gsearch db files.")
@click.option("-t", "--nthreads", default=1, help="Number of threads to use (default is 1).")
@click.option("-b", "--backup-dir", help="Directory for backing up the existing db files.")
# fmt: on
def update(indir, db, outdir, nthreads, backup_dir):
    update_db(db, indir, nthreads=nthreads, backup_dir=backup_dir)


if __name__ == "__main__":
    wipe()
