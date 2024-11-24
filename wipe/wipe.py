import click
from os.path import join
from wipe.modules.constants import MSG_WELCOME
from wipe.modules.linearize import linearize_genomes
from wipe.modules.metadata import generate_metadata
from wipe.modules.prodigal import run_prodigal_batch
from wipe.modules.checkm2 import run_checkm2_batch
from wipe.modules.collection import compile_results
from wipe.modules.gsearch import create_db, update_db

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
@click.option("-i", "--indir", type=click.Path(exists=True),
              help="Input directory containing fa or fa.gz files.")
@click.option("-l", "--logdir", required=True,
              help="Directory storing ./checkm2_summary.json.gz")
@click.option("-db", "--dbpath", default="/home/y1weng/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd",
              help="Path to the checkm2 database.")
@click.option("-t", "--threads", default=4)
@click.option("-m", "--metadata", type=click.Path(exists=True))
@click.option("-o", "--outdir")
# fmt: on
def qc(indir, logdir, dbpath, threads, metadata, outdir):
    run_checkm2_batch(
        indir, logdir, dbpath, threads, md=metadata, outdir=outdir
    )


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


# fmt: off
@wipe.command()
@click.option("-i", "--indir", required=True, type=click.Path(exists=True),
              help="Input directory containing all genome files.")
@click.option("-log", "--log_dir", required=True,
              help="E.g. ./tests/data/out.")
@click.option("-tmp", "--tmp_dir", required=True,
              help="Tmp dir storing decompressed genome data.")
# fmt: on
def annotate(indir, log_dir, tmp_dir):
    run_prodigal_batch(indir, log_dir, tmp_dir)


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
