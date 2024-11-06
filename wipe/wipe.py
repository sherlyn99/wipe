import click
from wipe.modules.constants import MSG_WELCOME
from wipe.modules.linearize import linearize_genomes

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
@click.option("-m", "--metadata", required=True, type=click.Path(exists=True),
              help="Tab-delimited mapping between genome filenames (without extension) and genome id (optional)")
@click.option("-e", "--ext", required=True,
              help="Filename extension following genome ID")
@click.option("-o", "--outdir", required=True,
              help="Output directory for multi-Fasta file")
@click.option("-g", "--gap", default=None,
              help='Fill sequence gaps with a string, use "*" to indicate repeats, e.g., "N*20"')
@click.option("-f", "--filt", default=None,
              help="Exclude sequences with any of the comma-delimited words in title, e.g., 'plasmid,phage'")
# fmt: on
def linearize(metadata, ext, outdir, gap, filt):
    linearize_genomes(metadata, ext, outdir, gap, filt)


## BASH
# wipe linearize -m <modified metadata>