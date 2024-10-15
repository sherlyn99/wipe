import click
from wipe.modules.constants import MSG_WELCOME

# takes in a directory of genomes
# linearize
# concatenate
# prodigal
# get coords.txt.gz
# get genome_length.gz
# get length.map.gz
# diamond
# get uniref.map.xz
# get uniref.names.xz


@click.group(help=MSG_WELCOME)
def wipe():
    pass
