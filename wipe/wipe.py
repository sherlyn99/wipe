import os
import click
import pandas as pd
from os.path import join
from wipe.modules.constants import MSG_WELCOME
from wipe.modules.checkm2 import run_checkm2_batch
from wipe.modules.collection import compile_results
from wipe.modules.linearization import linearization_batch

from wipe.modules.metadata import process_metadata, merge_qc
from wipe.modules.gsearch import create_db, update_db
from wipe.modules.annotate import annotate_multiple, annotate_single_kofamscan
from wipe.modules.barrnap import run_barrnap_single
from wipe.modules.recover_16s import extract_16s_from_tlp
from wipe.modules.utils import load_metadata, run_command, check_required_cols
from wipe.modules.functiondb import run_functional_annotation
from wipe.modules.uniref import process_uniref_xml, merge_uniref_maps, extract_uniref_names
from wipe.modules.bindash import bindash2_gen_dm


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
@click.option("--checkm")
@click.option("--ani")
@click.option("--fcs-gbk")
@click.option("--fcs-rfs")
@click.option("--checkm2")
@click.option("--linearization")
@click.option("--kofamscan")
@click.option("--barrnap")
# fmt: on
def metadata(
    metadata,
    outdir,
    checkm,
    ani,
    fcs_gbk,
    fcs_rfs,
    checkm2,
    linearization,
    kofamscan,
    barrnap,
):
    """
    Check genome path and add citation count and taxonomy to the metadata file.
    """
    if kofamscan and barrnap:
        df_md = pd.read_csv(metadata, sep="\t", low_memory=False)
        df_kofamscan = pd.read_csv(kofamscan, sep="\t", low_memory=False)
        df_barrnap = pd.read_csv(barrnap, sep="\t", low_memory=False)
        df_tmp = pd.merge(df_md, df_kofamscan, on="genome_id", how="left")
        df_tmp = pd.merge(df_tmp, df_barrnap, on="genome_id", how="left")
        df_tmp.to_csv(
            os.path.join(outdir, "metadata_final.tsv"),
            sep="\t",
            header=True,
            index=False,
        )
        return
    if linearization:
        df_md = pd.read_csv(metadata, sep="\t", low_memory=False)
        df_lin = pd.read_csv(linearization, sep="\t", low_memory=False)
        df_combined = pd.merge(df_md, df_lin, on="genome_id", how="left")
        df_combined.to_csv(
            os.path.join(outdir, "metadata_linearization.tsv"),
            sep="\t",
            header=True,
            index=False,
        )
        return
    if checkm and ani and fcs_gbk and fcs_rfs and checkm2:
        md = merge_qc(metadata, checkm, ani, fcs_gbk, fcs_rfs, checkm2)
        md.to_csv(
            os.path.join(outdir, "metadata_qc.tsv"),
            sep="\t",
            index=False,
            header=True,
        )
        return
    # Check genome path and add citation count and taxonomy to the metadata file.
    process_metadata(metadata, outdir)


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
@click.option("--linearization", is_flag=True, default=None, help="Compile linearization results")
@click.option("--kofamscan", is_flag=True, default=None, help="Compile kofamscan results")
@click.option("--barrnap", is_flag=True, default=None, help="Compile barrnap results")
@click.option("--proteins", is_flag=True, default=None, help="Compile prodigal results")
# fmt: on
def compile(indir, outdir, coords, checkm2, linearization, kofamscan, barrnap, proteins):
    print("hi world")
    compile_results(
        indir,
        outdir,
        checkm2=checkm2,
        linearization=linearization,
        kofamscan=kofamscan,
        barrnap=barrnap,
        proteins=proteins,
        coords=coords,
    )


# fmt: off
@wipe.command()
@click.option("-m", "--metadata", type=click.Path(exists=True),
              help="Tab-delimited mapping between genome filenames (without extension) and genome id (optional)")
@click.option("-o", "--outdir", required=True,
              help="Output directory for multi-Fasta file")
@click.option("-g", "--gap", default=None,
              help='Fill sequence gaps with a string, use "*" to indicate repeats, e.g., "N*20"')
@click.option("-f", "--filt", default=None,
              help="Exclude sequences with any of the comma-delimited words in title, e.g., 'plasmid,phage'")
# fmt: on
def linearize(metadata, outdir, gap, filt):
    linearization_batch(metadata, outdir, gap, filt)


# fmt: off
@wipe.command()
@click.option('--metadata', required=True, type=click.Path(exists=True),
              help="Path to the metadata file (e.g., metadata_fqc.tsv).")
@click.option("-o", "--outdir", required=True,
              help="Directory storing compiled annotation_summary.txt.")
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
    metadata, outdir, rrna_cutoff, ko_profiles, ko_list, tmp_dir, nthreads
):
    """
    Annotates genomes based on provided metadata and KO profiles.
    """
    annotate_multiple(
        metadata, outdir, rrna_cutoff, ko_profiles, ko_list, tmp_dir, nthreads
    )


# fmt: off
@wipe.command()
@click.option('--metadata', required=True, type=click.Path(exists=True),
              help="Path to the metadata file (e.g., metadata_fqc.tsv).")
@click.option("-o", "--outdir", required=True,
              help="Directory storing compiled annotation_summary.txt.")
@click.option('--ko-profiles', required=True, type=click.Path(exists=True),
              help="Path to KO profiles file.")
@click.option('--ko-list', required=True, type=click.Path(exists=True),
              help="Path to KO list file.")
@click.option('--tmp-dir', required=True, type=click.Path(),
              help="Temporary directory for intermediate files.")
@click.option('--nthreads', required=True, type=int,
              help="Number of threads to use.")
# fmt: on
def annotate_ko(
    metadata, outdir, ko_profiles, ko_list, tmp_dir, nthreads
):
    """
    Annotates genomes based on provided metadata and KO profiles.
    """
    df_md = pd.read_csv(metadata, sep="\t", low_memory=False)
    check_required_cols(df_md, ["genome_id", "lgenome_path"])
    summary_path = os.path.join(outdir, "annotation_summary.txt")

    for row in df_md.itertuples():
        genome_id = row.genome_id
        in_file_fna_gz = row.lgenome_path
        outdir_per_genome = os.path.join(
            outdir,
            genome_id[0],
            genome_id[1:4],
            genome_id[4:7],
            genome_id[7:10],
        )
        try:
            annotate_single_kofamscan(
                genome_id,
                in_file_fna_gz,
                outdir_per_genome,
                ko_profiles,
                ko_list,
                nthreads,
                tmp_dir,
            )
        except:
            with open(summary_path, "a") as log:
                log.write(f"{genome_id}\n")


# fmt: off
@wipe.command()
@click.option("-i", "--indir", required=True,
              help="Directory storing input fa(.gz) files.")
@click.option("-o", "--outdir", required=True,
              help="Directory storing compiled annotation_summary.txt.")
@click.option('--nthreads', required=True, type=int,
              help="Number of threads to use.")
# fmt: on
def gen_dm(indir, outdir, nthreads):
    """
    Run bindash2 and geneate distance matrix.
    """
    bindash2_gen_dm(indir, outdir, nthreads)


# protocol selection
# - completion
# - contamination
# - filter and select based on marker gene content
# - bindash2 and selection

# tree building
# - marker gene extraction
# - multiple sequence alignment
# - build tree


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
@wipe.command()
@click.option("-i", "--indir", required=True, type=click.Path(exists=True),
              help="Input directory containing protein files")
@click.option("-db", "--diamonddb", required=True, type=click.Path(exists=True),
              help="Path to DIAMOND database")
@click.option("-o", "--outdir", required=True,
              help="Output directory for functional annotation results")
@click.option("-t", "--threads", default=4,
              help="Number of threads to use (default: 4)")
# fmt: on
def functiondb(indir, diamonddb, outdir, threads):
    """
    Run functional annotation pipeline using DIAMOND and UniRef database.
    
    This command expects:
    1. Input directory with protein files (all.faa)
    2. Pre-built DIAMOND database
    3. Output directory for results
    
    Produces:
    - uniref.map.xz: Mapping between query proteins and UniRef hits
    - diamond_results.m8: Raw DIAMOND output
    """
    run_functional_annotation(indir, diamonddb, outdir, threads)


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


@wipe.group()
def uniref():
    """Commands for working with UniRef databases and annotations."""
    pass


@uniref.command()
@click.option("--level", type=click.Choice(['50', '90']), required=True,
              help="UniRef database level to download (50 or 90)")
@click.option("-o", "--outdir", required=True,
              help="Output directory for downloaded files")
def download(level, outdir):
    """Download UniRef database files from UniProt FTP server."""
    try:
        os.makedirs(outdir, exist_ok=True)
        
        base_url = f"ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref{level}/"
        wget_cmd = [
            "wget", "-r", "-np", "-nH",
            "--cut-dirs=6",
            "-R", "index.html*",
            "-P", outdir,
            base_url
        ]
        
        click.echo(f"Downloading UniRef{level} database files...")
        run_command(wget_cmd)
        click.echo(f"UniRef{level} download complete")
        
    except Exception as e:
        click.echo(f"Error downloading UniRef{level}: {str(e)}", err=True)
        raise


@uniref.command()
@click.option("-i", "--input-fasta", required=True, type=click.Path(exists=True),
              help="Input FASTA file (UniRef50 or UniRef90)")
@click.option("-o", "--output-db", required=True,
              help="Output DIAMOND database path")
@click.option("-t", "--threads", default=4,
              help="Number of threads to use (default: 4)")
def build(input_fasta, output_db, threads):
    """Create a DIAMOND database from a UniRef FASTA file."""
    try:
        os.makedirs(os.path.dirname(output_db), exist_ok=True)
        diamond_cmd = [
            "diamond", "makedb",
            "--in", input_fasta,
            "--db", output_db,
            "--threads", str(threads)
        ]
        run_command(diamond_cmd)
    except Exception as e:
        click.echo(f"Error building DIAMOND database: {str(e)}", err=True)
        raise


@uniref.command()
@click.option("-i", "--xml-file", required=True, type=click.Path(exists=True),
              help="Input UniRef XML file (can be gzipped)")
@click.option("-o", "--output-tsv", required=True,
              help="Output TSV file path")
def process(xml_file, output_tsv):
    """Process UniRef XML file to extract relevant information."""
    try:
        process_uniref_xml(xml_file, output_tsv)
    except Exception as e:
        click.echo(f"Error processing UniRef XML: {str(e)}", err=True)
        raise


@uniref.command()
@click.option("-i", "--indir", required=True, type=click.Path(exists=True),
              help="Input directory containing protein files (all.faa)")
@click.option("-db", "--diamonddb", required=True, type=click.Path(exists=True),
              help="Path to DIAMOND database")
@click.option("-o", "--outdir", required=True,
              help="Output directory for DIAMOND results")
@click.option("-t", "--threads", default=4,
              help="Number of threads to use (default: 4)")
def blastp(indir, diamonddb, outdir, threads):
    """Run DIAMOND BLASTP against UniRef database."""
    try:
        run_functional_annotation(indir, diamonddb, outdir, threads)
    except Exception as e:
        click.echo(f"Error running DIAMOND BLASTP: {str(e)}", err=True)
        raise


@uniref.command()
@click.option("--uniref90-map", required=True, type=click.Path(exists=True),
              help="Path to UniRef90 mapping file")
@click.option("--uniref50-map", required=True, type=click.Path(exists=True),
              help="Path to UniRef50 mapping file")
@click.option("-o", "--output-file", required=True,
              help="Output merged map file path")
@click.option("--simplify", is_flag=True, default=False,
              help="Simplify UniRef IDs by taking only part after '_'")
def merge_maps(uniref90_map, uniref50_map, output_file, simplify):
    """Merge UniRef90 and UniRef50 mapping files."""
    try:
        merge_uniref_maps(uniref90_map, uniref50_map, output_file, simplify)
    except Exception as e:
        click.echo(f"Error merging UniRef maps: {str(e)}", err=True)
        raise


@uniref.command()
@click.option("--map-file", required=True, type=click.Path(exists=True),
              help="Path to UniRef mapping file")
@click.option("--uniref90-names", required=True, type=click.Path(exists=True),
              help="Path to UniRef90 names TSV file")
@click.option("--uniref50-names", required=True, type=click.Path(exists=True),
              help="Path to UniRef50 names TSV file")
@click.option("-o", "--output-names", required=True,
              help="Output names file path")
def extract_names(map_file, uniref90_names, uniref50_names, output_names):
    """Extract UniRef names for the IDs in the mapping file."""
    try:
        extract_uniref_names(map_file, uniref90_names, uniref50_names, output_names)
    except Exception as e:
        click.echo(f"Error extracting UniRef names: {str(e)}", err=True)
        raise


if __name__ == "__main__":
    wipe()
