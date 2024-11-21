import os
import shutil
import subprocess
import pandas as pd
from pathlib import Path
from ete3 import NCBITaxa
from wipe.modules.utils import (
    check_required_cols,
    create_new_genomes_dir,
    run_command,
)

ncbi = NCBITaxa()


def get_taxonomy_ncbi(species_taxid):
    lineage = ncbi.get_lineage(species_taxid)
    names = ncbi.get_taxid_translator(lineage)
    ranks = ncbi.get_rank(lineage)

    # Define a dictionary to store the taxonomy ranks
    taxonomy_ranks = {
        "superkingdom": "d__",
        "phylum": "p__",
        "class": "c__",
        "order": "o__",
        "family": "f__",
        "genus": "g__",
        "species": "s__",
    }

    # Create a dictionary to store the taxonomy output with placeholders for missing ranks
    taxonomy_output = {
        rank: f"{taxonomy_ranks[rank]}" for rank in taxonomy_ranks
    }

    # Populate the taxonomy_output dictionary with the appropriate names
    for taxid in lineage:
        rank = ranks.get(taxid)
        if rank in taxonomy_output:
            taxonomy_output[rank] = f"{taxonomy_ranks[rank]}{names[taxid]}"

    # Generate the formatted taxonomy string
    formatted_taxonomy = ";".join(
        taxonomy_output[rank] for rank in taxonomy_ranks
    )

    return formatted_taxonomy


# unfinished
def recover_taxonomy_nonncbi(df, gsearch_db, nthreads, gsearch_outdir):
    """tmpdir is for storing the intermediate files generated from gsearch request."""
    check_required_cols(df, ["taxonomy", "genome_id", "genome_path"])

    glist = df["genome_path"].to_list()
    gdir = create_new_genomes_dir(glist)
    gsearch_neighbors = search_genomes(
        gsearch_db, 5, gdir, gsearch_outdir, nthreads=nthreads
    )
    gseach_anis = os.path.join(gsearch_outdir, "gsearch_results.tsv")
    gsearch_anis = reformat_search_results(gsearch_neighbors, gsearch_anis)
    # taxonomy_from_gsearch = process_gsearch_output(outdir)
    df_selected = process_gsearch_output(gseach_anis)
    # join with existing taxonomy
    df_selected["genome_id"] = ""

    taxonomy_from_gsearch = fill_taxonomy(df, df_tax)

    return taxonomy_from_gsearch


# unfinished
def recover_taxonomy(assembly, gsearch_db):
    check_required_cols(assembly, ["taxid"])

    assembly["taxonomy"] = assembly["species_taxid"].apply(get_taxonomy_ncbi)

    if assembly["taxonomy"].isna().sum() == 0:
        return assembly_with_tax

    if not gsearch_db:
        raise ValueError("Gsearch db needed.")

    assembly_with_tax = assembly.loc[assembly["taxonomy"].notna()]
    assembly_without_tax = assembly.loc[assembly["taxonomy"].isna()]

    assembly_without_tax_filled = recover_taxonomy_nonncbi(
        assembly_without_tax, reference_tax, gsearch_db
    )

    res = combine_tax_res(assembly_with_tax, assembly_without_tax_filled)

    return res


def create_db(indir, outdir, nthreads=1):
    """
    Output
    ------
    {outdir}/hnswdump.hnsw.data
    {outdir}/hnswdump.hnsw.graph
    {outdir}/parameters.json
    {outdir}/processing_state.json
    {outdir}/seqdict.json
    """
    Path(outdir).mkdir(parents=True, exist_ok=True)
    commands = [
        "gsearch",
        "--pio",
        "2000",
        "--nbthreads",
        str(nthreads),
        "tohnsw",
        "-d",
        indir,
        "-k",
        "16",
        "-s",
        "18000",
        "-n",
        "128",
        "--ef",
        "1600",
        "--algo",
        "optdens",
    ]
    run_command(commands, outdir)


def update_db(existing_db, indir, nthreads=1, backup_dir=None):
    if backup_dir:
        Path(backup_dir).mkdir(parents=True, exist_ok=True)
        shutil.copytree(existing_db, backup_dir, dirs_exist_ok=True)

    commands = [
        "gsearch",
        "--pio",
        "2000",
        "--nbthreads",
        str(nthreads),
        "add",
        "--hnsw",
        existing_db,
        "--new",
        indir,
    ]
    run_command(commands)


def search_genomes(db, n_neighbors, indir, outdir, nthreads=1):
    Path(outdir).mkdir(parents=True, exist_ok=True)
    commands = [
        "gsearch",
        "--pio",
        "2000",
        "--nbthreads",
        str(nthreads),
        "request",
        "--hnsw",
        db,
        "--nbanswers",
        str(n_neighbors),
        "-r",
        indir,
    ]
    run_command(commands, outdir)
    return f"{outdir}/gsearch.neighbors.txt"


def reformat_search_results(inpath, outpath, model=1):
    commands = ["reformat", "16", str(model), inpath, outpath]
    run_command(commands)


def process_gsearch_output(outdir):
    filepath = os.path.join(outdir, "gsearch_results.tsv")
    df = pd.read_csv(filepath, sep="\t", low_memory=True)
    df_selected = df.loc[df.groupby("Query_Name")["ANI"].idxmax()]


# unfinished
def process_search_results(inpath, assembly_with_tax, outpath):
    df_tax = pd.read_csv(assembly_with_tax, sep="\t", low_memory=True)
    check_required_cols(df_tax, ["genome_id", "taxonomy"])

    df = pd.read_csv(inpath, sep="\t", low_memory=True)
    df_result = df.loc[df.groupby("Query_Name")["ANI"].idxmax()]
    df_result["genome_id"] = df_result["Query_Name"].apply("some_function")
    df_result.merge(
        df_tax[["genome_id", "taxonomy"]],
        how="left",
        left_on="genome_id",
        right_on="genome_id",
    )

    # write df_result
    # write assembly_without_tax_filled

    print(df_result.head(5))


if __name__ == "__main__":
    pass
