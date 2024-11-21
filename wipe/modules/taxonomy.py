import pandas as pd
from ete3 import NCBITaxa
from wipe.modules.utils import (
    check_required_cols,
    create_new_genomes_dir,
)
from wipe.modules.gsearch import search_genomes, reformat_search_results

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


def recover_taxonomy_ncbi(df):
    check_required_cols(df, ["genome_id", "species_taxid"])
    df["taxonomy"] = df["species_taxid"].apply(get_taxonomy_ncbi)
    return df


def get_closest_neighbors(gdir, tmpdir, gsearch_db):
    """tmpdir is for storing the intermediate files generated from gsearch request."""
    neighbors_filepath = search_genomes(gsearch_db, 5, gdir, tmpdir)
    anis_filepath = reformat_search_results(neighbors_filepath)
    df = pd.read_csv(anis_filepath, sep="\t", low_memory=True)
    df_selected = df.loc[df.groupby("Query_Name")["ANI"].idxmax()]
    return df_selected


def get_taxonomy_nonncbi(df_selected, df_tax):
    df_selected_with_tax = df_selected.merge(
        df_tax[["genome_id", "taxonomy"]],
        how="left",
        left_on="genome_id",
        right_on="genome_id",
    )
    return df_selected_with_tax


def recover_taxonomy_nonncbi(df, df_tax, tmpdir, gsearch_db):
    check_required_cols(df, ["genome_id", "genome_path"])
    glist = df["genome_path"].to_list()
    gdir = create_new_genomes_dir(glist)
    df_selected = get_closest_neighbors(gdir, tmpdir, gsearch_db)
    df_selected_with_tax = get_taxonomy_nonncbi(df_selected, df_tax)


def recover_taxonomy(assembly, ref_tax, gsearch_db, tmpdir):
    """tmpdir is for storing the intermediate files generated from gsearch request."""
    taxid_col = "species_taxid"
    assembly_ncbi = assembly[assembly[taxid_col].notna()]
    assembly_nonncbi = assembly[assembly[taxid_col].isna()]

    assembly_ncbi_with_tax = recover_taxonomy_ncbi(assembly_ncbi)
    if ref_tax:
        assembly_nonncbi_with_tax = recover_taxonomy_nonncbi(
            assembly_nonncbi, ref_tax, tmpdir, gsearch_db
        )
    else:
        assembly_nonncbi_with_tax = recover_taxonomy_nonncbi(
            assembly_nonncbi, assembly_ncbi_with_tax, tmpdir, gsearch_db
        )

    res = pd.concat([assembly_ncbi_with_tax, assembly_nonncbi_with_tax])
    return res


if __name__ == "__main__":
    pass
