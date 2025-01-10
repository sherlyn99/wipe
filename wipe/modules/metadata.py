import os
import requests
import pandas as pd
from ete3 import NCBITaxa
from wipe.modules.utils import get_files, check_required_cols, create_outdir

ncbi = NCBITaxa()


def generate_gids(start_genome_id, num_ids):
    prefix = start_genome_id[0]
    number_part = int(start_genome_id[1:])

    genome_ids = [
        f"{prefix}{str(number_part + i).zfill(len(start_genome_id) - 1)}"
        for i in range(num_ids)
    ]

    return genome_ids


def generate_metadata(indir, ext, gid_start=None):
    files = get_files(indir, ext)
    ext = "." + ext.lstrip(".")
    files_without_ext = [f.replace(ext, "") for f in files]

    if not gid_start:
        return pd.DataFrame(files_without_ext, columns=["filepath"])

    n = len(files)
    gids = generate_gids(gid_start, n)
    return pd.DataFrame({"filepath": files_without_ext, "genome_id": gids})


def check_genome_paths(df, output_file):
    if "genome_path" not in df.columns:
        raise ValueError(
            "The DataFrame does not contain a column named 'genome_path'."
        )
    non_existent_paths = []
    for path in df["genome_path"]:
        if not os.path.exists(path):
            non_existent_paths.append(path)
    with open(output_file, "w") as f:
        for path in non_existent_paths:
            f.write(f"{path}\n")
    print(f"Non-existent paths written to {output_file}.")


def get_citation_count(pmid):
    url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:{pmid}&format=json"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        if (
            "resultList" in data
            and "result" in data["resultList"]
            and len(data["resultList"]["result"]) > 0
        ):
            citation_count = data["resultList"]["result"][0].get(
                "citedByCount", 0
            )
            return citation_count
    return 0


def process_citations(df, col_name):
    citation_sums = []
    for pubmed_ids in df[col_name]:
        if pd.isna(pubmed_ids) or pubmed_ids == "na":
            citation_sums.append(0)
            continue
        ids = pubmed_ids.strip(";").split(";")
        citation_ct = sum(
            get_citation_count(pmid.strip()) for pmid in ids if pmid
        )
        citation_sums.append(citation_ct)
    return citation_sums


def get_taxonomy_ncbi(species_taxid):
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

    if pd.isna(species_taxid):
        formatted_taxonomy = ";".join(
            taxonomy_ranks[rank] for rank in taxonomy_ranks
        )

    else:
        lineage = ncbi.get_lineage(species_taxid)
        names = ncbi.get_taxid_translator(lineage)
        ranks = ncbi.get_rank(lineage)

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


def process_metadata(inpath, outdir):
    md = pd.read_csv(inpath, sep="\t", low_memory=False)
    create_outdir(outdir)
    check_genome_paths(md, os.path.join(outdir, "nonexistent_path.txt"))
    md["citation"] = process_citations(md, "pubmed_id")
    md = recover_taxonomy_ncbi(md)
    md.to_csv(
        os.path.join(outdir, "metadata.tsv"),
        sep="\t",
        header=True,
        index=False,
    )


def merge_qc(metadata, checkm, ani, fcs_gbk, fcs_rfs, checkm2):
    # dm
    md = pd.read_csv(metadata, sep="\t", low_memory=False)

    # checkm
    df_checkm = pd.read_csv(checkm, sep="\t", low_memory=False).rename(
        columns={"#genbank-accession": "genbank-accession"}
    )
    df_checkm["assembly_accession"] = df_checkm["refseq-accession"].where(
        df_checkm["refseq-accession"].notna(), df_checkm["genbank-accession"]
    )
    df_checkm = df_checkm[
        ["assembly_accession", "checkm-completeness", "checkm-contamination"]
    ]

    # ani
    df_ani = pd.read_csv(ani, sep="\t", low_memory=False).rename(
        columns={"# genbank-accession": "genbank-accession"}
    )
    df_ani["assembly_accession"] = df_ani["refseq-accession"].where(
        df_ani["refseq-accession"] != "na", df_ani["genbank-accession"]
    )
    df_ani = df_ani[["assembly_accession", "best-match-status"]]

    # fcs
    df_fcs_gbk = pd.read_csv(
        fcs_gbk, sep="\t", low_memory=False, compression="gzip"
    ).rename(columns={"#assembly_accession": "assembly_accession"})
    df_fcs_rfs = pd.read_csv(
        fcs_rfs, sep="\t", low_memory=False, compression="gzip"
    ).rename(columns={"#assembly_accession": "assembly_accession"})
    df_fcs_rfs["genbank_assembly_accession"] = df_fcs_rfs[
        "assembly_accession"
    ].str.replace("GCF", "GCA")
    new_rows = df_fcs_gbk[
        ~df_fcs_gbk["assembly_accession"].isin(
            df_fcs_rfs["genbank_assembly_accession"]
        )
    ]
    df_fcs_rfs = df_fcs_rfs.drop("genbank_assembly_accession", axis=1)
    df_fcs = pd.concat([df_fcs_rfs, new_rows])
    df_fcs = df_fcs[
        [
            "assembly_accession",
            "total_contaminant_length",
            "fraction_contaminated",
        ]
    ]

    # checkm2
    df_checkm2 = pd.read_csv(checkm2, sep="\t", low_memory=False)

    # merge
    md_tmp = pd.merge(md, df_checkm, on="assembly_accession", how="left")
    md_tmp = pd.merge(md_tmp, df_checkm2, on="genome_id", how="left")
    md_tmp = pd.merge(md_tmp, df_ani, on="assembly_accession", how="left")
    md_tmp = pd.merge(md_tmp, df_fcs, on="assembly_accession", how="left")

    return md_tmp
