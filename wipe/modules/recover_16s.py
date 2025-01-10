import os
import pandas as pd
from wipe.modules.utils import run_command


def extract_header_list(species_name):
    inpath = "./wipe/database/LTP_10_2024.csv"
    df_ltp = pd.read_csv(
        inpath,
        sep="\t",
        encoding="cp1252",
        names=["header", "name", "taxonomy", "c_col", "type_strain"],
    )

    header_list = df_ltp[df_ltp["name"].str.contains(species_name, na=False)][
        "header"
    ].values
    return species_name, header_list


def run_seqkit(species_name, header_list, outdir):
    infile = "./wipe/database/LTP_10_2024_compressed.fasta"
    outfile = os.path.join(
        outdir, f"{'_'.join(species_name.split(' '))}_16s_from_tlp.fasta"
    )
    command = f"seqkit grep -p {','.join(header_list)} {infile} > {outfile}"
    outfile = run_command(command, use_shell=True)
    return outfile


def extract_16s_from_tlp(speices_name, outdir):
    species_name, header_list = extract_header_list(speices_name)
    if len(header_list) == 0:
        print(f"No headers found for species: {species_name}")
        return None
    outfile = run_seqkit(species_name, header_list, outdir)
    return outfile


if __name__ == "__main__":
    name = "Shigella flexneri"
    outdir = "tests/tmp/barrnap_test"
    outfile = extract_16s_from_tlp(name, outdir)
