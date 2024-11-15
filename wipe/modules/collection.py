import os
import lzma
import skbio
import pandas as pd
from glob import glob
from wipe.modules.utils import search_dirs


def compile_results_checkm2(indir):
    dirs = search_dirs(indir, "checkm2_out*")
    dfs = []

    for dir in dirs:
        quality_report_path = os.path.join(dir, "quality_report.tsv")
        if os.path.exists(quality_report_path):
            df = pd.read_csv(quality_report_path, sep="\t")
            dfs.append(df)

    combined_df = pd.concat(dfs, ignore_index=True)
    combined_df = combined_df.rename(columns={"Name": "genome_id"})
    return combined_df


def compile_results_prodigal(indir, coords=True):
    dirs = search_dirs(indir, "prodigal_out*")
    dfs = []

    for dir in dirs:
        proteins_fps = glob(f"{dir}/proteins_*.tsv.xz")
        if len(proteins_fps) > 1:
            raise ValueError(
                f"There is one more protein count file in the directory: {dir}"
            )
        proteins_fp = proteins_fps[0]
        if os.path.exists(proteins_fp):
            proteins_df = pd.read_csv(proteins_fp, compression="xz", sep="\t")
            dfs.append(proteins_df)

    combined_df = pd.concat(dfs, ignore_index=True)
    combined_df = combined_df.rename(columns={"Name": "genome_id"})

    return combined_df


def generate_coords(indir, outfile_fp):
    dirs = search_dirs(indir, "prodigal_out*")

    with lzma.open(outfile_fp, mode="wb") as coords_fo:
        for dir in dirs:
            coords_fps = glob(f"{dir}/coords_*.txt.xz")
            if len(coords_fps) > 1:
                raise ValueError(
                    f"There is one more protein count file in the directory: {dir}"
                )
            coords_fp = coords_fps[0]

            with lzma.open(coords_fp, mode="rb") as f:
                coords_fo.write(f.read())


def concatenate_xz_files(file_list, output_file):
    with lzma.open(output_file, mode="wb") as out_file:
        for file in file_list:
            with lzma.open(file, mode="rb") as f:
                out_file.write(f.read())
