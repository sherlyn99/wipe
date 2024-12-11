import os
import lzma
import gzip
import json
import pandas as pd
from glob import glob
from wipe.modules.utils import search_dirs


def compile_results_linearization(indir):
    dirs = search_dirs(indir)
    data = []
    dfs = []

    for dir in dirs:
        quality_report_path = os.path.join(dir, "linearization.tsv")
        if os.path.exists(quality_report_path):
            df = pd.read_csv(quality_report_path, sep="\t")
            dfs.append(df)

    combined_df = pd.concat(dfs, ignore_index=True)
    return combined_df


def compile_results_checkm2(indir):
    dirs = search_dirs(indir)
    dfs = []
    genome_ids = []

    for dir in dirs:
        quality_report_path = os.path.join(dir, "quality_report.tsv")
        if os.path.exists(quality_report_path):
            df = pd.read_csv(quality_report_path, sep="\t")
            dfs.append(df)

            gid = "".join(dir.split("/")[-6:-2])
            genome_ids.append(gid)

    combined_df = pd.concat(dfs, ignore_index=True)
    combined_df.insert(0, "genome_id", genome_ids)
    return combined_df


def compile_results_kofamscan(indir):
    dirs = search_dirs(indir, "kofamscan_out")
    dfs = []

    for dir in dirs: 
        quality_report_path = search_dirs(dir, "*marker_gene_ct.tsv.xz")[0]
        if os.path.exists(quality_report_path):
            df = pd.read_csv(
                quality_report_path, sep="\t", compression="xz"
            )
        else:
            gid = "".join(dir.split("/")[-5:-1])
            df = pd.DataFrame({"genome_id": gid, "marker_gene_ct": 0})
        dfs.append(df)

    combined_df = pd.concat(dfs, ignore_index=True)
    return combined_df


def compile_results_barrnap(indir):
    dirs = search_dirs(indir, "barrnap_out")
    dfs = []
    gids = []

    for dir in dirs:
        quality_report_path = search_dirs(dir, "*barrnap.tsv.xz")[0]
        if os.path.exists(quality_report_path):
            df_tmp = pd.read_csv(
                quality_report_path, sep="\t", compression="xz"
            )
            counts = df_tmp["rrna_name"].value_counts()
            desired_order = ["5S_rRNA", "16S_rRNA", "23S_rRNA"]
            counts = counts.reindex(desired_order).fillna(0)
            df = pd.DataFrame(counts, dtype=int).T.reset_index(drop=True)
            df.columns.name = None
        else:
            columns = ["5S_rRNA", "16S_rRNA", "23S_rRNA"]
            df = pd.DataFrame([[0] * len(columns)], columns=columns, dtype=int)

        dfs.append(df)
        gid = "".join(dir.split("/")[-5:-1])
        gids.append(gid)
    combined_df = pd.concat(dfs, ignore_index=True)
    combined_df.insert(0, "genome_id", gids)
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


def compile_results(
    indir,
    outdir,
    checkm2=None,
    linearization=None,
    kofamscan=None,
    proteins=None,
    barrnap=None,
    coords=False,
):
    outpath_res = os.path.join(outdir, "metadata_stats.tsv")
    if checkm2:
        df_checkm2 = compile_results_checkm2(indir)
        out_file = os.path.join(outdir, "results_checkm2.tsv")
        df_checkm2.to_csv(out_file, index=False, header=True, sep="\t")

    if linearization:
        df_lin = compile_results_linearization(indir)
        out_file = os.path.join(outdir, "results_linearization.tsv")
        df_lin.to_csv(out_file, index=False, header=True, sep="\t")

    if kofamscan:
        df_kofamscan = compile_results_kofamscan(indir)
        out_file = os.path.join(outdir, "results_kofamscan.tsv")
        df_kofamscan.to_csv(out_file, index=False, header=True, sep="\t")

    if barrnap:
        df_barrnap = compile_results_barrnap(indir)
        out_file = os.path.join(outdir, "results_barrnap.tsv")
        df_barrnap.to_csv(out_file, index=False, header=True, sep="\t")

    if proteins:
        df_proteins = compile_results_prodigal(indir)
        out_file = os.path.join(outdir, "results_proteins.tsv")
        df_proteins.to_csv(out_file, index=False, header=True, sep="\t")

        if coords:
            outpath_coords = os.path.join(outdir, "coords.txt.xz")
            generate_coords(indir, outpath_coords)

    # df = pd.merge(df_lin, df_qc, how="left", on="genome_id").merge(
    #     df_proteins, how="left", on="genome_id"
    # )

    # df.to_csv(outpath_res, index=False, header=True, sep="\t")
