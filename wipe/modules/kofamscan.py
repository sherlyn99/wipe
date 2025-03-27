import os
import shutil
import lzma
import pandas as pd
from pathlib import Path
from wipe.modules.utils import (
    run_command,
    load_metadata,
    check_required_cols,
    create_outdir,
    check_outputs,
    xz_compress_files,
)

# exec_annotation \
#   -f detail-tsv \
#   -p /projects/greengenes2/20231117_annotations_prelim/kofam_scan/profiles \
#   -k /projects/greengenes2/20231117_annotations_prelim/kofam_scan/ko_list \
#   -o kofamcan_out_G000965005/kofamscan_output_all.tsv \
#   prodigal_out_G000965005/G000965005.faa \
#   --no-report-unannotated \
#   --tmp-dir /panfs/y1weng/tmp_G000965005_all \
#   --cpu 64


def gen_command_kofamscan(
    input_file_faa,
    output_file_tsv,
    tmp_dir,
    ko_profiles="/projects/greengenes2/20231117_annotations_prelim/kofam_scan/profiles/",
    ko_list="/projects/greengenes2/20231117_annotations_prelim/kofam_scan/ko_list",
    nthreads=4,
):
    """
    Runs the exec_annotation command with the provided parameters.

    Args:
        outdir (str): Output directory path.
        filename (str): Name of the file to process.
        ko_number (str): KO number.
        condaenv (str): Conda environment name.

    Returns:
        None
    """
    # Construct the command
    command = [
        "exec_annotation",
        "-f",
        "detail-tsv",
        "-p",
        ko_profiles,
        "-k",
        ko_list,
        "-o",
        output_file_tsv,
        input_file_faa,
        "--no-report-unannotated",
        "--tmp-dir",
        tmp_dir,
        "--cpu",
        str(nthreads),
    ]
    return command


def filter_genome_output(output_file_tsv):
    headers = [
        "gene_name",
        "KO",
        "thrshld",
        "score",
        "E-value",
        "KO_definition",
    ]
    df = pd.read_csv(
        output_file_tsv, sep="\t", skiprows=2, header=None, names=headers
    ).reset_index(drop=True)

    # df = df[df["E-value"] < 1e-40].sort_values(by="KO")
    df = df.sort_values(by="KO") # comment out the e value filter
    df = df[df["score"] >= df["thrshld"]]
    df = df.loc[df.groupby(["KO"])["score"].idxmax()]

    return df


def process_genome_output(output_file_tsv, genome_id):
    headers = [
        "gene name",
        "KO",
        "thrshld",
        "score",
        "E-value",
        "KO definition",
    ]
    df = pd.read_csv(
        output_file_tsv, sep="\t", skiprows=1, header=None, names=headers
    )
    unique_ko_ct = df["KO"].nunique()

    genome_id = os.path.basename(output_file_tsv).split("_")[0]
    outdir = os.path.dirname(output_file_tsv)
    outpath = os.path.join(outdir, f"{genome_id}_marker_gene_ct.tsv.xz")
    with lzma.open(outpath, "wt") as f:
        f.write(f"genome_id\tmarker_gene_ct\n")
        f.write(f"{genome_id}\t{unique_ko_ct}\n")


def run_seqkit(infile, outfile, pattern, log_file):
    command = f"seqkit grep -p {pattern} {infile} > {outfile}"
    outfile = run_command(command, use_shell=True, logfile=log_file)
    return outfile


def extract_marker_genes(kofamscan_tsv_filtered, faa_xz, outdir, log_file):
    df = pd.read_csv(kofamscan_tsv_filtered, sep="\t", low_memory=False)
    for row in df.itertuples():
        ko = row.KO
        gene_name = row.gene_name
        outfile = os.path.join(outdir, f"{ko}.faa")
        run_seqkit(faa_xz, outfile, gene_name, log_file)


def run_kofamscan_single(
    genome_id,
    in_file_faa_xz,
    outdir,
    tmp_dir,
    ko_profiles,
    ko_list,
    nthreads,
):
    out_file_tsv = os.path.join(outdir, f"{genome_id}_kofamscan.tsv")
    out_file_tsv_xz = out_file_tsv + ".xz"
    out_file_tsv_filtered = out_file_tsv.replace(".tsv", "_filtered.tsv")
    out_file_tsv_filtered_xz = out_file_tsv_filtered + ".xz"
    log_file = os.path.join(outdir, f"{genome_id}_run_kofamscan.log")
    log_file_xz = log_file + ".xz"
    marker_gene_ct_xz = os.path.join(
        outdir, f"{genome_id}_marker_gene_ct.tsv.xz"
    )
    output_to_compress = [out_file_tsv, out_file_tsv_filtered, log_file]
    outputs_xz = [
        out_file_tsv_xz,
        out_file_tsv_filtered_xz,
        log_file_xz,
        marker_gene_ct_xz,
    ]

    # Skip if all outputs already exist
    if check_outputs(outputs_xz):
        return
    else:
        for o in outputs_xz:
            if os.path.exists(o):
                os.remove(o)

    # Create tmpdir for kofamscan
    tmp_dir = os.path.join(tmp_dir, f"{genome_id}_ko_tmp")
    Path(tmp_dir).mkdir(parents=True, exist_ok=True)

    # Decompress xz input file
    in_file_faa = os.path.join(
        tmp_dir, os.path.basename(in_file_faa_xz).replace(".xz", "")
    )
    with lzma.open(in_file_faa_xz, "rb") as compressed_file, open(
        in_file_faa, "wb"
    ) as decompressed_file:
        decompressed_file.write(compressed_file.read())

    if os.path.exists(in_file_faa):
        outdir = os.path.dirname(out_file_tsv)
        create_outdir(outdir)

        command = gen_command_kofamscan(
            in_file_faa,
            out_file_tsv,
            tmp_dir,
            ko_profiles,
            ko_list,
            nthreads,
        )
        run_command(command, logfile=log_file)

    if os.path.exists(out_file_tsv):
        df = filter_genome_output(out_file_tsv)
        df.to_csv(out_file_tsv_filtered, sep="\t", index=False, header=True)
        process_genome_output(out_file_tsv_filtered, genome_id)
        outdir_mg = os.path.join(outdir, "marker_genes")
        create_outdir(outdir_mg)
        extract_marker_genes(
            out_file_tsv_filtered, in_file_faa_xz, outdir_mg, log_file
        )
        xz_compress_files(output_to_compress)

    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)
