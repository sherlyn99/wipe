import os
import pandas as pd
from pathlib import Path
from wipe.modules.utils import run_command, check_outputs, xz_compress_files


def run_barrnap(
    in_file_fna_gz,
    out_file_gff,
    out_file_faa,
    reject_cutoff=0.67,
    log_file=None,
):
    """
    Example
    --------
    run_barrnap("SMGC_10.fna", "SMGC_10.gff", "SMGC_10.faa", reject_cutoff=0.67)
    """
    command = (
        f"zcat {in_file_fna_gz} | "
        f"barrnap --reject {reject_cutoff} -o {out_file_faa} > {out_file_gff}"
    )
    run_command(command, use_shell=True, logfile=log_file)


def process_gff(gff_file, out_file_tsv):
    """
    NB
    --
    barrnap only takes in decompressed file.
    --reject <N> N > 0

    Example
    -------
    process_gff("./SMGC_10.gff", "SMGC_10_barrnap_out.tsv")
    """
    gff_data = pd.read_csv(
        gff_file,
        sep="\t",
        skiprows=1,
        header=None,
        names=[
            "seqid",
            "source",
            "type",
            "start",
            "end",
            "score",
            "strand",
            "phase",
            "attributes",
        ],
    )
    gff_data["length"] = (gff_data["end"] - gff_data["start"]).abs()
    gff_data[["rrna_name", "product"]] = gff_data["attributes"].str.extract(
        r"Name=([^;]+);product=(.+)"
    )
    gff_data = gff_data.drop(columns=["attributes"])
    gff_data.to_csv(out_file_tsv, sep="\t", header=True, index=False)


def run_barrnap_single(
    in_file_fna_gz,
    outdir,
    reject_cutoff=0.67,
):
    """
    Example
    -------
    run_barrnap_single(
        "SMGC_10.fa.gz",
        "H/001/002/003/barrnap_out",
        .67
    )
    """
    file_stem = os.path.basename(in_file_fna_gz).split(".")[0]
    out_file_gff = os.path.join(outdir, f"{file_stem}_barrnap.gff")
    out_file_faa = os.path.join(outdir, f"{file_stem}_barrnap.faa")
    out_file_tsv = os.path.join(outdir, f"{file_stem}_barrnap.tsv")
    log_file = os.path.join(outdir, f"{file_stem}_run_barrnap.log")
    outputs = [out_file_gff, out_file_faa, out_file_tsv, log_file]
    outputs_xz = [o + ".xz" for o in outputs]

    # Skip if all outputs already exist
    if check_outputs(outputs_xz):
        return
    else:
        for o in outputs_xz:
            if os.path.exists(o):
                os.remove(o)

    Path(outdir).mkdir(parents=True, exist_ok=True)
    run_barrnap(
        in_file_fna_gz, out_file_gff, out_file_faa, reject_cutoff, log_file
    )
    process_gff(out_file_gff, out_file_tsv)
    xz_compress_files(outputs)

    return out_file_tsv


# def run_barrnap_batch(filelist, outdir, reject_cutoff):
#     Path(outdir).mkdir(parents=True, exist_ok=True)
#     for file in filelist:
#         file_name = os.path.basename(file)
#         file_stem = file_name.split(".")[0]
#         out_file_gff = os.path.join(outdir, f"{file_stem}_barrnap.gff")
#         out_file_faa = os.path.join(outdir, f"{file_stem}_barrnap.faa")
#         out_file_tsv = os.path.join(outdir, f"{file_stem}_barrnap.tsv")
#         run_barrnap_single(
#             file, out_file_gff, out_file_faa, out_file_tsv, reject_cutoff
#         )
