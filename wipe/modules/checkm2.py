import os
import gzip
import json
import subprocess
from datetime import datetime
from wipe.modules.utils import get_files


def gen_command_checkm2(
    infile,
    threads,
    outdir,
    dbpath,
    genes=None,
):
    commands = [
        "checkm2",
        "predict",
        "-i",
        infile,
        "-o",
        outdir,
        "-t",
        str(threads),
        "--force",
        "--database_path",
        dbpath,
        "--remove_intermediates",
    ]
    if infile.endswith(".gz"):
        commands += ["-x", "gz"]
    if genes:
        commands += ["--genes", genes]

    return commands


def get_files_checkm2(indir):
    filelist_fa_gz = get_files(indir, "fa.gz")
    filelist_fna_gz = get_files(indir, "fna.gz")
    filelist_fasta_gz = get_files(indir, "fasta.gz")
    filelist_fa = get_files(indir, ".fa")
    filelist_fna = get_files(indir, ".fna")
    filelist_fasta = get_files(indir, ".fasta")
    filelist = (
        filelist_fa_gz
        + filelist_fna_gz
        + filelist_fasta_gz
        + filelist_fa
        + filelist_fna
        + filelist_fasta
    )
    return filelist


def gen_summary_path_checkm2(logdir):
    return os.path.join(logdir, "checkm2_summary.json.gz")


def gen_summary_data_checkm2():
    summary_data = {
        "process": "checkm2_run",
        "start_time": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "end_time": None,
        "status": "in_progress",
        "error": "no error",
    }
    return summary_data


def gen_summary_checkm2(logdir):
    summary_path = gen_summary_path_checkm2(logdir)
    summary_data = gen_summary_data_checkm2()
    return summary_path, summary_data


def run_checkm2_single(infile, outdir, dbpath, threads, genes=None):
    """
    Example
    -------
    run_checkm2_single("/data/001/002/003/G001002003.fa.gz",
                       "/data/001/002/003/checkm2_out",
                       4)
    """
    commands = gen_command_checkm2(infile, threads, outdir, dbpath, genes)
    if not os.path.exists(outdir):
        p = subprocess.run(
            commands, capture_output=True, text=True, check=True
        )


def run_checkm2_batch(indir, logdir, dbpath, threads):
    filelist = get_files_checkm2(indir)
    summary_path, summary_data = gen_summary_checkm2(logdir)

    with gzip.open(summary_path, "wt") as file:
        json.dump(summary_data, file, indent=4, ensure_ascii=False)

    for file in filelist:
        outdir_single = os.path.join(os.path.dirname(file), "checkm2_out")
        try:
            run_checkm2_single(file, outdir_single, dbpath, threads)
        except Exception as e:
            if isinstance(summary_data["error"], list):
                summary_data["error"].append(f"{file}: {e}")
            else:
                summary_data["error"] = [f"{file}: {e}"]

    summary_data["status"] = "completed"
    summary_data["end_time"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with gzip.open(summary_path, "wt") as file:
        json.dump(summary_data, file, indent=4, ensure_ascii=False)
