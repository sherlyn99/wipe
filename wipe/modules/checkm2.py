import os
import gzip
import json
import subprocess
from pathlib import Path
from datetime import datetime
from wipe.modules.utils import (
    gen_output_paths,
    gen_path,
    gen_stats_data,
    gen_summary_data,
    check_outputs,
    check_log_and_retrieve_gid,
    get_files_all_fa,
)


def gen_command_checkm2(
    inpath,
    threads,
    outdir,
    dbpath,
    genes=None,
):
    """NB: --genes parameter not used currently"""
    commands = [
        "checkm2",
        "predict",
        "-i",
        inpath,
        "-o",
        outdir,
        "-t",
        str(threads),
        "--force",
        "--database_path",
        dbpath,
        "--remove_intermediates",
    ]
    if inpath.endswith(".gz"):
        commands += ["-x", "gz"]
    if genes:
        commands += ["--genes", genes]

    return commands


def gen_stats_data_checkm2(inpath, outdir_path):
    log_data = gen_stats_data("checkm2_run", inpath)
    log_data["details"]["output_files"][
        "checkm2_report"
    ] = f"{outdir_path}/quality_report.tsv"
    return log_data


def gen_summary_path_checkm2(logdir):
    """
    Example
    -------
    gen_stats_path_checkm2("/data/G/")

    Outputs:
    -------
    "/data/G/checkm2_summary.json.gz"
    """
    return gen_path(logdir, "checkm2_summary.json.gz")


def gen_summary_data_checkm2():
    summary_data = gen_summary_data("checkm2_run")
    return summary_data


def gen_summary_checkm2(logdir):
    summary_path = gen_summary_path_checkm2(logdir)
    summary_data = gen_summary_data_checkm2()
    return summary_path, summary_data


def check_outputs_checkm2(outdir_path, stats_path):
    """Check existence for key output files and log status"""
    output_fps = [
        outdir_path,
        os.path.join(outdir_path, "quality_report.tsv"),
        stats_path,
    ]
    file_existence = check_outputs(output_fps)
    status, gid = check_log_and_retrieve_gid(stats_path)
    return file_existence and status


def run_checkm2_single(inpath, outdir, dbpath, threads, genes=None):
    """
    NB: Filenames should be in the format of info1_info2.fna(.gz)

    Example
    -------
    run_checkm2_single("/data/001/002/003/G001002003.fa.gz",
                       "/data/001/002/003/",
                       4)

    Example Outputs
    ----------------
    "/data/001/002/003/checkm2_out_G001002003"
    "/data/001/002/003/checkm2_stats_G001002003.json.gz"
    """
    outdir_path, stats_path, _ = gen_output_paths("checkm2", inpath, outdir)

    # skip if already completed
    if check_outputs_checkm2(outdir_path, stats_path):
        return

    Path(outdir_path).mkdir(parents=True, exist_ok=True)
    commands = gen_command_checkm2(inpath, threads, outdir_path, dbpath, genes)
    stats_data = gen_stats_data_checkm2(inpath, outdir_path)

    with gzip.open(stats_path, "wt") as file:
        json.dump(stats_data, file, indent=4)

    try:
        p = subprocess.run(
            commands, capture_output=True, text=True, check=True
        )
        stats_data["status"] = "success"

    except subprocess.CalledProcessError as e:
        stats_data["status"] = "failure"
        stats_data["error"] = e.stderr  # Error message from Prodigal
    except FileNotFoundError as e:
        stats_data["status"] = "failure"
        stats_data["error"] = str(e)
    finally:
        stats_data["end_time"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        with gzip.open(stats_path, "wt") as file:
            json.dump(stats_data, file, indent=4)


def run_checkm2_batch(indir, logdir, dbpath, threads):
    """
    Example
    -------
    run_checkm2_batch("./tests/data/999/",
                      "./tests/data/",
                      "/home/y1weng/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd",
                      4)
    """
    filelist = get_files_all_fa(indir)
    summary_path, summary_data = gen_summary_checkm2(logdir)

    with gzip.open(summary_path, "wt") as file:
        json.dump(summary_data, file, indent=4)

    for file in filelist:
        outdir = os.path.dirname(file)
        try:
            run_checkm2_single(file, outdir, dbpath, threads)
        except Exception as e:
            if isinstance(summary_data["error"], list):
                summary_data["error"].append(f"{file}: {e}")
            else:
                summary_data["error"] = [f"{file}: {e}"]

    summary_data["status"] = "completed"
    summary_data["end_time"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with gzip.open(summary_path, "wt") as file:
        json.dump(summary_data, file, indent=4, ensure_ascii=False)
