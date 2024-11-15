import os
import gzip
import json
import subprocess
from datetime import datetime
from wipe.modules.utils import (
    gen_output_paths,
    extract_gid_from_inpath,
    decompress,
    xz_compress_files,
    check_log_and_retrieve_gid,
    get_files,
)


def gen_commands_prodigal(infna, gid, outdir):
    # NB: Prodiigal does not take .fna.gz as of 2024-11-11
    commands = [
        "prodigal",
        "-p",
        "single",
        "-g",
        "11",
        "-m",
        "-i",
        infna,
        "-f",
        "gff",
        "-o",
        f"{outdir}/{gid}.gff",
        "-a",
        f"{outdir}/{gid}.faa",
        "-d",
        f"{outdir}/{gid}.ffn",
    ]
    return commands


def get_files_prodigal(indir):
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


def gen_log_data_prodigal_all():
    log_data = {
        "process": "prodigal_run",
        "start_time": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "end_time": None,
        "status": "in_progress",
        "error": [],
    }
    return log_data


def gen_log_data_prodigal(inpath, infna, gid, outdir):
    log_data = {
        "genome_id": gid,
        "process": "prodigal_run",
        "start_time": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "end_time": None,
        "status": "in_progress",
        "details": {
            "input_file_compressed": inpath,
            "input_file": infna,
            "output_files": {
                "gff": f"{outdir}/{gid}.gff.xz",
                "faa": f"{outdir}/{gid}.faa.xz",
                "ffn": f"{outdir}/{gid}.ffn.xz",
            },
        },
        "error": None,
    }
    return log_data


def gen_log_filepath_prodigal(outdir, gid=None):
    if gid:
        return os.path.join(outdir, f"prodigal_stats_{gid}.json.gz")
    else:
        return os.path.join(outdir, "prodigal_stats.json.gz")


def gen_log_prodigal(inpath, infna, gid, outdir):
    log_filepath = gen_log_filepath_prodigal(outdir, gid)
    log_data = gen_log_data_prodigal(inpath, infna, gid, outdir)
    return log_filepath, log_data


def check_stats_prodigal(outdir, id=None):
    if id:
        logpath = os.path.join(outdir, f"prodigal_stats_{id}.json.gz")
    else:
        logpath = os.path.join(outdir, "prodigal_stats.json.gz")
    status, gid = check_log_and_retrieve_gid(logpath)
    return status, gid


def compress_outputs_prodigal(outdir_prodigal, gid):
    filelist = [
        os.path.join(outdir_prodigal, f"{gid}.gff"),
        os.path.join(outdir_prodigal, f"{gid}.faa"),
        os.path.join(outdir_prodigal, f"{gid}.ffn"),
    ]
    xz_compress_files(filelist)


def run_prodigal_single_genome(inpath, outdir, tmpdir, gid=None):
    """
    NB: tmpdir has to be an existing directory

    Example:
    run_prodigal_single_genome(
        "./linearized_genomes/M/001/002/003/M001002003.fna.gz",
        "./linearized_genomes/M/001/002/003/",
        "/tmp/run1"
    )

    Example Outputs:
        "./linearized_genomes/M/001/002/003/prodigal_stats.json.gz"
        "./linearized_genomes/M/001/002/003/prodigal_out/M001002003.gff.xz"
        "./linearized_genomes/M/001/002/003/prodigal_out/M001002003.faa.xz"
        "./linearized_genomes/M/001/002/003/prodigal_out/M001002003.ffn.xz"
    """
    outdir_path, stats_path, gid = gen_output_paths("prodigal", inpath, outdir)

    # Skip if the stats.json.gz says status == success
    status, retrieved_gid = check_stats_prodigal(outdir, gid)
    if status:
        return

    # Decompress input fna.gz
    infna = decompress(inpath, gid, tmpdir)

    # Create outdir/prodigal_out
    outdir_prodigal = os.path.join(outdir, "prodigal_out")
    os.makedirs(outdir_prodigal, exist_ok=True)

    # Generate prodigal commands and logs
    commands_prodigal = gen_commands_prodigal(infna, gid, outdir_prodigal)
    log_filepath, log_data = gen_log_prodigal(inpath, infna, gid, outdir)

    with gzip.open(log_filepath, "wt") as file:
        json.dump(log_data, file, indent=4, ensure_ascii=False)

    try:
        p = subprocess.run(
            commands_prodigal, capture_output=True, text=True, check=True
        )

        log_data["status"] = "success"

        # Compress the Prodigal output files using `xz`
        compress_outputs_prodigal(outdir_prodigal, gid)

    except subprocess.CalledProcessError as e:
        # Capture any errors encountered
        log_data["status"] = "failure"
        log_data["error"] = e.stderr  # Error message from Prodigal
    except FileNotFoundError as e:
        # Handle cases where Prodigal or the input file doesn't exist
        log_data["status"] = "failure"
        log_data["error"] = str(e)
    finally:
        # Update end time and write final log
        log_data["end_time"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        with gzip.open(log_filepath, "wt") as file:
            json.dump(log_data, file, indent=4, ensure_ascii=False)

        os.remove(infna)


def run_prodigal_multiple_genomes(indir, suffix, tmpdir, logdir):
    filelist = get_files(indir, suffix)
    log_data = gen_log_data_prodigal_all()
    log_filepath = os.path.join(logdir, "prodigal_summary.json.gz")

    with gzip.open(log_filepath, "wt") as file:
        json.dump(log_data, file, indent=4, ensure_ascii=False)

    for file in filelist:
        outdir_file = os.path.dirname(file)
        gid = extract_gid_from_inpath(file)
        run_prodigal_single_genome(file, outdir_file, tmpdir, gid)

        logpath = os.path.join(outdir_file, f"prodigal_stats_{gid}.json.gz")
        status, gid = check_log_and_retrieve_gid(logpath)
        if status == False:
            log_data["error"].append(gid)

        with gzip.open(log_filepath, "wt") as file:
            json.dump(log_data, file, indent=4, ensure_ascii=False)

    log_data["status"] = "completed"
    log_data["end_time"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with gzip.open(log_filepath, "wt") as file:
        json.dump(log_data, file, indent=4, ensure_ascii=False)
