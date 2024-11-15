import os
import bz2
import json
import gzip
import lzma
import shutil
import subprocess
import pandas as pd
from glob import glob
from pathlib import Path
from os.path import splitext
from datetime import datetime


def convert_bash_commands_to_subprocess_commands(bash_cmd):
    return bash_cmd.replace(" \\ ", " ").split(" ")


def read(fname):
    """Helper to read compressed files"""
    # compressed filename pattern
    zipdic = {".gz": gzip, ".bz2": bz2, ".xz": lzma, ".lz": lzma}
    ext = splitext(fname)[1]
    zipfunc = getattr(zipdic[ext], "open") if ext in zipdic else open
    return zipfunc(fname, "rt")


def load_md(md_path):
    md_df = pd.read_csv(md_path, sep="\t", header=None, low_memory=False)
    if md_df.shape[1] == 2:
        md_df.columns = ["filepath", "genome_id"]
    elif md_df.shape[1] == 1:
        md_df.columns = ["filepath"]
    else:
        raise ValueError(
            f"Unexpected number of columns in metadata {md_path}. Expected 1 or 2 columns."
        )
    return md_df


def check_duplicated_genome_ids(md_df):
    duplicaed_gids = md_df[md_df["genome_id"].duplicated()]
    if not duplicaed_gids.empty:
        raise ValueError(f"Duplicated genome ids found:\n{duplicaed_gids}")


def write_json_log(log_data, outdir, filename="log.json.gz", append=True):
    """
    Writes log data to a JSON file.

    Parameters:
    - log_data (dict): The data to log, typically in dictionary format.
    - outdir (str): Directory where the log file will be saved.
    - filename (str): The name of the log file (default is 'log.json').
    - append (bool): If True, appends to an existing file. If False, overwrites the file.
    """

    os.makedirs(outdir, exist_ok=True)
    filepath = os.path.join(outdir, filename)
    log_data["timestamp"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    if append and os.path.exists(filepath):
        with gzip.open(filepath, "rt") as file:
            try:
                existing_data = json.load(file)
                if not isinstance(existing_data, list):
                    existing_data = [existing_data]
                existing_data.append(log_data)
            except json.JSONDecodeError:
                existing_data = [log_data]
    else:
        existing_data = [log_data]

    with gzip.open(filepath, "wt") as file:
        json.dump(existing_data, file, indent=4, ensure_ascii=False)


def extract_gid_from_inpath(inpath):
    inpath = inpath.replace(".gz", "")
    inpath = inpath.replace(".xz", "")
    gid = splitext(os.path.split(inpath)[-1])[0]
    return gid


def decompress(inpath, gid, tmpdir):
    outpath = os.path.join(tmpdir, f"{gid}.fna")
    with gzip.open(inpath, "rb") as fi:
        with open(outpath, "wb") as fo:
            shutil.copyfileobj(fi, fo)
    return outpath


def check_outputs(filelist):
    return all(os.path.exists(file) for file in filelist)


def xz_compress_files(filelist):
    for file in filelist:
        subprocess.run(["xz", "-z", file], check=True)


def check_log_and_retrieve_gid(filepath):
    if os.path.isfile(filepath):
        with gzip.open(filepath, "rt") as file:
            log_data = json.load(file)
        return log_data.get("status") == "success", log_data.get("genome_id")
    else:
        return False, None


def get_files(indir, suffix):
    suffix = "." + suffix.lstrip(".")  # make sure suffix starts with '.'
    filepaths = Path(indir).rglob(f"*{suffix}")
    abs_filepaths = [
        os.path.abspath(f) for f in filepaths if "/.ipynb" not in str(f)
    ]
    return abs_filepaths


def get_files_all_fa(indir):
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


def search_dirs(start_dir, basename_pattern):
    """Does recursive search for directory/file.
    Example
    -------
    search_dirs("./tests/data/, "checkm2_out*)
    """
    dir_pattern = os.path.join(start_dir, "**", basename_pattern)
    dirs = glob(dir_pattern, recursive=True)
    return dirs


def gen_output_paths(process, inpath, outdir):
    gid = os.path.basename(inpath).split(".")[0]

    outdir_path = os.path.join(outdir, f"{process}_out_{gid}")
    stats_path = os.path.join(outdir, f"{process}_stats_{gid}.json.gz")
    return outdir_path, stats_path, gid


def gen_path(pdir, filename):
    return os.path.join(pdir, filename)


def gen_stats_data(process, inpath):
    gid = os.path.basename(inpath).split(".")[0]

    stats_data = {
        "genome_id": gid,
        "process": process,
        "start_time": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "end_time": None,
        "status": "in_progress",
        "details": {
            "input_file": inpath,
            "output_files": {},
        },
        "error": "no error",
    }
    return stats_data


def gen_summary_data(process):
    log_data = {
        "process": process,
        "start_time": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "end_time": None,
        "status": "in_progress",
        "error": "no error",
    }
    return log_data


# def setup_logger_with_stream(output_dir):
#     logger = logging.getLogger("wipe")
#     logger.setLevel(logging.INFO)

#     # Define formatter
#     formatter = logging.Formatter(
#         "%(name)s-%(asctime)s-%(levelname)s: %(message)s"
#     )

#     # File handler
#     file_handler = logging.FileHandler(
#         f"{output_dir}/mopp_workflow_{timestamp}.log"
#     )
#     file_handler.setFormatter(formatter)

#     # Stream handler
#     stream_handler = logging.StreamHandler()
#     stream_handler.setFormatter(formatter)

#     # Clear any existing handlers to avoid duplicate logging
#     if logger.hasHandlers():
#         logger.handlers.clear()

#     # Add handlers
#     logger.addHandler(file_handler)
#     logger.addHandler(stream_handler)
#     return logger


# def configure_logger():
#     logger = logging.getLogger("wipe")
#     logger.setLevel(logging.INFO)
#     if not logger.hasHandlers():
#         handler = logging.StreamHandler()
#         formatter = logging.Formatter(
#             "%(asctime)s - %(name)s - %(levelname)s: %(message)s"
#         )
#         handler.setFormatter(formatter)
#         logger.addHandler(handler)
#     return logger


# logger = configure_logger()
