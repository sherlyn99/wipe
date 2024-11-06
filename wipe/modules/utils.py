import os
import bz2
import json
import gzip
import lzma
import time
import logging
import pandas as pd
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
