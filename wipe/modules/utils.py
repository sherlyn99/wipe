import bz2
import gzip
import lzma
import logging
from os.path import splitext


def convert_bash_commands_to_subprocess_commands(bash_cmd):
    return bash_cmd.replace(" \\ ", " ").split(" ")


def read(fname):
    """Helper to read compressed files"""
    # compressed filename pattern
    zipdic = {".gz": gzip, ".bz2": bz2, ".xz": lzma, ".lz": lzma}
    ext = splitext(fname)[1]
    zipfunc = getattr(zipdic[ext], "open") if ext in zipdic else open
    return zipfunc(fname, "rt")


def configure_logger():
    logger = logging.getLogger("wipe")
    logger.setLevel(logging.INFO)
    if not logger.hasHandlers():
        handler = logging.StreamHandler()
        formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s: %(message)s"
        )
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    return logger


logger = configure_logger()
