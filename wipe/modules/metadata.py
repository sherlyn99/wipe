import os
import pandas as pd
from pathlib import Path
from wipe.modules.utils import get_files


def generate_gids(start_genome_id, num_ids):
    prefix = start_genome_id[0]
    number_part = int(start_genome_id[1:])

    genome_ids = [
        f"{prefix}{str(number_part + i).zfill(len(start_genome_id) - 1)}"
        for i in range(num_ids)
    ]

    return genome_ids


def generate_metadata(indir, ext, gid_start=None):
    files = get_files(indir, ext)
    ext = "." + ext.lstrip(".")
    files_without_ext = [f.replace(ext, "") for f in files]

    if not gid_start:
        return pd.DataFrame(files_without_ext, columns=["filepath"])

    n = len(files)
    gids = generate_gids(gid_start, n)
    return pd.DataFrame({"filepath": files_without_ext, "genome_id": gids})
