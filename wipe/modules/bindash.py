import os
import numpy as np
import pandas as pd
from pathlib import Path
from wipe.modules.utils import get_files_all_fa, run_command


def bindash_sketch(
    indir,
    outdir,
    nthreads=4,
    kmerlen=16,
    sketchsize64=200,
    dens=3,
    bbits=64,
    minhashtype=2,
):
    filelist = get_files_all_fa(indir)
    filelist_fp = os.path.join(outdir, "filelist_for_bindash.txt")
    Path(outdir).mkdir(parents=True, exist_ok=True)

    with open(filelist_fp, "w") as f:
        for fp in filelist:
            f.write(fp + "\n")

    outfname = os.path.join(outdir, "sketch.bindash")
    commands_sketch = [
        "bindash",
        "sketch",
        f"--kmerlen={str(kmerlen)}",
        f"--sketchsize64={str(sketchsize64)}",
        f"--dens={str(dens)}",
        f"--bbits={str(bbits)}",
        f"--minhashtype={str(minhashtype)}",
        f"--listfname={filelist_fp}",
        f"--outfname={outfname}",
        f"--nthreads={str(nthreads)}",
    ]
    run_command(commands_sketch)
    return outfname


def bindash_dist(outdir, dist_fp, nthreads=4):
    commands_dist = [
        "bindash",
        "dist",
        f"--nthreads={str(nthreads)}",
        f"--outfname={dist_fp}",
        os.path.join(outdir, "sketch.bindash"),
    ]
    run_command(commands_dist)
    return dist_fp


def gen_dm(dist_fp):
    pdir = os.path.dirname(dist_fp)
    dm_fp = os.path.join(pdir, "bindash.dm")

    df = pd.read_csv(
        dist_fp,
        sep="\t",
        names=[
            "genome_1",
            "genome_2",
            "distance",
            "significance",
            "shared_minhash",
        ],
        low_memory=False,
    )
    df["genome_1"] = (
        df["genome_1"].str.split("/").str[-1].str.split(".").str[0]
    )
    df["genome_2"] = (
        df["genome_2"].str.split("/").str[-1].str.split(".").str[0]
    )

    dist_matrix = df.pivot(
        index="genome_1", columns="genome_2", values="distance"
    )
    dist_matrix = dist_matrix.fillna(0)
    print("****")
    print(dist_matrix.head(5))
    # dist_matrix = (
    #     dist_matrix + dist_matrix.T - np.diag(dist_matrix.values.diagonal())
    # )
    return dist_matrix


def bindash2_gen_dm(indir, outdir, nthreads):
    pass


if __name__ == "__main__":
    bindash_sketch(
        "/home/y1weng/47_wipe/wipe/wol3_prototype/output2",
        "/home/y1weng/47_wipe/wipe/wol3_prototype/output2/bindash_out/",
    )
    bindash_dist(
        "/home/y1weng/47_wipe/wipe/wol3_prototype/output2/bindash_out/",
        "/home/y1weng/47_wipe/wipe/wol3_prototype/output2/bindash_out/dist.txt",
    )
