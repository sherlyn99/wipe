import os
import shutil
from pathlib import Path
from wipe.modules.utils import run_command


def create_db(indir, outdir, nthreads=1):
    """
    Output
    ------
    {outdir}/hnswdump.hnsw.data
    {outdir}/hnswdump.hnsw.graph
    {outdir}/parameters.json
    {outdir}/processing_state.json
    {outdir}/seqdict.json
    """
    Path(outdir).mkdir(parents=True, exist_ok=True)
    commands = [
        "gsearch",
        "--pio",
        "2000",
        "--nbthreads",
        str(nthreads),
        "tohnsw",
        "-d",
        indir,
        "-k",
        "16",
        "-s",
        "18000",
        "-n",
        "128",
        "--ef",
        "1600",
        "--algo",
        "optdens",
    ]
    run_command(commands, outdir)


def update_db(existing_db, indir, nthreads=1, backup_dir=None):
    if backup_dir:
        Path(backup_dir).mkdir(parents=True, exist_ok=True)
        shutil.copytree(existing_db, backup_dir, dirs_exist_ok=True)

    commands = [
        "gsearch",
        "--pio",
        "2000",
        "--nbthreads",
        str(nthreads),
        "add",
        "--hnsw",
        existing_db,
        "--new",
        indir,
    ]
    run_command(commands)


def search_genomes(db, n_neighbors, indir, outdir, nthreads=1):
    Path(outdir).mkdir(parents=True, exist_ok=True)
    commands = [
        "gsearch",
        "--pio",
        "2000",
        "--nbthreads",
        str(nthreads),
        "request",
        "--hnsw",
        db,
        "--nbanswers",
        str(n_neighbors),
        "-r",
        indir,
    ]
    run_command(commands, outdir)
    return f"{outdir}/gsearch.neighbors.txt"


def reformat_search_results(inpath, outpath, model=1):
    pdir = os.path.dirname(inpath)
    outpath = os.path.join(pdir, "gsearch_results.tsv")
    commands = ["reformat", "16", str(model), inpath, outpath]
    run_command(commands)
    return outpath
