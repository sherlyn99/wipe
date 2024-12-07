import os
from wipe.modules.utils import (
    check_outputs,
    check_required_cols,
    create_outdir,
    run_command,
    load_metadata,
    check_required_cols,
)


def gen_command_checkm2(
    inpath,
    threads,
    outdir,
    dbpath,
    genes=None,
):
    """NB: --genes parameter not used currently
    NB: checkm2 takes fna.gz directly
    """
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


def run_checkm2_single(
    infile_path, outdir_path, dbpath, threads, logfile_path
):
    """
    NB
    ---
    Takes about 2-3 minutes using 64 threads.

    Example
    -------
    run_checkm2_single(
        "./wol3_prototype/rawdata/external/SMGC_1.fa.gz",
        "./output2/checkm2_test_out",
        "/home/y1weng/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd",
        64,
        logfile_path="./output2/checkm2_test_out/run_checkm2.log",
    )
    """
    # skip if already completed
    if check_outputs(
        [
            os.path.join(outdir_path, "checkm2.log"),
            os.path.join(outdir_path, "quality_report.tsv"),
        ]
    ):
        return

    create_outdir(outdir_path)
    commands = gen_command_checkm2(
        infile_path,
        threads,
        outdir_path,
        dbpath,
    )
    run_command(commands, logfile=logfile_path)


def run_checkm2_batch(metadata, outdir, dbpath, threads):
    """
    Example
    -------
    run_checkm2_batch(
        "./wol3_prototype/md_prototype_test.txt",
        "./output2/checkm2_test_out",
        "/home/y1weng/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd",
        64,
    )
    """
    df_md = load_metadata(metadata)
    check_required_cols(df_md, ["genome_id", "genome_path"])
    summary_path = os.path.join(outdir, "checkm2_summary.txt")
    for row in df_md.itertuples():
        gid = row.genome_id
        infile_path = row.genome_path
        outdir_path = os.path.join(
            outdir, gid[0], gid[1:4], gid[4:7], gid[7:10]
        )
        logfile_path = os.path.join(outdir_path, f"{gid}_run_checkm2.log")
        try:
            run_checkm2_single(
                infile_path, outdir_path, dbpath, threads, logfile_path
            )
        except:
            with open(summary_path, "a") as log:
                log.write(f"{gid}\n")
