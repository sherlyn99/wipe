import os
import gzip
import lzma
import json
import skbio
import shutil
import subprocess
from pathlib import Path
from datetime import datetime
from wipe.modules.utils import (
    gen_output_paths,
    check_outputs,
    check_log_and_retrieve_gid,
    gen_stats_data,
    xz_compress_files,
    get_files_all_fa,
    gen_path,
    gen_summary_data,
    run_command,
)


def check_outputs_prodigal(outdir_path, stats_path, gid):
    output_fps = [
        outdir_path,
        stats_path,
        os.path.join(outdir_path, f"{gid}.faa.xz"),
        os.path.join(outdir_path, f"{gid}.ffn.xz"),
        os.path.join(outdir_path, f"{gid}.gff.xz"),
        os.path.join(outdir_path, f"{gid}_coords.txt.xz"),
        os.path.join(outdir_path, f"{gid}_proteins_ct.tsv.xz"),
    ]
    file_existence = check_outputs(output_fps)
    status, gid = check_log_and_retrieve_gid(stats_path)
    return file_existence and status


def decompress(inpath, tmpdir):
    """Decompress the genome file into tmpdir/{gid}.fna if it
    is gz compressed. If not, return inpath.

    Example
    -------
    decompress("/G/999/G999.fna.gz", "/tmp")
    >> "/tmp/G999.fna", TRUE

    decompress("/G/999/G999.fna", "/tmp")
    >> "/G/999/G999.fna", FALSE
    """
    gid = os.path.basename(inpath).split(".")[0]
    if inpath.endswith(".gz"):
        outpath = os.path.join(tmpdir, f"{gid}.fna")
        with gzip.open(inpath, "rb") as fi:
            with open(outpath, "wb") as fo:
                shutil.copyfileobj(fi, fo)
        return outpath, True
    else:
        return inpath, False


def gen_commands_prodigal(infna, gid, outdir):
    # NB: Prodiigal does not take .fna.gz and one thread only as of 2024-11-11
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


def gen_stats_data_prodigal(inpath, outdir_path):
    stats_data = gen_stats_data("prodigal_run", inpath)

    gid = os.path.basename(inpath).split(".")[0]
    stats_data["details"]["output_files"]["faa"] = os.path.join(
        outdir_path, f"{gid}.faa.xz"
    )
    stats_data["details"]["output_files"]["ffn"] = os.path.join(
        outdir_path, f"{gid}.ffn.xz"
    )
    stats_data["details"]["output_files"]["gff"] = os.path.join(
        outdir_path, f"{gid}.gff.xz"
    )
    stats_data["details"]["output_files"]["coords"] = os.path.join(
        outdir_path, f"{gid}_coords.txt.xz"
    )
    stats_data["details"]["output_files"]["proteins"] = os.path.join(
        outdir_path, f"{gid}_proteins_ct.tsv.xz"
    )

    return stats_data


def compress_outputs_prodigal(outdir_prodigal, gid):
    filelist = [
        os.path.join(outdir_prodigal, f"{gid}.gff"),
        os.path.join(outdir_prodigal, f"{gid}.faa"),
        os.path.join(outdir_prodigal, f"{gid}.ffn"),
    ]
    xz_compress_files(filelist)


def get_coords_and_proteins(faa_fp, out_coords_fp, out_proteins_fp):
    proteins_ct = 0
    with lzma.open(out_coords_fp, "wt") as coords_fo:
        current_nucl = None
        for seq in skbio.read(faa_fp, format="fasta"):
            header = seq.metadata["id"] + " " + seq.metadata["description"]
            name, pos5, pos3, strand, _ = header.split(" # ", 4)
            nucl, idx = name.rsplit("_", 1)
            beg, end = (pos5, pos3) if strand == "1" else (pos3, pos5)
            if nucl != current_nucl:
                current_nucl = nucl
                coords_fo.write(f">{nucl}\n")

            coords_fo.write(f"{idx}\t{beg}\t{end}\n")
            proteins_ct += 1

    with lzma.open(out_proteins_fp, "wt") as proteins_fo:
        proteins_fo.write("genome_id\tproteins_ct\n")
        proteins_fo.write(f"{current_nucl}\t{str(proteins_ct)}\n")


def gen_summary_path_prodigal(logdir):
    return gen_path(logdir, "prodigal_summary.json.gz")


def gen_summary_data_prodigal():
    return gen_summary_data("prodigal_run")


def gen_summary_prodigal(logdir):
    summary_path = gen_summary_path_prodigal(logdir)
    summary_data = gen_summary_data_prodigal()
    return summary_path, summary_data


def run_prodigal_single(gid, infile_path, outdir_path, tmpdir):
    """
    NB: tmpdir has to be an existing directory
    infile can be fna or fna.gz

    Example:
    run_prodigal_single_genome("M001002003",
                               "./linearized_genomes/M/001/002/003/M001002003.fna.gz",
                               "./linearized_genomes/M/001/002/003/prodigal_out",
                               "/tmp/run1")

    Example Outputs:
        "./linearized_genomes/M/001/002/003/prodigal_out/M001002003_prodigal_log.json.gz"
        "./linearized_genomes/M/001/002/003/prodigal_out/M001002003.gff.xz"
        "./linearized_genomes/M/001/002/003/prodigal_out/M001002003.faa.xz"
        "./linearized_genomes/M/001/002/003/prodigal_out/M001002003.ffn.xz"
        "./linearized_genomes/M/001/002/003/prodigal_out/M001002003_coords.txt.xz"
        "./linearized_genomes/M/001/002/003/prodigal_out/M001002003_proteins_ct.tsv.xz"
    """
    gff = os.path.join(outdir_path, f"{gid}.gff")
    gff_xz = gff + ".xz"
    faa = os.path.join(outdir_path, f"{gid}.faa")
    faa_xz = faa + ".xz"
    ffn = os.path.join(outdir_path, f"{gid}.ffn")
    ffn_xz = ffn + ".xz"
    coords_xz = os.path.join(outdir_path, f"{gid}_coords.txt.xz")
    proteins_ct_xz = os.path.join(outdir_path, f"{gid}_proteins_ct.tsv.xz")
    log_file = os.path.join(outdir_path, f"{gid}_run_prodigal.log")
    log_file_xz = log_file + ".xz"
    outputs = [
        gff_xz,
        faa_xz,
        ffn_xz,
        coords_xz,
        proteins_ct_xz,
        log_file_xz,
    ]

    if check_outputs(outputs):
        return faa_xz
    else:
        for o in outputs:
            if os.path.exists(o):
                os.remove(o)

    # decompress if the inpath is gz compressed
    infna, to_delete = decompress(infile_path, tmpdir)

    Path(outdir_path).mkdir(parents=True, exist_ok=True)
    commands = gen_commands_prodigal(infna, gid, outdir_path)
    try:
        run_command(commands, logfile=log_file)
        get_coords_and_proteins(faa, coords_xz, proteins_ct_xz)
        xz_compress_files([gff, faa, ffn, log_file])
        if to_delete:
            os.remove(infna)
    finally:
        return os.path.join(outdir_path, f"{gid}.faa.xz")


# def run_prodigal_batch(indir, logdir, tmpdir):
#     """
#     Example
#     -------
#     """
#     filelist = get_files_all_fa(indir)
#     summary_path, summary_data = gen_summary_prodigal(logdir)

#     with gzip.open(summary_path, "wt") as file:
#         json.dump(summary_data, file, indent=4)

#     for file in filelist:
#         outdir_path = os.path.dirname(file)
#         try:
#             run_prodigal_single(file, outdir_path, tmpdir)
#         except Exception as e:
#             if isinstance(summary_data["error"], list):
#                 summary_data["error"].append(f"{file}: {e}")
#             else:
#                 summary_data["error"] = [f"{file}: {e}"]

#     summary_data["status"] = "completed"
#     summary_data["end_time"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
#     with gzip.open(summary_path, "wt") as file:
#         json.dump(summary_data, file, indent=4, ensure_ascii=False)
