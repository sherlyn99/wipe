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
)


def check_outputs_prodigal(outdir_path, stats_path, gid):
    output_fps = [
        outdir_path,
        stats_path,
        os.path.join(outdir_path, f"{gid}.faa.xz"),
        os.path.join(outdir_path, f"{gid}.ffn.xz"),
        os.path.join(outdir_path, f"{gid}.gff.xz"),
        os.path.join(outdir_path, f"coords_{gid}.txt.xz"),
        os.path.join(outdir_path, f"proteins_{gid}.tsv.xz"),
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
        outdir_path, f"coords_{gid}.txt.xz"
    )
    stats_data["details"]["output_files"]["proteins"] = os.path.join(
        outdir_path, f"proteins_{gid}.tsv.xz"
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


def run_prodigal_single(inpath, outdir, tmpdir):
    """
    NB: tmpdir has to be an existing directory

    Example:
    run_prodigal_single_genome("./linearized_genomes/M/001/002/003/M001002003.fna.gz",
                               "./linearized_genomes/M/001/002/003/",
                               "/tmp/run1")

    Example Outputs:
        "./linearized_genomes/M/001/002/003/prodigal_stats_M001002003.json.gz"
        "./linearized_genomes/M/001/002/003/prodigal_out/M001002003.gff.xz"
        "./linearized_genomes/M/001/002/003/prodigal_out/M001002003.faa.xz"
        "./linearized_genomes/M/001/002/003/prodigal_out/M001002003.ffn.xz"
        "./linearized_genomes/M/001/002/003/prodigal_out/coords_M001002003.txt.xz"
        "./linearized_genomes/M/001/002/003/prodigal_out/proteins_{gid}.tsv.xz"
    """
    # outdir_path = "./linearized_genomes/M/001/002/003/prodigal_out/"
    # stats_path = "./linearized_genomes/M/001/002/003/prodigal_stats_M001002003.json.gz"
    # gid = "M001002003"
    outdir_path, stats_path, gid = gen_output_paths("prodigal", inpath, outdir)

    # skip if already completed
    if check_outputs_prodigal(outdir_path, stats_path, gid):
        return

    # decompress if the inpath is gz compressed
    infna, to_delete = decompress(inpath, tmpdir)

    Path(outdir_path).mkdir(parents=True, exist_ok=True)
    commands = gen_commands_prodigal(infna, gid, outdir_path)
    stats_data = gen_stats_data_prodigal(infna, outdir_path)

    with gzip.open(stats_path, "wt") as file:
        json.dump(stats_data, file, indent=4)

    try:
        p = subprocess.run(
            commands, capture_output=True, text=True, check=True
        )

        stats_data["status"] = "success"

        # Get coords and proteins
        faa_fp = os.path.join(outdir_path, f"{gid}.faa")
        out_coords_fp = stats_data["details"]["output_files"]["coords"]
        out_proteins_fp = stats_data["details"]["output_files"]["proteins"]
        get_coords_and_proteins(faa_fp, out_coords_fp, out_proteins_fp)

        # Compress the Prodigal output files using `xz`
        compress_outputs_prodigal(outdir_path, gid)

    except subprocess.CalledProcessError as e:
        # Capture any errors encountered
        stats_data["status"] = "failure"
        stats_data["error"] = e.stderr  # Error message from Prodigal
    except FileNotFoundError as e:
        # Handle cases where Prodigal or the input file doesn't exist
        stats_data["status"] = "failure"
        stats_data["error"] = str(e)
    except Exception as e:
        stats_data["status"] = "failure"
        stats_data["error"] = str(e)
    finally:
        # Update end time and write final log
        stats_data["end_time"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        with gzip.open(stats_path, "wt") as file:
            json.dump(stats_data, file, indent=4, ensure_ascii=False)

        if to_delete:
            os.remove(infna)


def run_prodigal_batch(indir, logdir, tmpdir):
    """
    Example
    -------
    """
    filelist = get_files_all_fa(indir)
    summary_path, summary_data = gen_summary_prodigal(logdir)

    with gzip.open(summary_path, "wt") as file:
        json.dump(summary_data, file, indent=4)

    for file in filelist:
        outdir_path = os.path.dirname(file)
        try:
            run_prodigal_single(file, outdir_path, tmpdir)
        except Exception as e:
            if isinstance(summary_data["error"], list):
                summary_data["error"].append(f"{file}: {e}")
            else:
                summary_data["error"] = [f"{file}: {e}"]

    summary_data["status"] = "completed"
    summary_data["end_time"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with gzip.open(summary_path, "wt") as file:
        json.dump(summary_data, file, indent=4, ensure_ascii=False)
