import os
import re
import sys
import gzip
import skbio
from os.path import join, exists
from wipe.modules.utils import (
    write_json_log,
    load_md,
    check_duplicated_genome_ids,
    infer_gid_ncbi,
)


def generate_gap_string(gap):
    """If the input is erroneous. gap_string will default to ''."""
    gap_string = ""
    if "*" in gap:
        str_, n_ = gap.rsplit("*", 1)
        if n_.isdigit():
            gap_string = str_ * int(n_)
    return gap_string


def should_filter_contig_name(filt, contig_name):
    """Returns True if _filt_ exists in the contig name"""
    if filt:
        ptn = re.compile(
            r"[^a-zA-Z0-9](%s)[^a-zA-Z0-9]" % filt.replace(",", "|"),
            re.IGNORECASE,
        )
        return ptn.search(contig_name) is not None
    else:
        return False


def generate_inpath_outpath(inpath, ext, outdir, gid):
    gid = gid or infer_gid_ncbi(inpath)
    ext = ext.lstrip(".")
    inpath = inpath + "." + ext
    outdir = join(outdir, gid[0], gid[1:4], gid[4:7], gid[7:])
    out_ext = ext.rstrip(".gz").rstrip(".bz2").rstrip(".xz").rstrip(".lz")
    outpath = join(outdir, f"{gid}.{out_ext}.gz")
    return gid, inpath, outdir, outpath


def generate_log_entries(gid, n_written, n_char, n_filtered, outpath):
    return {
        "genome_id": gid,
        "genome_path": os.path.abspath(outpath),
        "contigs_written": n_written,
        "chars_written": n_char,
        "contigs_filtered": n_filtered,
    }


def read_fasta(inpath):
    """This works for .fna as well .fna.gz but not .fna.xz"""
    records = {}
    for seq in skbio.read(inpath, format="fasta"):
        header = (
            f">{seq.metadata['id']} {seq.metadata.get('description', '')}"
        ).strip()
        records[header] = seq
    return records


def check_outputs(outdir, outpath, gid):
    """Check if expected outputs already exist."""
    genome_file_exists = False
    log_file_exists = False
    if exists(outdir):
        if exists(outpath):
            genome_file_exists = True
        if exists(join(outdir, f"linearization_stats_{gid}.json.gz")):
            log_file_exists = True
    return genome_file_exists and log_file_exists


def generate_outdir(outdir, gid):
    return join(outdir, gid[0], gid[1:4], gid[4:7], gid[7:])


def linearize_single_genome(
    inpath, ext, outdir, gid=None, gap=None, filt=None
):
    gid, inpath, outdir, outpath = generate_inpath_outpath(
        inpath, ext, outdir, gid
    )

    # Skip if expected outputs already exist
    if check_outputs(outdir, outpath, gid):
        return

    # Error out if the input path does not exist
    if not exists(inpath):
        raise FileNotFoundError(f"Input file {inpath} does not exist.")

    os.makedirs(outdir, exist_ok=True)

    linearized = read_fasta(inpath)
    n_written, n_char, n_filtered = 0, 0, 0

    with gzip.open(outpath, "wt") as output_file:
        if gap:
            output_file.write(f">{gid}\n")
            gap_string = generate_gap_string(gap)
            seqs = []
            for key, val in linearized.items():
                if not should_filter_contig_name(filt, key):
                    seqs.append(str(val))
                    n_written += 1
                    n_char += len(val)
                else:
                    n_filtered += 1
            full_seq = gap_string.join(seqs)
            n_char += (len(seqs) - 1) * len(gap_string)
            output_file.write(f"{full_seq}\n")
        else:
            for key, val in linearized.items():
                if not should_filter_contig_name(filt, key):
                    output_file.write(f">{gid}_{key[1:]}\n")
                    output_file.write(f"{val}\n")
                    n_char += len(val)
                    n_written += 1
                else:
                    n_filtered += 1

        log_data = generate_log_entries(
            gid, n_written, n_char, n_filtered, outpath
        )
        write_json_log(
            log_data, outdir, f"linearization_stats_{gid}.json.gz", append=True
        )


def linearize_genomes(metadata, ext, outdir, gap, filt):
    md_df = load_md(metadata)

    if "genome_id" not in md_df.columns:
        md_df["genome_id"] = md_df["filepath"].apply(infer_gid_ncbi)

    try:
        check_duplicated_genome_ids(md_df)
    except Exception as e:
        log_data = {"input_path": "NA", "error": f"{e}"}
        write_json_log(
            log_data, outdir, "linearization_summary.gz", append=True
        )
        sys.exit(1)

    for row in md_df.itertuples(index=False):
        inpath = row.filepath.strip()
        gid = row.genome_id.strip()

        try:
            linearize_single_genome(inpath, ext, outdir, gid, gap, filt)
        except Exception as e:
            log_data = {"input_path": inpath, "error": f"{e}"}
            write_json_log(
                log_data, outdir, "linearization_summary.json.gz", append=True
            )

    log_data = {"status": "complete"}
    write_json_log(
        log_data, outdir, "linearization_summary.json.gz", append=True
    )
