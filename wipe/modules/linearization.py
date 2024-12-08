import os
import re
import gzip
import skbio
import pandas as pd
from wipe.modules.utils import (
    check_outputs,
    create_outdir,
    load_metadata,
    check_required_cols,
    check_duplicated_genome_ids,
)


def read_fasta(infile_path):
    """This works for .fna as well .fna.gz but not .fna.xz"""
    records = {}
    for seq in skbio.read(infile_path, format="fasta"):
        header = (
            f">{seq.metadata['id']} {seq.metadata.get('description', '')}"
        ).strip()
        records[header] = seq
    return records


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


def generate_log_entries(gid, n_written, n_char, n_filtered, outfile_path):
    return {
        "genome_id": gid,
        "lgenome_path": outfile_path,
        "contigs_written": n_written,
        "chars_written": n_char,
        "contigs_filtered": n_filtered,
    }


def linearization_single(infile_path, gid, outdir_path, gap, filt=None):
    """
    Example
    -------
    linearization_single(
    "wol3_prototype/rawdata/GCA/000/965/005/GCA_000965005.1_ASM96500v1/GCA_000965005.1_ASM96500v1_genomic.fna.gz",
    "G000965005",
    "./output2",
    "N*20")
    """
    outfile_path = os.path.join(outdir_path, f"{gid}.fna.gz")
    linearization_result_path = os.path.join(outdir_path, "linearization.tsv")

    if not os.path.exists(infile_path):
        raise FileNotFoundError(f"Input file {infile_path} does not exist.")

    # Skip if expected outputs already exist
    if check_outputs([outfile_path, linearization_result_path]):
        return

    create_outdir(outdir_path)
    linearized = read_fasta(infile_path)
    n_written, n_char, n_filtered = 0, 0, 0

    with gzip.open(outfile_path, "wt") as output_file:
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
            gid, n_written, n_char, n_filtered, outfile_path
        )
        pd.DataFrame([log_data]).to_csv(
            linearization_result_path,
            sep="\t",
            index=False,
            header=True,
        )


def linearization_batch(metadata, outdir, gap, filt):
    df_md = load_metadata(metadata)
    check_required_cols(df_md, ["genome_id", "genome_path"])
    check_duplicated_genome_ids(df_md)
    summary_path = os.path.join(outdir, "linearization_summary.txt")
    for row in df_md.itertuples():
        gid = row.genome_id
        infile_path = row.genome_path
        outdir_path = os.path.join(
            outdir, gid[0], gid[1:4], gid[4:7], gid[7:10]
        )
        try:
            linearization_single(infile_path, gid, outdir_path, gap, filt)
        except:
            with open(summary_path, "a") as log:
                log.write(f"{gid}\n")
