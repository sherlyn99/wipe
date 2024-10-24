import os
import re
import sys
import skbio
from os.path import join, exists
from wipe.modules.utils import logger, read


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


def infer_gid_ncbi(inpath):
    """Infer genome id (GXXXXX) from ncbi-stype genome file name (GCF_XXXXXX_XX_genomic.fna)"""
    fname = inpath.split("/")[-1]
    ptn = re.compile(r"GC[FA]_([0-9]{9})\.[0-9]+_.*")
    m = ptn.match(fname)
    if m:
        gid = "G" + m.group(1)
        return gid
    else:
        raise ValueError("No valid genome ID provided or extracted.")


def generate_inpath_outpath(inpath, ext, outdir, gid):
    gid = gid or infer_gid_ncbi(inpath)
    ext = ext.lstrip(".")
    inpath = inpath + "." + ext
    outdir = join(outdir, gid)
    out_ext = ext.rstrip(".gz").rstrip(".bz2").rstrip(".xz").rstrip(".lz")
    outpath = join(outdir, f"{gid}.{out_ext}")
    return gid, inpath, outdir, outpath


def generate_log_msg(n_written, n_char, n_filtered, outpath):
    log_message = (
        f"Wrote {n_written} contigs ({n_char} characters in total) into {outpath}.\n"
        + (f"{n_filtered} contigs were filtered.\n")
    )
    return log_message


def read_fasta(inpath):
    try:
        records = {}
        for seq in skbio.read(inpath, format="fasta"):
            header = (
                f">{seq.metadata['id']} {seq.metadata.get('description', '')}"
            ).strip()
            records[header] = seq
    except Exception as e:
        logger.error(f"Reading {inpath} failed due to the error {str(e)}.")
        sys.exit(1)
    return records


def write_log(outdir, filename, message):
    """Writes a log message to the linearization log."""
    with open(f"{outdir}/{filename}", "w") as fo_log:
        fo_log.write(message)


def linearize_single_genome(
    inpath, ext, outdir, gid=None, gap=None, filt=None
):
    gid, inpath, outdir, outpath = generate_inpath_outpath(
        inpath, ext, outdir, gid
    )
    os.makedirs(outdir, exist_ok=True)

    if not exists(inpath):
        logger.error(f"Input file {inpath} does not exist.")
        sys.exit(1)

    linearized = read_fasta(inpath)
    n_written, n_char, n_filtered = 0, 0, 0

    with open(outpath, "w") as output_file:
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

        log_msg = generate_log_msg(n_written, n_char, n_filtered, outpath)
        write_log(outdir, "linearization.log", log_msg)

    return


def linearize_genomes(metadata, ext, outdir, gap, filt):
    gid_set = set()
    with read(metadata) as fi:
        for line in fi:
            line = line.strip()
            if not line:
                continue

            fields = list(map(str.strip, line.split("\t")))
            inpath = fields[0]
            gid = fields[1] if len(fields) == 2 else infer_gid_ncbi(inpath)

            if gid in gid_set:
                write_log(
                    outdir,
                    "linearization_all.log",
                    f"Dupliacted genome ID: {gid} with input path {inpath}\n",
                )
            gid_set.add(gid)

            linearize_single_genome(inpath, ext, outdir, gid, gap, filt)
