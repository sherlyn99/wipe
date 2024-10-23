import os
import re
import gzip
import bz2
import lzma
from os.path import join, exists, splitext


def generate_gap_string(gap):
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


def read(fname):
    """Helper to read compressed files"""
    # compressed filename pattern
    zipdic = {".gz": gzip, ".bz2": bz2, ".xz": lzma, ".lz": lzma}
    ext = splitext(fname)[1]
    zipfunc = getattr(zipdic[ext], "open") if ext in zipdic else open
    return zipfunc(fname, "rt")


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


def generate_log_msg(ncontig, nchar, contigs_empty, contigs_filtered, outpath):
    log_message = (
        f"Wrote {ncontig} contigs ({nchar} characters in total) into {outpath}\n"
        + (
            f"\n{len(contigs_empty)} contigs were empty:\n{contigs_empty}\n"
            if contigs_empty
            else ""
        )
        + (
            f"\n{len(contigs_filtered)} contigs were filtered.\n{contigs_filtered}\n"
            if contigs_filtered
            else ""
        )
    )
    return log_message


def update_ncontig_stats(contig, ncontig, nchar):
    ncontig += 1
    nchar += len(contig)
    return ncontig, nchar


def write_contig(gid, curr_header, contig, contigs, gap, fo):
    if gap:
        contigs.append(contig)
    else:
        fo.write(f">{gid}_{curr_header[1:]}\n{contig}\n")


def write_log(outdir, message):
    """Writes a log message to the linearization log."""
    with open(f"{outdir}/linearization.log", "w") as fo_log:
        fo_log.write(message)


def check_file_exists(inpath, outdir):
    """Check if the input file exists, if not, log the error."""
    if not exists(inpath):
        write_log(outdir, "Input file does not exist")
        return False
    return True


def handle_contig(
    filt,
    gap,
    gid,
    curr_header,
    contig,
    contigs,
    contigs_filtered,
    contigs_empty,
    fo,
    ncontig,
    nchar,
):
    if should_filter_contig_name(filt, curr_header):
        contigs_filtered.append(curr_header)
    else:
        if contig:
            write_contig(gid, curr_header, contig, contigs, gap, fo)
            ncontig, nchar = update_ncontig_stats(contig, ncontig, nchar)
        else:
            contigs_empty.append(curr_header)
    return ncontig, nchar


def generate_inpath_outpath(inpath, ext, outdir, gid):
    ext = ext.lstrip(".")
    inpath = inpath + "." + ext
    outdir = join(outdir, gid)
    out_ext = ext.rstrip(".gz").rstrip(".bz2").rstrip(".xz").rstrip(".lz")
    outpath = join(outdir, f"{gid}.{out_ext}")
    return inpath, outdir, outpath


def linearize_single_genome(
    inpath, ext, outdir, gid=None, gap=None, filt=None
):
    gid = gid or infer_gid_ncbi(inpath)
    inpath, outdir, outpath = generate_inpath_outpath(inpath, ext, outdir, gid)

    if not check_file_exists(inpath, outdir):
        return

    os.makedirs(outdir, exist_ok=True)

    with read(inpath) as fi, open(outpath, "w") as fo:
        if gap:
            fo.write(
                f">{gid}\n"
            )  # Write gid as a header line when in concat mode
        contigs, contigs_filtered, contigs_empty, contig = [], [], [], ""
        ncontig, nchar = 0, 0

        curr_header = fi.readline().rstrip("\r\n")

        for line in fi:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                # fmt: off
                ncontig, nchar = handle_contig(filt, gap, gid, curr_header, contig, contigs, contigs_filtered, contigs_empty, fo, ncontig, nchar)
                # fmt: on
                contig = ""
                curr_header = line
            else:
                contig += line

        # fmt: off
        ncontig, nchar = handle_contig(filt, gap, gid, curr_header, contig, contigs, contigs_filtered, contigs_empty, fo, ncontig, nchar)
        # fmt: on

        if gap:
            gap_string = generate_gap_string(gap)
            fo.write(f"{gap_string.join(contigs)}\n")

    log_msg = generate_log_msg(
        ncontig, nchar, contigs_empty, contigs_filtered, outpath
    )
    write_log(outdir, log_msg)

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
                with open(join(outdir, "linearization_all.log"), "w") as fo:
                    fo.write(
                        f"Dupliacted genome ID: {gid} with input path {inpath}\n"
                    )
            gid_set.add(gid)

            linearize_single_genome(inpath, ext, outdir, gid, gap, filt)
