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


def filter_contig_name(filt, contig_name):
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


def linearize_single_genome(
    inpath, ext, outdir, gid=None, gap=None, filt=None
):
    if not gid:
        gid = infer_gid_ncbi(inpath)

    ext = ext.lstrip(".")
    inpath = inpath + "." + ext

    outdir = join(outdir, gid)
    if not exists(inpath):
        with open(f"{outdir}/linearization.log", "w") as fo_log:
            fo_log.write("Input file does not exist")
        return

    outpath = join(outdir, f"{gid}.{ext}")
    os.makedirs(outdir, exist_ok=True)

    with read(inpath) as fi, open(outpath, "w") as fo:
        # Write gid as a header line when in concat mode
        if gap:
            fo.write(f">{gid}\n")
        contigs, contig = [], ""
        ncontig, nchar = 0, 0
        contigs_filtered, contigs_empty = [], []
        curr_header = fi.readline().rstrip("\r\n")

        def linearize_helper():
            nonlocal ncontig, nchar
            todel = filter_contig_name(filt, curr_header)
            if not todel:
                if contig:
                    if gap:
                        contigs.append(contig)
                    else:
                        fo.write(f">{gid}_{curr_header[1:]}\n")
                        fo.write(f"{contig}\n")
                    ncontig += 1
                    nchar += len(contig)

                else:
                    contigs_empty.append(curr_header)
            else:
                contigs_filtered.append(curr_header)

        for line in fi:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                linearize_helper()
                contig = ""
                curr_header = line
            else:
                contig += line

        linearize_helper()

        if gap:
            gap_string = generate_gap_string(gap)
            fo.write(f"{gap_string.join(contigs)}\n")

    with open(f"{outdir}/linearization.log", "w") as fo_log:
        fo_log.write(
            f"Wrote {ncontig} contigs ({nchar} characters in total) into {outpath}\n"
        )
        if contigs_empty:
            fo_log.write(
                f"\n{len(contigs_empty)} contigs were empty:\n{contigs_empty}\n"
            )
        if contigs_filtered:
            fo_log.write(
                f"\n{len(contigs_filtered)} contigs were filtered.\n{contigs_filtered}\n"
            )
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
