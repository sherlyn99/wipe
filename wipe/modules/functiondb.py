import os
import lzma
import pandas as pd
from pathlib import Path
from wipe.modules.utils import check_required_cols


def get_coords(
    metadata_path,
    outdir,
):
    """
    NB: metadata_path has to have the two columns 'genome_id' and 'lgenome_dir'.
    """
    with open(metadata_path, "r") as f:
        lines = f.read().splitlines()

    md_df = pd.read_csv(metadata_path, sep="\t", low_memory=True)
    check_required_cols(md_df, ["genome_id", "lgenome_dir"])

    with lzma.open(f"{outdir}/coords.txt.xz", "wt") as fo:
        for row in md_df.itertuples():
            g = row.genome_id
            lgdir = row.lgenome_dir
            prodigal_outdir = os.path.join(lgdir, "prodigal_out")
            with lzma.open(f"{prodigal_outdir}/{g}.faa.xz", "rb") as fi:
                lines = fi.read().decode("utf-8").splitlines()
            cnuc = None
            for line in lines:
                if line.startswith(">"):
                    name, pos5, pos3, strand, _ = line[1:].split(" # ", 4)
                    nucl, idx = name.rsplit("_", 1)
                    if nucl != cnuc:
                        cnuc = nucl
                        fo.write(f">{nucl}\n")
                    beg, end = (pos5, pos3) if strand == "1" else (pos3, pos5)
                    fo.write(f"{idx}\t{beg}\t{end}\n")
            fi.close()
            print(g)
    fo.close()


def get_uniref_map():
    from pathlib import Path


def merge_uniref(uniref90, uniref50, outfile, simplify=False):
    """
    Merge UniRef90 (preferred) and UniRef50 maps.

    Example
    -------
    merge_uniref('./output/proteins/uniref90.m8',
                 './output/proteins/uniref50.m8',
                 './output/functions/uniref/uniref_map.txt',
                 True)
    """
    res = {}
    res_add = res.setdefault

    # Ensure the parent directory of outfile exists
    Path(outfile).parent.mkdir(parents=True, exist_ok=True)

    # Define how to process each entry based on the simplify flag
    if simplify:
        process_entry = lambda entry: entry.split("_")[
            1
        ]  # Simplified: Take only the part after '_'
    else:
        process_entry = lambda entry: entry  # No simplification

    # Read UniRef90 file
    with open(uniref90, "r") as f:
        for line in f:
            orf, entry = line.rstrip().split("\t")[:2]
            g, i = orf.split("_")
            res_add(g, {})[int(i)] = process_entry(entry)

    # Read UniRef50 file
    with open(uniref50, "r") as f:
        for line in f:
            orf, entry = line.rstrip().split("\t")[:2]
            g, i = orf.split("_")
            if g not in res or int(i) not in res[g]:
                res_add(g, {})[int(i)] = process_entry(entry)

    # Write output to the specified outfile
    with open(outfile, "w") as out_f:
        for g, entries in sorted(res.items()):
            for i, entry in sorted(entries.items()):
                out_f.write(f"{g}_{i}\t{entry}\n")


def get_uniref_names():
    # cd /panfs/y1weng/uniprot_db/uniref50
    # zcat uniref50.xml.gz | python ~/31_hcom2/hcom2_subset/extract_uniref_xml.py > uniref50.tsv
    # cd /panfs/y1weng/uniprot_db/uniref90
    # zcat uniref90.xml.gz | python ~/31_hcom2/hcom2_subset/extract_uniref_xml.py > uniref90.tsv
    # do this for a new dataset every time
    # get_uniref_name.py in sketch.ipynb

    umap = "/home/y1weng/31_hcom2/hcom2_subset/output/functions/uniref/uniref_map.txt"
    uniref50_names = "/panfs/y1weng/uniprot_db/uniref50/uniref50.tsv"
    uniref90_names = "/panfs/y1weng/uniprot_db/uniref90/uniref90.tsv"
    output_uniref_names = "/home/y1weng/31_hcom2/hcom2_subset/output/functions/uniref/uniref_names.txt"

    urs = set()
    with open(umap, "r") as f:
        lines = f.readlines()
        for line in lines:
            uid = line.strip().split("\t")[1]
            urs.add(uid)

    # with open('uniref.lst', 'r') as f:
    #     urs = set(f.read().splitlines())

    with open(output_uniref_names, "w") as fo:

        seen = []
        with open(uniref90_names, "r") as fi:
            for line in fi:
                if not (line.startswith("#") or line.startswith("Error: ")):
                    ur, name, _ = line.split("\t", 2)
                    if ur in urs:
                        print(ur, name, sep="\t", file=fo)
                        seen.append(ur)

        seen = set(seen)
        with open(uniref50_names, "r") as fi:
            for line in fi:
                if not (line.startswith("#") or line.startswith("Error: ")):
                    ur, name, _ = line.split("\t", 2)
                    if ur in urs and ur not in seen:
                        print(ur, name, sep="\t", file=fo)

        fo.close()
