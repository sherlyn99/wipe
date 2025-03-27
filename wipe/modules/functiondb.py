import os
import lzma
import pandas as pd
from pathlib import Path
from wipe.modules.utils import check_required_cols, write_json_log, run_command

def run_functional_annotation(indir, diamonddb, outdir, threads):
    """
    Run functional annotation using DIAMOND against UniRef database.
    
    Args:
        indir (str): Input directory containing protein files (all.faa)
        diamonddb (str): Path to DIAMOND database
        outdir (str): Output directory for results
        threads (int): Number of threads to use
    """
    # Create output directory
    os.makedirs(outdir, exist_ok=True)
    
    # Input validation
    protein_file = os.path.join(indir, "all.faa")
    if not os.path.exists(protein_file):
        raise FileNotFoundError(f"Protein file not found: {protein_file}")
    
    if not os.path.exists(diamonddb):
        raise FileNotFoundError(f"DIAMOND database not found: {diamonddb}")

    # Output files
    diamond_out = os.path.join(outdir, "diamond_results.m8")
    uniref_map = os.path.join(outdir, "uniref.map")
    
    try:
        # Run DIAMOND blastp
        diamond_cmd = [
            "diamond", "blastp",
            "--threads", str(threads),
            "--db", diamonddb,
            "--query", protein_file,
            "--out", diamond_out,
            "--id", "90",
            "--subject-cover", "80",
            "--query-cover", "80",
            "--index-chunks", "1",
            "--max-target-seqs", "1"
        ]
        
        run_command(diamond_cmd)
        
        # Process DIAMOND results to create uniref.map
        with open(diamond_out) as f_in, open(uniref_map, 'w') as f_out:
            for line in f_in:
                query, target = line.split('\t')[:2]
                f_out.write(f"{query}\t{target}\n")
        
        # Compress uniref.map to uniref.map.xz
        run_command(["xz", "-z", uniref_map])
        
        # Log success
        write_json_log(
            {
                "status": "success",
                "process": "functiondb",
                "input_file": protein_file,
                "output_files": {
                    "diamond_results": diamond_out,
                    "uniref_map": uniref_map + ".xz"
                }
            },
            outdir
        )
        
    except Exception as e:
        write_json_log(
            {
                "status": "error",
                "process": "functiondb",
                "error": str(e)
            },
            outdir
        )
        raise

def get_coords(metadata_path, outdir):
    """
    Extract coordinates from protein files.
    
    Args:
        metadata_path (str): Path to metadata file with genome_id and lgenome_dir columns
        outdir (str): Output directory for coordinates file
    """
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

def merge_uniref(uniref90, uniref50, outfile, simplify=False):
    """
    Merge UniRef90 (preferred) and UniRef50 maps.

    Args:
        uniref90 (str): Path to UniRef90 mapping file
        uniref50 (str): Path to UniRef50 mapping file
        outfile (str): Path to output merged map file
        simplify (bool): If True, simplify UniRef IDs by taking only part after '_'
    """
    res = {}
    res_add = res.setdefault
    Path(outfile).parent.mkdir(parents=True, exist_ok=True)

    process_entry = lambda entry: entry.split("_")[1] if simplify else entry

    # Read UniRef90 file (preferred mappings)
    with open(uniref90, "r") as f:
        for line in f:
            orf, entry = line.rstrip().split("\t")[:2]
            g, i = orf.split("_")
            res_add(g, {})[int(i)] = process_entry(entry)

    # Read UniRef50 file (supplementary mappings)
    with open(uniref50, "r") as f:
        for line in f:
            orf, entry = line.rstrip().split("\t")[:2]
            g, i = orf.split("_")
            if g not in res or int(i) not in res[g]:
                res_add(g, {})[int(i)] = process_entry(entry)

    # Write merged output
    with open(outfile, "w") as out_f:
        for g, entries in sorted(res.items()):
            for i, entry in sorted(entries.items()):
                out_f.write(f"{g}_{i}\t{entry}\n")
