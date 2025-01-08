import os
import subprocess
from pathlib import Path
from wipe.modules.utils import write_json_log, run_command

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
        # This assumes the diamond output is in standard blast tabular format
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
        # Log error
        write_json_log(
            {
                "status": "error",
                "process": "functiondb",
                "error": str(e)
            },
            outdir
        )
        raise 