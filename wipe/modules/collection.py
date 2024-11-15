import os
import lzma
import skbio
import pandas as pd
from glob import glob
from wipe.modules.utils import search_dirs


def compile_results_checkm2(indir):
    dirs = search_dirs(indir, "checkm2_out*")
    dfs = []

    for dir in dirs:
        quality_report_path = os.path.join(dir, "quality_report.tsv")
        if os.path.exists(quality_report_path):
            df = pd.read_csv(quality_report_path, sep="\t")
            dfs.append(df)

    combined_df = pd.concat(dfs, ignore_index=True)
    combined_df = combined_df.rename(columns={"Name": "genome_id"})
    return combined_df


def compile_results_prodigal(indir, coords=False, outdir=None):
    dirs = search_dirs(indir, "prodigal_out*")
    results = {"genome_id": [], "protein_ct": {}}

    if coords:
        if outdir is None:
            raise ValueError(
                "Output directory must be specified when coords=True."
            )
        coords_fp = os.path.join(outdir, "coords.txt.xz")
        with lzma.open(coords_fp, "wt") as coords_fo:
            for dir in dirs:
                genome_id = os.path.basename(dir).split("_")[-1]
                protein_ct = 0
                faa_xz_fp = os.path.join(dir, f"{genome_id}.faa.xz")
                with lzma.open(faa_xz_fp, "rb") as fi:
                    current_nucl = None
                    for seq in skbio.read(fi, format="fasta"):

                        header = seq.metadata["id"]
                        name, pos5, pos3, strand, _ = header.split(" # ", 4)
                        nucl, idx = name.rsplit("_", 1)
                        beg, end = (
                            (pos5, pos3) if strand == "1" else (pos3, pos5)
                        )
                        if nucl != current_nucl:
                            current_nucl = nucl
                            coords_fo.write(f">{nucl}\n")

                        coords_fo.write(f"{idx}\t{beg}\t{end}\n")
                        protein_ct += 1

                results["genome_id"].append(genome_id)
                results["protein_ct"].append(protein_ct)

    protein_ct_df = pd.DataFrame(results)
    return protein_ct_df


if __name__ == "__main__":
    pass
