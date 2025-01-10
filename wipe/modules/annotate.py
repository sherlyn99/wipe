import os
import pandas as pd
from wipe.modules.utils import check_required_cols
from wipe.modules.barrnap import run_barrnap_single
from wipe.modules.prodigal import run_prodigal_single
from wipe.modules.kofamscan import run_kofamscan_single

# design
# take single file

# take batch
#   metadata provided


def annotate_single(
    genome_id,
    in_file_fna_gz,
    outdir,
    rrna_cutoff,
    ko_profiles,
    ko_list,
    nthreads,
    tmp_dir,
):
    outdir_barrnap = os.path.join(outdir, "barrnap_out")
    run_barrnap_single(in_file_fna_gz, outdir_barrnap, rrna_cutoff)

    outdir_prodigal = os.path.join(outdir, "prodigal_out")
    in_file_faa_xz = run_prodigal_single(
        genome_id, in_file_fna_gz, outdir_prodigal, tmpdir=tmp_dir
    )

    outdir_kofamscan = os.path.join(outdir, "kofamscan_out")
    if in_file_faa_xz and os.path.exists(in_file_faa_xz):
        run_kofamscan_single(
            genome_id,
            in_file_faa_xz,
            outdir_kofamscan,
            tmp_dir,
            ko_profiles,
            ko_list,
            nthreads,
        )


def annotate_multiple(
    metadata, outdir, rrna_cutoff, ko_profiles, ko_list, tmp_dir, nthreads
):
    df_md = pd.read_csv(metadata, sep="\t", low_memory=False)
    check_required_cols(df_md, ["genome_id", "lgenome_path"])
    summary_path = os.path.join(outdir, "annotation_summary.txt")

    for row in df_md.itertuples():
        genome_id = row.genome_id
        in_file_fna_gz = row.lgenome_path
        outdir_per_genome = os.path.join(
            outdir,
            genome_id[0],
            genome_id[1:4],
            genome_id[4:7],
            genome_id[7:10],
        )
        try:
            annotate_single(
                genome_id,
                in_file_fna_gz,
                outdir_per_genome,
                rrna_cutoff,
                ko_profiles,
                ko_list,
                nthreads,
                tmp_dir,
            )
        except:
            with open(summary_path, "a") as log:
                log.write(f"{genome_id}\n")


# if __name__ == "__main__":
#     annotate_single(
#         "G000964075",
#         "wol3_prototype/output2/G/000/964/075/G000964075.fna.gz",
#         "wol3_prototype/output2/G/000/964/075/",
#         0.67,
#         "/projects/greengenes2/20231117_annotations_prelim/kofam_scan/profiles/test_10.hal",
#         "/projects/greengenes2/20231117_annotations_prelim/kofam_scan/ko_list_10",
#         64,
#         "/ddn_scratch/y1weng",
#     )

# metadata = "./wol3_prototype/output/metadata_fqc.tsv"
# lin_dir = "./wol3_prototype/output"
# rrna_cutoff = 0.67
# tmp_dir = "/panfs/y1weng"
# ko_profiles = "/projects/greengenes2/20231117_annotations_prelim/kofam_scan/profiles/test_10.hal"
# ko_list = "/projects/greengenes2/20231117_annotations_prelim/kofam_scan/ko_list_10"
# nthreads = 64

# annotate_multiple(
#     metadata, lin_dir, rrna_cutoff, ko_profiles, ko_list, tmp_dir, nthreads
# )
