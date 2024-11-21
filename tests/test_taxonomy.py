import os
import unittest
import tempfile
import pandas as pd
from wipe.modules.utils import load_assembly
from wipe.modules.taxonomy import (
    create_db,
    update_db,
    search_genomes,
    reformat_search_results,
    process_search_results,
)


class TaxonomyTests(unittest.TestCase):
    # # commented out because of the thread creating issue
    # def test_create_db(self):
    #     indir = "/projects/greengenes2/gg2_genomes/linearized_test/G/964/307"
    #     # outdir = os.path.join(tempfile.mkdtemp(), "gsearch_db"
    #     outdir = "./tests/out/gsearch_db"
    #     create_db(indir, outdir, 1)

    # def test_update_db(self):
    #     existing_db = "/home/y1weng/47_wipe/wipe/tests/out/gsearch_db"
    #     indir = "/projects/greengenes2/gg2_genomes/ncbi/GCF/964/306"
    #     outdir = "./tests/out/gsearch_db_updated_2"
    #     nbthreads = 1
    #     update_db(existing_db, indir, outdir, nbthreads)

    # need to swicth output dir to tempdir
    # def test_search_genomes(self):
    #     db = "/home/y1weng/47_wipe/wipe/tests/out/gsearch_db"
    #     n_neighbors = 5
    #     indir = "/projects/greengenes2/gg2_genomes/ncbi/GCF/964/306"
    #     outdir = "./tests/out/gsearch_results"
    #     outpath = search_genomes(db, n_neighbors, indir, outdir)
    #     print(outpath)

    # def test_reformat_search_results(self):
    #     inpath = "./tests/out/gsearch_results/gsearch.neighbors.txt"
    #     outpath = "./tests/out/gsearch_results/gsearch_results.tsv"
    #     reformat_search_results(inpath, outpath)

    def test_process_search_results(self):
        inpath = "tests/out/gsearch_results/gsearch_results.tsv"
        assembly_with_tax = ""
        outpath = ""

        process_search_results(inpath, assembly_with_tax, outpath)

    # def test_decompress(self):
    #     assembly_df_subset = pd.read_csv("./tests/data/assembly.tsv", sep="\t")

    #     assembly_df_subset["taxonomy"] = assembly_df_subset[
    #         "species_taxid"
    #     ].apply(get_taxonomy_ncbi)
    #     assembly_df_subset.to_csv(
    #         "./tests/tmp/tax.tsv", sep="\t", header=True, index=False
    #     )
    #     # print(assembly_df_subset.shape)

    # def test_create_new_genomes_dir(md, dir):
    #     create_new_genomes_dir()

    # def test_create_db(self):
    #     # errors out due to thread-creating issue in login node, must use srun
    #     indir = "/projects/greengenes2/gg2_genomes/linearized_test/G"
    #     outdir = "./tests/tmp/gsearch_out/"
    #     db_dir = create_db(indir, outdir)
    #     print(db_dir)  # ./tests/tmp/gsearch_out/

    # def test_search_genomes(self):
    #     db_dir = "./tests/tmp/gsearch_out/"
    #     n = 5
    #     indir = "/projects/greengenes2/gg2_genomes/linearized_test/G"
    #     search_genomes(db_dir, n, indir)


if __name__ == "__main__":
    unittest.main()
