import os
import unittest
import tempfile
import pandas as pd
from wipe.modules.taxonomy import (
    process_search_results,
)


class TaxonomyTests(unittest.TestCase):

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
