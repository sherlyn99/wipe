import unittest
from wipe.modules.gsearch import (
    create_db,
    update_db,
    search_genomes,
)


class TaxonomyTests(unittest.TestCase):
    pass
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
    #     outpath = reformat_search_results(inpath)


if __name__ == "__main__":
    unittest.main()
