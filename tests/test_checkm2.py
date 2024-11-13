import unittest
from wipe.modules.checkm2 import gen_command_checkm2, run_checkm2_single


class LinearizeTests(unittest.TestCase):
    def test_gen_command_checkm2(self):
        test_infile = "/test/infile"
        test_threads = 4
        test_outdir = "/test/out"
        test_dbpath = "/test/db/path"
        obs = gen_command_checkm2(
            test_infile, test_threads, test_outdir, test_dbpath
        )

        exp = [
            "checkm2",
            "predict",
            "-i",
            "/test/infile",
            "-o",
            "/test/out",
            "-t",
            "4",
            "--force",
            "--database_path",
            "/test/db/path",
            "--remove_intermediates",
        ]
        self.assertEqual(obs, exp)

    def test_gen_command_checkm2_genes(self):
        test_infile = "/test/infile.gz"
        test_threads = 4
        test_outdir = "/test/out"
        test_dbpath = "/test/db/path"
        test_genes = "/test/genes"
        obs = gen_command_checkm2(
            test_infile, test_threads, test_outdir, test_dbpath, test_genes
        )
        exp = [
            "checkm2",
            "predict",
            "-i",
            "/test/infile.gz",
            "-o",
            "/test/out",
            "-t",
            "4",
            "--force",
            "--database_path",
            "/test/db/path",
            "--remove_intermediates",
            "-x",
            "gz",
            "--genes",
            "/test/genes",
        ]
        self.assertEqual(obs, exp)

    # # comment out because of thread-creating issue on barnacle2
    # def test_run_checkm2_single(self):
    #     test_infile = "./tests/data/999/M000000999.fa.gz"
    #     test_outdir = "./tests/data/999/checkm2_out2"
    #     test_dbpath = (
    #         "/home/y1weng/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd"
    #     )
    #     test_threads = 4
    #     run_checkm2_single(test_infile, test_outdir, test_dbpath, test_threads)


if __name__ == "__main__":
    unittest.main()
