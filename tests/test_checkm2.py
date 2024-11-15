import unittest
from datetime import datetime
from unittest.mock import patch
from wipe.modules.checkm2 import (
    gen_command_checkm2,
    gen_stats_data_checkm2,
    run_checkm2_single,
    run_checkm2_batch,
)


class LinearizeTests(unittest.TestCase):
    def test_gen_command_checkm2(self):
        test_infile = "/test/infile"
        test_threads = 4
        test_outdir = "/test/checkm2_out"
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
            "/test/checkm2_out",
            "-t",
            "4",
            "--force",
            "--database_path",
            "/test/db/path",
            "--remove_intermediates",
        ]
        self.assertEqual(obs, exp)

    def test_gen_command_checkm2_genes(self):
        test_infile = "/test/infile"
        test_threads = 4
        test_outdir = "/test/checkm2_out"
        test_dbpath = "/test/db/path"
        test_genes = "/test/genes"
        obs = gen_command_checkm2(
            test_infile, test_threads, test_outdir, test_dbpath, test_genes
        )
        exp = [
            "checkm2",
            "predict",
            "-i",
            "/test/infile",
            "-o",
            "/test/checkm2_out",
            "-t",
            "4",
            "--force",
            "--database_path",
            "/test/db/path",
            "--remove_intermediates",
            "--genes",
            "/test/genes",
        ]
        self.assertEqual(obs, exp)

    def test_gen_command_checkm2_gz(self):
        test_infile = "/test/infile.gz"
        test_threads = 4
        test_outdir = "/test/checkm2_out"
        test_dbpath = "/test/db/path"
        test_genes = None
        obs = gen_command_checkm2(
            test_infile, test_threads, test_outdir, test_dbpath, test_genes
        )
        exp = [
            "checkm2",
            "predict",
            "-i",
            "/test/infile.gz",
            "-o",
            "/test/checkm2_out",
            "-t",
            "4",
            "--force",
            "--database_path",
            "/test/db/path",
            "--remove_intermediates",
            "-x",
            "gz",
        ]
        self.assertEqual(obs, exp)

    @patch("wipe.modules.utils.datetime")
    def test_gen_stats_data(self, mock_datetime):
        mock_datetime.now.return_value = datetime(2024, 11, 13, 15, 30, 0)
        inpath = "/path/to/input/file.fna.gz"
        outdir = "/path/to/input/checkm2"

        obs = gen_stats_data_checkm2(inpath, outdir)
        exp = {
            "genome_id": "file",
            "process": "checkm2_run",
            "start_time": "2024-11-13 15:30:00",
            "end_time": None,
            "status": "in_progress",
            "details": {
                "input_file": "/path/to/input/file.fna.gz",
                "output_files": {
                    "checkm2_report": "/path/to/input/checkm2/quality_report.tsv"
                },
            },
            "error": "no error",
        }
        self.assertTrue(obs, exp)

    # comment out because of thread-creating issue on barnacle2 and abs dbpath
    # def test_run_checkm2_single(self):
    #     test_infile = "./tests/data/999/M000000999.fa.gz"
    #     test_outdir = "./tests/out/999/checkm2_out"
    #     test_dbpath = (
    #         "/home/y1weng/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd"
    #     )
    #     test_threads = 4
    #     run_checkm2_single(test_infile, test_outdir, test_dbpath, test_threads)

    # # comment out because of thread-creating issue on barnacle2 and abs dbpath
    # def test_run_checkm2_batch(self):
    #     test_indir = "./tests/data/999/"
    #     test_logdir = "./tests/out/"
    #     test_dbpath = (
    #         "/home/y1weng/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd"
    #     )
    #     test_threads = 4
    #     run_checkm2_batch(test_indir, test_logdir, test_dbpath, test_threads)


if __name__ == "__main__":
    unittest.main()
