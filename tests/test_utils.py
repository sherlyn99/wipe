import os
import gzip
import json
import shutil
import tempfile
import unittest
import pandas as pd
import pandas.testing as pdt
from pathlib import Path
from datetime import datetime
from unittest.mock import patch
from wipe.modules.utils import (
    load_md,
    check_duplicated_genome_ids,
    write_json_log,
    extract_gid_from_inpath,
    decompress,
    check_outputs,
    check_log_and_retrieve_gid,
    get_files,
    gen_path,
    gen_stats_data,
    gen_summary_data,
    gen_output_paths,
)


class UtilsTests(unittest.TestCase):
    def test_load_md_with_gid(self):
        test_md_path = "./tests/data/metadata.tsv"
        obs = load_md(test_md_path)
        exp = pd.DataFrame(
            {
                "filepath": [
                    "./tests/data/GCF_000981955.1_ASM98195v1_genomic",
                    "./tests/data/GCF_000981955.1_ASM98195v1_genomic2",
                ],
                "genome_id": ["G000000001", "G000000002"],
            }
        )
        pdt.assert_frame_equal(obs, exp)

    def test_load_md_without_gid(self):
        test_md_path = "./tests/data/metadata_nogid.tsv"
        obs = load_md(test_md_path)
        exp = pd.DataFrame(
            {
                "filepath": [
                    "./tests/data/GCF_000981955.1_ASM98195v1_genomic",
                ],
            }
        )
        pdt.assert_frame_equal(obs, exp)

    def test_check_duplicated_genome_ids(self):
        test_md_df = pd.DataFrame(
            {
                "filepath": [
                    "file/path/1",
                    "file/path/2",
                    "file/path/3",
                    "file/path/4",
                    "file/path/5",
                    "file/path/6",
                ],
                "genome_id": [
                    "G001",
                    "G002",
                    "G002",
                    "G002",
                    "G003",
                    "G003",
                ],
            }
        )
        with self.assertRaises(ValueError) as context:
            obs = check_duplicated_genome_ids(test_md_df)
        self.assertEqual(
            str(context.exception),
            f"Duplicated genome ids found:\n{test_md_df.iloc[[2, 3, 5]]}",
        )

    def test_check_duplicated_genome_ids_noduplication(self):
        test_md_df = pd.DataFrame(
            {
                "filepath": [
                    "file/path/1",
                    "file/path/2",
                    "file/path/3",
                    "file/path/4",
                    "file/path/5",
                    "file/path/6",
                ],
                "genome_id": [
                    "G001",
                    "G002",
                    "G003",
                    "G004",
                    "G005",
                    "G006",
                ],
            }
        )

        check_duplicated_genome_ids(test_md_df)

    def test_write_json_log(self):
        test_log_data = {"field1": "value1", "field2": "value2\nmore stuff"}
        test_outdir = tempfile.mkdtemp()
        test_filename = "linearization.err.gz"
        test_append = False

        try:
            write_json_log(
                test_log_data, test_outdir, test_filename, test_append
            )
            with gzip.open(
                f"{test_outdir}/linearization.err.gz", "rt"
            ) as file:
                obs = json.load(file)
            self.assertEqual(len(obs), 1)
            self.assertEqual(len(obs[0]), 3)
            self.assertIn("field1", obs[0])
            self.assertEqual(obs[0]["field1"], "value1")

            self.assertIn("field2", obs[0])
            self.assertEqual(obs[0]["field2"], "value2\nmore stuff")
        finally:
            shutil.rmtree(test_outdir)

    def test_extract_gid_from_inpath(self):
        test_inpath = "/panfs/y1weng/41_marine_data/db/24195.GOMC_genomes_linearized/M/000/002/999/M000002999.fa.gz"
        obs = extract_gid_from_inpath(test_inpath)
        exp = "M000002999"
        self.assertEqual(obs, exp)

    def test_extract_gid_from_inpath_xz(self):
        test_inpath = "M000002999.fa.xz"
        obs = extract_gid_from_inpath(test_inpath)
        exp = "M000002999"
        self.assertEqual(obs, exp)

    def test_decompress(self):
        test_inpath = "./tests/data/GCF_000981955.1_ASM98195v1_genomic.fna.gz"
        test_gid = "G000981955"
        test_tmpdir = tempfile.mkdtemp()
        try:
            obs = decompress(test_inpath, test_gid, test_tmpdir)
            exp = f"{test_tmpdir}/G000981955.fna"
            self.assertEqual(obs, exp)
            self.assertTrue(os.path.exists(exp))
        finally:
            shutil.rmtree(test_tmpdir)

    def test_check_outputs(self):
        test_filelist = [
            "./tests/data/GCF_000981955.1_ASM98195v1_genomic.fna.gz",
            "./tests/data/metadata.tsv",
        ]
        self.assertTrue(check_outputs(test_filelist))

    def test_check_outputs(self):
        test_filelist = [
            "./tests/data/GCF_000981955.1_ASM98195v1_genomic.fna.gz",
            "./tests/data/metadata2.tsv",
        ]
        self.assertFalse(check_outputs(test_filelist))

    def test_check_log_and_retrieve_gid(self):
        test_filepath = "./tests/data/checkm2_stats_M000000999.json.gz"
        obs_bool, obs_gid = check_log_and_retrieve_gid(test_filepath)
        self.assertTrue(obs_bool)
        self.assertEqual(obs_gid, "M000000999")

    def test_check_log_and_retrieve_gid_failed(self):
        test_filepath = "./tests/data/prodigal_stats_fail.json.gz"
        obs_bool, obs_gid = check_log_and_retrieve_gid(test_filepath)
        self.assertFalse(obs_bool)
        self.assertEqual(obs_gid, "GCF_000981955.1_ASM98195v1_genomic")

    def test_get_files(self):
        suffix = ".fa"
        indir = tempfile.mkdtemp()

        try:
            files = [
                Path(indir) / "file1.fa",
                Path(indir) / "subdir" / "file2.fa",
                Path(indir) / "file3.md",
                Path(indir) / "subdir" / ".ipynb_checkpoints" / "file4.fa",
            ]
            for file in files:
                file.parent.mkdir(parents=True, exist_ok=True)
                file.touch()

            obs = get_files(indir, suffix)

            exp = [
                f"{indir}/file1.fa",
                f"{indir}/subdir/file2.fa",
            ]
            self.assertEqual(obs, exp)
        finally:
            shutil.rmtree(indir)

    def test_gen_path(self):
        logdir = "/path/to/log"
        filename = "stats_G001.json.gz"
        obs = gen_path(logdir, filename)
        exp = "/path/to/log/stats_G001.json.gz"
        self.assertEqual(obs, exp)

    @patch("wipe.modules.utils.datetime")
    def test_gen_stats_data(self, mock_datetime):
        mock_datetime.now.return_value = datetime(2024, 11, 13, 15, 30, 0)
        process = "prodigal_run"
        inpath = "/path/to/input/file.fna"

        obs = gen_stats_data(process, inpath)
        exp = {
            "genome_id": "file",
            "process": process,
            "start_time": "2024-11-13 15:30:00",
            "end_time": None,
            "status": "in_progress",
            "details": {
                "input_file": inpath,
                "output_files": {},
            },
            "error": "no error",
        }
        self.assertTrue(obs, exp)

    @patch("wipe.modules.utils.datetime")
    def test_gen_summary_data(self, mock_datetime):
        mock_datetime.now.return_value = datetime(2024, 11, 13, 15, 30, 0)
        process = "prodigal_run"

        obs = gen_summary_data(process)
        exp = {
            "process": process,
            "start_time": "2024-11-13 15:30:00",
            "end_time": None,
            "status": "in_progress",
            "error": "no error",
        }
        self.assertEqual(obs, exp)

    def test_gen_output_paths(self):
        process = "checkm2"
        inpath = "/G/001/002/003/G001002003.fna.gz"
        outdir = "/G/001/002/003/"
        obs_outdir_path, obs_stats_path, obs_gid = gen_output_paths(
            process, inpath, outdir
        )
        self.assertEqual(
            obs_outdir_path, "/G/001/002/003/checkm2_out_G001002003"
        )
        self.assertEqual(
            obs_stats_path, "/G/001/002/003/checkm2_stats_G001002003.json.gz"
        )
        self.assertEqual(obs_gid, "G001002003")

    def test_gen_output_paths_nonncbi(self):
        process = "checkm2"
        inpath = "/data/my_genome.fna.gz"
        outdir = "/out/"
        obs_outdir_path, obs_stats_path, obs_gid = gen_output_paths(
            process, inpath, outdir
        )
        self.assertEqual(obs_outdir_path, "/out/checkm2_out_my_genome")
        self.assertEqual(
            obs_stats_path, "/out/checkm2_stats_my_genome.json.gz"
        )
        self.assertEqual(obs_gid, "my_genome")


if __name__ == "__main__":
    unittest.main()
