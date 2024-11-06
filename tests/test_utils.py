import gzip
import json
import shutil
import tempfile
import unittest
import pandas as pd
import pandas.testing as pdt
from wipe.modules.utils import (
    load_md,
    check_duplicated_genome_ids,
    write_json_log,
)


class LinearizeTests(unittest.TestCase):
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


if __name__ == "__main__":
    unittest.main()
