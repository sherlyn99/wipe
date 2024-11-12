import unittest
import pandas as pd
import pandas.testing as pdt
from wipe.modules.metadata import generate_gids, generate_metadata


class MetadataTests(unittest.TestCase):
    def test_generate_gids(self):
        test_start_gid = "M000001"
        test_n = 5
        obs = generate_gids(test_start_gid, test_n)
        exp = ["M000001", "M000002", "M000003", "M000004", "M000005"]
        self.assertEqual(obs, exp)

    def test_generate_metadata(self):
        test_indir = "./tests/data"
        test_ext = "fna"
        test_gid_start = "M00000001"
        obs = generate_metadata(test_indir, test_ext, test_gid_start)
        exp = pd.DataFrame(
            {
                "filepath": [
                    "/home/y1weng/47_wipe/wipe/tests/data/GCF_000981955.1_ASM98195v1_genomic2",
                    "/home/y1weng/47_wipe/wipe/tests/data/GCF_000981955.1_ASM98195v1_genomic",
                    "/home/y1weng/47_wipe/wipe/tests/data/GCF_000981956.1_ASM98195v1_genomic_empty",
                ],
                "genome_id": ["M00000001", "M00000002", "M00000003"],
            }
        )
        pdt.assert_frame_equal(obs, exp)

    def test_generate_metadata_nogid(self):
        test_indir = "./tests/data"
        test_ext = "fna"
        obs = generate_metadata(test_indir, test_ext)
        exp = pd.DataFrame(
            {
                "filepath": [
                    "/home/y1weng/47_wipe/wipe/tests/data/GCF_000981955.1_ASM98195v1_genomic2",
                    "/home/y1weng/47_wipe/wipe/tests/data/GCF_000981955.1_ASM98195v1_genomic",
                    "/home/y1weng/47_wipe/wipe/tests/data/GCF_000981956.1_ASM98195v1_genomic_empty",
                ],
            }
        )
        pdt.assert_frame_equal(obs, exp)


if __name__ == "__main__":
    unittest.main()
