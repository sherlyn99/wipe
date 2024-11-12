import os
import shutil
import unittest
import tempfile
from wipe.modules.prodigal import (
    gen_commands_prodigal,
    gen_log_data_prodigal,
    gen_log_filepath_prodigal,
    run_prodigal_single_genome,
)
from wipe.modules.utils import check_outputs


class ProdigalTests(unittest.TestCase):
    def test_gen_commands_prodigal(self):
        test_infna = "./tests/tmp/GCF_000981955.1_ASM98195v1_genomic.fna"  # decompressed file in /tmp
        test_gid = "G000981955"
        test_outdir = "./tests/out/G/000/981/955"  # distributed output
        obs = gen_commands_prodigal(test_infna, test_gid, test_outdir)

        exp = [
            "prodigal",
            "-p",
            "single",
            "-g",
            "11",
            "-m",
            "-i",
            "./tests/tmp/GCF_000981955.1_ASM98195v1_genomic.fna",
            "-f",
            "gff",
            "-o",
            "./tests/out/G/000/981/955/G000981955.gff",
            "-a",
            "./tests/out/G/000/981/955/G000981955.faa",
            "-d",
            "./tests/out/G/000/981/955/G000981955.ffn",
        ]
        self.assertEqual(obs, exp)

    def test_gen_log_data_prodigal(self):
        test_inpath = "./tests/data/GCF_000981955.1_ASM98195v1_genomic.fna.gz"
        test_infna = "./tests/tmp/GCF_000981955.1_ASM98195v1_genomic.fna"
        test_gid = "G000981955"
        test_outdir = "./tests/out/G/000/981/955"
        obs = gen_log_data_prodigal(
            test_inpath, test_infna, test_gid, test_outdir
        )
        exp = {
            "genome_id": "G000981955",
            "process": "prodigal_run",
            "status": "in_progress",
            "details": {
                "input_file_compressed": "./tests/data/GCF_000981955.1_ASM98195v1_genomic.fna.gz",
                "input_file": "./tests/tmp/GCF_000981955.1_ASM98195v1_genomic.fna",
                "output_files": {
                    "gff": "./tests/out/G/000/981/955/G000981955.gff.xz",
                    "faa": "./tests/out/G/000/981/955/G000981955.faa.xz",
                    "ffn": "./tests/out/G/000/981/955/G000981955.ffn.xz",
                },
            },
            "error": None,
        }
        for key, val in exp.items():
            self.assertEqual(obs[key], val)
        self.assertTrue("start_time" in obs.keys())
        self.assertTrue("end_time" in obs.keys())

    def test_gen_log_filepath_prodigal(self):
        obs = gen_log_filepath_prodigal("./G/000/000/001/")
        exp = "./G/000/000/001/prodigal_stats.json.gz"
        self.assertEqual(obs, exp)

    def test_gen_log_filepath_prodigal(self):
        obs = gen_log_filepath_prodigal("./G/000/000/001/", "G000000001")
        exp = "./G/000/000/001/prodigal_stats_G000000001.json.gz"
        self.assertEqual(obs, exp)

    def test_run_prodigal_single_genome(self):
        test_inpath = "./tests/data/GCF_000981955.1_ASM98195v1_genomic.fna.gz"
        test_outdir = tempfile.mkdtemp()
        test_tmp = "./tests/tmp"
        filelist = [
            os.path.join(test_outdir, "prodigal_out/"),
            os.path.join(test_outdir, "prodigal_stats.json.gz"),
        ]
        try:
            obs = run_prodigal_single_genome(test_inpath, test_outdir, test_tmp)
            self.assertTrue(check_outputs(filelist))
        finally:
            shutil.rmtree(filelist[0])
            os.remove(filelist[1])


if __name__ == "__main__":
    unittest.main()
