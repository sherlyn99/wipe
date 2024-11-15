import os
import lzma
import tempfile
import unittest
from datetime import datetime
from unittest.mock import patch
from wipe.modules.prodigal import (
    decompress,
    gen_commands_prodigal,
    gen_stats_data_prodigal,
    get_coords_and_proteins,
    gen_summary_path_prodigal,
    gen_summary_data_prodigal,
    run_prodigal_single,
    run_prodigal_batch,
)


class ProdigalTests(unittest.TestCase):
    def test_decompress(self):
        inpath = "./tests/data/999/M000000999.fa.gz"
        tmpdir = tempfile.mkdtemp()
        try:
            outpath, to_del = decompress(inpath, tmpdir)
            self.assertEqual(f"{tmpdir}/M000000999.fna", outpath)
            self.assertTrue(to_del)
        finally:
            os.remove(outpath)

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

    @patch("wipe.modules.utils.datetime")
    def test_gen_stats_data_prodigal(self, mock_datetime):
        inpath = "/path/to/my/genome.fna.gz"
        outdir_path = "/path/to/outdir/prodigal_out"
        mock_datetime.now.return_value = datetime(2024, 11, 13, 15, 30, 0)
        obs = gen_stats_data_prodigal(inpath, outdir_path)
        exp = {
            "genome_id": "genome",
            "process": "prodigal_run",
            "start_time": "2024-11-13 15:30:00",
            "end_time": None,
            "status": "in_progress",
            "details": {
                "input_file": "/path/to/my/genome.fna.gz",
                "output_files": {
                    "faa": "/path/to/outdir/prodigal_out/genome.faa.xz",
                    "ffn": "/path/to/outdir/prodigal_out/genome.ffn.xz",
                    "gff": "/path/to/outdir/prodigal_out/genome.gff.xz",
                    "coords": "/path/to/outdir/prodigal_out/coords_genome.txt.xz",
                    "proteins": "/path/to/outdir/prodigal_out/proteins_genome.tsv.xz",
                },
            },
            "error": "no error",
        }
        self.assertEqual(obs, exp)

    def test_get_coords_and_proteins(self):
        faa_fp = "./tests/data/M000000999_subset.faa"
        outdir = tempfile.mkdtemp()
        out_coords_fp = f"{outdir}/coords.txt.xz"
        out_proteins_fp = f"{outdir}/proteins.tsv.xz"
        exp_coords = ">M000000999\n1\t1\t1809\n2\t3187\t1799\n"
        exp_proteins = "genome_id\tproteins_ct\nM000000999\t2\n"
        try:
            get_coords_and_proteins(faa_fp, out_coords_fp, out_proteins_fp)
            with lzma.open(out_coords_fp, "rt") as file:
                obs_coords = file.read()
            with lzma.open(out_proteins_fp, "rt") as file:
                obs_proteins = file.read()

            self.assertEqual(obs_coords, exp_coords)
            self.assertEqual(obs_proteins, exp_proteins)
        finally:
            os.remove(out_coords_fp)
            os.remove(out_proteins_fp)

    def test_gen_summary_path_prodigal(self):
        obs = gen_summary_path_prodigal("/path/to/logdir")
        self.assertEqual(obs, "/path/to/logdir/prodigal_summary.json.gz")

    @patch("wipe.modules.utils.datetime")
    def test_gen_summary_data_prodigal(self, mock_datetime):
        mock_datetime.now.return_value = datetime(2024, 11, 13, 15, 30, 0)
        obs = gen_summary_data_prodigal()
        exp = {
            "process": "prodigal_run",
            "start_time": "2024-11-13 15:30:00",
            "end_time": None,
            "status": "in_progress",
            "error": "no error",
        }
        self.assertEqual(obs, exp)

    # # def test_run_prodigal_single(self):
    # inpath = "/home/y1weng/47_wipe/wipe/tests/data/999/M000000999.fa.gz"
    # outdir = "/home/y1weng/47_wipe/wipe/tests/data/999/"
    # tmpdir = "/home/y1weng/47_wipe/wipe/tests/tmp"
    # run_prodigal_single(inpath, outdir, tmpdir)

    # def test_run_prodigal_batch(self):
    #     indir = "/home/y1weng/47_wipe/wipe/tests/data/999/"
    #     logdir = "/home/y1weng/47_wipe/wipe/tests/tmp"
    #     tmpdir = "/home/y1weng/47_wipe/wipe/tests/tmp"
    #     run_prodigal_batch(indir, logdir, tmpdir)


if __name__ == "__main__":
    unittest.main()
