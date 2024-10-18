import os
import shutil
import unittest
import tempfile
from wipe.modules.linearize import (
    generate_gap_string,
    filter_contig_name,
    linearize_single_genome,
    linearize_genomes,
)


class LinearizeTests(unittest.TestCase):
    def test_generate_gap_string(self):
        exp = "NNNNNNNNNNNNNNNNNNNN"
        obs = generate_gap_string("N*20")
        self.assertEqual(exp, obs)

    def test_generate_gap_string_erroneous_input(self):
        exp = ""
        obs = generate_gap_string("20N")
        self.assertEqual(exp, obs)

    def test_filter_contig_name_plasmid(self):
        exp = True
        obs = filter_contig_name("plasmid,phage", ">contig1_plasmid_C1")
        self.assertEqual(exp, obs)

    def test_filter_contig_name_phage(self):
        exp = True
        obs = filter_contig_name("plasmid,phage", ">contig2_phage_C2")
        self.assertEqual(exp, obs)

    def test_filter_contig_name_ignore_case(self):
        exp = True
        obs = filter_contig_name("plasmid,phage", ">contig2_phAge_C2")
        self.assertEqual(exp, obs)

    def test_linearize_single_genome_ncbi_noconcat_nofilt(self):
        test_inpath = "./tests/data/GCF_000981955.1_ASM98195v1_genomic"
        test_ext = ".fna"
        test_outdir = tempfile.mkdtemp()
        test_gid = None
        test_gap = None
        test_filt = None

        try:
            linearize_single_genome(
                test_inpath,
                test_ext,
                test_outdir,
                test_gid,
                test_gap,
                test_filt,
            )

            outpath = os.path.join(test_outdir, "G000981955", "G000981955.fna")
            self.assertTrue(os.path.exists(outpath), "Output file not found.")

            with open(outpath, "r") as f:
                obs = f.read()
                exp = (
                    ">G000981955_contig_1_plasmid_at_the_start\n"
                    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
                    ">G000981955_contig_2\n"
                    "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
                    ">G000981955_contig_3_phage_in_the_middle\n"
                    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
                    ">G000981955_contig_3\n"
                    "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
                    ">G000981955_contig_4_plasmid_at_the_end\n"
                    "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
                )
                self.assertEqual(
                    obs,
                    exp,
                    "Output content does not match expected content.",
                )

            logpath = os.path.join(
                test_outdir, "G000981955", "linearization.log"
            )
            with open(logpath, "r") as f:
                obs = f.read()
                exp = f"Wrote 5 contigs (240 characters in total) into {test_outdir}/G000981955/G000981955.fna\n"
                self.assertAlmostEqual(
                    obs, exp, "Log does not match expected log content."
                )
        finally:
            shutil.rmtree(test_outdir)

    def test_linearize_single_genome_ncbi_noconcat_nofilt_empty_contig(self):
        test_inpath = "./tests/data/GCF_000981955.1_ASM98195v1_genomic_empty"
        test_ext = ".fna"
        test_outdir = tempfile.mkdtemp()
        test_gid = None
        test_gap = None
        test_filt = None

        try:
            linearize_single_genome(
                test_inpath,
                test_ext,
                test_outdir,
                test_gid,
                test_gap,
                test_filt,
            )

            outpath = os.path.join(test_outdir, "G000981955", "G000981955.fna")
            self.assertTrue(os.path.exists(outpath), "Output file not found.")

            with open(outpath, "r") as f:
                obs = f.read()
                exp = (
                    ">G000981955_contig_1_plasmid_at_the_start\n"
                    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
                    ">G000981955_contig_2\n"
                    "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
                    ">G000981955_contig_3_phage_in_the_middle\n"
                    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
                    ">G000981955_contig_3\n"
                    "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
                    ">G000981955_contig_4_plasmid_at_the_end\n"
                    "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
                )
                self.assertEqual(
                    obs,
                    exp,
                    "Output content does not match expected content.",
                )

            logpath = os.path.join(
                test_outdir, "G000981955", "linearization.log"
            )
            with open(logpath, "r") as f:
                obs = f.read()
                exp = (
                    f"Wrote 5 contigs (240 characters in total) into {test_outdir}/G000981955/G000981955.fna\n"
                    "\n3 contigs were empty:\n['>contig_empty_1', '>contig_empty_2', '>contig_empty_5']\n"
                )
                self.assertAlmostEqual(
                    obs, exp, "Log does not match expected log content."
                )
        finally:
            shutil.rmtree(test_outdir)

    def test_linearize_single_genome_ncbi_noconcat_filt_empty_contig(self):
        test_inpath = "./tests/data/GCF_000981955.1_ASM98195v1_genomic_empty"
        test_ext = ".fna"
        test_outdir = tempfile.mkdtemp()
        test_gid = None
        test_gap = None
        test_filt = "plasmid,phage"

        try:
            linearize_single_genome(
                test_inpath,
                test_ext,
                test_outdir,
                test_gid,
                test_gap,
                test_filt,
            )

            outpath = os.path.join(test_outdir, "G000981955", "G000981955.fna")
            self.assertTrue(os.path.exists(outpath), "Output file not found.")

            with open(outpath, "r") as f:
                obs = f.read()
                exp = (
                    ">G000981955_contig_2\n"
                    "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
                    ">G000981955_contig_3\n"
                    "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
                )
                self.assertEqual(
                    obs,
                    exp,
                    "Output content does not match expected content.",
                )

            logpath = os.path.join(
                test_outdir, "G000981955", "linearization.log"
            )
            with open(logpath, "r") as f:
                obs = f.read()
                exp = (
                    f"Wrote 2 contigs (90 characters in total) into {test_outdir}/G000981955/G000981955.fna\n"
                    "\n3 contigs were empty:\n['>contig_empty_1', '>contig_empty_2', '>contig_empty_5']\n"
                    "\n3 contigs were filtered.\n['>contig_1_plasmid_at_the_start', '>contig_3_phage_in_the_middle', '>contig_4_plasmid_at_the_end']\n"
                )
                self.assertAlmostEqual(
                    obs, exp, "Log does not match expected log content."
                )
        finally:
            shutil.rmtree(test_outdir)

    def test_linearize_single_genome_ncbi_concat_filt_empty_contig(self):
        test_inpath = "./tests/data/GCF_000981955.1_ASM98195v1_genomic_empty"
        test_ext = ".fna"
        test_outdir = tempfile.mkdtemp()
        test_gid = None
        test_gap = "N*20"
        test_filt = "plasmid,phage"

        try:
            linearize_single_genome(
                test_inpath,
                test_ext,
                test_outdir,
                test_gid,
                test_gap,
                test_filt,
            )

            outpath = os.path.join(test_outdir, "G000981955", "G000981955.fna")
            self.assertTrue(os.path.exists(outpath), "Output file not found.")

            with open(outpath, "r") as f:
                obs = f.read()
                exp = (
                    ">G000981955\n"
                    "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
                    "NNNNNNNNNNNNNNNNNNNN"
                    "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
                )
                self.assertEqual(
                    obs,
                    exp,
                    "Output content does not match expected content.",
                )

            logpath = os.path.join(
                test_outdir, "G000981955", "linearization.log"
            )
            with open(logpath, "r") as f:
                obs = f.read()
                exp = (
                    f"Wrote 2 contigs (90 characters in total) into {test_outdir}/G000981955/G000981955.fna\n"
                    "\n3 contigs were empty:\n['>contig_empty_1', '>contig_empty_2', '>contig_empty_5']\n"
                    "\n3 contigs were filtered.\n['>contig_1_plasmid_at_the_start', '>contig_3_phage_in_the_middle', '>contig_4_plasmid_at_the_end']\n"
                )
                self.assertAlmostEqual(
                    obs, exp, "Log does not match expected log content."
                )
        finally:
            shutil.rmtree(test_outdir)

    def test_linearize_single_genome_nonncbi_concat_filt_empty_contig(self):
        test_inpath = "./tests/data/GCF_000981955.1_ASM98195v1_genomic_empty"
        test_ext = ".fna"
        test_outdir = tempfile.mkdtemp()
        test_gid = "G000000001"
        test_gap = "N*20"
        test_filt = "plasmid,phage"

        try:
            linearize_single_genome(
                test_inpath,
                test_ext,
                test_outdir,
                test_gid,
                test_gap,
                test_filt,
            )

            outpath = os.path.join(test_outdir, "G000000001", "G000000001.fna")
            self.assertTrue(os.path.exists(outpath), "Output file not found.")

            with open(outpath, "r") as f:
                obs = f.read()
                exp = (
                    ">G000000001\n"
                    "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
                    "NNNNNNNNNNNNNNNNNNNN"
                    "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
                )
                self.assertEqual(
                    obs,
                    exp,
                    "Output content does not match expected content.",
                )

            logpath = os.path.join(
                test_outdir, "G000000001", "linearization.log"
            )
            with open(logpath, "r") as f:
                obs = f.read()
                exp = (
                    f"Wrote 2 contigs (90 characters in total) into {test_outdir}/G000000001/G000000001.fna\n"
                    "\n3 contigs were empty:\n['>contig_empty_1', '>contig_empty_2', '>contig_empty_5']\n"
                    "\n3 contigs were filtered.\n['>contig_1_plasmid_at_the_start', '>contig_3_phage_in_the_middle', '>contig_4_plasmid_at_the_end']\n"
                )
                self.assertAlmostEqual(
                    obs, exp, "Log does not match expected log content."
                )
        finally:
            shutil.rmtree(test_outdir)

    def test_linearize_genomes(self):
        test_metadata = "./tests/data/metadata.tsv"
        test_ext = "fna"
        test_outdir = tempfile.mkdtemp()
        test_gap = "N*20"
        test_filt = None
        try:
            linearize_genomes(
                test_metadata, test_ext, test_outdir, test_gap, test_filt
            )

            paths = [
                os.path.join(test_outdir, "G000000001", "G000000001.fna"),
                os.path.join(test_outdir, "G000000001", "linearization.log"),
                os.path.join(test_outdir, "G000000002", "G000000002.fna"),
                os.path.join(test_outdir, "G000000001", "linearization.log"),
            ]

            for p in paths:
                self.assertTrue(os.path.exists(p), "Output file not found.")

        finally:
            shutil.rmtree(test_outdir)

    def test_linearize_genomes_nogid(self):
        test_metadata = "./tests/data/metadata_nogid.tsv"
        test_ext = "fna"
        test_outdir = tempfile.mkdtemp()
        test_gap = "N*20"
        test_filt = None
        try:
            linearize_genomes(
                test_metadata, test_ext, test_outdir, test_gap, test_filt
            )

            paths = [
                os.path.join(test_outdir, "G000981955", "G000981955.fna"),
                os.path.join(test_outdir, "G000981955", "linearization.log"),
            ]

            for p in paths:
                self.assertTrue(os.path.exists(p), "Output file not found.")

        finally:
            shutil.rmtree(test_outdir)

    def test_linearize_genomes_repeated_gid(self):
        test_metadata = "./tests/data/metadata_repeated_gid.tsv"
        test_ext = "fna"
        test_outdir = tempfile.mkdtemp()
        test_gap = "N*20"
        test_filt = None

        try:
            linearize_genomes(
                test_metadata, test_ext, test_outdir, test_gap, test_filt
            )

            paths = [
                os.path.join(test_outdir, "G000000001", "G000000001.fna"),
                os.path.join(test_outdir, "G000000001", "linearization.log"),
                os.path.join(test_outdir, "linearization_all.log"),
            ]
            for p in paths:
                self.assertTrue(os.path.exists(p), "Output file not found.")

            with open(
                os.path.join(test_outdir, "linearization_all.log"), "r"
            ) as f:
                obs = f.read()
                exp = "Dupliacted genome ID: G000000001 with input path ./tests/data/GCF_000981955.1_ASM98195v1_genomic\n"
                self.assertAlmostEqual(
                    obs, exp, "Log does not match expected log content."
                )
        finally:
            shutil.rmtree(test_outdir)


if __name__ == "__main__":
    unittest.main()
