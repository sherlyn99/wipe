import os
import shutil
import unittest
import tempfile
import warnings
from skbio.io import FormatIdentificationWarning
from wipe.modules.linearize import (
    generate_gap_string,
    should_filter_contig_name,
    infer_gid_ncbi,
    generate_inpath_outpath,
    generate_log_msg,
    read_fasta,
    linearize_single_genome,
    linearize_genomes,
)
from wipe.modules.utils import logger


class LinearizeTests(unittest.TestCase):
    def test_generate_gap_string(self):
        obs = generate_gap_string("N*20")
        exp = "NNNNNNNNNNNNNNNNNNNN"
        self.assertEqual(obs, exp)

    def test_generate_gap_string_erroneous_input(self):
        obs = generate_gap_string("20N")
        exp = ""
        self.assertEqual(obs, exp)

    def test_should_filter_contig_name_plasmid(self):
        obs = should_filter_contig_name("plasmid,phage", ">contig1_plasmid_C1")
        exp = True
        self.assertEqual(obs, exp)

    def test_should_filter_contig_name_phage(self):
        obs = should_filter_contig_name("plasmid,phage", ">contig2_phage_C2")
        exp = True
        self.assertEqual(obs, exp)

    def test_should_filter_contig_name_at_the_end_of_contig_name(self):
        obs = should_filter_contig_name("plasmid,phage", ">contig2_phage")
        exp = False
        self.assertEqual(obs, exp)

    def test_should_filter_contig_name_ignore_case(self):
        obs = should_filter_contig_name("plasmid,phage", ">contig2_phAge_C2")
        exp = True
        self.assertEqual(obs, exp)

    def test_infer_gid_ncbi_fna(self):
        obs = infer_gid_ncbi("/path/to/GCF_000981955.1_ASM98195v1_genomic.fna")
        exp = "G000981955"
        self.assertEqual(obs, exp)

    def test_infer_gid_ncbi_gz(self):
        obs = infer_gid_ncbi(
            "/path/to/GCF_000981955.1_ASM98195v1_genomic.fa.gz"
        )
        exp = "G000981955"
        self.assertEqual(obs, exp)

    def test_infer_gid_ncbi_xz(self):
        obs = infer_gid_ncbi(
            "/path/to/GCA_000981955.1_ASM98195v1_genomic.fasta.xz"
        )
        exp = "G000981955"
        self.assertEqual(obs, exp)

    def test_infer_gid_ncbi_error(self):
        with self.assertRaises(ValueError) as context:
            obs = infer_gid_ncbi(
                "/path/to/GCF_000981955.1ASM98195v1_genomic.fasta.xz"
            )
        self.assertEqual(
            str(context.exception), "No valid genome ID provided or extracted."
        )

    def test_generate_inpath_outpath_fna_gid(self):
        test_inpath = "/path/to/GCA_000981955.1_ASM98195v1_genomic"
        test_ext = "fna"
        test_outdir = "/path/to/outdir"
        test_gid = "G000000001"
        obs_gid, obs_inpath, obs_outdir, obs_outpath = generate_inpath_outpath(
            test_inpath, test_ext, test_outdir, test_gid
        )
        exp_gid = "G000000001"
        exp_inpath = "/path/to/GCA_000981955.1_ASM98195v1_genomic.fna"
        exp_outdir = "/path/to/outdir/G000000001"
        exp_outpath = "/path/to/outdir/G000000001/G000000001.fna"
        self.assertEqual(obs_gid, exp_gid)
        self.assertEqual(obs_inpath, exp_inpath)
        self.assertEqual(obs_outdir, exp_outdir)
        self.assertEqual(obs_outpath, exp_outpath)

    def test_generate_inpath_outpath_gz_gid(self):
        test_inpath = "/path/to/GCA_000981955.1_ASM98195v1_genomic"
        test_ext = "fa.gz"
        test_outdir = "/path/to/outdir"
        test_gid = "G000000001"
        obs_gid, obs_inpath, obs_outdir, obs_outpath = generate_inpath_outpath(
            test_inpath, test_ext, test_outdir, test_gid
        )
        exp_gid = "G000000001"
        exp_inpath = "/path/to/GCA_000981955.1_ASM98195v1_genomic.fa.gz"
        exp_outdir = "/path/to/outdir/G000000001"
        exp_outpath = "/path/to/outdir/G000000001/G000000001.fa"
        self.assertEqual(obs_gid, exp_gid)
        self.assertEqual(obs_inpath, exp_inpath)
        self.assertEqual(obs_outdir, exp_outdir)
        self.assertEqual(obs_outpath, exp_outpath)

    def test_generate_inpath_outpath_xz_gid(self):
        test_inpath = "/path/to/GCA_000981955.1_ASM98195v1_genomic"
        test_ext = ".fa.xz"
        test_outdir = "/path/to/outdir"
        test_gid = "G000000001"
        obs_gid, obs_inpath, obs_outdir, obs_outpath = generate_inpath_outpath(
            test_inpath, test_ext, test_outdir, test_gid
        )
        exp_gid = "G000000001"
        exp_inpath = "/path/to/GCA_000981955.1_ASM98195v1_genomic.fa.xz"
        exp_outdir = "/path/to/outdir/G000000001"
        exp_outpath = "/path/to/outdir/G000000001/G000000001.fa"
        self.assertEqual(obs_gid, exp_gid)
        self.assertEqual(obs_inpath, exp_inpath)
        self.assertEqual(obs_outdir, exp_outdir)
        self.assertEqual(obs_outpath, exp_outpath)

    def test_generate_inpath_outpath_xz_nogid(self):
        test_inpath = "/path/to/GCA_000981955.1_ASM98195v1_genomic"
        test_ext = ".fa.xz"
        test_outdir = "/path/to/outdir"
        test_gid = None
        obs_gid, obs_inpath, obs_outdir, obs_outpath = generate_inpath_outpath(
            test_inpath, test_ext, test_outdir, test_gid
        )
        exp_gid = "G000981955"
        exp_inpath = "/path/to/GCA_000981955.1_ASM98195v1_genomic.fa.xz"
        exp_outdir = "/path/to/outdir/G000981955"
        exp_outpath = "/path/to/outdir/G000981955/G000981955.fa"
        self.assertEqual(obs_gid, exp_gid)
        self.assertEqual(obs_inpath, exp_inpath)
        self.assertEqual(obs_outdir, exp_outdir)
        self.assertEqual(obs_outpath, exp_outpath)

    def test_generate_log_msg(self):
        test_n_written = 100
        test_n_char = 10000
        test_n_filtered = 50
        test_outpath = "/path/to/outpath"
        obs = generate_log_msg(
            test_n_written, test_n_char, test_n_filtered, test_outpath
        )
        exp = "Wrote 100 contigs (10000 characters in total) into /path/to/outpath.\n50 contigs were filtered.\n"
        self.assertEqual(obs, exp)

    def test_read_fasta(self):
        test_inpath = "./tests/data/GCF_000981955.1_ASM98195v1_genomic.fna"
        obs_seqs = [seq for seq in read_fasta(test_inpath)]
        exp_seqs = [
            ">contig_1_plasmid_at_the_start",
            ">contig_2",
            ">contig_3_phage_in_the_middle",
            ">contig_3",
            ">contig_4_plasmid_at_the_end",
        ]
        self.assertEqual(obs_seqs, exp_seqs)

    def test_read_fasta_gz(self):
        test_inpath = "./tests/data/GCF_000981955.1_ASM98195v1_genomic.fna.gz"
        obs_seqs = [seq for seq in read_fasta(test_inpath)]
        exp_seqs = [
            ">contig_1_plasmid_at_the_start",
            ">contig_2",
            ">contig_3_phage_in_the_middle",
            ">contig_3",
            ">contig_4_plasmid_at_the_end",
        ]
        self.assertEqual(obs_seqs, exp_seqs)

    def test_read_fasta_empty(self):
        test_inpath = (
            "./tests/data/GCF_000981955.1_ASM98195v1_genomic_empty.fna"
        )
        with self.assertRaises(SystemExit) as context:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                with self.assertLogs(logger, level="ERROR") as log:
                    obs_seqs = [seq for seq in read_fasta(test_inpath)]
                    # verify skbio warning
                    self.assertTrue(
                        any(
                            issubclass(
                                warn.category, FormatIdentificationWarning
                            )
                            for warn in w
                        )
                    )
                    # verify wipe logger
                    self.assertIn(
                        f"Reading {test_inpath} failed due to the error",
                        log.output[0],
                    )
        # verify sys.exit(1)
        self.assertEqual(context.exception.code, 1)

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
                exp = (
                    f"Wrote 5 contigs (240 characters in total) into {test_outdir}/G000981955/G000981955.fna.\n"
                    "0 contigs were filtered.\n"
                )
                self.assertAlmostEqual(
                    obs, exp, "Log does not match expected log content."
                )
        finally:
            shutil.rmtree(test_outdir)

    def test_linearize_single_genome_ncbi_noconcat_filt(self):
        test_inpath = "./tests/data/GCF_000981955.1_ASM98195v1_genomic"
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
                    f"Wrote 2 contigs (90 characters in total) into {test_outdir}/G000981955/G000981955.fna.\n"
                    "3 contigs were filtered.\n"
                )
                self.assertAlmostEqual(
                    obs, exp, "Log does not match expected log content."
                )
        finally:
            shutil.rmtree(test_outdir)

    def test_linearize_single_genome_ncbi_concat_filt(self):
        test_inpath = "./tests/data/GCF_000981955.1_ASM98195v1_genomic"
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
                    f"Wrote 2 contigs (110 characters in total) into {test_outdir}/G000981955/G000981955.fna.\n"
                    "3 contigs were filtered.\n"
                )
                self.assertAlmostEqual(
                    obs, exp, "Log does not match expected log content."
                )
        finally:
            shutil.rmtree(test_outdir)

    def test_linearize_single_genome_ncbi_concat_nofilt(self):
        test_inpath = "./tests/data/GCF_000981955.1_ASM98195v1_genomic"
        test_ext = ".fna"
        test_outdir = tempfile.mkdtemp()
        test_gid = None
        test_gap = "N*20"
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
                    ">G000981955\n"
                    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"
                    "NNNNNNNNNNNNNNNNNNNN"
                    "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
                    "NNNNNNNNNNNNNNNNNNNN"
                    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                    "NNNNNNNNNNNNNNNNNNNN"
                    "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"
                    "NNNNNNNNNNNNNNNNNNNN"
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
                    f"Wrote 5 contigs (320 characters in total) into {test_outdir}/G000981955/G000981955.fna.\n"
                    "0 contigs were filtered.\n"
                )
                self.assertAlmostEqual(
                    obs, exp, "Log does not match expected log content."
                )
        finally:
            shutil.rmtree(test_outdir)

    def test_linearize_single_genome_noncbi_concat_filt(self):
        test_inpath = "./tests/data/GCF_000981955.1_ASM98195v1_genomic"
        test_ext = ".fna"
        test_outdir = tempfile.mkdtemp()
        test_gid = "H000000001"
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

            outpath = os.path.join(test_outdir, "H000000001", "H000000001.fna")
            self.assertTrue(os.path.exists(outpath), "Output file not found.")

            with open(outpath, "r") as f:
                obs = f.read()
                exp = (
                    ">H000000001\n"
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
                test_outdir, "H000000001", "linearization.log"
            )
            with open(logpath, "r") as f:
                obs = f.read()
                exp = (
                    f"Wrote 2 contigs (110 characters in total) into {test_outdir}/H000000001/H000000001.fna.\n"
                    "3 contigs were filtered.\n"
                )
                self.assertAlmostEqual(
                    obs, exp, "Log does not match expected log content."
                )
        finally:
            shutil.rmtree(test_outdir)

    def test_linearize_single_genome_noncbi_concat_filt_zip(self):
        test_inpath = "./tests/data/GCF_000981955.1_ASM98195v1_genomic"
        test_ext = ".fna.gz"
        test_outdir = tempfile.mkdtemp()
        test_gid = "H000000001"
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

            outpath = os.path.join(test_outdir, "H000000001", "H000000001.fna")
            self.assertTrue(os.path.exists(outpath), "Output file not found.")

            with open(outpath, "r") as f:
                obs = f.read()
                exp = (
                    ">H000000001\n"
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
                test_outdir, "H000000001", "linearization.log"
            )
            with open(logpath, "r") as f:
                obs = f.read()
                exp = (
                    f"Wrote 2 contigs (110 characters in total) into {test_outdir}/H000000001/H000000001.fna.\n"
                    "3 contigs were filtered.\n"
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

    def test_linearize_genomes_repeated_gid_nogid(self):
        test_metadata = "./tests/data/metadata_repeated_gid_nogid.tsv"
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
                os.path.join(test_outdir, "linearization_all.log"),
            ]
            for p in paths:
                self.assertTrue(os.path.exists(p), "Output file not found.")

            with open(
                os.path.join(test_outdir, "linearization_all.log"), "r"
            ) as f:
                obs = f.read()
                exp = "Dupliacted genome ID: G000981955 with input path ./tests/data/GCF_000981955.1_ASM98195v1_genomic\n"
                self.assertAlmostEqual(
                    obs, exp, "Log does not match expected log content."
                )
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
