import os
import gzip
import shutil
import unittest
import tempfile
import warnings
from skbio.io._exception import FASTAFormatError
from skbio.io import FormatIdentificationWarning
from wipe.modules.linearization import (
    generate_gap_string,
    should_filter_contig_name,
    generate_log_entries,
    read_fasta,
    linearization_single,
    linearization_batch,
)

warnings.filterwarnings("ignore", category=DeprecationWarning)


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

    def test_generate_log_entries(self):
        test_gid = "G001"
        test_n_written = 100
        test_n_char = 10000
        test_n_filtered = 50
        test_outpath = "/path/to/outpath"
        obs = generate_log_entries(
            test_gid,
            test_n_written,
            test_n_char,
            test_n_filtered,
            test_outpath,
        )
        exp = {
            "genome_id": "G001",
            "lgenome_path": "/path/to/outpath",
            "contigs_written": 100,
            "chars_written": 10000,
            "contigs_filtered": 50,
        }
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
            "./tests/data/GCF_000981956.1_ASM98195v1_genomic_empty.fna"
        )
        with self.assertRaises(FASTAFormatError) as context:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                obs_seqs = [seq for seq in read_fasta(test_inpath)]
                self.assertTrue(
                    any(
                        issubclass(warn.category, FormatIdentificationWarning)
                        for warn in w
                    )
                )

    def test_linearize_single_genome_ncbi_noconcat_nofilt(self):
        test_inpath = "./tests/data/GCF_000981955.1_ASM98195v1_genomic.fna.gz"
        test_gid = "G000981955"
        test_outdir = tempfile.mkdtemp()
        test_gap = None
        test_filt = None

        try:
            linearization_single(
                test_inpath,
                test_gid,
                test_outdir,
                test_gap,
                test_filt,
            )

            outpath = os.path.join(test_outdir, "G000981955.fna.gz")
            self.assertTrue(
                os.path.exists(outpath), f"Output file not found: {outpath}."
            )

            with gzip.open(outpath, "rt") as f:
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
                test_outdir,
                "linearization.tsv",
            )
            self.assertTrue(
                os.path.exists(logpath),
                f"Linearization log file not found: {logpath}.",
            )
        finally:
            shutil.rmtree(test_outdir)

    def test_linearize_single_genome_ncbi_noconcat_filt(self):
        test_inpath = "./tests/data/GCF_000981955.1_ASM98195v1_genomic.fna.gz"
        test_gid = "G000981955"
        test_outdir = tempfile.mkdtemp()
        test_gap = None
        test_filt = "plasmid,phage"

        try:
            linearization_single(
                test_inpath,
                test_gid,
                test_outdir,
                test_gap,
                test_filt,
            )

            outpath = os.path.join(test_outdir, "G000981955.fna.gz")
            self.assertTrue(
                os.path.exists(outpath), f"Output file not found: {outpath}."
            )

            with gzip.open(outpath, "rt") as f:
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
                test_outdir,
                "linearization.tsv",
            )
            self.assertTrue(
                os.path.exists(logpath),
                f"Linearization log file not found: {logpath}.",
            )
        finally:
            shutil.rmtree(test_outdir)

    def test_linearize_single_genome_ncbi_concat_filt(self):
        test_inpath = "./tests/data/GCF_000981955.1_ASM98195v1_genomic.fna.gz"
        test_gid = "G000981955"
        test_outdir = tempfile.mkdtemp()
        test_gap = "N*20"
        test_filt = "plasmid,phage"

        try:
            linearization_single(
                test_inpath,
                test_gid,
                test_outdir,
                test_gap,
                test_filt,
            )

            outpath = os.path.join(test_outdir, "G000981955.fna.gz")
            self.assertTrue(
                os.path.exists(outpath), f"Output file not found: {outpath}."
            )

            with gzip.open(outpath, "rt") as f:
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
                test_outdir,
                "linearization.tsv",
            )
            self.assertTrue(
                os.path.exists(logpath),
                f"Linearization log log file not found: {logpath}.",
            )
        finally:
            shutil.rmtree(test_outdir)

    def test_linearize_single_genome_ncbi_concat_nofilt(self):
        test_inpath = "./tests/data/GCF_000981955.1_ASM98195v1_genomic.fna.gz"
        test_gid = "G000981955"
        test_outdir = tempfile.mkdtemp()
        test_gap = "N*20"
        test_filt = None

        try:
            linearization_single(
                test_inpath,
                test_gid,
                test_outdir,
                test_gap,
                test_filt,
            )

            outpath = os.path.join(test_outdir, "G000981955.fna.gz")
            self.assertTrue(os.path.exists(outpath), "Output file not found.")

            with gzip.open(outpath, "rt") as f:
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
                test_outdir,
                "linearization.tsv",
            )
            self.assertTrue(
                os.path.exists(logpath),
                f"Linearization log file not found: {logpath}.",
            )
        finally:
            shutil.rmtree(test_outdir)

    def test_linearize_single_genome_nonncbi_concat_filt(self):
        test_inpath = "./tests/data/GCF_000981955.1_ASM98195v1_genomic.fna.gz"
        test_gid = "H000000001"
        test_outdir = tempfile.mkdtemp()
        test_gap = "N*20"
        test_filt = "plasmid,phage"

        try:
            linearization_single(
                test_inpath,
                test_gid,
                test_outdir,
                test_gap,
                test_filt,
            )

            outpath = os.path.join(test_outdir, "H000000001.fna.gz")
            self.assertTrue(
                os.path.exists(outpath), f"Output file not found: {outpath}."
            )

            with gzip.open(outpath, "rt") as f:
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
                test_outdir,
                "linearization.tsv",
            )
            self.assertTrue(
                os.path.exists(logpath),
                f"Linearization log file not found: {logpath}.",
            )
        finally:
            shutil.rmtree(test_outdir)

    def test_linearize_single_genome_nonncbi_concat_filt_nonzip(self):
        test_inpath = "./tests/data/GCF_000981955.1_ASM98195v1_genomic.fna"
        test_gid = "H000000001"
        test_outdir = tempfile.mkdtemp()
        test_gap = "N*20"
        test_filt = "plasmid,phage"

        try:
            linearization_single(
                test_inpath,
                test_gid,
                test_outdir,
                test_gap,
                test_filt,
            )

            outpath = os.path.join(test_outdir, "H000000001.fna.gz")
            self.assertTrue(
                os.path.exists(outpath), f"Output file not found: {outpath}."
            )

            with gzip.open(outpath, "rt") as f:
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
                test_outdir,
                "linearization.tsv",
            )
            self.assertTrue(
                os.path.exists(logpath),
                f"Linearization log file not found: {logpath}.",
            )
        finally:
            shutil.rmtree(test_outdir)

    def test_linearize_genomes(self):
        test_metadata = "./tests/data/metadata.tsv"
        test_outdir = tempfile.mkdtemp()
        test_gap = "N*20"
        test_filt = None
        try:
            linearization_batch(
                test_metadata, test_outdir, test_gap, test_filt
            )

            paths = [
                os.path.join(
                    test_outdir, "G", "000", "000", "001", "G000000001.fna.gz"
                ),
                os.path.join(
                    test_outdir,
                    "G",
                    "000",
                    "000",
                    "001",
                    "linearization.tsv",
                ),
                os.path.join(
                    test_outdir, "G", "000", "000", "002", "G000000002.fna.gz"
                ),
                os.path.join(
                    test_outdir,
                    "G",
                    "000",
                    "000",
                    "001",
                    "linearization.tsv",
                ),
            ]

            for p in paths:
                self.assertTrue(
                    os.path.exists(p), f"Output file not found: {p}."
                )

        finally:
            shutil.rmtree(test_outdir)

    def test_linearize_genomes_repeated_gid(self):
        test_metadata = "./tests/data/metadata_repeated_gid.tsv"
        test_outdir = tempfile.mkdtemp()
        test_gap = "N*20"
        test_filt = None

        try:
            with self.assertRaises(ValueError) as context:
                linearization_batch(
                    test_metadata, test_outdir, test_gap, test_filt
                )

                paths = [
                    os.path.join(test_outdir, "linearization_summary.txt"),
                ]
                for p in paths:
                    self.assertTrue(
                        os.path.exists(p), f"Output file not found: {p}."
                    )
            self.assertIn(
                "Duplicated genome ids found:", str(context.exception)
            )
        finally:
            shutil.rmtree(test_outdir)

    def test_linearize_genomes_empty(self):
        test_metadata = "./tests/data/metadata_empty.tsv"
        test_outdir = tempfile.mkdtemp()
        test_gap = "N*20"
        test_filt = None

        try:
            linearization_batch(
                test_metadata, test_outdir, test_gap, test_filt
            )

            paths = [
                os.path.join(test_outdir, "linearization_summary.txt"),
            ]
            for p in paths:
                self.assertTrue(
                    os.path.exists(p), f"Output file not found: {p}."
                )
        finally:
            shutil.rmtree(test_outdir)


if __name__ == "__main__":
    unittest.main()
