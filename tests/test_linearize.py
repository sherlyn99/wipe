import os
import gzip
import shutil
import unittest
import tempfile
import warnings
from skbio.io._exception import FASTAFormatError
from skbio.io import FormatIdentificationWarning
from wipe.modules.linearize import (
    generate_gap_string,
    should_filter_contig_name,
    generate_inpath_outpath,
    generate_log_entries,
    read_fasta,
    generate_outdir,
    linearize_single_genome,
    linearize_genomes,
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
        exp_outdir = "/path/to/outdir/G/000/000/001"
        exp_outpath = "/path/to/outdir/G/000/000/001/G000000001.fna.gz"
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
        exp_outdir = "/path/to/outdir/G/000/000/001"
        exp_outpath = "/path/to/outdir/G/000/000/001/G000000001.fa.gz"
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
        exp_outdir = "/path/to/outdir/G/000/000/001"
        exp_outpath = "/path/to/outdir/G/000/000/001/G000000001.fa.gz"
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
        exp_outdir = "/path/to/outdir/G/000/981/955"
        exp_outpath = "/path/to/outdir/G/000/981/955/G000981955.fa.gz"
        self.assertEqual(obs_gid, exp_gid)
        self.assertEqual(obs_inpath, exp_inpath)
        self.assertEqual(obs_outdir, exp_outdir)
        self.assertEqual(obs_outpath, exp_outpath)

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
            "genome_path": "/path/to/outpath",
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

    def test_generate_outdir(self):
        test_outdir = "/parent/path"
        test_gid = "G000981955"
        obs = generate_outdir(test_outdir, test_gid)
        exp = "/parent/path/G/000/981/955"
        self.assertEqual(obs, exp)

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

            outpath = os.path.join(
                test_outdir, "G", "000", "981", "955", "G000981955.fna.gz"
            )
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
                "G",
                "000",
                "981",
                "955",
                "linearization_stats_G000981955.json.gz",
            )
            self.assertTrue(
                os.path.exists(logpath),
                f"Stats log file not found: {logpath}.",
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

            outpath = os.path.join(
                test_outdir, "G", "000", "981", "955", "G000981955.fna.gz"
            )
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
                "G",
                "000",
                "981",
                "955",
                "linearization_stats_G000981955.json.gz",
            )
            self.assertTrue(
                os.path.exists(logpath),
                f"Stats log file not found: {logpath}.",
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

            outpath = os.path.join(
                test_outdir, "G", "000", "981", "955", "G000981955.fna.gz"
            )
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
                "G",
                "000",
                "981",
                "955",
                "linearization_stats_G000981955.json.gz",
            )
            self.assertTrue(
                os.path.exists(logpath),
                f"Stats log file not found: {logpath}.",
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

            outpath = os.path.join(
                test_outdir, "G", "000", "981", "955", "G000981955.fna.gz"
            )
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
                "G",
                "000",
                "981",
                "955",
                "linearization_stats_G000981955.json.gz",
            )
            self.assertTrue(
                os.path.exists(logpath),
                f"Stats log file not found: {logpath}.",
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

            outpath = os.path.join(
                test_outdir, "H", "000", "000", "001", "H000000001.fna.gz"
            )
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
                "H",
                "000",
                "000",
                "001",
                "linearization_stats_H000000001.json.gz",
            )
            self.assertTrue(
                os.path.exists(logpath),
                f"Stats log file not found: {logpath}.",
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

            outpath = os.path.join(
                test_outdir, "H", "000", "000", "001", "H000000001.fna.gz"
            )
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
                "H",
                "000",
                "000",
                "001",
                "linearization_stats_H000000001.json.gz",
            )
            self.assertTrue(
                os.path.exists(logpath),
                f"Stats log file not found: {logpath}.",
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
                os.path.join(
                    test_outdir, "G", "000", "000", "001", "G000000001.fna.gz"
                ),
                os.path.join(
                    test_outdir,
                    "G",
                    "000",
                    "000",
                    "001",
                    "linearization_stats_G000000001.json.gz",
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
                    "linearization_stats_G000000001.json.gz",
                ),
            ]

            for p in paths:
                self.assertTrue(
                    os.path.exists(p), f"Output file not found: {p}."
                )

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
                os.path.join(
                    test_outdir, "G", "000", "981", "955", "G000981955.fna.gz"
                ),
                os.path.join(
                    test_outdir,
                    "G",
                    "000",
                    "981",
                    "955",
                    "linearization_stats_G000981955.json.gz",
                ),
            ]

            for p in paths:
                self.assertTrue(
                    os.path.exists(p), f"Output file not found: {p}."
                )

        finally:
            shutil.rmtree(test_outdir)

    def test_linearize_genomes_nogid_error(self):
        test_metadata = "./tests/data/metadata_nogid_error.tsv"
        test_ext = "fna"
        test_outdir = tempfile.mkdtemp()
        test_gap = "N*20"
        test_filt = None
        try:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                linearize_genomes(
                    test_metadata, test_ext, test_outdir, test_gap, test_filt
                )

                paths = [
                    os.path.join(
                        test_outdir,
                        "G",
                        "000",
                        "981",
                        "955",
                        "G000981955.fna.gz",
                    ),
                    os.path.join(
                        test_outdir,
                        "G",
                        "000",
                        "981",
                        "955",
                        "linearization_stats_G000981955.json.gz",
                    ),
                    os.path.join(test_outdir, "linearization_summary.json.gz"),
                ]

                for p in paths:
                    self.assertTrue(
                        os.path.exists(p), f"Output file not found: {p}."
                    )
                self.assertTrue(
                    any(
                        issubclass(warn.category, FormatIdentificationWarning)
                        for warn in w
                    )
                )

        finally:
            shutil.rmtree(test_outdir)

    def test_linearize_genomes_repeated_gid_nogid(self):
        test_metadata = "./tests/data/metadata_repeated_gid_nogid.tsv"
        test_ext = "fna"
        test_outdir = tempfile.mkdtemp()
        test_gap = "N*20"
        test_filt = None

        try:
            with self.assertRaises(SystemExit) as context:
                linearize_genomes(
                    test_metadata, test_ext, test_outdir, test_gap, test_filt
                )

                paths = [
                    os.path.join(test_outdir, "linearization_summary.json.gz"),
                ]
                for p in paths:
                    self.assertTrue(
                        os.path.exists(p), f"Output file not found: {p}."
                    )
            self.assertEqual(context.exception.code, 1)
        finally:
            shutil.rmtree(test_outdir)

    def test_linearize_genomes_repeated_gid(self):
        test_metadata = "./tests/data/metadata_repeated_gid.tsv"
        test_ext = "fna"
        test_outdir = tempfile.mkdtemp()
        test_gap = "N*20"
        test_filt = None

        try:
            with self.assertRaises(SystemExit) as context:
                linearize_genomes(
                    test_metadata, test_ext, test_outdir, test_gap, test_filt
                )

                paths = [
                    os.path.join(test_outdir, "linearization_summary.json.gz"),
                ]
                for p in paths:
                    self.assertTrue(
                        os.path.exists(p), f"Output file not found: {p}."
                    )
            self.assertEqual(context.exception.code, 1)
        finally:
            shutil.rmtree(test_outdir)


if __name__ == "__main__":
    unittest.main()
