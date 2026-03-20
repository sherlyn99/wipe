import os
import shutil
import tempfile
import unittest
import subprocess
from pathlib import Path
from unittest.mock import patch
from click.testing import CliRunner
from wipe.wipe import functional_db as functional_db_cli
from wipe.modules.functional_db import make_eggnog_mappings


class TestFunctionalDbDownload(unittest.TestCase):
    def setUp(self):
        self.runner = CliRunner()
        self.test_outdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.test_outdir)

    @patch("subprocess.run")
    def test_download_uniref_only(self, mock_run):
        """--uniref downloads uniref90/50 via wget and builds diamond databases."""
        mock_run.return_value = subprocess.CompletedProcess([], 0, "", "")
        result = self.runner.invoke(
            functional_db_cli,
            ["download", "--uniref", "--no-eggnog", "-o", self.test_outdir],
        )
        assert result.exit_code == 0, result.output
        calls = mock_run.call_args_list
        assert any("uniref90.fasta.gz" in str(c) for c in calls)
        assert any("uniref50.fasta.gz" in str(c) for c in calls)
        assert any("uniref90.dmnd" in str(c) for c in calls)
        assert any("uniref50.dmnd" in str(c) for c in calls)
        diamond_calls = [c for c in calls if "diamond" in str(c) and "makedb" in str(c)]
        assert len(diamond_calls) == 2
        assert all("download_eggnog_data.py" not in str(c) for c in calls)

    @patch("subprocess.run")
    def test_download_eggnog_only(self, mock_run):
        """--eggnog downloads all five EggNOG files via wget and extracts each."""
        mock_run.return_value = subprocess.CompletedProcess([], 0, "", "")
        result = self.runner.invoke(
            functional_db_cli,
            ["download", "--no-uniref", "--eggnog", "-o", self.test_outdir],
        )
        assert result.exit_code == 0, result.output
        calls = mock_run.call_args_list
        expected_files = [
            "eggnog.db.gz",
            "eggnog.taxa.tar.gz",
            "eggnog_proteins.dmnd.gz",
            "mmseqs.tar.gz",
            "pfam.tar.gz",
        ]
        for filename in expected_files:
            assert any(filename in str(c) for c in calls), f"{filename} not downloaded"
        assert all("uniref" not in str(c) for c in calls)

    @patch("subprocess.run")
    def test_download_both_by_default(self, mock_run):
        """Both databases are downloaded when no flags are given (default True)."""
        mock_run.return_value = subprocess.CompletedProcess([], 0, "", "")
        result = self.runner.invoke(
            functional_db_cli,
            ["download", "-o", self.test_outdir],
        )
        assert result.exit_code == 0, result.output
        calls = mock_run.call_args_list
        assert any("uniref90.fasta.gz" in str(c) for c in calls)
        assert any("uniref50.fasta.gz" in str(c) for c in calls)
        assert any("eggnog.db.gz" in str(c) for c in calls)

    @patch("subprocess.run")
    def test_download_eggnog_skips_existing_file(self, mock_run):
        """A file that already exists and is non-empty is not re-downloaded."""
        mock_run.return_value = subprocess.CompletedProcess([], 0, "", "")
        eggnog_dir = os.path.join(self.test_outdir, "eggnog")
        os.makedirs(eggnog_dir, exist_ok=True)
        existing = os.path.join(eggnog_dir, "eggnog.db.gz")
        with open(existing, "w") as f:
            f.write("placeholder")

        result = self.runner.invoke(
            functional_db_cli,
            ["download", "--no-uniref", "--eggnog", "-o", self.test_outdir],
        )
        assert result.exit_code == 0, result.output
        wget_calls = [c for c in mock_run.call_args_list if "wget" in str(c)]
        assert all("eggnog.db.gz" not in str(c) for c in wget_calls)

    @patch("subprocess.run")
    def test_download_eggnog_gz_uses_gunzip(self, mock_run):
        """Plain .gz files (non-tar) are extracted with gunzip."""
        mock_run.return_value = subprocess.CompletedProcess([], 0, "", "")
        result = self.runner.invoke(
            functional_db_cli,
            ["download", "--no-uniref", "--eggnog", "-o", self.test_outdir],
        )
        assert result.exit_code == 0, result.output
        gunzip_calls = [c for c in mock_run.call_args_list if "gunzip" in str(c)]
        gunzip_targets = " ".join(str(c) for c in gunzip_calls)
        assert "eggnog.db.gz" in gunzip_targets
        assert "eggnog_proteins.dmnd.gz" in gunzip_targets

    @patch("subprocess.run")
    def test_download_eggnog_tar_uses_tar(self, mock_run):
        """.tar.gz files are extracted with tar -xzf."""
        mock_run.return_value = subprocess.CompletedProcess([], 0, "", "")
        result = self.runner.invoke(
            functional_db_cli,
            ["download", "--no-uniref", "--eggnog", "-o", self.test_outdir],
        )
        assert result.exit_code == 0, result.output
        tar_calls = [c for c in mock_run.call_args_list if "tar" in str(c) and "-xzf" in str(c)]
        tar_targets = " ".join(str(c) for c in tar_calls)
        assert "eggnog.taxa.tar.gz" in tar_targets
        assert "mmseqs.tar.gz" in tar_targets
        assert "pfam.tar.gz" in tar_targets

    def test_download_neither_fails(self):
        """--no-uniref --no-eggnog should exit with an error."""
        result = self.runner.invoke(
            functional_db_cli,
            ["download", "--no-uniref", "--no-eggnog", "-o", self.test_outdir],
        )
        assert result.exit_code != 0

    @patch("subprocess.run")
    def test_download_uniref_creates_subdirectory(self, mock_run):
        """download_uniref creates dbs/uniref/ subdirectory."""
        mock_run.return_value = subprocess.CompletedProcess([], 0, "", "")
        result = self.runner.invoke(
            functional_db_cli,
            ["download", "--uniref", "--no-eggnog", "-o", self.test_outdir],
        )
        assert result.exit_code == 0, result.output
        assert os.path.isdir(os.path.join(self.test_outdir, "uniref"))

    @patch("subprocess.run")
    def test_download_eggnog_creates_subdirectory(self, mock_run):
        """download_eggnog creates dbs/eggnog/ subdirectory."""
        mock_run.return_value = subprocess.CompletedProcess([], 0, "", "")
        result = self.runner.invoke(
            functional_db_cli,
            ["download", "--no-uniref", "--eggnog", "-o", self.test_outdir],
        )
        assert result.exit_code == 0, result.output
        assert os.path.isdir(os.path.join(self.test_outdir, "eggnog"))


class TestFunctionalDbBuild(unittest.TestCase):
    def setUp(self):
        self.runner = CliRunner()
        self.test_outdir = tempfile.mkdtemp()
        self.faa = os.path.join(self.test_outdir, "test.faa")
        with open(self.faa, "w") as f:
            f.write(">G000000001_1\nMVKVKYFGCELNGYTDVVEVKDGK\n")

        self.uniref_db_dir = tempfile.mkdtemp()
        Path(os.path.join(self.uniref_db_dir, "uniref90.dmnd")).touch()
        Path(os.path.join(self.uniref_db_dir, "uniref50.dmnd")).touch()

        self.eggnog_db_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.test_outdir)
        shutil.rmtree(self.uniref_db_dir)
        shutil.rmtree(self.eggnog_db_dir)

    def test_build_no_flags_fails(self):
        """Invoking build without --uniref or --eggnog should fail."""
        result = self.runner.invoke(
            functional_db_cli,
            ["build", "-i", self.faa, "-o", self.test_outdir],
        )
        assert result.exit_code != 0
        assert "At least one of --uniref or --eggnog" in result.output

    @patch("wipe.modules.functional_db.merge_uniref")
    @patch("subprocess.run")
    def test_build_uniref_calls_diamond_twice(self, mock_run, mock_merge):
        """build --uniref runs diamond blastp against uniref90 and uniref50."""
        mock_run.return_value = subprocess.CompletedProcess([], 0, "", "")
        merged = os.path.join(self.test_outdir, "uniref_map.txt")
        mock_merge.side_effect = lambda *a, **kw: Path(merged).touch()

        result = self.runner.invoke(
            functional_db_cli,
            [
                "build", "-i", self.faa, "-o", self.test_outdir,
                "--uniref", "--uniref-db", self.uniref_db_dir, "-t", "8",
            ],
        )
        assert result.exit_code == 0, result.output

        diamond_calls = [c for c in mock_run.call_args_list if "diamond" in str(c)]
        assert len(diamond_calls) == 2
        assert any("uniref90.dmnd" in str(c) for c in diamond_calls)
        assert any("uniref50.dmnd" in str(c) for c in diamond_calls)
        for c in diamond_calls:
            assert "--threads" in str(c) and "8" in str(c)
            assert "--id" in str(c)

    @patch("wipe.modules.functional_db.merge_uniref")
    @patch("subprocess.run")
    def test_build_uniref_merges_with_simplify(self, mock_run, mock_merge):
        """build --uniref calls merge_uniref with simplify=True."""
        mock_run.return_value = subprocess.CompletedProcess([], 0, "", "")
        merged = os.path.join(self.test_outdir, "uniref_map.txt")
        mock_merge.side_effect = lambda *a, **kw: Path(merged).touch()

        self.runner.invoke(
            functional_db_cli,
            [
                "build", "-i", self.faa, "-o", self.test_outdir,
                "--uniref", "--uniref-db", self.uniref_db_dir,
            ],
        )
        mock_merge.assert_called_once()
        assert mock_merge.call_args[1].get("simplify") is True

    @patch("wipe.modules.functional_db.merge_uniref")
    @patch("subprocess.run")
    def test_build_uniref_compresses_output(self, mock_run, mock_merge):
        """build --uniref calls xz on the merged map file."""
        mock_run.return_value = subprocess.CompletedProcess([], 0, "", "")
        merged = os.path.join(self.test_outdir, "uniref_map.txt")
        mock_merge.side_effect = lambda *a, **kw: Path(merged).touch()

        self.runner.invoke(
            functional_db_cli,
            [
                "build", "-i", self.faa, "-o", self.test_outdir,
                "--uniref", "--uniref-db", self.uniref_db_dir,
            ],
        )
        xz_calls = [c for c in mock_run.call_args_list if "xz" in str(c)]
        assert len(xz_calls) == 1
        assert merged in str(xz_calls[0])

    @patch("subprocess.run")
    def test_build_eggnog_calls_emapper(self, mock_run):
        """build --eggnog runs emapper.py with correct arguments."""
        mock_run.return_value = subprocess.CompletedProcess([], 0, "", "")

        result = self.runner.invoke(
            functional_db_cli,
            [
                "build", "-i", self.faa, "-o", self.test_outdir,
                "--eggnog", "--eggnog-db", self.eggnog_db_dir, "-t", "8",
            ],
        )
        assert result.exit_code == 0, result.output

        emapper_calls = [c for c in mock_run.call_args_list if "emapper.py" in str(c)]
        assert len(emapper_calls) == 1
        call_str = str(emapper_calls[0])
        assert self.faa in call_str
        assert self.eggnog_db_dir in call_str
        assert "--cpu" in call_str and "8" in call_str

    @patch("subprocess.run")
    def test_build_eggnog_renames_annotations(self, mock_run):
        """eggnog.emapper.annotations is renamed to eggnog_map.tsv and mapping files are created."""
        mock_run.return_value = subprocess.CompletedProcess([], 0, "", "")
        annotations = os.path.join(self.test_outdir, "eggnog.emapper.annotations")
        header = "#query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs\n"
        row = "G000000001_1\tX\t1e-10\t100\tX\tX\tJ\tdesc\tname\tGO:0001\t1.1.1.1\tko:K00001\t-\t-\t-\t-\t-\t-\tGH1\t-\tPF00001\n"
        with open(annotations, "w") as f:
            f.write(header)
            f.write(row)

        result = self.runner.invoke(
            functional_db_cli,
            [
                "build", "-i", self.faa, "-o", self.test_outdir,
                "--eggnog", "--eggnog-db", self.eggnog_db_dir,
            ],
        )
        assert result.exit_code == 0, result.output
        assert os.path.exists(os.path.join(self.test_outdir, "eggnog_map.tsv"))
        assert not os.path.exists(annotations)
        for fname in ("orf_to_go.tsv", "orf_to_ec.tsv", "orf_to_ko.tsv",
                      "orf_to_cazy.tsv", "orf_to_cog.tsv", "orf_to_pfam.tsv"):
            assert os.path.exists(os.path.join(self.test_outdir, fname)), f"{fname} missing"

    @patch("wipe.modules.functional_db.merge_uniref")
    @patch("subprocess.run")
    def test_build_custom_uniref_db_path(self, mock_run, mock_merge):
        """Custom --uniref-db path is passed to diamond."""
        mock_run.return_value = subprocess.CompletedProcess([], 0, "", "")
        merged = os.path.join(self.test_outdir, "uniref_map.txt")
        mock_merge.side_effect = lambda *a, **kw: Path(merged).touch()

        custom_db = tempfile.mkdtemp()
        try:
            Path(os.path.join(custom_db, "uniref90.dmnd")).touch()
            Path(os.path.join(custom_db, "uniref50.dmnd")).touch()

            result = self.runner.invoke(
                functional_db_cli,
                [
                    "build", "-i", self.faa, "-o", self.test_outdir,
                    "--uniref", "--uniref-db", custom_db,
                ],
            )
            assert result.exit_code == 0, result.output
            diamond_calls = [c for c in mock_run.call_args_list if "diamond" in str(c)]
            assert all(custom_db in str(c) for c in diamond_calls)
        finally:
            shutil.rmtree(custom_db)

    @patch("wipe.modules.functional_db.merge_uniref")
    @patch("subprocess.run")
    def test_build_both_flags(self, mock_run, mock_merge):
        """build with --uniref and --eggnog runs both annotation pipelines."""
        mock_run.return_value = subprocess.CompletedProcess([], 0, "", "")
        merged = os.path.join(self.test_outdir, "uniref_map.txt")
        mock_merge.side_effect = lambda *a, **kw: Path(merged).touch()

        result = self.runner.invoke(
            functional_db_cli,
            [
                "build", "-i", self.faa, "-o", self.test_outdir,
                "--uniref", "--uniref-db", self.uniref_db_dir,
                "--eggnog", "--eggnog-db", self.eggnog_db_dir,
            ],
        )
        assert result.exit_code == 0, result.output
        diamond_calls = [c for c in mock_run.call_args_list if "diamond" in str(c)]
        emapper_calls = [c for c in mock_run.call_args_list if "emapper.py" in str(c)]
        assert len(diamond_calls) == 2
        assert len(emapper_calls) == 1


class TestMakeEggnogMappings(unittest.TestCase):
    """Unit tests for make_eggnog_mappings using a synthetic annotations file."""

    # Minimal annotations file with one row covering all six annotation types.
    # Columns (0-based): 0=query, 6=COG_category, 9=GOs, 10=EC, 11=KEGG_ko,
    #                    18=CAZy, 20=PFAMs
    HEADER = "#query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs\n"
    ROW = "ORF1\tX\t1e-10\t100\tX\tX\tJK\tdesc\tname\tGO:0001,GO:0002\t1.2.3.4\tko:K00001,ko:K00002\t-\t-\t-\t-\t-\t-\tGH13,GT2\t-\tPF00001,PF00002\n"
    ROW_EMPTY = "ORF2\tX\t1e-5\t50\tX\tX\t-\tdesc\tname\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\n"

    def setUp(self):
        self.outdir = tempfile.mkdtemp()
        self.annotations = os.path.join(self.outdir, "test.emapper.annotations")
        with open(self.annotations, "w") as f:
            f.write(self.HEADER)
            f.write(self.ROW)
            f.write(self.ROW_EMPTY)

    def tearDown(self):
        shutil.rmtree(self.outdir)

    def _read(self, fname):
        path = os.path.join(self.outdir, fname)
        with open(path) as f:
            return [line.rstrip("\n").split("\t") for line in f]

    def test_go_mapping(self):
        make_eggnog_mappings(self.annotations, self.outdir)
        rows = self._read("orf_to_go.tsv")
        assert ["ORF1", "GO:0001"] in rows
        assert ["ORF1", "GO:0002"] in rows
        assert all(r[0] != "ORF2" for r in rows)

    def test_ec_mapping(self):
        make_eggnog_mappings(self.annotations, self.outdir)
        rows = self._read("orf_to_ec.tsv")
        assert ["ORF1", "1.2.3.4"] in rows

    def test_ko_mapping(self):
        make_eggnog_mappings(self.annotations, self.outdir)
        rows = self._read("orf_to_ko.tsv")
        assert ["ORF1", "ko:K00001"] in rows
        assert ["ORF1", "ko:K00002"] in rows

    def test_cazy_mapping(self):
        make_eggnog_mappings(self.annotations, self.outdir)
        rows = self._read("orf_to_cazy.tsv")
        assert ["ORF1", "GH13"] in rows
        assert ["ORF1", "GT2"] in rows

    def test_cog_mapping_splits_by_character(self):
        """COG_category 'JK' should produce two rows, one per letter."""
        make_eggnog_mappings(self.annotations, self.outdir)
        rows = self._read("orf_to_cog.tsv")
        assert ["ORF1", "J"] in rows
        assert ["ORF1", "K"] in rows

    def test_pfam_mapping(self):
        make_eggnog_mappings(self.annotations, self.outdir)
        rows = self._read("orf_to_pfam.tsv")
        assert ["ORF1", "PF00001"] in rows
        assert ["ORF1", "PF00002"] in rows

    def test_dash_values_skipped(self):
        """ORF2 has '-' for all annotations and should appear in no mapping file."""
        make_eggnog_mappings(self.annotations, self.outdir)
        for fname in ("orf_to_go.tsv", "orf_to_ec.tsv", "orf_to_ko.tsv",
                      "orf_to_cazy.tsv", "orf_to_cog.tsv", "orf_to_pfam.tsv"):
            rows = self._read(fname)
            assert all(r[0] != "ORF2" for r in rows), f"ORF2 should not appear in {fname}"

    def test_all_output_files_created(self):
        make_eggnog_mappings(self.annotations, self.outdir)
        for fname in ("orf_to_go.tsv", "orf_to_ec.tsv", "orf_to_ko.tsv",
                      "orf_to_cazy.tsv", "orf_to_cog.tsv", "orf_to_pfam.tsv"):
            assert os.path.exists(os.path.join(self.outdir, fname)), f"{fname} missing"


class TestFunctionalDbBuildOutputFiles(unittest.TestCase):
    """Integration-style tests using the real example .faa from tests/data/.

    External tools (diamond, emapper, xz) are mocked with side_effects that
    create realistic output files so the full Python pipeline runs end-to-end.
    """

    FAA = "tests/data/M000000999_subset.faa"


    def setUp(self):
        self.runner = CliRunner()
        self.test_outdir = tempfile.mkdtemp()
        self.uniref_db_dir = tempfile.mkdtemp()
        Path(os.path.join(self.uniref_db_dir, "uniref90.dmnd")).touch()
        Path(os.path.join(self.uniref_db_dir, "uniref50.dmnd")).touch()
        self.eggnog_db_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.test_outdir)
        shutil.rmtree(self.uniref_db_dir)
        shutil.rmtree(self.eggnog_db_dir)

    _ANNOT_HEADER = (
        "#query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\t"
        "COG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\t"
        "KEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\t"
        "KEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs\n"
    )

    def _fake_run(self, cmd, **kwargs):
        """Simulate external tools by creating their expected output files."""
        if cmd[0] == "diamond":
            out_idx = cmd.index("--out")
            m8_path = cmd[out_idx + 1]
            with open(m8_path, "w") as f:
                f.write("M000000999_1\tUniRef90_Q7VK91\n")
                f.write("M000000999_2\tUniRef90_Q7VK90\n")
        elif cmd[0] == "xz":
            src = cmd[-1]
            if os.path.exists(src):
                os.rename(src, src + ".xz")
        elif cmd[0] == "emapper.py":
            out_idx = cmd.index("--output_dir")
            outdir = cmd[out_idx + 1]
            with open(os.path.join(outdir, "eggnog.emapper.annotations"), "w") as f:
                f.write(self._ANNOT_HEADER)
                f.write("M000000999_1\tX\t1e-10\t100\tX\tX\tJ\tdesc\tname\tGO:0001\t1.1.1.1\tko:K00001\t-\t-\t-\t-\t-\t-\tGH1\t-\tPF00001\n")
                f.write("M000000999_2\tX\t1e-5\t50\tX\tX\tK\tdesc\tname\tGO:0002\t2.2.2.2\tko:K00002\t-\t-\t-\t-\t-\t-\tGT2\t-\tPF00002\n")
        return subprocess.CompletedProcess(cmd, 0, "", "")

    @patch("subprocess.run")
    def test_build_uniref_creates_uniref_map(self, mock_run):
        """uniref_map.txt.xz is created in outdir when --uniref is used."""
        mock_run.side_effect = self._fake_run
        result = self.runner.invoke(
            functional_db_cli,
            [
                "build", "-i", self.FAA, "-o", self.test_outdir,
                "--uniref", "--uniref-db", self.uniref_db_dir,
            ],
        )
        assert result.exit_code == 0, result.output
        assert os.path.exists(os.path.join(self.test_outdir, "uniref_map.txt.xz"))

    @patch("subprocess.run")
    def test_build_uniref_no_extra_files(self, mock_run):
        """Running --uniref does not create eggnog_map.tsv."""
        mock_run.side_effect = self._fake_run
        self.runner.invoke(
            functional_db_cli,
            [
                "build", "-i", self.FAA, "-o", self.test_outdir,
                "--uniref", "--uniref-db", self.uniref_db_dir,
            ],
        )
        assert not os.path.exists(os.path.join(self.test_outdir, "eggnog_map.tsv"))

    @patch("subprocess.run")
    def test_build_eggnog_creates_eggnog_map(self, mock_run):
        """eggnog_map.tsv and all six mapping files are created when --eggnog is used."""
        mock_run.side_effect = self._fake_run
        result = self.runner.invoke(
            functional_db_cli,
            [
                "build", "-i", self.FAA, "-o", self.test_outdir,
                "--eggnog", "--eggnog-db", self.eggnog_db_dir,
            ],
        )
        assert result.exit_code == 0, result.output
        assert os.path.exists(os.path.join(self.test_outdir, "eggnog_map.tsv"))
        for fname in ("orf_to_go.tsv", "orf_to_ec.tsv", "orf_to_ko.tsv",
                      "orf_to_cazy.tsv", "orf_to_cog.tsv", "orf_to_pfam.tsv"):
            assert os.path.exists(os.path.join(self.test_outdir, fname)), f"{fname} missing"

    @patch("subprocess.run")
    def test_build_eggnog_no_extra_files(self, mock_run):
        """Running --eggnog does not create uniref_map.txt.xz."""
        mock_run.side_effect = self._fake_run
        self.runner.invoke(
            functional_db_cli,
            [
                "build", "-i", self.FAA, "-o", self.test_outdir,
                "--eggnog", "--eggnog-db", self.eggnog_db_dir,
            ],
        )
        assert not os.path.exists(os.path.join(self.test_outdir, "uniref_map.txt.xz"))

    @patch("subprocess.run")
    def test_build_both_creates_both_output_files(self, mock_run):
        """Both output files are created when --uniref and --eggnog are used."""
        mock_run.side_effect = self._fake_run
        result = self.runner.invoke(
            functional_db_cli,
            [
                "build", "-i", self.FAA, "-o", self.test_outdir,
                "--uniref", "--uniref-db", self.uniref_db_dir,
                "--eggnog", "--eggnog-db", self.eggnog_db_dir,
            ],
        )
        assert result.exit_code == 0, result.output
        assert os.path.exists(os.path.join(self.test_outdir, "uniref_map.txt.xz"))
        assert os.path.exists(os.path.join(self.test_outdir, "eggnog_map.tsv"))


_UNIREF_DB = "dbs/uniref"
_EGGNOG_DB = "dbs/eggnog"
_UNIREF_READY = (
    os.path.exists(os.path.join(_UNIREF_DB, "uniref90.dmnd"))
    and os.path.exists(os.path.join(_UNIREF_DB, "uniref50.dmnd"))
)
_EGGNOG_READY = os.path.exists(os.path.join(_EGGNOG_DB, "eggnog.db"))


class TestFunctionalDbBuildReal(unittest.TestCase):
    """End-to-end tests that run real tools against real databases.

    Skipped automatically when the required databases are not present.
    Run 'wipe functional-db download' followed by 'wipe uniref build' first.
    """

    FAA = "tests/data/M000000999_subset.faa"

    def setUp(self):
        self.runner = CliRunner()
        self.test_outdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.test_outdir)

    @unittest.skipUnless(_UNIREF_READY, "UniRef diamond databases not found in dbs/uniref/")
    def test_build_uniref_real(self):
        """build --uniref produces a non-empty uniref_map.txt.xz."""
        result = self.runner.invoke(
            functional_db_cli,
            [
                "build", "-i", self.FAA, "-o", self.test_outdir,
                "--uniref", "--uniref-db", _UNIREF_DB,
            ],
        )
        assert result.exit_code == 0, result.output
        out = os.path.join(self.test_outdir, "uniref_map.txt.xz")
        assert os.path.exists(out), "uniref_map.txt.xz was not created"
        assert os.path.getsize(out) > 0, "uniref_map.txt.xz is empty"

    @unittest.skipUnless(_EGGNOG_READY, "EggNOG database not found in dbs/eggnog/")
    def test_build_eggnog_real(self):
        """build --eggnog produces a non-empty eggnog_map.tsv."""
        result = self.runner.invoke(
            functional_db_cli,
            [
                "build", "-i", self.FAA, "-o", self.test_outdir,
                "--eggnog", "--eggnog-db", _EGGNOG_DB,
            ],
        )
        assert result.exit_code == 0, result.output
        out = os.path.join(self.test_outdir, "eggnog_map.tsv")
        assert os.path.exists(out), "eggnog_map.tsv was not created"
        assert os.path.getsize(out) > 0, "eggnog_map.tsv is empty"

    @unittest.skipUnless(
        _UNIREF_READY and _EGGNOG_READY,
        "UniRef and/or EggNOG databases not found",
    )
    def test_build_both_real(self):
        """build --uniref --eggnog produces both output files."""
        result = self.runner.invoke(
            functional_db_cli,
            [
                "build", "-i", self.FAA, "-o", self.test_outdir,
                "--uniref", "--uniref-db", _UNIREF_DB,
                "--eggnog", "--eggnog-db", _EGGNOG_DB,
            ],
        )
        assert result.exit_code == 0, result.output
        assert os.path.exists(os.path.join(self.test_outdir, "uniref_map.txt.xz"))
        assert os.path.exists(os.path.join(self.test_outdir, "eggnog_map.tsv"))


if __name__ == "__main__":
    unittest.main()
