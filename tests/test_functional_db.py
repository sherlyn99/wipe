import os
import shutil
import tempfile
import unittest
import subprocess
from pathlib import Path
from unittest.mock import patch
from click.testing import CliRunner
from wipe.wipe import functional_db as functional_db_cli


class TestFunctionalDbDownload(unittest.TestCase):
    def setUp(self):
        self.runner = CliRunner()
        self.test_outdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.test_outdir)

    @patch("subprocess.run")
    def test_download_uniref_only(self, mock_run):
        """--uniref flag downloads uniref90 and uniref50 via wget."""
        mock_run.return_value = subprocess.CompletedProcess([], 0, "", "")
        result = self.runner.invoke(
            functional_db_cli,
            ["download", "--uniref", "--no-eggnog", "-o", self.test_outdir],
        )
        assert result.exit_code == 0, result.output
        calls = mock_run.call_args_list
        assert any("uniref90.fasta.gz" in str(c) for c in calls)
        assert any("uniref50.fasta.gz" in str(c) for c in calls)
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
        """eggnog.emapper.annotations is renamed to eggnog_map.tsv."""
        mock_run.return_value = subprocess.CompletedProcess([], 0, "", "")
        annotations = os.path.join(self.test_outdir, "eggnog.emapper.annotations")
        with open(annotations, "w") as f:
            f.write("G000000001_1\tCOG001\n")

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


if __name__ == "__main__":
    unittest.main()
