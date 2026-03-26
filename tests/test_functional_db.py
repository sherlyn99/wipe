import os
import shutil
import tempfile
import unittest
import subprocess
from pathlib import Path
from unittest.mock import patch
from click.testing import CliRunner
from wipe.wipe import functional_db as functional_db_cli
from wipe.modules.functional_db import make_eggnog_mappings, annotate_uniref, annotate_eggnog

_HERE = os.path.dirname(__file__)


class TestFunctionalDbDownload(unittest.TestCase):
    def setUp(self):
        self.runner = CliRunner()
        self.test_outdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.test_outdir)

    @patch("wipe.modules.functional_db.process_uniref_xml")
    @patch("subprocess.run")
    def test_download_uniref_only(self, mock_run, mock_process_xml):
        """--uniref downloads uniref90/50 FASTAs and XMLs, extracts names, and builds diamond databases."""
        mock_run.return_value = subprocess.CompletedProcess([], 0, "", "")
        result = self.runner.invoke(
            functional_db_cli,
            ["download", "--uniref", "--no-eggnog", "-o", self.test_outdir],
        )
        assert result.exit_code == 0, result.output
        calls = mock_run.call_args_list
        assert any("uniref90.fasta.gz" in str(c) for c in calls)
        assert any("uniref50.fasta.gz" in str(c) for c in calls)
        assert any("uniref90.xml.gz" in str(c) for c in calls)
        assert any("uniref50.xml.gz" in str(c) for c in calls)
        assert any("uniref90.dmnd" in str(c) for c in calls)
        assert any("uniref50.dmnd" in str(c) for c in calls)
        diamond_calls = [c for c in calls if "diamond" in str(c) and "makedb" in str(c)]
        assert len(diamond_calls) == 2
        assert mock_process_xml.call_count == 2
        xml_calls = mock_process_xml.call_args_list
        assert any("uniref90.xml.gz" in str(c) and "uniref90_names.tsv" in str(c) for c in xml_calls)
        assert any("uniref50.xml.gz" in str(c) and "uniref50_names.tsv" in str(c) for c in xml_calls)
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

    @patch("wipe.modules.functional_db.process_uniref_xml")
    @patch("subprocess.run")
    def test_download_both_by_default(self, mock_run, mock_process_xml):
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

    @patch("wipe.modules.functional_db.process_uniref_xml")
    @patch("subprocess.run")
    def test_download_uniref_creates_subdirectory(self, mock_run, mock_process_xml):
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

    @patch("wipe.modules.functional_db.extract_uniref_names")
    @patch("wipe.modules.functional_db.merge_uniref")
    @patch("subprocess.run")
    def test_build_uniref_extracts_names_when_present(self, mock_run, mock_merge, mock_extract):
        """extract_uniref_names is called when name TSVs exist in uniref_db_dir."""
        mock_run.return_value = subprocess.CompletedProcess([], 0, "", "")
        merged = os.path.join(self.test_outdir, "uniref_map.txt")
        mock_merge.side_effect = lambda *a, **kw: Path(merged).touch()
        Path(os.path.join(self.uniref_db_dir, "uniref90_names.tsv")).touch()
        Path(os.path.join(self.uniref_db_dir, "uniref50_names.tsv")).touch()

        self.runner.invoke(
            functional_db_cli,
            [
                "build", "-i", self.faa, "-o", self.test_outdir,
                "--uniref", "--uniref-db", self.uniref_db_dir,
            ],
        )
        mock_extract.assert_called_once()
        args = mock_extract.call_args[0]
        assert args[0] == merged
        assert "uniref90_names.tsv" in args[1]
        assert "uniref50_names.tsv" in args[2]
        assert args[3] == os.path.join(self.test_outdir, "uniref_names.txt")

    @patch("wipe.modules.functional_db.extract_uniref_names")
    @patch("wipe.modules.functional_db.merge_uniref")
    @patch("subprocess.run")
    def test_build_uniref_skips_names_when_absent(self, mock_run, mock_merge, mock_extract):
        """extract_uniref_names is not called when name TSVs are missing."""
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
        mock_extract.assert_not_called()

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

    @patch("wipe.modules.functional_db.extract_uniref_names")
    @patch("subprocess.run")
    def test_build_uniref_creates_names_file_when_tsv_present(self, mock_run, mock_extract):
        """uniref_names.txt is created when name TSVs exist in uniref_db_dir."""
        mock_run.side_effect = self._fake_run
        Path(os.path.join(self.uniref_db_dir, "uniref90_names.tsv")).touch()
        Path(os.path.join(self.uniref_db_dir, "uniref50_names.tsv")).touch()

        result = self.runner.invoke(
            functional_db_cli,
            [
                "build", "-i", self.FAA, "-o", self.test_outdir,
                "--uniref", "--uniref-db", self.uniref_db_dir,
            ],
        )
        assert result.exit_code == 0, result.output
        mock_extract.assert_called_once()
        assert mock_extract.call_args[0][3] == os.path.join(self.test_outdir, "uniref_names.txt")

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


class TestFunctionalDbBuildSubset(unittest.TestCase):
    """Integration tests using real (tiny) databases stored in tests/data/.

    These tests run real DIAMOND blastp against a 2-sequence subset UniRef
    database so no mocking of external tools is needed for the UniRef path.
    For EggNOG, emapper.py is replaced by a side-effect that copies the
    pre-baked annotations fixture, letting all Python-side logic run for real.
    """

    FAA = os.path.join(_HERE, "data", "M000000999_subset.faa")
    UNIREF_DB_DIR = os.path.join(_HERE, "data", "uniref")
    EGGNOG_ANNOT_FIXTURE = os.path.join(_HERE, "data", "eggnog", "eggnog_test.emapper.annotations")

    def setUp(self):
        self.outdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.outdir)

    # ------------------------------------------------------------------
    # UniRef tests — fully real (diamond + merge_uniref + xz)
    # ------------------------------------------------------------------

    def test_uniref_produces_compressed_map(self):
        """annotate_uniref creates a non-empty uniref_map.txt.xz."""
        annotate_uniref(self.FAA, self.UNIREF_DB_DIR, self.outdir, threads=1)
        out = os.path.join(self.outdir, "uniref_map.txt.xz")
        self.assertTrue(os.path.exists(out), "uniref_map.txt.xz not created")
        self.assertGreater(os.path.getsize(out), 0, "uniref_map.txt.xz is empty")

    def test_uniref_map_contains_both_orfs(self):
        """Both ORFs in the test FAA appear in the merged map."""
        annotate_uniref(self.FAA, self.UNIREF_DB_DIR, self.outdir, threads=1)
        import lzma
        with lzma.open(os.path.join(self.outdir, "uniref_map.txt.xz"), "rt") as f:
            content = f.read()
        self.assertIn("M000000999_1", content)
        self.assertIn("M000000999_2", content)

    def test_uniref_names_extracted_when_tsv_present(self):
        """uniref_names.txt is created because the test DB has names TSVs."""
        annotate_uniref(self.FAA, self.UNIREF_DB_DIR, self.outdir, threads=1)
        names_file = os.path.join(self.outdir, "uniref_names.txt")
        self.assertTrue(os.path.exists(names_file), "uniref_names.txt not created")
        with open(names_file) as f:
            content = f.read()
        self.assertTrue(len(content.strip()) > 0, "uniref_names.txt is empty")

    def test_uniref_simplified_ids(self):
        """Merged map uses simplified IDs (part after '_' in UniRef header)."""
        annotate_uniref(self.FAA, self.UNIREF_DB_DIR, self.outdir, threads=1)
        import lzma
        with lzma.open(os.path.join(self.outdir, "uniref_map.txt.xz"), "rt") as f:
            ids = [line.strip().split("\t")[1] for line in f if line.strip()]
        # Simplified IDs should not contain "UniRef" prefix
        for uid in ids:
            self.assertNotIn("UniRef", uid, f"ID not simplified: {uid}")

    # ------------------------------------------------------------------
    # EggNOG tests — real Python pipeline, mocked emapper.py call
    # ------------------------------------------------------------------

    def _fake_emapper(self, cmd, **kwargs):
        """Copy pre-baked annotations fixture to emapper's expected output path."""
        out_idx = cmd.index("--output_dir")
        outdir = cmd[out_idx + 1]
        dest = os.path.join(outdir, "eggnog.emapper.annotations")
        shutil.copy(self.EGGNOG_ANNOT_FIXTURE, dest)
        return subprocess.CompletedProcess(cmd, 0, "", "")

    @patch("subprocess.run")
    def test_eggnog_creates_map_and_per_annotation_files(self, mock_run):
        """annotate_eggnog creates eggnog_map.tsv and all six orf_to_*.tsv files."""
        mock_run.side_effect = self._fake_emapper
        annotate_eggnog(self.FAA, os.path.join(_HERE, "data", "eggnog"), self.outdir, threads=1)
        self.assertTrue(os.path.exists(os.path.join(self.outdir, "eggnog_map.tsv")))
        for fname in ("orf_to_go.tsv", "orf_to_ec.tsv", "orf_to_ko.tsv",
                      "orf_to_cazy.tsv", "orf_to_cog.tsv", "orf_to_pfam.tsv"):
            self.assertTrue(
                os.path.exists(os.path.join(self.outdir, fname)),
                f"{fname} was not created",
            )

    @patch("subprocess.run")
    def test_eggnog_map_content(self, mock_run):
        """eggnog_map.tsv contains the two ORF rows from the fixture."""
        mock_run.side_effect = self._fake_emapper
        annotate_eggnog(self.FAA, os.path.join(_HERE, "data", "eggnog"), self.outdir, threads=1)
        with open(os.path.join(self.outdir, "eggnog_map.tsv")) as f:
            content = f.read()
        self.assertIn("M000000999_1", content)
        self.assertIn("M000000999_2", content)

    @patch("subprocess.run")
    def test_eggnog_mapping_file_content(self, mock_run):
        """orf_to_go.tsv and orf_to_cog.tsv have correct entries from fixture."""
        mock_run.side_effect = self._fake_emapper
        annotate_eggnog(self.FAA, os.path.join(_HERE, "data", "eggnog"), self.outdir, threads=1)

        with open(os.path.join(self.outdir, "orf_to_go.tsv")) as f:
            go_rows = [line.strip().split("\t") for line in f]
        go_ids = {row[1] for row in go_rows}
        self.assertIn("GO:0008150", go_ids)
        self.assertIn("GO:0055114", go_ids)

        with open(os.path.join(self.outdir, "orf_to_cog.tsv")) as f:
            cog_rows = [line.strip().split("\t") for line in f]
        cog_ids = {row[1] for row in cog_rows}
        self.assertIn("S", cog_ids)
        self.assertIn("C", cog_ids)
        self.assertIn("P", cog_ids)

    @patch("subprocess.run")
    def test_eggnog_original_annotations_removed(self, mock_run):
        """eggnog.emapper.annotations is moved to eggnog_map.tsv (not left behind)."""
        mock_run.side_effect = self._fake_emapper
        annotate_eggnog(self.FAA, os.path.join(_HERE, "data", "eggnog"), self.outdir, threads=1)
        self.assertFalse(
            os.path.exists(os.path.join(self.outdir, "eggnog.emapper.annotations")),
            "Original annotations file should have been renamed",
        )


if __name__ == "__main__":
    unittest.main()
