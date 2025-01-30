import os
import unittest
from click.testing import CliRunner
from unittest.mock import patch, MagicMock
from wipe.wipe import uniref as uniref_cli
import subprocess
import gzip
from pathlib import Path

class TestUniRef(unittest.TestCase):
    def setUp(self):
        self.test_outdir = os.path.join("tests", "out", "uniref")
        self.runner = CliRunner()
        os.makedirs(self.test_outdir, exist_ok=True)
        
        # create test data directory and file
        self.test_data_dir = os.path.join("tests", "data", "uniref")
        os.makedirs(self.test_data_dir, exist_ok=True)
        
        # create a small test FASTA file
        self.test_fasta = os.path.join(self.test_data_dir, "mini_uniref90.fasta")
        with open(self.test_fasta, 'w') as f:
            f.write(">UniRef90_A0A000|Test protein 1\n")
            f.write("MVKVKYFGCELNGYTDVVEVKDGKVIGESRKILQTAIDRLKENGVKTIVFAR\n")
            f.write(">UniRef90_A0A001|Test protein 2\n")
            f.write("MEEAKQKVVDFLNSKSGSKSKGFHPDTFVGQIAAMVDKLAAQGKEVAVFEAR\n")
        
    def tearDown(self):
        # Clean up test directories
        if os.path.exists(self.test_outdir):
            import shutil
            shutil.rmtree(self.test_outdir)
        if os.path.exists(self.test_data_dir):
            import shutil
            shutil.rmtree(self.test_data_dir)

    @patch('subprocess.run')
    def test_download_uniref50(self, mock_subprocess_run):
        """Test downloading UniRef50 database"""
        mock_subprocess_run.return_value = subprocess.CompletedProcess(
            args=[],
            returncode=0,
            stdout="",
            stderr=""
        )
        
        result = self.runner.invoke(
            uniref_cli,
            ['download', '--level', '50', '--outdir', self.test_outdir]
        )
        
        # verify request was successful
        assert result.exit_code == 0
        
        expected_url = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref50/"
        expected_cmd = [
            "wget", "-r", "-np", "-nH", 
            "--cut-dirs=6", "-R", "index.html*",
            "-P", self.test_outdir,
            expected_url
        ]
        
        # verify request was correct
        mock_subprocess_run.assert_called_once_with(
            expected_cmd,
            cwd=None,
            capture_output=True,
            text=True,
            check=True,
            shell=False
        )

    def test_download_invalid_level(self):
        """Test that invalid UniRef level raises error"""
        result = self.runner.invoke(
            uniref_cli,
            ['download', '--level', '100', '--outdir', self.test_outdir]
        )
        
        # verify request failed with correct error
        assert result.exit_code != 0
        assert "Invalid value for '--level'" in result.output

    @patch('subprocess.run')
    def test_build_diamond_db(self, mock_subprocess_run):
        """Test building DIAMOND database from UniRef FASTA"""
        mock_subprocess_run.return_value = subprocess.CompletedProcess(
            args=[],
            returncode=0,
            stdout="",
            stderr=""
        )
        
        result = self.runner.invoke(
            uniref_cli,
            ['build',
             '--input-fasta', self.test_fasta,
             '--output-db', self.test_outdir,
             '--threads', '4']
        )
        
        # verify request was successful
        assert result.exit_code == 0
        
        expected_cmd = [
            "diamond", "makedb",
            "--in", self.test_fasta,
            "--db", self.test_outdir,
            "--threads", "4"
        ]
        
        # verify request was correct
        mock_subprocess_run.assert_called_once_with(
            expected_cmd,
            cwd=None,
            capture_output=True,
            text=True,
            check=True,
            shell=False
        )

    def test_build_missing_input(self):
        """Test building database with missing input file"""
        input_fasta = os.path.join(self.test_data_dir, "nonexistent.fasta")
        
        result = self.runner.invoke(
            uniref_cli,
            ['build',
             '--input-fasta', input_fasta,
             '--output-db', self.test_outdir]
        )
        
        # verify request failed with correct error
        assert result.exit_code != 0
        assert "does not exist" in result.output

    @patch('subprocess.run')
    def test_process_uniref_xml(self, mock_subprocess_run):
        """Test processing UniRef XML file"""
        test_xml = os.path.join(self.test_data_dir, "test_uniref.xml.gz")
        with gzip.open(test_xml, 'wt') as f:
            f.write('<entry id="UniRef90_A0A000" updated="2024-01-01">\n')
            f.write('<name>Cluster: Test protein 1</name>\n')
            f.write('<property type="member count" value="1"/>\n')
            f.write('<property type="common taxon" value="Bacteria"/>\n')
            f.write('<property type="common taxon ID" value="2"/>\n')
            f.write('<property type="other property" value="test"/>\n')
            f.write('</entry>\n')
        
        output_tsv = os.path.join(self.test_outdir, "test_output.tsv")
        
        result = self.runner.invoke(
            uniref_cli,
            ['process', '--xml-file', test_xml, '--output-tsv', output_tsv]
        )
        
        assert result.exit_code == 0
        
        # verify output file exists and contains correct content
        assert os.path.exists(output_tsv)
        with open(output_tsv, 'r') as f:
            lines = f.readlines()
            assert len(lines) >= 2  # at least header and one data line
            assert lines[0].strip() == "A0A000\tTest protein 1\t1\tBacteria\t2"
            assert lines[1].strip() == "# other property\ttest"

    def test_process_invalid_xml(self):
        """Test processing invalid UniRef XML file"""
        test_xml = os.path.join(self.test_data_dir, "invalid_uniref.xml.gz")
        with gzip.open(test_xml, 'wt') as f:
            # write entry line that will fail the regex validation
            f.write('<entry id="INVALID_FORMAT" updated="2024-01-01">\n')
        
        output_tsv = os.path.join(self.test_outdir, "invalid_output.tsv")
        
        result = self.runner.invoke(
            uniref_cli,
            ['process', '--xml-file', test_xml, '--output-tsv', output_tsv]
        )
        
        print(f"\nExit code: {result.exit_code}")
        print(f"Exception: {result.exception}")
        print(f"Output: {result.output}")
        
        # verify command failed appropriately
        assert result.exit_code != 0
        assert "Error processing UniRef XML" in result.output

    def test_process_missing_xml(self):
        """Test processing non-existent XML file"""
        nonexistent_xml = os.path.join(self.test_data_dir, "nonexistent.xml.gz")
        output_tsv = os.path.join(self.test_outdir, "output.tsv")
        
        result = self.runner.invoke(
            uniref_cli,
            ['process', '--xml-file', nonexistent_xml, '--output-tsv', output_tsv]
        )
        
        # verify command failed appropriately
        assert result.exit_code != 0
        assert "does not exist" in result.output

    @patch('subprocess.run')
    def test_blastp_success(self, mock_subprocess_run):
        """Test running DIAMOND BLASTP against UniRef database"""
        mock_subprocess_run.return_value = subprocess.CompletedProcess(
            args=[],
            returncode=0,
            stdout="",
            stderr=""
        )
        
        test_indir = os.path.join(self.test_data_dir, "proteins")
        os.makedirs(test_indir, exist_ok=True)
        
        protein_file = os.path.join(test_indir, "all.faa")
        with open(protein_file, 'w') as f:
            f.write(">protein1\n")
            f.write("MVKVKYFGCELNGYTDVVEVKDGKVIGESRKILQTAIDRLKENGVKTIVFAR\n")
        
        test_db = os.path.join(self.test_data_dir, "uniref90.dmnd")
        Path(test_db).touch()
        
        os.makedirs(self.test_outdir, exist_ok=True)
        diamond_results = os.path.join(self.test_outdir, "diamond_results.m8")
        with open(diamond_results, 'w') as f:
            f.write("protein1\tUniRef90_TEST\t100.0\n")
        
        result = self.runner.invoke(
            uniref_cli,
            ['blastp', 
             '--indir', test_indir,
             '--diamonddb', test_db,
             '--outdir', self.test_outdir,
             '--threads', '4']
        )
        
        assert result.exit_code == 0
        
        expected_diamond_cmd = [
            "diamond", "blastp",
            "--threads", "4",
            "--db", test_db,
            "--query", protein_file,
            "--out", diamond_results,
            "--id", "90",
            "--subject-cover", "80",
            "--query-cover", "80",
            "--index-chunks", "1",
            "--max-target-seqs", "1"
        ]
        
        expected_xz_cmd = [
            "xz", "-z",
            os.path.join(self.test_outdir, "uniref.map")
        ]
        
        # verify both commands were called in correct order
        assert mock_subprocess_run.call_count == 2
        mock_subprocess_run.assert_has_calls([
            unittest.mock.call(
                expected_diamond_cmd,
                cwd=None,
                capture_output=True,
                text=True,
                check=True,
                shell=False
            ),
            unittest.mock.call(
                expected_xz_cmd,
                cwd=None,
                capture_output=True,
                text=True,
                check=True,
                shell=False
            )
        ])

    def test_blastp_missing_input(self):
        """Test BLASTP with missing input directory"""
        missing_dir = os.path.join(self.test_data_dir, "nonexistent")
        test_db = os.path.join(self.test_data_dir, "uniref90.dmnd")
        
        result = self.runner.invoke(
            uniref_cli,
            ['blastp',
             '--indir', missing_dir,
             '--diamonddb', test_db,
             '--outdir', self.test_outdir]
        )
        
        # verify command failed appropriately
        assert result.exit_code != 0
        assert "does not exist" in result.output

    def test_blastp_missing_db(self):
        """Test BLASTP with missing database file"""
        test_indir = os.path.join(self.test_data_dir, "proteins")
        os.makedirs(test_indir, exist_ok=True)
        missing_db = os.path.join(self.test_data_dir, "nonexistent.dmnd")
        
        result = self.runner.invoke(
            uniref_cli,
            ['blastp',
             '--indir', test_indir,
             '--diamonddb', missing_db,
             '--outdir', self.test_outdir]
        )
        
        # verify command failed appropriately
        assert result.exit_code != 0
        assert "does not exist" in result.output

    def test_merge_maps_success(self):
        """Test merging UniRef90 and UniRef50 mapping files"""
        # Create test mapping files
        uniref90_map = os.path.join(self.test_data_dir, "uniref90.map")
        with open(uniref90_map, 'w') as f:
            f.write("genome1_1\tUniRef90_A0A000\n")
            f.write("genome1_2\tUniRef90_A0A001\n")
        
        uniref50_map = os.path.join(self.test_data_dir, "uniref50.map")
        with open(uniref50_map, 'w') as f:
            f.write("genome1_2\tUniRef50_B0B000\n")  # should be ignored in favor of UniRef90
            f.write("genome1_3\tUniRef50_B0B001\n")
        
        output_file = os.path.join(self.test_outdir, "merged.map")
        
        result = self.runner.invoke(
            uniref_cli,
            ['merge-maps',
             '--uniref90-map', uniref90_map,
             '--uniref50-map', uniref50_map,
             '--output-file', output_file]
        )
        
        # verify command was successful
        assert result.exit_code == 0
        
        # verify output file exists and contains correct content
        assert os.path.exists(output_file)
        with open(output_file, 'r') as f:
            lines = set(f.readlines())
            assert "genome1_1\tUniRef90_A0A000\n" in lines
            assert "genome1_2\tUniRef90_A0A001\n" in lines 
            assert "genome1_3\tUniRef50_B0B001\n" in lines

    def test_merge_maps_with_simplify(self):
        """Test merging UniRef maps with simplify option"""
        uniref90_map = os.path.join(self.test_data_dir, "uniref90.map")
        with open(uniref90_map, 'w') as f:
            f.write("genome1_1\tUniRef90_A0A000\n")
        
        uniref50_map = os.path.join(self.test_data_dir, "uniref50.map")
        with open(uniref50_map, 'w') as f:
            f.write("genome1_2\tUniRef50_B0B000\n")
        
        output_file = os.path.join(self.test_outdir, "merged_simple.map")
        
        result = self.runner.invoke(
            uniref_cli,
            ['merge-maps',
             '--uniref90-map', uniref90_map,
             '--uniref50-map', uniref50_map,
             '--output-file', output_file,
             '--simplify']
        )
        
        assert result.exit_code == 0
        
        # verify simplified IDs in output
        with open(output_file, 'r') as f:
            lines = set(f.readlines())
            assert "genome1_1\tA0A000\n" in lines  # UniRef90_ prefix removed
            assert "genome1_2\tB0B000\n" in lines  # UniRef50_ prefix removed

    def test_merge_maps_missing_input(self):
        """Test merge-maps with missing input file"""
        nonexistent = os.path.join(self.test_data_dir, "nonexistent.map")
        output_file = os.path.join(self.test_outdir, "merged.map")
        
        result = self.runner.invoke(
            uniref_cli,
            ['merge-maps',
             '--uniref90-map', nonexistent,
             '--uniref50-map', nonexistent,
             '--output-file', output_file]
        )
        
        # verify command failed appropriately
        assert result.exit_code != 0
        assert "does not exist" in result.output 

    def test_extract_names_success(self):
        """Test extracting UniRef names from mapping file"""
        map_file = os.path.join(self.test_data_dir, "test.map")
        with open(map_file, 'w') as f:
            f.write("genome1_1\tUniRef90_A0A000\n")
            f.write("genome1_2\tUniRef50_B0B000\n")
        
        uniref90_names = os.path.join(self.test_data_dir, "uniref90.names")
        with open(uniref90_names, 'w') as f:
            f.write("UniRef90_A0A000\tProtein A\tExtra info\n")
        
        uniref50_names = os.path.join(self.test_data_dir, "uniref50.names")
        with open(uniref50_names, 'w') as f:
            f.write("UniRef50_B0B000\tProtein B\tExtra info\n")
        
        output_names = os.path.join(self.test_outdir, "extracted.names")
        
        result = self.runner.invoke(
            uniref_cli,
            ['extract-names',
             '--map-file', map_file,
             '--uniref90-names', uniref90_names,
             '--uniref50-names', uniref50_names,
             '--output-names', output_names]
        )
        
        assert result.exit_code == 0
        
        # verify output file exists and contains correct content
        assert os.path.exists(output_names)
        with open(output_names, 'r') as f:
            lines = set(f.readlines())
            assert "UniRef90_A0A000\tProtein A\n" in lines
            assert "UniRef50_B0B000\tProtein B\n" in lines

    def test_extract_names_missing_id(self):
        """Test extracting names with missing IDs"""
        # create test mapping file with ID that doesn't exist in names files
        map_file = os.path.join(self.test_data_dir, "test.map")
        with open(map_file, 'w') as f:
            f.write("genome1_1\tUniRef90_MISSING\n")
        

        uniref90_names = os.path.join(self.test_data_dir, "uniref90.names")
        with open(uniref90_names, 'w') as f:
            f.write("UniRef90_A0A000\tProtein A\tExtra info\n")
        
        uniref50_names = os.path.join(self.test_data_dir, "uniref50.names")
        with open(uniref50_names, 'w') as f:
            f.write("UniRef50_B0B000\tProtein B\tExtra info\n")
        
        output_names = os.path.join(self.test_outdir, "extracted.names")
        
        result = self.runner.invoke(
            uniref_cli,
            ['extract-names',
             '--map-file', map_file,
             '--uniref90-names', uniref90_names,
             '--uniref50-names', uniref50_names,
             '--output-names', output_names]
        )
        
        # verify command was successful (missing IDs should be skipped)
        assert result.exit_code == 0
        
        # verify output file exists but doesn't contain missing ID
        assert os.path.exists(output_names)
        with open(output_names, 'r') as f:
            content = f.read()
            assert "UniRef90_MISSING" not in content

    def test_extract_names_missing_input(self):
        """Test extract-names with missing input file"""
        nonexistent = os.path.join(self.test_data_dir, "nonexistent.map")
        output_names = os.path.join(self.test_outdir, "extracted.names")
        
        result = self.runner.invoke(
            uniref_cli,
            ['extract-names',
             '--map-file', nonexistent,
             '--uniref90-names', nonexistent,
             '--uniref50-names', nonexistent,
             '--output-names', output_names]
        )
        
        # verify command failed appropriately
        assert result.exit_code != 0
        assert "does not exist" in result.output 

        