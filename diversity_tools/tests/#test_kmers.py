import unittest

from pathlib import Path
from shutil import rmtree as remove_folder

from src.kmers import run_meryl, count_meryl_kmers, group_kmers_by_species
from src.plots import venn_diagram_of_kmers


class TestKmers(unittest.TestCase):
    def setUp(self):
        self.test_path = Path(__file__).parent.absolute() / "test_data"
        self.counts_fpath = self.test_path / "counts"

    def test01_run_meryl(self):
        sequences_fpath = self.test_path / "sequences_to_kmer.fa"
        meryl_results = run_meryl(sequences_fpath, self.counts_fpath, threads=1, kmer_size=3)
        assert meryl_results["return_code"] == 0

    def test02_count_kmers(self):
        kmers = count_meryl_kmers(self.counts_fpath)
        assert kmers["AGC"] == 1
        assert kmers["CTC"] == 1
        remove_folder(self.counts_fpath)

    def test_03_group_kmers(self):
        kmers = [{"AAAAA": 1, "AAAAT": 1, "ATCAA": 1},
                 {"AAAAA": 2, "AAAAT": 1, "GACTT": 1}]
        sps_names = ["sp1", "sp2"]
        grouped_kmers = group_kmers_by_species(kmers, sps_names)
        assert grouped_kmers["sp1"]["AAAAA"] == 1
        assert grouped_kmers["sp2"]["AAAAA"] == 2
        assert grouped_kmers["sp1"]["ATCAA"] == 1
        assert grouped_kmers["sp2"]["ATCAA"] == 0
        assert grouped_kmers["sp1"]["GACTT"] == 0
        assert grouped_kmers["sp2"]["GACTT"] == 1

    def test_04_venn_diagram(self):
        kmers = [{"AAAAA": 1, "AAAAT": 1, "ATCAA": 1},
                 {"AAAAA": 2, "AAAAT": 1, "GACTT": 1},
                 {"AAAAA": 1, "AAAAT": 1, "ATCGA": 1}]
        sps_names = ["sp1", "sp2", "sp3"]
        grouped_kmers = group_kmers_by_species(kmers, sps_names)

