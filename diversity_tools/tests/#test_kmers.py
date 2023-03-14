import unittest

from pathlib import Path
from shutil import rmtree as remove_folder

from src.kmers import run_meryl, count_meryl_kmers, group_kmers


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

    
    # def test_03_group_kmers(self):
    #     kmers = {"sp1": {"AAAAA": 1, "AAAAT": 1, "ATCAA": 1},
    #              "sp2": {"AAAAA": 2, "AAAAT": 1, "GACTT": 1}}
    #     grouped_kmers = group_kmers(kmers)
    #     assert group_kmers[""]