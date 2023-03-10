import unittest

from pathlib import Path
from shutil import rmtree as remove_folder
from time import sleep

from src.kmers import run_meryl, count_meryl_kmers


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


    