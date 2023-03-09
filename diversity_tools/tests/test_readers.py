import unittest

from pathlib import Path

from src.matrix_readers import read_orthovenn2_composition_output, read_matrix_from_file

class  TestReaders(unittest.TestCase):

    def setUp(self):
        self.test_path = Path(__file__).parent.absolute() / "test_data"

    def test_read_orthovenn2_matrix(self):
        fpath = self.test_path / "orthovenn_out.txt"
        with open(fpath) as orthovenn2_input_test_fhand:
            family_matrix_test = read_orthovenn2_composition_output(orthovenn2_input_test_fhand)
            assert "Niund100Scf103505g0000010.1" in family_matrix_test["CL00001"]["Nicotiana_undulata"]
            assert "NipanScf047676g0002.1" in family_matrix_test["CL00002"]["Nicotiana_paniculata"]
            assert "Nideb015S394759g0000010.1" in family_matrix_test["CL00001"]["Nicotiana_debneyi"]
            assert "Niafr015S145251g004.1" in family_matrix_test["CL00001"]["Nicotiana_africana"]
            assert "NipanScf149893g0003.1" not in family_matrix_test["CL00001"]["Nicotiana_undulata"]
            assert "Niund100Scf057563g0000010.1" not in family_matrix_test["CL00002"]["Nicotiana_paniculata"]
            assert len(family_matrix_test["CL00001"]["Nicotiana_undulata"]) == 420
    
    def test_read_matrix_from_file(self):
        fpath = self.test_path / "geneFamilies_ov2_out.csv"
        with open(fpath) as geneFamilies_ov2_input_fhand:
            geneFamilies_matrix_test = read_matrix_from_file(geneFamilies_ov2_input_fhand)
            assert geneFamilies_matrix_test["CL00001"]["Nicotiana_accuminata"] == 0
            assert len(list(geneFamilies_matrix_test.keys())) == 88902