import unittest
import csv

from pathlib import Path

from src.matrix_writers import write_file_from_matrix

class  TestReaders(unittest.TestCase):

    def setUp(self):
        self.test_path = Path(__file__).parent.absolute() / "test_data"
    
    def test_write_file_from_matrix(self):
        fpath = self.test_path / "geneFamiliesMatrix_to_File_input.csv"
        family_matrix = {"CL0001": {"SP1": 2, 
                                    "SP2": 3, 
                                    "SP3": 1},
                        "CL0002": {"SP1": 1, 
                                    "SP2": 2, 
                                    "SP3": 3, 
                                    "SP4": 2, 
                                    "SP5": 2}}
        write_file_from_matrix(family_matrix, fpath)
        assert family_matrix["CL0001"]["SP2"] == 3
        for line in csv.DictReader(self):