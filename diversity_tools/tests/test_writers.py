import unittest
import csv

from os import remove as remove_file
from pathlib import Path

from src.matrix_writers import write_file_from_matrix

class  TestReaders(unittest.TestCase):

    def setUp(self):
        self.test_path = Path(__file__).parent.absolute() / "test_data"

    def tearDown(self):
        remove_file(self.test_path/ "geneFamiliesMatrix_to_File_input.csv")
    
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
        with open(fpath) as fhand:
            for line in csv.DictReader(fhand):
                print(line)
                assert line == {'#ID':'CL0001', 'SP1': '2', 'SP2': '3', 'SP3': '1', 'SP4': '0', 'SP5': '0'}
                assert line["#ID"] == "CL0001"
                assert int(line["SP1"]) == 2
                assert int(line["SP4"]) == 0
                break