import unittest
import pandas as pd
import os.path

from pathlib import Path
from os import remove as remove_file

from src.matrix_operations import calculate_shannon_diversity_index
from src.plots import convert_diversity_matrix_to_graph

class  TestReaders(unittest.TestCase):

    def setUp(self):
        self.test_path = Path(__file__).parent.absolute() / "test_data"

    def tearDown(self):
        remove_file(self.test_path/ "Families_Diversity_plot.svg")

    def test_seaborn(self):
        fpath = self.test_path / "Families_Diversity_plot.svg"
        example_matrix = {"Community1" : {"black": 12, "purple": 21, "striped": 5,
                                          "green": 25, "brown": 2, "lblue": 17,
                                          "sblue": 9}, 
                          "Community2" : {"black": 10, "purple": 0, "striped": 15,
                                          "green": 15, "brown": 2, "lblue": 30,
                                          "sblue": 0}}
        df = pd.DataFrame(example_matrix)
        example_diversity_df = calculate_shannon_diversity_index(df)
        convert_diversity_matrix_to_graph(example_diversity_df, fpath)