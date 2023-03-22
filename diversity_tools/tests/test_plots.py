import unittest
import pandas as pd

from pathlib import Path

from src.matrix_operations import calculate_shannon_diversity_index
from src.plots import convert_diversity_matrix_to_graph

class  TestReaders(unittest.TestCase):

    def setUp(self):
        self.test_path = Path(__file__).parent.absolute() / "test_data"

    def seaborn_test(self):
        example_matrix = {"Community1" : {"black": 12, "purple": 21, "striped": 5,
                                          "green": 25, "brown": 2, "lblue": 17,
                                          "sblue": 9}, 
                          "Community2" : {"black": 10, "purple": 0, "striped": 15,
                                          "green": 15, "brown": 2, "lblue": 30,
                                          "sblue": 0}}
        df = pd.DataFrame(example_matrix)
        example_diversity_df = calculate_shannon_diversity_index(df)
        example_diversity_sns = convert_diversity_matrix_to_graph(example_diversity_df)