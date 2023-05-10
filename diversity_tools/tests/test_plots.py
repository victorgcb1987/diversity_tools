import unittest
import pandas as pd
import os.path

from pathlib import Path
from os import remove as remove_file

from src.matrix_operations import (calculate_shannon_diversity_index, calculate_shannon_specialization_index,
                                   calculate_dataframe_frecuencies_col, calculate_dataframe_frecuencies_row, 
                                   calculate_shannon_specificity_index)
from src.plots import convert_diversity_matrix_to_graph, convert_dataframe_to_scatter

class  TestReaders(unittest.TestCase):

    def setUp(self):
        self.test_path = Path(__file__).parent.absolute() / "test_data"

    #def tearDown(self):
    #    remove_file(self.test_path/ "Families_Diversity_plot.svg")

    def test_seaborn(self):
        fpath = self.test_path / "Families_Diversity_plot.svg"
        example_matrix = {"Community1" : {"black": 0.13186, "purple": 0.23077, "striped": 0.05495,
                                          "green": 0.27473, "brown": 0.02198, "lblue": 0.18681,
                                          "sblue": 0.09890}, 
                          "Community2" : {"black": 0.13889, "purple": 0, "striped": 0.20833,
                                          "green": 0.20833, "brown": 0.02778, "lblue": 0.41667,
                                          "sblue": 0}}
        df = pd.DataFrame(example_matrix)
        example_diversity_df = calculate_shannon_diversity_index(df)
        #convert_diversity_matrix_to_graph(example_diversity_df, fpath)

    def test_scatter(self):
        fpath = self.test_path / "Scatter_plot.svg"
        example_matrix = {"Community1" : {"black": 0.13186, "purple": 0.23077, "striped": 0.05495,
                                          "green": 0.27473, "brown": 0.02198, "lblue": 0.18681,
                                          "sblue": 0.09890}, 
                          "Community2" : {"black": 0.13889, "purple": 0, "striped": 0.20833,
                                          "green": 0.20833, "brown": 0.02778, "lblue": 0.41667,
                                          "sblue": 0}}
        df = pd.DataFrame(example_matrix)
        diversity_df = calculate_shannon_diversity_index(df)
        frecuency_rows_df = calculate_dataframe_frecuencies_row(df)
        specificity_df = calculate_shannon_specificity_index(frecuency_rows_df)
        frecuency_Cols_df = calculate_dataframe_frecuencies_col(df)
        Shannon_index_df = calculate_shannon_specialization_index(frecuency_Cols_df, specificity_df)
        convert_dataframe_to_scatter(diversity_df, Shannon_index_df, fpath)