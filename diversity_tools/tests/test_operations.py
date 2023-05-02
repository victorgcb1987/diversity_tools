import unittest
import pandas as pd

from pathlib import Path

from src.matrix_operations import (convert_list_to_numbers, convert_into_dataframe, 
                                   calculate_shannon_diversity_index, get_dataframe_with_limited_families, 
                                   filter_dataframe_cols_by_value_occurrence, 
                                   select_families_with_highest_number_of_genes, filter_dataframe_by_cols_name,
                                   filter_dataframe_by_rows_name, calculate_shannon_specificity_index)                            

class  TestReaders(unittest.TestCase):

    def setUp(self):
        self.test_path = Path(__file__).parent.absolute() / "test_data"

    def test_convert_to_numbers(self):
        gene_families = {"CL0001": {"SP1": ["gen1", "gen2"], 
                                    "SP2": ["gen1", "gen2", "gen4"], 
                                    "SP3": ["gen3"]},
                        "CL0002": {"SP1": ["gen7"], 
                                    "SP2": ["gen6", "gen9"], 
                                    "SP3": ["gen1", "gen2", "gen4"], 
                                    "SP4": ["gen3"], 
                                    "SP5": ["gen7", "gen9"]}}
        gene_families_into_numbers = convert_list_to_numbers(gene_families)
        assert gene_families_into_numbers == {"CL0001": {"SP1": 2, 
                                                        "SP2": 3, 
                                                        "SP3": 1},
                                              "CL0002": {"SP1": 1, 
                                                        "SP2": 2, 
                                                        "SP3": 3, 
                                                        "SP4": 1, 
                                                        "SP5": 2}} 
        
    def test_matrix_dataframe(self):
        families_matrix = {"CL0001": {"SP1": 2, 
                                      "SP2": 3, 
                                      "SP3": 1},
                            "CL0002": {"SP1": 1, 
                                       "SP2": 2, 
                                       "SP3": 3, 
                                       "SP4": 1, 
                                       "SP5": 2}} 
        df = convert_into_dataframe(families_matrix)
    #    assert df.at["CL0001", "SP1"] == 2
    #    assert df.at["CL0002", "SP4"] == 1

    def test_get_dataframe_with_limited_families(self):
        families_matrix = {"CL0001": {"SP1": 2, 
                                      "SP2": 3, 
                                      "SP3": 1},
                            "CL0002": {"SP1": 1, 
                                       "SP2": 2, 
                                       "SP3": 3, 
                                       "SP4": 1, 
                                       "SP5": 2},
                            "CL0003": {"SP1": 1, 
                                       "SP2": 2}}
        df = convert_into_dataframe(families_matrix)
        df = get_dataframe_with_limited_families(df, value= 1)

    
    def test_calculate_shannon_index(self):
        example_matrix = {"Community1" : {"black": 0.13186, "purple": 0.23077, "striped": 0.05495,
                                          "green": 0.27473, "brown": 0.02198, "lblue": 0.18681,
                                          "sblue": 0.09890}, 
                          "Community2" : {"black": 0.13889, "purple": 0, "striped": 0.20833,
                                          "green": 0.20833, "brown": 0.02778, "lblue": 0.41667,
                                          "sblue": 0}}
        #Dos comunidades mas con diferentes diversidad
        df = pd.DataFrame(example_matrix)
        example_diversity_df = calculate_shannon_diversity_index(df)
        assert round(example_diversity_df['Community1'], 3) == 1.746
        assert round(example_diversity_df['Community2'], 3) == 1.392

    def test_calculate_specificity_index(self):
        example_matrix = {"Community1" : {"P1": 1, "P2": 0, "P3": 0,
                                          "P4": 0}, 
                          "Community2" : {"P1": 0, "P2": 0.5, "P3": 0,
                                          "P4": 0.5}, 
                          "Community3" : {"P1": 0.25, "P2": 0.25, "P3": 0.25,
                                          "P4": 0.25}, 
                          "Community4" : {"P1": 0.1, "P2": 0.3, "P3": 0.2,
                                          "P4": 0.4}}
        df = pd.DataFrame(example_matrix)
        example_specificity_df = calculate_shannon_specificity_index(df)

    def test_filter_dataframe_cols_by_value_occurrence(self):
        families_matrix = {"CL0001": {"SP1": 2, 
                                      "SP2": 3, 
                                      "SP3": 2,
                                      "SP4": 4, 
                                       "SP5": 2},
                            "CL0002": {"SP1": 1, 
                                       "SP2": 2, 
                                       "SP3": 3, 
                                       "SP4": 2, 
                                       "SP5": 2},
                            "CL0003": {"SP1": 1, 
                                       "SP2": 1,
                                       "SP3": 3, 
                                       "SP4": 5, 
                                       "SP5": 2}}
        df = convert_into_dataframe(families_matrix)
        df = filter_dataframe_cols_by_value_occurrence(df, value= 2,ignore_zeros=False, threshold=0.3, mode="equal")

    def test_families_with_higher_number_of_genes(self):
        families_matrix = {"CL0001": {"SP1": 3, 
                                      "SP2": 1, 
                                      "SP3": 3},
                            "CL0002": {"SP1": 3, 
                                       "SP2": 2, 
                                       "SP3": 2, 
                                       "SP4": 2, 
                                       "SP5": 2},
                            "CL0003": {"SP1": 1, 
                                       "SP2": 2}}
        df = convert_into_dataframe(families_matrix)
        df = select_families_with_highest_number_of_genes(df, 1)

    def test_filter_dataframe_by_cols_name(self):
        families_matrix = {"CL0001": {"SP1": 3, 
                                      "SP2": 1, 
                                      "SP3": 3},
                            "CL0002": {"SP1": 3, 
                                       "SP2": 2, 
                                       "SP3": 2, 
                                       "SP4": 2, 
                                       "SP5": 2},
                            "CL0003": {"SP1": 1, 
                                       "SP2": 2}}
        df = convert_into_dataframe(families_matrix)
        column_list = "CL0002"
        df = filter_dataframe_by_cols_name(df, column_list, keep_columns=True)
    
    def test_filter_dataframe_by_rows_name(self):
        families_matrix = {"CL0001": {"SP1": 3, 
                                      "SP2": 1, 
                                      "SP3": 3},
                            "CL0002": {"SP1": 3, 
                                       "SP2": 2, 
                                       "SP3": 2, 
                                       "SP4": 2, 
                                       "SP5": 2},
                            "CL0003": {"SP1": 1, 
                                       "SP2": 2}}
        df = convert_into_dataframe(families_matrix)
        row_list = ["SP1", "SP2"]
        df = filter_dataframe_by_rows_name(df, row_list, keep_row=True)