import unittest
import pandas as pd

from pathlib import Path

from src.matrix_operations import convert_list_to_numbers, convert_into_dataframe, calculate_shannon_diversity_index, convert_into_limited_dataframe, select_families_with_x_genes_from_dataframe                               

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

    def test_matrix_limited_dataframe(self):
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
        df = convert_into_limited_dataframe(families_matrix, row_numbers= 2)

    
    def test_calculate_shannon_index(self):
        example_matrix = {"Community1" : {"black": 12, "purple": 21, "striped": 5,
                                          "green": 25, "brown": 2, "lblue": 17,
                                          "sblue": 9}, 
                          "Community2" : {"black": 10, "purple": 0, "striped": 15,
                                          "green": 15, "brown": 2, "lblue": 30,
                                          "sblue": 0}}
        #Dos comunidades mas con diferentes diversidad
        df = pd.DataFrame(example_matrix)
        example_diversity_df = calculate_shannon_diversity_index(df)
        assert round(example_diversity_df['Community1'], 3) == 1.746
        assert round(example_diversity_df['Community2'], 3) == 1.392

    def test_select_families_with_x_genes_from_dataframe(self):
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
        df = select_families_with_x_genes_from_dataframe(families_matrix, gene_number = 1)