import ast
import unittest
from pathlib import Path

import pandas as pd
from pandas.testing import assert_frame_equal

from src.matrix_operations import filter_df_by_domain

class FilterDom(unittest.TestCase):

    def setUp(self):
        test_path = Path(__file__).parent.absolute()
        test_path = test_path / "test_data"
        self.test_path = test_path / "test_filter_df_by_domain.csv"

    def test_filter_dom(self):
        with open(self.test_path) as input_fhand:
            input_fhand = pd.read_csv(input_fhand, converters={"domains": ast.literal_eval})
            filtered_df = filter_df_by_domain(input_fhand, [], [], [{"RT": "LINE"}])

        test_df = pd.DataFrame([{"sw": "452", "per div": "30.6", "per del": "1.4",
              "per ins": "1.4", "seqid": "Peame105C00", "start": "10027969",
              "end": "10028180", "q left": "45826121", "match": "+",
              "repeat": "rnd-5_family-987", "class": "LINE",
              "superfamily": "L1", "r start": "659", "r end": "870",
              "r left": "467", "id": "7765", "domains": [{"RT": "LINE"}],
              "tes order": "LINE", "tes superfamily": "unknown",
              "complete": "unknown", "strand": "+", "clade": "LINE", "length": 211}])
        
        convert_dict = {
        "sw": "int32", "per div": "float32", "per del": "float32",
        "per ins": "float32","start": "int64", "end": "int64",
        "q left": "int64", "r start": "int32", "r end": "int32",
        "r left": "int32", "id": "int32", "length": "int32"
        }
        test_df = test_df.astype(convert_dict)
        filtered_df = filtered_df.astype(convert_dict)

        assert_frame_equal(filtered_df.reset_index(drop=True), test_df.reset_index(drop=True)) 
        
if __name__ == "__main__":
    unittest.main()