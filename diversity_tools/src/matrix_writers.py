import pandas as pd

def write_file_from_matrix(families_matrix, fpath):
    df = pd.DataFrame(families_matrix)
    df = df.fillna(value=0)
    df.to_csv(fpath, float_format='%.2f')