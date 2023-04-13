import pandas as pd

def write_file_from_matrix(families_matrix, fpath):
    df = pd.DataFrame(families_matrix)
    df = df.fillna(value=0)
    df = df.astype(int)
    df = df.T
    df.to_csv(fpath, index_label="#ID")

def write_csv_from_matrix(families_matrix, fpath):
    families_matrix.to_csv(fpath, index_label="#ID")