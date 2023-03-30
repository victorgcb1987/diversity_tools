import pandas as pd
import numpy as np

def convert_list_to_numbers(gene_families):
    gene_families_counts = {}
    for family_names, species in gene_families.items():
            for species_names, gene_names in species.items():
                if family_names not in gene_families_counts:
                    gene_families_counts[family_names] = {}
                if species_names not in gene_families_counts[family_names]:
                    gene_families_counts[family_names][species_names] = len(gene_names)
    return gene_families_counts

def convert_into_dataframe(gene_families):
    df = pd.DataFrame(gene_families)
    df = df.fillna(value=0)
    df = df.astype(int)
    #df = df.T
    return df

def get_dataframe_with_limited_families(df, family_number = 3):
    family_number = int(family_number)
    limited_df = df.iloc[:, :family_number]
    limited_df = limited_df.fillna(value=0)
    limited_df = limited_df.astype(int)
    return limited_df

def calculate_shannon_diversity_index(df_matrix):
    cols = list(df_matrix.columns)
    df_frequencies = pd.DataFrame(df_matrix[cols].div(df_matrix[cols].sum(axis=0), axis=1))
    df_diversity = pd.DataFrame(df_frequencies[cols].transform(lambda x: -(x*np.log(x))))
    df_shannon = df_diversity.sum(axis=0)
    return df_shannon

def filter_dataframe_cols_by_value_occurrence(df_matrix, value=1, ignore_zeros=True, threshold=1, mode = "equal"):
    df = df_matrix.fillna(value=0)
    df = df.astype(int)
    if ignore_zeros:
        threshold_cutoff = threshold * (df != 0).sum()
    else:
        threshold_cutoff = threshold * df.sum()
    if mode == "greater_than":
        df = df.loc[:, (df >= value).sum() >= threshold_cutoff]
    elif mode == "equal":
        df = df.loc[:, (df == value).sum() >= threshold_cutoff]
    elif mode == "less_than":
        df = df.loc[:, (df <= value).sum() <= threshold_cutoff]
    return df

def select_families_with_highest_number_of_genes(df_matrix, num_values):
    selected_df = df_matrix.apply(lambda x: x.nlargest(num_values, keep = 'all'))
    selected_df = selected_df.fillna(value=0)
    selected_df = selected_df.astype(int)
    return selected_df

def filter_dataframe_by_cols_name(df_matrix, column_name_list, keep_columns=True):
    df = df_matrix.fillna(value=0)
    df = df.astype(int)
    if keep_columns:
        df = df.loc[:, column_name_list]
    else:
        df = df.drop(column_name_list, axis = 1)
    return df

def filter_dataframe_by_rows_name(df_matrix, row_name_list, keep_row=True):
    df = df_matrix.fillna(value=0)
    df = df.astype(int)
    if keep_row:
        df = df.loc[row_name_list]
    else:
        df = df.drop(row_name_list, axis = 0)
    return df