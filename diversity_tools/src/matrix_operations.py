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

def get_dataframe_with_limited_families(df, value = 3):
    value = int(value)
    limited_df = df.iloc[:, :value]
    return limited_df

def calculate_shannon_diversity_index(df_matrix):
    cols = list(df_matrix.columns)
    df_frequencies = pd.DataFrame(df_matrix[cols].div(df_matrix[cols].sum(axis=0), axis=1))
    np.seterr(divide = 'ignore')
    df_diversity = pd.DataFrame(df_frequencies[cols].transform(lambda x: -(x*np.log(x))))
    df_shannon = df_diversity.sum(axis=0)
    return df_shannon

def calculate_shannon_specificity_index(df_matrix):
    cols = list(df_matrix.columns)
    #df_frequencies = pd.DataFrame(df_matrix[cols].div(df_matrix[cols].sum(axis=0), axis=1))
    np.seterr(divide = 'ignore')
    print(df_matrix[cols])
    t = df_matrix[cols].shape[0]
    print(t)
    # df_average_frequency = (1/t)*((df_frequencies.sum(axis=0)))
    # df_copy = df_frequencies[cols].divide(df_average_frequency)
    df_average_frequency = (1/t)*((df_matrix.sum(axis=0)))
    df_copy = df_matrix[cols].divide(df_average_frequency)
    df_family_specificity = pd.DataFrame(df_copy[cols].transform(lambda x: (x)*np.log(x)))
    df_family_specificity = df_family_specificity.fillna(value=0)
    df_family_specificity = df_family_specificity / t
    #df_family_specificity = pd.DataFrame(df_frequencies[cols].transform(lambda x: (1/t)*((x/df_average_frequency)*np.log(x/df_average_frequency))))
    df_shannon = df_family_specificity.sum(axis=0)
    print(df_shannon)
    return df_shannon

def filter_dataframe_cols_by_value_occurrence(df_matrix, value=1, ignore_zeros=False, threshold=1, mode = "equal"):
    df = df_matrix
    value = int(value)
    if ignore_zeros:
        threshold_cutoff = threshold * len(df != 0)
    else:
        threshold_cutoff = threshold * len(df)
    threshold_cutoff = round(threshold_cutoff, 1)
    if mode == "greater_than":
        mask = df.apply(lambda col: (col >= value).sum() > threshold_cutoff)    
        df = df.loc[:, mask]
    elif mode == "equal":
        mask = df.apply(lambda col: (col == value).sum() > threshold_cutoff)    
        df = df.loc[:, mask]
    elif mode == "less_than":
        mask = df.apply(lambda col: (col <= value).sum() < threshold_cutoff)    
        df = df.loc[:, mask]
    return df

def select_families_with_highest_number_of_genes(df_matrix, value):
    selected_df = df_matrix.apply(lambda x: x.nlargest(value, keep = 'all'))
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