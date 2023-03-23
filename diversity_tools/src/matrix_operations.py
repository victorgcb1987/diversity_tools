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

def convert_into_limited_dataframe(gene_families, row_numbers = 10):
    df = pd.DataFrame(gene_families)
    limited_df = df.iloc[:row_numbers]
    limited_df = limited_df.fillna(value=0)
    limited_df = limited_df.astype(int)
    return limited_df

def calculate_shannon_diversity_index(df_matrix):
    cols = list(df_matrix.columns)
    df_frequencies = pd.DataFrame(df_matrix[cols].div(df_matrix[cols].sum(axis=0), axis=1))
    df_diversity = pd.DataFrame(df_frequencies[cols].transform(lambda x: -(x*np.log(x))))
    df_shannon = df_diversity.sum(axis=0)
    return df_shannon

def select_families_with_x_genes_from_dataframe(df_matrix, gene_number = 1):
     df = pd.DataFrame(df_matrix)
     #df = df.loc[df == gene_number]
     print(df)
     return df