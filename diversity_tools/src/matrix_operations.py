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

def calculate_shannon_diversity_index(df_matrix):
    cols = list(df_matrix.columns)
    df_frequencies = pd.DataFrame(df_matrix[cols].div(df_matrix[cols].sum(axis=0), axis=1))
    df_diversity = pd.DataFrame(df_frequencies[cols].transform(lambda x: -(x*np.log(x))))
    df_shannon = df_diversity.sum(axis=0)
    return df_shannon
    