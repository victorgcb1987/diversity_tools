def convert_list_to_numbers(gene_families):
    gene_families_counts = {}
    for family_names, species in gene_families.items():
            for species_names, gene_names in species.items():
                if family_names not in gene_families_counts:
                    gene_families_counts[family_names] = {}
                if species_names not in gene_families_counts[family_names]:
                    gene_families_counts[family_names][species_names] = len(gene_names)
    return gene_families_counts


def calculate_shannon_diversity_index(gene_families):
    g_sum = df.groupby('Name_Receive')['Amount'].transform('sum')
    values = df['Amount']/g_sum
    df['Entropy'] = -(values*np.log(values))

    df1 = df.groupby('Name_Receive',as_index=False,sort=False)['Entropy'].sum()