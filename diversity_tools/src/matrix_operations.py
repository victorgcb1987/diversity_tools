def convert_list_to_numbers(gene_families):
    gene_families_counts = {}
    for family_names, species in gene_families.items():
            for species_names, gene_names in species.items():
                if family_names not in gene_families_counts:
                    gene_families_counts[family_names] = {}
                if species_names not in gene_families_counts[family_names]:
                    gene_families_counts[family_names][species_names] = len(gene_names)
    return gene_families_counts