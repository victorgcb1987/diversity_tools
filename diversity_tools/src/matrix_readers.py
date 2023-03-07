def read_orthovenn2_composition_output(orthovenn2_input_fhand):
    orthovenn2_matrix = {}
    for line in orthovenn2_input_fhand:
        line = line.rstrip()
        line = line.split()
        family_name = line[0]
        species_gene_names = line[1]
        for species_gene_name in species_gene_names.split(','):
            species_name, gene_name = species_gene_name.split('|')
            if family_name not in orthovenn2_matrix:
                orthovenn2_matrix[family_name] = {}
            if species_name not in orthovenn2_matrix[family_name]:
                orthovenn2_matrix[family_name][species_name] = [gene_name]
            else:
                orthovenn2_matrix[family_name][species_name].append(gene_name)
    return orthovenn2_matrix