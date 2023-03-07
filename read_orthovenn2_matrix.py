import argparse

def parse_arguments():
    desc = "Get gene families with txt file"
    parser = argparse.ArgumentParser(description=desc)

    help_input_fpath = "input fpath"
    parser.add_argument("--input_fpath" ,
                        "-i", type=str,
                        help=help_input_fpath,
                        required=True)
    
    #Meter el test en options ademas de escribir 
    help_test_program = "Test program and quit"
    parser.add_argument("--test", 
                        "-t", action = "store_true", 
                        default=False, 
                        help=help_test_program)
    return parser

def get_options():
    parser = parse_arguments()  
    options = parser.parse_args()
    input_fpath = options.input_fpath
    test = options.test
    return {"input_fpath" : input_fpath, "test" : test}

def get_gene_families_with_selected_species(orthovenn2_input_fhand):
    family_matrix = {}
    for line in orthovenn2_input_fhand:
        line = line.rstrip()
        line = line.split()
        family_name = line[0]
        species_gene_names = line[1]
        for species_gene_name in species_gene_names.split(','):
            species_name, gene_name = species_gene_name.split('|')
            if family_name not in family_matrix:
                family_matrix[family_name] = {}
            if species_name not in family_matrix[family_name]:
                family_matrix[family_name][species_name] = [gene_name]
            else:
                family_matrix[family_name][species_name].append(gene_name)
    return family_matrix

def test(orthovenn2_input_test_fhand):
    family_matrix_test = get_gene_families_with_selected_species(orthovenn2_input_test_fhand)
    assert "Niund100Scf103505g0000010.1" in family_matrix_test["CL00001"]["Nicotiana_undulata"]
    assert "NipanScf047676g0002.1" in family_matrix_test["CL00002"]["Nicotiana_paniculata"]
    assert "Nideb015S394759g0000010.1" in family_matrix_test["CL00001"]["Nicotiana_debneyi"]
    assert "Niafr015S145251g004.1" in family_matrix_test["CL00001"]["Nicotiana_africana"]
    assert "NipanScf149893g0003.1" not in family_matrix_test["CL00001"]["Nicotiana_undulata"]
    assert "Niund100Scf057563g0000010.1" not in family_matrix_test["CL00002"]["Nicotiana_paniculata"]
    assert len(family_matrix_test["CL00001"]["Nicotiana_undulata"]) == 420
    

def main():
    options = get_options()
    if test:
        input_fpath_test = options["input_fpath"]
        with open(input_fpath_test) as input_fhand:
            test(input_fhand)
    else:
        input_fpath = options["input_fpath"]
        with open(input_fpath) as input_fhand:
            family_matrix = get_gene_families_with_selected_species(input_fhand)
        print(family_matrix)

if __name__ == "__main__":
    main()