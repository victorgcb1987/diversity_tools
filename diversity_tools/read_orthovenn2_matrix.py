import argparse

from src.matrix_readers import read_orthovenn2_composition_output

def parse_arguments():
    desc = "Get gene family matrix from orthovenn2 composition output"
    parser = argparse.ArgumentParser(description=desc)

    help_input_fpath = "input fpath"
    parser.add_argument("--input_fpath" ,
                        "-i", type=str,
                        help=help_input_fpath,
                        required=True)
    return parser

def get_options():
    parser = parse_arguments()  
    options = parser.parse_args()
    input_fpath = options.input_fpath
    return {"input_fpath" : input_fpath}


def main():
    options = get_options()
    input_fpath = options["input_fpath"]
    with open(input_fpath) as input_fhand:
        family_matrix = read_orthovenn2_composition_output(input_fhand)
        print(family_matrix)


if __name__ == "__main__":
    main()