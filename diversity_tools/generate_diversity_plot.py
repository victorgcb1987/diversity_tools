import argparse

from pathlib import Path

from src.matrix_readers import read_orthovenn2_composition_output, read_matrix_from_file
from src.matrix_operations import convert_list_to_numbers, convert_into_dataframe, calculate_shannon_diversity_index
from src.matrix_writers import write_file_from_matrix
from src.plots import convert_diversity_matrix_to_graph

def parse_arguments():
    desc = "Extract biological features from OrthoVenn2 (.txt file) and generate diversity plot"
    parser = argparse.ArgumentParser(description=desc)

    help_OrthoVenn2_fpath = "OrthoVenn2 file"
    parser.add_argument("--ovFile" ,
                        "-f", type=str,
                        help=help_OrthoVenn2_fpath,
                        required=True)
    
    help_Final_Name = "Name of svg File to Generate"
    parser.add_argument("--svgName" ,
                        "-s", type=str,
                        help=help_Final_Name,
                        required=True)
    
    help_csv_Name = "Name of csv File to Generate"
    parser.add_argument("--csvName" ,
                        "-c", type=str,
                        help=help_csv_Name,
                        required=True)
    
    help_out_fdir = "Output folder"
    parser.add_argument("--out" ,
                        "-o", type=str,
                        help=help_out_fdir,
                        required=True)
    return parser

def get_options():
    parser = parse_arguments()  
    options = parser.parse_args()
    OrthoVenn2_fpath = Path(options.ovFile)
    svgName = options.svgName
    csvName = options.csvName
    output_fpath = Path(options.out)
    return {"OrthoVenn2_fpath" : OrthoVenn2_fpath,
            "svgName": svgName,
            "csvName": csvName,
            "out_fpath": output_fpath}

def main():
    options = get_options()
    out_fpath = options["out_fpath"]
    if not out_fpath.exists():
        out_fpath.mkdir()
    orthoVenn_output = read_orthovenn2_composition_output(options["OrthoVenn2_fpath"])
    List_to_numbers = convert_list_to_numbers(orthoVenn_output)
    generate_csvFile = write_file_from_matrix(List_to_numbers, options["csvName"])
    read_csvFile = read_matrix_from_file(generate_csvFile, family_field="CLNAME")
    dataframe = convert_into_dataframe(read_csvFile)
    Diversity_df = calculate_shannon_diversity_index(dataframe)
    convert_diversity_matrix_to_graph(Diversity_df, options["svgName"])
    plot_fpath = out_fpath / "diversity.svg"

if __name__ == "__main__":
    main()