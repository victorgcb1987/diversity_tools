import argparse

from pathlib import Path
from sys import argv

from src.matrix_readers import read_matrix_from_file
from src.matrix_operations import (convert_into_dataframe, calculate_shannon_diversity_index, 
                                   calculate_shannon_specificity_index, calculate_dataframe_frecuencies)
from src.matrix_writers import write_csv_from_matrix
from src.plots import convert_diversity_matrix_to_graph


def parse_arguments():
    desc = "Extract biological features from 'ov' -OrthoVenn2- or 'csv' file and generate diversity plot"
    parser = argparse.ArgumentParser(description=desc)

    help_Input_fpath = "Input file"
    parser.add_argument("--iFile" ,
                        "-i", type=str,
                        help=help_Input_fpath,
                        required=True)
    
    help_Operations = "Operations on the DataFrame: <diversity> or <specificity>"
    parser.add_argument("--Operations" ,
                        "-o", nargs='+',
                        help=help_Operations,
                        required=True)
    
    help_out_fdir = "Output folder"
    parser.add_argument("--out" ,
                        "-f", type=str,
                        help=help_out_fdir,
                        required=True)
    return parser


def get_options():
    parser = parse_arguments()  
    options = parser.parse_args()
    Input_fpath = Path(options.iFile)
    operations = options.Operations
    output_fdir = Path(options.out)
    return {"Input_fpath" : Input_fpath,
            "operations": operations,
            "out_fdir": output_fdir}


def main():
    options = get_options()
    operations = options["operations"]
    out_fdir = options["out_fdir"]

    main_fdir = Path("C:/Users/aleja/OneDrive - UPV/Documentos/Trabajo_Fin_de_Grado/results_folder")
    dir_path = main_fdir /  out_fdir

    if not dir_path.exists():
        dir_path.mkdir()

    with open(Path(options["Input_fpath"])) as fhand:
            List_to_numbers = read_matrix_from_file(fhand, family_field="#ID")

    df_matrix = convert_into_dataframe(List_to_numbers)
    df_matrix = df_matrix.T
    print(df_matrix)
    message = ""

    for operation in operations:
        if operation == "diversity":
            frecuency_df = calculate_dataframe_frecuencies(df_matrix)
            Diversity_df = calculate_shannon_diversity_index(frecuency_df)
            message += f"Has been applied the shannon diversity index to the dataframe"

        if operation == "specificity":
            frecuency_df = calculate_dataframe_frecuencies(df_matrix)
            Diversity_df = calculate_shannon_specificity_index(frecuency_df)
            message += f"Has been applied the shannon diversity index to the dataframe"

    out_fpath = dir_path / "Shannon_index_DataFrame.csv"
    with open(out_fpath, "w") as out_fhand:
        write_csv_from_matrix(Diversity_df, out_fhand)

    run_log_fpath = dir_path / "run.log.txt"
    with open(run_log_fpath, "w") as log_fhand:
        log_fhand.write("\t".join(argv))
        log_fhand.write(message)


if __name__ == "__main__":
    main()