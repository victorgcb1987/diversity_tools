import argparse

from pathlib import Path

from src.matrix_readers import read_orthovenn2_composition_output, read_matrix_from_file
from src.matrix_operations import convert_list_to_numbers, convert_into_dataframe, calculate_shannon_diversity_index
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
    
    help_file_format = "Input file format"
    parser.add_argument("--Format" ,
                        "-f", type=str,
                        help=help_file_format,
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
    Input_fpath = Path(options.iFile)
    fileFormat = options.Format
    output_fdir = Path(options.out)
    return {"Input_fpath" : Input_fpath,
            "fileFormat": fileFormat,
            "out_fdir": output_fdir}

def main():
    options = get_options()
    format = options["fileFormat"]
    out_fdir = options["out_fdir"]

    main_fdir = Path("C:/Users/aleja/OneDrive - UPV/Documentos/Trabajo_Fin_de_Grado/results_folder")
    dir_path = main_fdir /  out_fdir

    if not dir_path.exists():
        dir_path.mkdir()

    if format == "ov":
        with open(Path(options["Input_fpath"])) as fhand:
            file_output = read_orthovenn2_composition_output(fhand)
            List_to_numbers = convert_list_to_numbers(file_output)

    if format == "csv":
        with open(Path(options["Input_fpath"])) as fhand:
            List_to_numbers = read_matrix_from_file(fhand, family_field="#ID")

    df_matrix = convert_into_dataframe(List_to_numbers)
    df_matrix = df_matrix.T
    Diversity_df = calculate_shannon_diversity_index(df_matrix)

    out_plot_fpath = dir_path / "Diversity_DataFrame.svg"
    convert_diversity_matrix_to_graph(Diversity_df, out_plot_fpath)

    out_fpath = dir_path / "Diversity_DataFrame.csv"
    with open(out_fpath, "w") as out_fhand:
        write_csv_from_matrix(Diversity_df, out_fhand)


if __name__ == "__main__":
    main()