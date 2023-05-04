import argparse

from pathlib import Path

from src.matrix_readers import read_orthovenn2_composition_output, read_matrix_from_file
from src.matrix_operations import convert_list_to_numbers, convert_into_dataframe, calculate_shannon_diversity_index
from src.plots import convert_diversity_matrix_to_graph, convert_specialization_matrix_to_graph

def parse_arguments():
    desc = "Extract biological features from 'ov' -OrthoVenn2- or 'csv' file and generate diversity plot"
    parser = argparse.ArgumentParser(description=desc)

    help_Input_fpath = "Input file"
    parser.add_argument("--iFile" ,
                        "-i", type=str,
                        help=help_Input_fpath,
                        required=True)
    
    help_Operations = "Operations on the DataFrame: <diversity> / <specialization> / <both>"
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
            dataframe = read_matrix_from_file(fhand, family_field="#ID")
            df_matrix = convert_into_dataframe(dataframe)
            df_matrix = df_matrix.T

    for operation in operations:
        if operation == "diversity":
            out_plot_fpath = dir_path / "Diversity_Analyzed_DataFrame.svg"
            convert_diversity_matrix_to_graph(df_matrix, out_plot_fpath)
            message += "\n" + f"Has been applied the shannon diversity index to the dataframe"

        if operation == "specialization":
            out_plot_fpath = dir_path / "Specificity_Specialization_Analyzed_DataFrameDataFrame.svg"
            convert_specialization_matrix_to_graph(df_matrix, out_plot_fpath)
            message += "\n" + f"Has been applied the shannon diversity index to the dataframe"


if __name__ == "__main__":
    main()