import argparse

from pathlib import Path
from sys import argv

from src.matrix_readers import read_matrix_from_file
from src.matrix_operations import convert_into_dataframe
from src.plots import convert_diversity_matrix_to_graph, convert_specialization_matrix_to_graph, convert_dataframe_to_scatter

def parse_arguments():
    desc = "Extract biological features from 'ov' -OrthoVenn2- or 'csv' file and generate diversity plot"
    parser = argparse.ArgumentParser(description=desc)

    help_Input_fpath = "Input file"
    parser.add_argument("--input_file" ,
                        "-i", type=str,
                        help=help_Input_fpath,
                        required=True)
    
    help_Input_fpath = "Input file2"
    parser.add_argument("--input_file2" ,
                        "-i2", type=str,
                        help=help_Input_fpath,
                        required=False, default=False)
    
    help_Operations = "Operations on the DataFrame: <diversity> / <specialization> / <both>"
    parser.add_argument("--operations" ,
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
    input_fpath = Path(options.input_file)
    if options.input_file2:
        input_fpath2 = Path(options.input_file2)
    else:
        input_fpath2 = False
    operations = options.operations
    output_fdir = Path(options.out)
    return {"input_fpath" : input_fpath,
            "input_fpath2": input_fpath2,
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

    with open(Path(options["input_fpath"])) as fhand:
            dataframe = read_matrix_from_file(fhand, family_field="#ID")
            df_matrix = convert_into_dataframe(dataframe)
            df_matrix = df_matrix.T
            if options["input_fpath2"]:
                with open(Path(options["input_fpath2"])) as fhand2:
                    dataframe2 = read_matrix_from_file(fhand2, family_field="#ID")
                    df_matrix2 = convert_into_dataframe(dataframe2)
                    df_matrix2 = df_matrix2.T

    message = ""

    for operation in operations:
        if operation == "diversity":
            out_plot_fpath = dir_path / "Diversity_Analyzed_DataFrame.svg"
            convert_diversity_matrix_to_graph(df_matrix, out_plot_fpath)
            message += "\n" + f"Has been applied a Barplot conversion from the shannon diversity index."

        if operation == "specialization":
            out_plot_fpath = dir_path / "Specialization_Analyzed_DataFrame.svg"
            convert_specialization_matrix_to_graph(df_matrix, out_plot_fpath)
            message += "\n" + f"Has been applied a Barplot conversion from the shannon specialization index."

        if operation == "scatter":
            out_plot_fpath = dir_path / "Diversity_Specialization_DataFrame.svg"
            convert_dataframe_to_scatter(df_matrix, df_matrix2, out_plot_fpath)
            message += "\n" + f"Has been applied a scatter_plot conversion of the shannon diversity and specialization index."


    run_log_fpath = dir_path / "run.log.txt"
    with open(run_log_fpath, "w") as log_fhand:
        log_fhand.write("\t".join(argv))
        log_fhand.write(message)

if __name__ == "__main__":
    main()