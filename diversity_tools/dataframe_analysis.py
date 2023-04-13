import argparse

from pathlib import Path
from sys import argv

from src.matrix_readers import read_orthovenn2_composition_output
from src.matrix_operations import (convert_list_to_numbers, convert_into_dataframe, get_dataframe_with_limited_families,
                                   filter_dataframe_cols_by_value_occurrence, select_families_with_highest_number_of_genes,
                                   filter_dataframe_by_cols_name, filter_dataframe_by_rows_name, calculate_shannon_diversity_index)
from src.matrix_writers import write_csv_from_matrix #, write_file_from_matrix
from src.plots import convert_diversity_matrix_to_graph

def parse_arguments():
    desc = "Analysis biological features from OrthoVenn2 (.txt file) DataFrame"
    parser = argparse.ArgumentParser(description=desc)

    help_Input_fpath = "Input file"
    parser.add_argument("--File" ,
                        "-f", type=str,
                        help=help_Input_fpath,
                        required=True)

    help_Operations = "Operations on the DataFrame"
    parser.add_argument("--Operations" ,
                        "-o", nargs='+',
                        help=help_Operations,
                        required=True)
    
    help_out_fdir = "Output folder"
    parser.add_argument("--out" ,
                        "-u", type=str,
                        help=help_out_fdir,
                        required=True)
    return parser

def get_options():
    parser = parse_arguments()  
    options = parser.parse_args()
    Input_fpath = Path(options.File)
    operations = options.Operations
    output_fdir = Path(options.out)
    return {"Input_fpath" : Input_fpath,
            "operations": operations,
            "out_fdir": output_fdir}

def main():
    options = get_options()
    operations = options["operations"]
    out_fdir = Path(options["out_fdir"])

    if not out_fdir.exists():
        out_fdir.mkdir()

    with open(Path(options["Input_fpath"])) as fhand:
        input_list = read_orthovenn2_composition_output(fhand)

    list_to_numbers = convert_list_to_numbers(input_list)
    df_matrix = convert_into_dataframe(list_to_numbers) 

    for operation in operations:
        operation_command, arguments_to_parse = operation.split("=")

        if operation_command == "limitation":
            arguments = {}
            argument, value = arguments_to_parse.split(",")
            arguments[argument] = value
            arguments["df"] = df_matrix
            df_matrix = get_dataframe_with_limited_families(**arguments)

        if operation_command == "filter":
            arguments = {}
            for argument in arguments_to_parse.split(";"):
                argument, value = argument.split(",")
                arguments[argument] = value
            arguments["df_matrix"] = df_matrix
            if "threshold" in arguments:
                arguments["threshold"] = float(arguments["threshold"])
            if "ignore_zeros" in arguments:
                if arguments["ignore_zeros"] == "False":
                    arguments["ignore_zeros"] = False
                elif arguments["ignore_zeros"] == "True":
                    arguments["ignore_zeros"] = True
                else:
                    raise ValueError("ignore zeros should be True or False")
            df_matrix = filter_dataframe_cols_by_value_occurrence(**arguments)

        if operation_command == "Highest_value":
            arguments = {}
            argument, value = arguments_to_parse.split(",")
            arguments[argument] = int(value)
            arguments["df_matrix"] = df_matrix
            df_matrix = select_families_with_highest_number_of_genes(**arguments)
            print(df_matrix)

        if operation_command == "Cols":
            arguments = {}
            for argument in arguments_to_parse.split(";"):
                argument, value = argument.split(',', 1)
                if argument == "column_name_list":
                    value = value.split(",")
                arguments[argument] = value
            arguments["df_matrix"] = df_matrix
            if "keep_columns" in arguments:
                if arguments["keep_columns"] == "False":
                    arguments["keep_columns"] = False
                elif arguments["keep_columns"] == "True":
                    arguments["keep_columns"] = True
                else:
                    raise ValueError("keep_columns should be True or False")
            df_matrix = filter_dataframe_by_cols_name(**arguments)

        if operation_command == "Rows":
            arguments = {}
            for argument in arguments_to_parse.split(";"):
                argument, value = argument.split(',', 1)
                if argument == "row_name_list":
                    value = value.split(",")
                arguments[argument] = value
            arguments["df_matrix"] = df_matrix
            if "keep_row" in arguments:
                if arguments["keep_row"] == "False":
                    arguments["keep_row"] = False
                elif arguments["keep_row"] == "True":
                    arguments["keep_row"] = True
                else:
                    raise ValueError("keep_row should be True or False")
            df_matrix = filter_dataframe_by_rows_name(**arguments)
    
    out_fpath = out_fdir / "Analyzed_DataFrame.csv"
    with open(out_fpath, "w") as out_fhand:
        write_csv_from_matrix(df_matrix, out_fhand)
    
    out_plot_fpath = out_fdir / "Families_Diversity_plot.svg"
    df_to_plot = calculate_shannon_diversity_index(df_matrix)
    convert_diversity_matrix_to_graph(df_to_plot, out_plot_fpath)

    run_log_fpath = out_fdir / "run.log.txt"
    with open(run_log_fpath, "w") as log_fhand:
        log_fhand.write("\t".join(argv)) 

    #Añadir al Log una explicación de las operaciones que se ha hecho sobre el data frame
    

if __name__ == "__main__":
    main()