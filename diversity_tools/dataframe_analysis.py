import argparse

from pathlib import Path
from sys import argv

from src.matrix_readers import read_orthovenn2_composition_output
from src.matrix_operations import (convert_list_to_numbers, convert_into_dataframe, get_dataframe_with_limited_families,
                                   filter_dataframe_cols_by_value_occurrence, select_families_with_highest_number_of_genes,
                                   filter_dataframe_by_cols_name, filter_dataframe_by_rows_name, calculate_shannon_diversity_index)
from src.matrix_writers import write_csv_from_matrix

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

    main_fdir = Path("C:/Users/aleja/OneDrive - UPV/Documentos/Trabajo_Fin_de_Grado/results_folder")
    dir_path = main_fdir /  out_fdir

    if not dir_path.exists():
        dir_path.mkdir()

    with open(Path(options["Input_fpath"])) as fhand:
        input_list = read_orthovenn2_composition_output(fhand)

    list_to_numbers = convert_list_to_numbers(input_list)
    df_matrix = convert_into_dataframe(list_to_numbers)
    message = ""

    for operation in operations:
        operation_command, arguments_to_parse = operation.split("=")

        if operation_command == "limitation":
            arguments = {}
            argument, value = arguments_to_parse.split(",")
            arguments[argument] = value
            arguments["df"] = df_matrix
            message += "\n" + f"The DataFrame has a original length of {df_matrix.shape[1]} columns and"
            df_matrix = get_dataframe_with_limited_families(**arguments)
            message += f" has been reduced to a length of {df_matrix.shape[1]} columns"

        if operation_command == "filter":
            arguments = {}
            for argument in arguments_to_parse.split(";"):
                argument, value = argument.split(",")
                arguments[argument] = value
            arguments["df_matrix"] = df_matrix
            if "value" in arguments:
                message_value = arguments["value"]
                message += "\n" + f"The DataFrame has been applied a filter where {message_value} is the value to filter on,"
            if "threshold" in arguments:
                arguments["threshold"] = float(arguments["threshold"])
                message_threshold = arguments["threshold"]
                message += f" with a threshold cutoff of {message_threshold}/1." + "\n"
            if "mode" in arguments:
                message_mode = arguments["mode"]
                message += f"The mode applied is: {message_mode}, which determines whether to retain uppercase (greater_than), lowercase (lower_than), or same (equal) types of values, "
            if "ignore_zeros" in arguments:
                if arguments["ignore_zeros"] == "False":
                    arguments["ignore_zeros"] = False
                    message += f"without ingoring 0."
                elif arguments["ignore_zeros"] == "True":
                    arguments["ignore_zeros"] = True
                    message += f"ingoring 0."
                else:
                    raise ValueError("ignore zeros should be True or False")
            df_matrix = filter_dataframe_cols_by_value_occurrence(**arguments)
            
        if operation_command == "Highest_value":
            arguments = {}
            argument, value = arguments_to_parse.split(",")
            arguments[argument] = int(value)
            arguments["df_matrix"] = df_matrix
            df_matrix = select_families_with_highest_number_of_genes(**arguments)
            message += "\n" + f"Has been selected the {value} number of families with the highest number of genes to keep in the DataFrame."

        if operation_command == "Cols":
            arguments = {}
            for argument in arguments_to_parse.split(";"):
                argument, value = argument.split(',', 1)
                if argument == "column_name_list":
                    value = value.split(",")
                arguments[argument] = value
            arguments["df_matrix"] = df_matrix
            if "column_name_list" in arguments:
                message_cols_list = arguments["column_name_list"]
                message += "\n" + f"{message_cols_list} is a list of column names to "
            if "keep_columns" in arguments:
                if arguments["keep_columns"] == "False":
                    arguments["keep_columns"] = False
                    message += f"exclude  from DataFrame"
                elif arguments["keep_columns"] == "True":
                    arguments["keep_columns"] = True
                    message += f"keep from DataFrame"
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
            if "row_name_list" in arguments:
                message_list = arguments["row_name_list"]
                message += "\n" + f"{message_list} is a list of rows names to "
            if "keep_row" in arguments:
                if arguments["keep_row"] == "False":
                    arguments["keep_row"] = False
                    message += f"exclude from DataFrame"
                elif arguments["keep_row"] == "True":
                    arguments["keep_row"] = True
                    message += f"keep from DataFrame"
                else:
                    raise ValueError("keep_row should be True or False")
            df_matrix = filter_dataframe_by_rows_name(**arguments)
    
    out_fpath = dir_path / "Analyzed_DataFrame.csv"
    with open(out_fpath, "w") as out_fhand:
        write_csv_from_matrix(df_matrix, out_fhand)

    run_log_fpath = dir_path / "run.log.txt"
    with open(run_log_fpath, "w") as log_fhand:
        log_fhand.write("\t".join(argv))
        log_fhand.write(message)
    print(df_matrix)
if __name__ == "__main__":
    main()