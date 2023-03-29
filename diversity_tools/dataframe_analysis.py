import argparse

from pathlib import Path

from src.matrix_readers import read_orthovenn2_composition_output
from src.matrix_operations import (convert_list_to_numbers, convert_into_dataframe, get_dataframe_with_limited_families,
                                   filter_dataframe_cols_by_value_occurrence, select_families_with_highest_number_of_genes,
                                   filter_dataframe_by_cols_name, filter_dataframe_by_rows_name)

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
                        "-o", type=str,
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
    output_fpath = Path(options.out)
    return {"Input_fpath" : Input_fpath,
            "operations": operations,
            "out_fpath": output_fpath}

def main():
    options = get_options()
    operations = options["operations"]
    out_fpath = options["out_fpath"]

    if not out_fpath.exists():
        out_fpath.mkdir()
    input_list = read_orthovenn2_composition_output(options["Input_fpath"])
    list_to_numbers = convert_list_to_numbers(input_list)
    df_matrix = convert_into_dataframe(list_to_numbers)

    for operation in operations:
        operation_command, arguments_to_parse = operation.split("=")

        if operation_command == "limitation":
            arguments = {}
            for argument in arguments_to_parse:
                argument, value = argument.split(",")
                arguments[argument] = value
            arguments["df_matrix"] = df_matrix
            get_dataframe_with_limited_families(**arguments)

        if operation_command == "filter":
            arguments = {}
            for argument in arguments_to_parse.split(";"):
                argument, value = argument.split(",")
                arguments[argument] = value
            arguments["df_matrix"] = df_matrix
            filter_dataframe_cols_by_value_occurrence(**arguments)
    
        if operation_command == "Highest_value":
            arguments = {}
            for argument in arguments_to_parse:
                argument, value = argument.split(",")
                arguments[argument] = value
            arguments["df_matrix"] = df_matrix
            select_families_with_highest_number_of_genes(**arguments)
        
        if operation_command == "Cols":
            arguments = {}
            for argument in arguments_to_parse.split(";"):
                argument, value = argument.split(",")
                arguments[argument] = value
            arguments["df_matrix"] = df_matrix
            filter_dataframe_by_cols_name(**arguments)
        
        if operation_command == "Rows":
            arguments = {}
            for argument in arguments_to_parse.split(";"):
                argument, value = argument.split(",")
                arguments[argument] = value
            arguments["df_matrix"] = df_matrix
            filter_dataframe_by_rows_name(**arguments)
        

if __name__ == "__main__":
    main()






















# python scripy.py blao --operations filter=value,1;threshold,1;ignore_zeros,True;mode,equal  
#operations = options.operations



#main 
#operations = options["operations"]
# for operation in operations:
#     operation_command, arguments_to_parse = operation.split("=")
#     if operation_command == "filter":
#         arguments = {}
#         for argument in arguments_to_parse.split(";"):
#             argument, value = argument.split(",")
#             arguments[argument] = value
#         arguments["df_matrix"] = df_matrix
#         funcion(**arguments)
