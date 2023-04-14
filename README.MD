diversity_tools --description in progress--

"""
dataframe_analysis.py
Analysis biological features from OrthoVenn2 (.txt file) DataFrame.

This script reads an input file in OrthoVenn2 (.txt file) format, performs various operations on it
such as limiting the number of families, filtering columns by value occurrence, selecting families with the highest 
number of genes, and filtering rows by their names, and finally writes the processed DataFrame to a CSV file.

Usage:
    python main.py --File <input_file> --Operations <operations> --out <output_folder>

Arguments:
    --File (-f): Path to the input file in OrthoVenn2 (.txt file) format.
    --Operations (-o): Operations to perform on the DataFrame. Multiple operations should be separated by spaces.
        - limitation=<value,n>, where n is the maximum number of families to keep in the DataFrame.
        - filter=<value,v;ignore_zeros,y;threshold,t;mode,m>, where value is the value to filter on. Multiple filters should be separated by semicolons.
          Additionally, the following optional arguments can be used:
            - threshold=<threshold_value>: Specifies the minimum value (0-1) a column must have to be kept.
            - ignore_zeros=<True/False>: Specifies whether to ignore columns with all zero values.
            - mode=<greater_than/equal/less_than>: Specifies whether to compare the value with >= values, == values or <= values.
        - Highest_value=<n>, where n is the number of families with the highest number of genes to keep in the DataFrame.
        - Cols=<column_name_list,[col_name_1,col_name_2,...];keep_columns=True/False>, where column_name_list is a list of column names to keep,
          and keep_columns specifies whether to keep the specified columns or remove them.
        - Rows=<row_name_list,[row_name_1,row_name_2,...];keep_row=True/False>, where row_name_list is a list of row names to keep, and keep_row specifies whether to keep the specified rows or remove them.
    --out (-u): Path to the output folder.

Output:
    - Analyzed_DataFrame.csv: Processed DataFrame in CSV format.
    - run.log.txt: Log file containing the command line arguments used to run the script.

Examples:
    python main.py --File input.txt --Operations limitation=value,5 filter=value,1;threshold,0.5;ignore_zeros,False;mode,equal --out output_folder
    
    python main.py --File input.txt --Operations limitation=value,3 Cols=column_name_list,[COG0532,COG0052];keep_columns=True Rows=row_name_list,[Gene_1,Gene_2];keep_row,True --out output_folder
""" 