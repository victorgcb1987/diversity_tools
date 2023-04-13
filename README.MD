diversity_tools --description in progress--

"""
dataframe_analysis.py
Analysis biological features from OrthoVenn2 (.txt file) DataFrame.

This script reads an input file in OrthoVenn2 (.txt file) format, performs various operations on it
such as limiting the number of families, filtering columns by value occurrence, selecting families with the highest 
number of genes, and filtering rows by their names, calculates the Shannon diversity index, and finally writes the
processed DataFrame to a CSV file and generates a plot.

Usage:
    python main.py --File <input_file> --Operations <operations> --out <output_folder>

Arguments:
    --File (-f): Path to the input file in OrthoVenn2 (.txt file) format.
    --Operations (-o): Operations to perform on the DataFrame. Multiple operations should be separated by spaces.
        - limitation=<n>, where n is the maximum number of families to keep in the DataFrame.
        - filter=<column_name_1,value_1;column_name_2,value_2;...>, where column_name_i is the name of the column to filter
          and value_i is the value to filter on. Multiple filters should be separated by semicolons.
          Additionally, the following optional arguments can be used:
            - threshold=<threshold_value>: Specifies the minimum value a column must have to be kept.
            - ignore_zeros=<True/False>: Specifies whether to ignore columns with all zero values.
        - Highest_value=<n>, where n is the number of families with the highest number of genes to keep in the DataFrame.
        - Cols=<column_name_list=[col_1,col_2,...];keep_columns=True/False>, where column_name_list is a list of column names to keep,
          and keep_columns specifies whether to keep the specified columns or remove them.
        - Rows=<row_name_list=[row_1,row_2,...];keep_row=True/False>, where row_name_list is a list of row names to keep, and 
          keep_row specifies whether to keep the specified rows or remove them.
    --out (-u): Path to the output folder.

Output:
    - Analyzed_DataFrame.csv: Processed DataFrame in CSV format.
    - Families_Diversity_plot.svg: Plot showing the Shannon diversity index of the families.
    - run.log.txt: Log file containing the command line arguments used to run the script.

Example:
    python main.py --File input.txt --Operations limitation=5 filter=Species,Escherichia coli;threshold=0.5 Highest_value=3 
    Cols=column_name_list=[COG0532,COG0052];keep_columns=True Rows=row_name_list=[Gene_1,Gene_2];keep_row=True 
    --out output_folder
""" 