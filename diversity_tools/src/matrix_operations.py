import pandas as pd
import numpy as np

def convert_list_to_numbers(gene_families):
    gene_families_counts = {}
    for family_names, species in gene_families.items():
            for species_names, gene_names in species.items():
                if family_names not in gene_families_counts:
                    gene_families_counts[family_names] = {}
                if species_names not in gene_families_counts[family_names]:
                    gene_families_counts[family_names][species_names] = len(gene_names)
    return gene_families_counts

def convert_into_dataframe(gene_families):
    df = pd.DataFrame(gene_families)
    df = df.fillna(value=0)
    df = df.astype(float)
    #df = df.T
    return df

def get_dataframe_with_limited_families(df, value = 3):
    value = int(value)
    limited_df = df.iloc[:, :value]
    return limited_df

def calculate_dataframe_frecuencies_col(df_matrix):
    cols = list(df_matrix.columns)
    df_frequencies = pd.DataFrame(df_matrix[cols].div(df_matrix[cols].sum(axis=0), axis=1))
    return df_frequencies

def calculate_dataframe_frecuencies_row(df_matrix):
    cols = list(df_matrix.columns)
    df_frequencies = pd.DataFrame(df_matrix[cols].div(df_matrix[cols].sum(axis=1), axis=0))
    return df_frequencies

def calculate_shannon_diversity_index(df_frecuency_matrix):
    cols = list(df_frecuency_matrix.columns)
    np.seterr(divide = 'ignore')
    df_diversity = pd.DataFrame(df_frecuency_matrix[cols].transform(lambda x: -(x*np.log2(x))))
    df_shannon = df_diversity.sum(axis=0)
    return df_shannon

def calculate_shannon_specificity_index(df_frecuency_matrix):
    cols = list(df_frecuency_matrix.columns)
    np.seterr(divide = 'ignore')
    t = df_frecuency_matrix.shape[1]
    df_average_frequency = (df_frecuency_matrix.sum(axis=1)) / t
    df_divide_cal = df_frecuency_matrix.div(df_average_frequency[0], axis='columns')
    df_family_specificity = pd.DataFrame(df_divide_cal[cols].transform(lambda x: (x)*np.log2(x), axis = 1))
    df_family_specificity = df_family_specificity.fillna(value=0).sum(axis=1)
    df_shannon = df_family_specificity / t
    return df_shannon

def calculate_shannon_specialization_index(df_frecuency_matrix, df_specificity_matrix):
    np.seterr(divide = 'ignore')
    df_shannon_specialization = df_frecuency_matrix.mul(df_specificity_matrix, axis = 0)
    df_shannon_specialization = df_shannon_specialization.sum(axis=0)   
    return df_shannon_specialization

def filter_dataframe_cols_by_value_occurrence(df_matrix, value=1, ignore_zeros=False, threshold=1, mode = "equal"):
    df = df_matrix
    value = int(value)
    if ignore_zeros:
        threshold_cutoff = threshold * len(df != 0)
    else:
        threshold_cutoff = threshold * len(df)
    threshold_cutoff = round(threshold_cutoff, 1)
    if mode == "greater_than":
        mask = df.apply(lambda col: (col >= value).sum() > threshold_cutoff)    
        df = df.loc[:, mask]
    elif mode == "equal":
        mask = df.apply(lambda col: (col == value).sum() > threshold_cutoff)    
        df = df.loc[:, mask]
    elif mode == "less_than":
        mask = df.apply(lambda col: (col <= value).sum() < threshold_cutoff)    
        df = df.loc[:, mask]
    return df

def select_families_with_highest_number_of_genes(df_matrix, value):
    selected_df = df_matrix.apply(lambda x: x.nlargest(value, keep = 'all'))
    selected_df = selected_df.fillna(value=0)
    selected_df = selected_df.astype(int)
    return selected_df

def filter_dataframe_by_cols_name(df_matrix, column_name_list, keep_columns=True):
    df = df_matrix.fillna(value=0)
    df = df.astype(int)
    if keep_columns:
        df = df.loc[:, column_name_list]
    else:
        df = df.drop(column_name_list, axis = 1)
    return df

def filter_dataframe_by_rows_name(df_matrix, row_name_list, keep_row=True):
    df = df_matrix.fillna(value=0)
    df = df.astype(int)
    if keep_row:
        df = df.loc[row_name_list]
    else:
        df = df.drop(row_name_list, axis = 0)
    return df

def count_tes(input_df, species_name, col="superfamily"):
    """Counts each element of the selected column.

    Parameters
    ----------
    input_df : `pandas.DataFrame`
        Dataframe in which rows correspond to repeats and some
        of the columns consist of the available categories.

    species_name : str
        Name of the species to name the resulting Series.

    col : str, default: 'superfamily'
        Column that will be selected and counted.
        
    Returns
    -------
    counted_tes : `pandas.Series`
        Indexes are named after counted elements and the Series
        is named after the species it came from, so that it
        can be used as the name of the column when combining
        different Series.
    """

    if col == "domains":
        input_df["domains"] = input_df['domains'].apply(lambda x: ','.join(map(str, x)))

    counted_tes = input_df.value_counts(col).rename(species_name).astype("int32")
    return counted_tes

def create_te_count_matrix(list_of_inputs):
    """Combines the series from count_tes() into a dataframe.
    

    Parameters
    ----------
    list_of_inputs : list
        List containing all the `pandas.Series` to concatenate.

    Returns
    -------
    te_count_matrix : `pandas.DataFrame`
        Columns: name of the species; indexes: element. In case
        some element was not present in a species, it is filled
        with a 0. All numbers are converted into integers,
        as some Series can be added as floats.
    """
    te_count_matrix = pd.DataFrame()
    for input in list_of_inputs:
        te_count_matrix = pd.concat([te_count_matrix, input], axis=1)
    te_count_matrix = te_count_matrix.fillna(0).astype("int32")
    return te_count_matrix

def filter_df_by_domain(df_to_filter, doms, clades, special_features):
    """Filters a dataframe to remove repeats that do not contain data on their domains.
    
    It can also filter by specific domains, clades and domain:clade features.

    Parameters
    ----------
    df_to_filter : `pandas.DataFrame`
        Dataframe containing columns for clades (named clades)
        and domains (named domains).

    doms : list, optional
        List of domains to filter (include).

    clades : list, optional
        List of domains to filter (include).

    special_features : list of dictionaries, optional
        List of dictionaries in which each dictionary specifies
        a different domain:clade feature.

    Returns
    -------
    filtered_df : `pandas.DataFrame`
        Dataframe containing the specified data.
    """
    #Filtering is based on True/False series,
    #checking if there is at least one of the given elements in any row
    filtered_df = df_to_filter[df_to_filter.domains.apply(lambda x: x != [{"none": "none"}])]

    if doms:
        filtered_df = filtered_df[filtered_df.domains.apply(lambda x: any(any(feat.get(dom) for dom in doms) for feat in x))]

    if clades:
        filtered_df = filtered_df[filtered_df.clade.isin(clades)]

    if special_features:
        filtered_df = filtered_df[filtered_df.domains.apply(lambda x: any(any(feat.items() == dom.items() for feat in special_features) for dom in x))]

    return filtered_df

def filter_df_by_length(df_to_filter, length):
    """Filters the dataframe by eliminating rows whose length is lower than the given.

    Parameters
    ----------
    df_to_filter : `pandas.DataFrame`
        Dataframe containing columns for clades (named clades)
        and domains (named domains).

    length : int
        Minimum length to remain in the dataframe.

    Returns
    -------
    filtered_df : `pandas.DataFrame`
        Dataframe containing the specified data.    
    """
    filtered_df = df_to_filter[df_to_filter["length"] >= length]

    return filtered_df

def filter_df_by_percentages(df_to_filter, percentage="div=20.0", mode="lower_than"):
    """Filters the dataframe by percentages and mode

    Parameters
    ----------
    df_to_filter : `pandas.DataFrame`
        Dataframe containing columns for clades (named clades)
        and domains (named domains).

    percentage : str, default: 'div=20.0', optional
        Percentage to filter the dataframe. The number must
        be preceded by either 'div', 'del' or 'ins', and joined
        by a an equal sign (=).
    
    mode : str, default: 'lower_than', optional
        Mode to filter data to only include repeats 
        lower, higher or equal to the threshold
    
    Returns
    -------
    filtered_df : `pandas.DataFrame`
        Dataframe containing the specified data.    
    """
    name_and_perc = percentage.split("=")
    name = "per " + name_and_perc[0]
    perc = float(name_and_perc[1])
    names = ["per div", "per del", "per ins"]      
    modes = ["higher_than", "lower_than", "equal"]

    if name in names:
        if mode in modes:
            if mode == "higher_than":
                mode_filter = df_to_filter[name] >= perc
                filtered_df = df_to_filter[mode_filter]
            elif mode == "lower_than":
                mode_filter = df_to_filter[name] <= perc
                filtered_df = df_to_filter[mode_filter]
            elif mode == "equal":
                mode_filter = df_to_filter[name] == perc
                filtered_df = df_to_filter[mode_filter]

            return filtered_df
