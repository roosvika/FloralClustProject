import pandas as pd
import xlrd


def get_normilezed_xcms_df(data_path):
    """

    Args:
        data_path: Path to .xlsx xcms file

    Returns:
        DataFrame contaning only normalized samples columns. Each column
        normalized by mean devision.
    """
    data = pd.read_excel(data_path)
    data = data.sort_values(by=['featureidx'])
    data = data[[col for col in data.columns if 'sample' in col]]
    data = data / data.mean()
    return data

def get_xcms_df(data_path):
    """

    Args:
        data_path: Path to .xlsx xcms file

    Returns:
        DataFrame contaning only original samples columns (not normalized).
    """
    data = pd.read_excel(data_path)
    data = data.sort_values(by=['featureidx'])
    data = data[[col for col in data.columns if 'sample' in col]]
    return data

def get_flower_samples_xcms_df(xscm_df):
    """

    Args:
        xcms_df: DataFrame contaning only normalized samples columns (usually the output of get_normilezed_xcms_df)

    Returns:
        DataFrame containing only the flower samples intencities columns
    """
    headers = xscm_df.columns
    irrelevent_headers = [
        header for header in headers
        if ('blank' in header.lower() in header.lower() or 'stem' in header.lower())
    ]
    return xscm_df.copy().drop(irrelevent_headers, axis='columns')


def get_clumn_interval(col_name):
    return col_name.split('_')[1]


def get_time_intervals(samples_xcms_df):
    time_intervals = list(set([x.split('_')[1] for x in get_flower_samples_xcms_df(samples_xcms_df).columns]))
    time_intervals.sort(key=lambda a:int(a.split(':')[0]))
    return time_intervals
