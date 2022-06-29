"""This module contains common functions for the whole project"""
import os

import pandas as pd
import data_loading.load_xcms as load_xcms

def intersection(lst1, lst2):
    """Returns intersection of 2 lists
    Prtameters:
    lst1 (list): First list
    lst2 (list): Second list
    Returns:
    list
    """
    lst3 = [value for value in lst1 if value in lst2]
    return lst3


def ensure_dir_exists(path):
    """Creates path if it doesn't exist

        Parameters:
        path (string): The path to be create

       """
    if path.endswith('/'):
        path = path[0:-1]
    if not os.path.exists(path):
        os.makedirs(path)


def invert_dict(dic):
    """Returns inverted dictionary
    Parameters:
    dic (dict): Dictionary to be inverted
    Returns:
    dict
    """
    inverted_dict = dict()
    for key, value in dic.items():
        inverted_dict.setdefault(value, list()).append(key)
    return inverted_dict


def unite_repeated_columns(data_frame):
    res_df = pd.DataFrame()
    flower_only = load_xcms.get_flower_samples_xcms_df(data_frame)
    data_dict = {}
    time_intervals = set([x.split('_')[1] for x in flower_only.columns])
    not_flower_types = ['stem', 'blank']
    for interval in time_intervals:
        interval_headers = [item for item in flower_only.columns if (interval in item.replace(' ', ''))]
        interval_headers.sort()
        data_dict[interval] = interval_headers
        relevant_columns = flower_only[interval_headers]
        mean_column = relevant_columns.mean(axis=1)
        new_header = interval
        res_df[new_header] = mean_column
        for type in not_flower_types:
            relevant_columns = [col for col in data_frame.columns if (type in col and interval in col)]
            if relevant_columns:
                mean_column = data_frame[relevant_columns].mean(axis=1)
                res_df[interval+'-'+type] = mean_column
    return res_df
