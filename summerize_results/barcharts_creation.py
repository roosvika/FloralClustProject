import data_loading
from data_loading.load_xcms import get_clumn_interval
from floralclust_utils.utils import unite_repeated_columns, ensure_dir_exists
import matplotlib.pyplot as plt
import numpy as np


def get_new_column_name(header):
    time_interval = get_clumn_interval(header)
    if 'blank' in header.lower():
        return time_interval + '-blank'
    return time_interval


def autolabel(rects, ax):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        ax.annotate("{:.2f}".format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')


def create_barchart_from_partition(data, dest, cluster_id,
                                   part_data):
    """

    Args:
        data:
        dest:
        cluster_id:
        partition_title:
        partition_path:
        cluster_id_to_metabolites_mapping:

    Returns:

    """
    time_intervals = data_loading.load_xcms.get_time_intervals(data)
    data = unite_repeated_columns(data)
    cluster_id_to_metabolites_mapping = part_data.cluster_id_to_metabolites_mapping
    cluster_metabolites_list = cluster_id_to_metabolites_mapping[cluster_id]
    metabolites_volume_df = data.loc[cluster_metabolites_list].sum()
    headers = metabolites_volume_df.axes[0].values
    data_dict = {}
    for interval in time_intervals:
        interval_headers = [item for item in headers if (interval in item.replace(' ', ''))]
        interval_headers.sort()
        data_dict[interval] = interval_headers
    blank_volumes = []
    flower_volumes = []
    for interval in time_intervals:
        blank_header = [header for header in data_dict[interval] if ('blank' in header or 'control' in header)][0]
        blank_volumes.append(metabolites_volume_df[blank_header])
        flower_header = [header for header in data_dict[interval] if not (header in [blank_header])][0]
        flower_volumes.append(metabolites_volume_df[flower_header])
    y_pos = np.arange(len(time_intervals))
    plt.clf()
    x = np.arange(len(time_intervals))  # the label locations
    width = 0.9  # the width of the bars
    fig, ax = plt.subplots()
    rects1 = ax.bar(x - width / 3, flower_volumes, color='b', width=width / 3, alpha=0.5, label='Flower Volume')
    rects2 = ax.bar(x, blank_volumes, color='r', width=width / 3, alpha=0.25, label='Blank Volume')
    ax.set_xticks(x)
    ax.set_xticklabels(time_intervals)
    ax.set_ylabel('Metabolites Volume')
    ax.set_xlabel('Time Slot')
    title = f'Intensity Of Cluster {cluster_id} During The Day'
    ax.set_title(title)
    ax.legend()
    autolabel(rects1, ax)
    autolabel(rects2, ax)
    # autolabel(rects3, ax)
    fig.set_figwidth(10)
    fig.set_figheight(8)
    fig.tight_layout()
    plt.savefig(dest + '/' + str(cluster_id) + '.png')


def create_flower_barcharts_from_partition(data, dest,
                                           part_data):
    """

    Args:
        data:
        dest:
        partition_title:
        partition_path:
        cluster_id_to_metabolites_mapping:

    Returns:

    """
    # new_headers = [get_new_column_name(header) for header in data]
    # data.columns = new_headers
    # data = unite_repeated_columns(data)
    ensure_dir_exists(dest)
    cluster_ids = part_data.cluster_id_to_metabolites_mapping.keys()
    for cluster_id in cluster_ids:
        create_barchart_from_partition(data, dest, cluster_id,
                                   part_data)
