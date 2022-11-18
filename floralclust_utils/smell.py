import numpy
import pandas as pd
from data_loading.load_xcms import get_time_intervals
from floralclust_utils.flower_enums import MetabolitesType
from floralclust_utils.utils import unite_repeated_columns


def calc_daily_pattern_value(volumes_list):
    """Evaluates the smelling factor of a flower by daily pattern.

    Args:
        volumes_list: array of volumes during the day

    Returns:
        A number indicating the smelling factor of the flower,
        0 - extremely smelling
        1 - smelling
        >=2 - not smelling
    """
    first_derivative_deltas_list = []
    for idx in range(len(volumes_list) - 1):
        first_derivative_deltas_list.append(volumes_list[idx + 1] - volumes_list[idx])
    smelling_value = 0
    for idx in range(len(first_derivative_deltas_list) - 1):
        smelling_value += int(max(numpy.sign(first_derivative_deltas_list[idx]), 0)) ^ \
                          int(max(numpy.sign(first_derivative_deltas_list[idx + 1]), 0))
    return smelling_value


def calc_flower_blank_avg_dist(flower_volumes, blank_volumes):
    """Returns a ration indicating the smelling factor of a flower according
    to difference between flower and blank during the day

    Args:
        flower_volumes: array of flower volumes during the day
        blank_volumes: array of blank volumes during the day

    Returns: A float number indicating the how much the cluster smells.
     0 - not smelling
     > 0 - smelling, the higher the result - the more smelling is the cluster

    """
    flower_blank_ratio = 0
    for volume_f, volume_b in list(zip(flower_volumes, blank_volumes)):
        volume_f += 0.00000001
        volume_b -= 0.00000001
        if volume_f < volume_b:
            return 0
        flower_blank_ratio += volume_b / volume_f
    flower_blank_ratio /= len(flower_volumes)
    return 1 - flower_blank_ratio


def get_cluster_data(data, cluster_id, partition_data, metabolites_volume_type):
    """

    Args:
        data:
        cluster_id:
        partition_data:
        metabolites_volume_type:

    Returns:

    """
    time_intervals = get_time_intervals(data)
    # data = data.copy()
    data = unite_repeated_columns(data)
    cluster_metabolites_list = partition_data.cluster_id_to_metabolites_mapping[cluster_id]
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
    if metabolites_volume_type == MetabolitesType.FLOWER:
        return flower_volumes
    if metabolites_volume_type == MetabolitesType.BLANK:
        return blank_volumes
    raise TypeError('metabolites_volume_type should be MetabolitesType.FLOWER or MetabolitesType.BLANK')


def is_attractor(data, cluster_id, partition_data):
    """Returns weather a certain cluster is a smelling cluster

    Args:
        data: A DataFrame of xcms data (blank and flower samples)
        cluster_id: integer of cluster id
        partition_data: partition data object

    Returns: Boolean - weather the cluster is a smelling

    """
    cluster_metabolites_list = partition_data.cluster_id_to_metabolites_mapping[cluster_id]
    metabolites_volume_df = data.loc[cluster_metabolites_list]
    flower_volumes = get_cluster_data(data, cluster_id, partition_data,
                                      MetabolitesType.FLOWER)
    if len(metabolites_volume_df)<2:
        return False
    blank_volumes = get_cluster_data(data, cluster_id, partition_data, MetabolitesType.BLANK)
    is_cluster_attractor = calc_flower_blank_avg_dist(flower_volumes, blank_volumes) > 0
    is_daily_pattern = calc_daily_pattern_value(flower_volumes) < 2
    return is_cluster_attractor and is_daily_pattern


def get_attractors_list(data, partition_data):
    attractors_list = []
    return [cluster_id for cluster_id in partition_data.cluster_id_to_metabolites_mapping.keys() if
            is_attractor(data, int(cluster_id), partition_data)]


def get_total_smelling_factor(daily_pattern_value, flower_blank_ratio):
    """

    Args:
        daily_pattern_value:
        flower_blank_ratio:

    Returns:

    """
    return (max(2 - daily_pattern_value, 0) + flower_blank_ratio) / 2


def get_clusters_smelling_factors(data, partition_data):
    """

    Args:
        data:
        partition_data:

    Returns:

    """
    clusters_ids = partition_data.cluster_id_to_metabolites_mapping.keys()
    attractors_smelling_factors_df = pd.DataFrame(
        columns=['cluster_id', 'daily_pattern_value', 'flower_blank_ratio', 'smelling_factor'])
    for cluster_id in clusters_ids:
        flower_volumes = get_cluster_data(data, cluster_id, partition_data,
                                          MetabolitesType.FLOWER)
        blank_volumes = get_cluster_data(data, cluster_id, partition_data,
                                         MetabolitesType.BLANK)
        daily_pattern_value = calc_daily_pattern_value(flower_volumes)
        flower_blank_ratio = calc_flower_blank_avg_dist(flower_volumes, blank_volumes)
        total_smelling_factor = get_total_smelling_factor(daily_pattern_value, flower_blank_ratio)
        attractors_smelling_factors_df = attractors_smelling_factors_df.append(
            {'cluster_id': cluster_id, 'daily_pattern_value': daily_pattern_value,
             'flower_blank_ratio': flower_blank_ratio,
             'smelling_factor': total_smelling_factor}, ignore_index=True)
    return attractors_smelling_factors_df
