import itertools
import math
import numpy as np

from clustering.partition import PartitionData
from sklearn.metrics.cluster import normalized_mutual_info_score

import numpy as np
from sklearn import metrics

def culculate_purity(first_partition, second_partition):
    first_part_order = list(first_partition.id_to_cluster_mapping.items())
    second_part_order = list(second_partition.id_to_cluster_mapping.items())
    for x in [x[0] for x in first_part_order if not x[0] in (x[0] for x in second_part_order)]:
        second_part_order.append((x, -1))
    for z in [y[0] for y in second_part_order if not (y[0] in (y[0] for y in first_part_order))]:
        first_part_order.append((z, -1))
    first_part_order.sort(key=lambda a: a[0])
    second_part_order.sort(key=lambda a: a[0])
    second_part_order = [x[1] for x in second_part_order]
    first_part_order = [x[1] for x in first_part_order]
    # compute contingency matrix (also called confusion matrix)
    contingency_matrix = metrics.cluster.contingency_matrix(first_part_order, second_part_order)
    # return purity
    return np.sum(np.amax(contingency_matrix, axis=0)) / np.sum(contingency_matrix)

def binom(n, k):
    if n == 0 or k > n:
        return 0
    if k == 0:
        return 1
    return math.factorial(n) // math.factorial(k) // math.factorial(n - k)


def culculate_RI(first_partition, second_partition):
    first_mapping = first_partition.id_to_cluster_mapping
    second_mapping = second_partition.id_to_cluster_mapping
    first_keys = first_mapping.keys()
    second_keys = second_mapping.keys()
    common_ids = [value for value in first_keys if value in second_keys]
    pairs = list(itertools.combinations(common_ids, 2))
    always_together = 0
    always_apart = 0
    for pair in pairs:
        first_id, second_id = pair
        if first_mapping[first_id] == first_mapping[second_id] and second_mapping[first_id] == second_mapping[
            second_id]:
            always_together += 1
        if first_mapping[first_id] != first_mapping[second_id] and second_mapping[first_id] != second_mapping[
            second_id]:
            always_apart += 1
    denumerator = binom(len(common_ids), 2)
    if denumerator == 0:
        return 0
    return (always_together + always_apart) / binom(len(common_ids), 2)


def calculate_ARI(first_partition, second_partition):
    ids = list(first_partition.id_to_cluster_mapping.keys())
    ids.extend(list(second_partition.id_to_cluster_mapping.keys()))
    ids = [value for value in list(first_partition.id_to_cluster_mapping.keys()) if
           value in second_partition.id_to_cluster_mapping.keys()]
    ids = list(map(int, ids))
    first_clusters_mappings = first_partition.cluster_id_to_metabolites_mapping
    first_clusters_mappings = dict(
        zip([x for x in range(0, len(first_clusters_mappings.keys()))], first_clusters_mappings.values()))
    second_clusters_mappings = second_partition.cluster_id_to_metabolites_mapping
    second_clusters_mappings = dict(
        zip([x for x in range(0, len(second_clusters_mappings.keys()))], second_clusters_mappings.values()))
    sigma_n_ij = 0
    ai_array = np.zeros(max(first_clusters_mappings.keys()) + 1)
    bj_array = np.zeros(max(second_clusters_mappings.keys()) + 1)
    for i in first_clusters_mappings.keys():
        for j in second_clusters_mappings.keys():
            firts_cluster = first_clusters_mappings[i]
            second_cluster = second_clusters_mappings[j]
            n_ij = len([value for value in firts_cluster if value in second_cluster])
            ai_array[i] += n_ij
            bj_array[j] += n_ij
            sigma_n_ij += (binom(n_ij, 2))
    sigma_ai = sum(list(map(lambda x: binom(x, 2), list(ai_array))))
    sigma_bj = sum(list(map(lambda x: binom(x, 2), list(bj_array))))
    try:
        ari_numerator = sigma_n_ij - sigma_ai * sigma_bj / binom(len(ids), 2)
        ari_denominator = 0.5 * (sigma_ai + sigma_bj) - (sigma_ai * sigma_bj) / binom(len(ids), 2)
        return ari_numerator / ari_denominator
    except ZeroDivisionError:
        return 0

def calculate_nmi(first_partition, second_partition):
    first_part_order = list(first_partition.id_to_cluster_mapping.items())
    second_part_order = list(second_partition.id_to_cluster_mapping.items())
    for x in [x[0] for x in first_part_order if not x[0] in (x[0] for x in second_part_order)]:
        second_part_order.append((x,-1))
    for z in [y[0] for y in second_part_order if not (y[0] in (y[0] for y in first_part_order))]:
        first_part_order.append((z,-1))
    first_part_order.sort(key=lambda a: a[0])
    second_part_order.sort(key=lambda a: a[0])
    second_part_order = [x[1] for x in second_part_order]
    first_part_order = [x[1] for x in first_part_order]
    return normalized_mutual_info_score(first_part_order, second_part_order)

def calculate_ARI_from_path(first_partition_path, second_partition_path):
    return calculate_ARI(PartitionData(first_partition_path), PartitionData(second_partition_path))


def calculate_RI_from_path(first_partition_path, second_partition_path):
    return culculate_RI(PartitionData(first_partition_path), PartitionData(second_partition_path))
