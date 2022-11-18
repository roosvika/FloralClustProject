import itertools
import math
import numpy as np

from clustering.partition import PartitionData
from sklearn.metrics.cluster import normalized_mutual_info_score
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.cluster import rand_score

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



def culculate_RI(first_partition, second_partition):
    first_part_order = list(first_partition.id_to_cluster_mapping.items())
    second_part_order = list(second_partition.id_to_cluster_mapping.items())
    first_part_order = [x for x in first_part_order if (x[0] in [y[0] for y in second_part_order])]
    second_part_order = [x for x in second_part_order if (x[0] in [y[0] for y in first_part_order])]
    first_part_order.sort(key=lambda a: a[0])
    second_part_order.sort(key=lambda a: a[0])
    second_part_order = [x[1] for x in second_part_order]
    first_part_order = [x[1] for x in first_part_order]
    return rand_score(first_part_order, second_part_order)


def calculate_ARI(first_partition, second_partition):
    first_part_order = list(first_partition.id_to_cluster_mapping.items())
    second_part_order = list(second_partition.id_to_cluster_mapping.items())
    first_part_order = [x for x in first_part_order if (x[0] in [y[0] for y in second_part_order])]
    second_part_order = [x for x in second_part_order if (x[0] in [y[0] for y in first_part_order])]
    first_part_order.sort(key=lambda a: a[0])
    second_part_order.sort(key=lambda a: a[0])
    second_part_order = [x[1] for x in second_part_order]
    first_part_order = [x[1] for x in first_part_order]
    return adjusted_rand_score(first_part_order, second_part_order)


def calculate_nmi(first_partition, second_partition):
    first_part_order = list(first_partition.id_to_cluster_mapping.items())
    second_part_order = list(second_partition.id_to_cluster_mapping.items())
    first_part_order = [x for x in first_part_order if (x[0] in [y[0] for y in second_part_order])]
    second_part_order = [x for x in second_part_order if (x[0] in [y[0] for y in first_part_order])]
    first_part_order.sort(key=lambda a: a[0])
    second_part_order.sort(key=lambda a: a[0])
    second_part_order = [x[1] for x in second_part_order]
    first_part_order = [x[1] for x in first_part_order]
    return normalized_mutual_info_score(first_part_order, second_part_order)


def calculate_ARI_from_path(first_partition_path, second_partition_path):
    return calculate_ARI(PartitionData(first_partition_path), PartitionData(second_partition_path))


def calculate_RI_from_path(first_partition_path, second_partition_path):
    return culculate_RI(PartitionData(first_partition_path), PartitionData(second_partition_path))
