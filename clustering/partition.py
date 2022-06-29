import json
import utils


class PartitionData:
    def __init__(self, cluster_id_to_metabolites_mapping=[], id_to_cluster_mapping=[]):
        if id_to_cluster_mapping != [] and cluster_id_to_metabolites_mapping == []:
            inverted_partition = utils.invert_dict(id_to_cluster_mapping)
            for key in inverted_partition.keys():
                inverted_partition[key] = list(map(lambda x: int(x), inverted_partition[key]))
            cluster_id_to_metabolites_mapping = inverted_partition
        elif cluster_id_to_metabolites_mapping != [] and id_to_cluster_mapping:
            id_to_cluster_mapping = dict()
            for cluster_id in set(cluster_id_to_metabolites_mapping.keys()):
                metab_list = cluster_id_to_metabolites_mapping[cluster_id]
                for id in metab_list:
                    id_to_cluster_mapping[id] = cluster_id
        self.cluster_id_to_metabolites_mapping = cluster_id_to_metabolites_mapping
        self.id_to_cluster_mapping = id_to_cluster_mapping


def create_from_path(partition_path):
    with open(partition_path) as f:
        id_to_cluster_mapping = json.load(f)
        id_to_cluster_mapping_new = {}
        for key in id_to_cluster_mapping:
            id_to_cluster_mapping_new[int(key)] = int(id_to_cluster_mapping[key])
        id_to_cluster_mapping = id_to_cluster_mapping_new
        inverted_partition = utils.invert_dict(id_to_cluster_mapping)
        for key in inverted_partition.keys():
            inverted_partition[key] = list(map(lambda x: int(x), inverted_partition[key]))
        cluster_id_to_metabolites_mapping = inverted_partition
        return PartitionData(cluster_id_to_metabolites_mapping=cluster_id_to_metabolites_mapping,
                             id_to_cluster_mapping=id_to_cluster_mapping)


def get_sub_part_data(partition, clusters_list):
    cluster_id_to_metabolites_mapping = {}
    id_to_cluster_mapping = {}
    for cluster in clusters_list:
        if cluster in partition.cluster_id_to_metabolites_mapping:
            cluster_id_to_metabolites_mapping[cluster] = partition.cluster_id_to_metabolites_mapping[cluster]
            for id in partition.cluster_id_to_metabolites_mapping[cluster]:
                id_to_cluster_mapping[str(id)] = cluster
    return PartitionData(cluster_id_to_metabolites_mapping=cluster_id_to_metabolites_mapping,
                         id_to_cluster_mapping=id_to_cluster_mapping)
