import networkx as nx

from floralclust_utils import utils
from clustering.louvain import set_graph_communities, create_correlaton_graph
from clustering.partition import PartitionData

def unite_partitions(cluster_to_metabolites_mapping, single_recluster_mapping, original_cluster_id):
    """

    Args:
        cluster_to_metabolites_mapping:
        single_recluster_mapping:
        original_cluster_id:

    Returns:

    """
    new_cluster_id_to_metabolites_mapping = {}
    for cluster_id in single_recluster_mapping.keys():
        new_cluster_id_to_metabolites_mapping[1000 + 1000 * original_cluster_id + cluster_id] = \
            single_recluster_mapping[cluster_id]

    cluster_to_metabolites_mapping.pop(original_cluster_id)
    cluster_to_metabolites_mapping = {**cluster_to_metabolites_mapping,
                                      **new_cluster_id_to_metabolites_mapping}
    return cluster_to_metabolites_mapping


def recluster(original_graph_data, original_patrition_path, original_threshold, clusters_list, threshold):
    """

    Args:
        original_graph_data:
        original_patrition_path:
        original_threshold:
        clusters_list:
        threshold:

    Returns:

    """
    original_partition_data = PartitionData(original_patrition_path)
    original_mapping = original_partition_data.cluster_id_to_metabolites_mapping
    for cluster_id in clusters_list:
        cluster_metabolites = original_partition_data.cluster_id_to_metabolites_mapping[cluster_id]
        cluster_data = original_graph_data.iloc[cluster_metabolites]
        correlation_graph = create_correlaton_graph(cluster_data, threshold=threshold)
        correlation_graph, part = set_graph_communities(correlation_graph)
        original_mapping = unite_partitions(original_mapping, utils.invert_dict(part), cluster_id)
    correlation_graph = create_correlaton_graph(original_graph_data, threshold=original_threshold)
    inverted_partition = {}
    for key, value in original_mapping.items():
        for id in value:
            inverted_partition[id] = key
    nx.set_node_attributes(correlation_graph, original_mapping, 'cluster_id')
    return correlation_graph, inverted_partition


def recluster_from_graph(original_graph_data, part_data, original_graph, clusters_list, threshold):
    """

    Args:
        original_graph_data:
        original_patrition_path:
        original_threshold:
        clusters_list:
        threshold:

    Returns:

    """
    original_partition_data = part_data
    original_mapping = original_partition_data.cluster_id_to_metabolites_mapping.copy()
    for cluster_id in clusters_list:
        cluster_metabolites = original_partition_data.cluster_id_to_metabolites_mapping[cluster_id]
        cluster_data = original_graph_data.iloc[cluster_metabolites]
        correlation_graph = create_correlaton_graph(cluster_data, threshold=threshold)
        correlation_graph, part = set_graph_communities(correlation_graph)
        original_mapping = unite_partitions(original_mapping, utils.invert_dict(part), cluster_id)
    correlation_graph = original_graph.copy()
    inverted_partition = {}
    for key, value in original_mapping.items():
        for id in value:
            inverted_partition[id] = key
    for node in list(correlation_graph.nodes.keys()):
        if int(node) in inverted_partition.keys():
            correlation_graph.nodes[node]['cluster_id'] = inverted_partition[int(node)]
        else:
            correlation_graph.remove_node(node)
    return correlation_graph, inverted_partition
