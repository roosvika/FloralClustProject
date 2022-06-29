import data_loading
import pandas as pd
import numpy as np
import networkx as nx
import os
import json
import networkx as nx
import community.community_louvain as community_louvain
import data_loading.load_xcms as data_loading
from clustering import partition
from clustering.partition import PartitionData


def create_links(metabolites_data):
    """Returns links data according to correlations
    Parameters:
    metabolites_data (DataFrame): Flower metabolites data
    Returns:
    DataFrame: A data frame containing node id's and their correlation value

    """
    corr = metabolites_data.astype(float).corr()
    links = corr.stack()
    index_list = links.index.to_list()
    first = list(map(lambda t: t[0], index_list))
    second = list(map(lambda t: t[1], index_list))
    values = list(links.to_frame('c')['c'])
    links = pd.DataFrame(list(zip(first, second, values)),
                         columns=['var1', 'var2', 'value'])
    links = links.astype({'var1': int, 'var2': int, 'value': float})
    return links


def get_edge_list(links, threshold):
    """Returns links data according to threshold without self correlation
    Parameters:
    links (DataFrame): Links data
    threshold (float): correlation threshold
    Returns:
    DataFrame: Links data without rows who's correlation value below threshold
                and without self correlations

    """
    links_filtered = links.loc[(links['value'] > threshold)
                               & (links['var1'] != links['var2'])]
    return links_filtered


def create_correlaton_graph(data, threshold, axis=0):
    """Returns a graph from Pandas DataFrame containing an edge list.
    Parameters:
        data (String): Flower data path
        threshold (float): correlation threshold
        axis (int): 0 for correlation according to columns and 1 for correlation
                    according to rows


    """
    if axis == 0:
        metabolites_data = data.T
    else:
        assert (axis == 1, "Axis should be 0 or 1 only")
    links = create_links(metabolites_data)
    edge_list = get_edge_list(links, threshold)
    g = nx.from_pandas_edgelist(edge_list, 'var1', 'var2')
    # nodes_in_graph = set(list(edge_list['var1']) + list(edge_list['var2']))
    # missing_nodes = [node for node in data.index if not (node in nodes_in_graph)]
    # for node in missing_nodes:
    #     g.add_node(node)
    return g


def set_graph_communities(g):
    """Creates best partition based on Louvain clustering and sets cluster id for
    graph nodes
    Parameters:
    g (Graph): graph from Pandas DataFrame containing an edge list
    Returns:
    tuple(g, part): where g is the given graph with cluster_id attribute added to the nodes
                    and part is a dictionary representing the partition

    """
    part = community_louvain.best_partition(g)
    nx.set_node_attributes(g, part, 'cluster_id')
    return g, part


def get_graph_partition(g):
    """

    Args:
        g:

    Returns:

    """
    return community_louvain.best_partition(g)


def create_flower_louvain_clustered_graph(data_path, threshold):
    """Returns Flower only metabolites correlation graph, with cluster_id as attributes and partition

    Args:
        data: Total metabolites DataFrame
        threshold: Correlation threshold for edge creation between metabolites
        dest_dir: (Optional) dest directory to save graph and partition
        file_name: partition ang graph files name

    Returns:
        Flower only metabolites correlation graph, with cluster_id as attributes and partition
    """
    data = data_loading.get_normilezed_xcms_df(data_path)
    data = data_loading.get_flower_samples_xcms_df(data)
    correlation_graph = create_correlaton_graph(data, threshold=threshold)
    correlation_graph, part = set_graph_communities(correlation_graph)
    nodes_in_graph = list(correlation_graph.nodes)
    missing_nodes = [node for node in data.index if not (node in nodes_in_graph)]
    for node in missing_nodes:
        correlation_graph.add_node(node)
    mapping = part
    for id in missing_nodes:
        mapping[id] = -1
    part = PartitionData(id_to_cluster_mapping=mapping)
    data = pd.read_excel(data_path)
    data = data[[col for col in data.columns if not ('sample' in col)]]
    nodes_attributes = dict(zip(data.index, data.to_dict('records')))
    nx.set_node_attributes(correlation_graph, nodes_attributes)
    return correlation_graph, part


def create_flower_louvain_clustered_graph_only_nodes_with_edges(data_path, threshold):
    """Returns Flower only metabolites correlation graph, with cluster_id as attributes and partition

    Args:
        data: Total metabolites DataFrame
        threshold: Correlation threshold for edge creation between metabolites
        dest_dir: (Optional) dest directory to save graph and partition
        file_name: partition ang graph files name

    Returns:
        Flower only metabolites correlation graph, with cluster_id as attributes and partition
    """
    data = data_loading.get_normilezed_xcms_df(data_path)
    data = data_loading.get_flower_samples_xcms_df(data)
    correlation_graph = create_correlaton_graph(data, threshold=threshold)
    correlation_graph, part = set_graph_communities(correlation_graph)
    part = PartitionData(id_to_cluster_mapping=part)
    data = pd.read_excel(data_path)
    data = data[[col for col in data.columns if not ('sample' in col)]]
    nodes_attributes = dict(zip(data.index, data.to_dict('records')))
    nx.set_node_attributes(correlation_graph, nodes_attributes)
    return correlation_graph, part


def get_partition(data_path, threshold, dest_dir='', file_name=''):
    correlation_graph, part = create_flower_louvain_clustered_graph(data_path, threshold)
    if dest_dir != '':
        nx.write_graphml(correlation_graph, os.path.join(dest_dir, file_name + '_graph' + '.graphml'))
        with open(os.path.join(dest_dir, file_name + '_partition.json'), 'w') as outfile:
            json.dump(part, outfile)
    return partition.PartitionData(id_to_cluster_mapping=part.id_to_cluster_mapping)
