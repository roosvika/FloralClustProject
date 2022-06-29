import data_loading.load_xcms as data_loading
from clustering.louvain import create_correlaton_graph
import networkx as nx
import infomap
import os
import json
import pandas as pd
from clustering.partition import PartitionData


def find_communities(G):
    """
    Partition network with the Infomap algorithm.
    Annotates nodes with 'community' id and return number of communities found.
    """
    infomapWrapper = infomap.Infomap("--two-level --silent")

    for e in list(G.edges):
        infomapWrapper.addLink(*e)

    infomapWrapper.run();

    tree = infomapWrapper.tree
    communities = infomapWrapper.get_modules()

    nx.set_node_attributes(G, communities, 'cluster_id')
    return G, communities



def create_flower_infomap_clustered_graph(data_path, threshold):
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
    g, part = find_communities(correlation_graph)
    nodes_in_graph = list(g.nodes)
    missing_nodes = [node for node in data.index if not (node in nodes_in_graph)]
    for node in missing_nodes:
        g.add_node(node)
    mapping = part
    for id in missing_nodes:
        mapping[id] = -1
    part = PartitionData(id_to_cluster_mapping=mapping)
    return g, part


def create_flower_infomap_clustered_graph_only_nodes_with_edges(data_path, threshold):
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
    correlation_graph, part = find_communities(correlation_graph)
    part = PartitionData(id_to_cluster_mapping=part)
    data = pd.read_excel(data_path)
    data = data[[col for col in data.columns if not ('sample' in col)]]
    nodes_attributes = dict(zip(data.index, data.to_dict('records')))
    nx.set_node_attributes(correlation_graph, nodes_attributes)
    return correlation_graph, part


def get_partition(data_path, graph_threshold, dest_dir='', file_name=''):
    graph, part_data = create_flower_infomap_clustered_graph(data_path, graph_threshold)
    if dest_dir != '':
        nx.write_graphml(graph, os.path.join(dest_dir, file_name + '_graph' + '.graphml'))
        with open(os.path.join(dest_dir, file_name + '_partition.json'), 'w') as outfile:
            json.dump(part_data.id_to_cluster_mapping, outfile)
    return part_data

