import json
import os

import networkx as nx

import data_loading.load_xcms as data_loading
from clustering.FloralClust import reclustering
from clustering.louvain import create_flower_louvain_clustered_graph_only_nodes_with_edges
from clustering.partition import PartitionData
from floralclust_utils.smell import is_attractor, get_attractors_list
from floralclust_utils.utils import invert_dict, ensure_dir_exists


def improve_atractors(data, part_data, graph, relevant_attractors=[], threshold=0.95):
    if relevant_attractors == []:
        attractors = (get_attractors_list(data, part_data))
    else:
        attractors = relevant_attractors
    for cluster in attractors:
        temp_graph, temp_part = reclustering.recluster_from_graph(original_graph_data=data,
                                                                  part_data=PartitionData(
                                                                      cluster_id_to_metabolites_mapping=part_data.cluster_id_to_metabolites_mapping.copy(),
                                                                      id_to_cluster_mapping=part_data.id_to_cluster_mapping.copy()),
                                                                  original_graph=graph,
                                                                  clusters_list=[cluster],
                                                                  threshold=threshold)
        temp_part_data = PartitionData(cluster_id_to_metabolites_mapping=invert_dict(temp_part),
                                       id_to_cluster_mapping=temp_part)
        sub_clusters_ids = [id for id in temp_part_data.cluster_id_to_metabolites_mapping.keys() if
                            int(id / 1000) == 1 + cluster]
        sub_attractors = [id for id in sub_clusters_ids if
                          is_attractor(data, id, temp_part_data)]
        sub_non_attractors = [id for id in sub_clusters_ids if (not id in sub_attractors)]
        if len(sub_attractors) != len(sub_clusters_ids) and len(sub_attractors) > 0:
            graph, part_data = improve_atractors(data, temp_part_data, temp_graph,
                                                 sub_attractors)
            graph, part_data = improve_non_atracctors(data, part_data, graph, sub_non_attractors)
    return graph, part_data


def improve_non_atracctors(data, part_data, graph, relevant_clusters=[], threshold=0.95):
    if relevant_clusters == []:
        non_attractors = [id for id in part_data.cluster_id_to_metabolites_mapping.keys() if
                          not (id in (get_attractors_list(data, part_data)))]

    else:
        non_attractors = relevant_clusters
    copy_of_part_data = PartitionData(
        cluster_id_to_metabolites_mapping=part_data.cluster_id_to_metabolites_mapping.copy(),
        id_to_cluster_mapping=part_data.id_to_cluster_mapping.copy())
    for cluster in non_attractors:
        temp_graph, temp_part = reclustering.recluster_from_graph(original_graph_data=data,
                                                                  part_data=copy_of_part_data,
                                                                  original_graph=graph,
                                                                  clusters_list=[cluster],
                                                                  threshold=threshold)
        temp_part_data = PartitionData(cluster_id_to_metabolites_mapping=invert_dict(temp_part),
                                       id_to_cluster_mapping=temp_part)
        sub_clusters_ids = [id for id in temp_part_data.cluster_id_to_metabolites_mapping.keys() if
                            int(id / 1000) == 1 + cluster]
        sub_non_attractors = [id for id in sub_clusters_ids if
                              not (is_attractor(data, id, temp_part_data))]
        sub_attractors = [id for id in sub_clusters_ids if
                          is_attractor(data, id, temp_part_data)]
        if len(sub_non_attractors) != len(sub_clusters_ids) and len(sub_non_attractors) > 0:
            graph, part_data = improve_non_atracctors(data, temp_part_data, temp_graph,
                                                      sub_non_attractors)
    #           graph, part_data = improve_atractors(data, part_data, graph, sub_attractors)
    return graph, part_data


def create_flower_FloralClust_graph(data_path, graph_threshold, addaptive_threshold):
    correlation_graph, part = create_flower_louvain_clustered_graph_only_nodes_with_edges(data_path, graph_threshold)
    # part_data = PartitionData(id_to_cluster_mapping=part)
    part_data = part
    data = data_loading.get_normilezed_xcms_df(data_path)
    graph, part_data = improve_atractors(data, part_data, correlation_graph, threshold=addaptive_threshold)
    graph, part_data = improve_non_atracctors(data, part_data, graph, threshold=addaptive_threshold)
    nodes_in_graph = list(graph.nodes)
    missing_nodes = [node for node in data.index if not (node in nodes_in_graph)]
    for node in missing_nodes:
        graph.add_node(node)
    mapping = part_data.id_to_cluster_mapping
    none_edge_cluster = max(list(mapping.values()))+1
    for id in missing_nodes:
        mapping[id] = none_edge_cluster
    part_data = PartitionData(id_to_cluster_mapping=mapping)
    return graph, part_data


def get_partition(data_path, graph_threshold, addaptive_threshold, dest_dir='', file_name=''):
    graph, part_data = create_flower_FloralClust_graph(data_path, graph_threshold, addaptive_threshold)
    if dest_dir != '':
        ensure_dir_exists(dest_dir)
        nx.write_graphml(graph, os.path.join(dest_dir, file_name + '_graph' + '.graphml'))
        with open(os.path.join(dest_dir, file_name + '_partition.json'), 'w') as outfile:
            json.dump(part_data.id_to_cluster_mapping, outfile)
    return part_data
