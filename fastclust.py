import fastcluster

import data_loading.load_xcms as data_loading
from clustering.louvain import create_links
from scipy.spatial.distance import pdist
from dynamicTreeCut import cutreeHybrid


def some_func(data_path):
    data = data_loading.get_normilezed_xcms_df(data_path)
    data = data_loading.get_flower_samples_xcms_df(data)
    metabolites_data = data.T
    links = create_links(metabolites_data)
    a=5


def get_partition(data_path, dest_dir='', file_name=''):
    correlation_graph, part = some_func(data_path)
    # if dest_dir != '':
    #     nx.write_graphml(correlation_graph, os.path.join(dest_dir, file_name + '_graph' + '.graphml'))
    #     with open(os.path.join(dest_dir, file_name + '_partition.json'), 'w') as outfile:
    #         json.dump(part, outfile)
    # return partition.PartitionData(id_to_cluster_mapping=part)

get_partition('Data/combined_data.xls')
