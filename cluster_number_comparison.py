import data_loading.load_xcms as load_xcms
from clustering.partition import create_from_path
from compound_separation.partition_by_compounds import get_partition_after_separation
import clustering.FloralClust.louvain as voc_louvain
import clustering.FloralClust.corr_infomap as voc_infomap

from floralclust_utils.smell import get_attractors_list

print('RAMClust')
data_path = 'Data/combined_data.xls'
data = load_xcms.get_normilezed_xcms_df(data_path=data_path)
part = create_from_path('ramclust_results/combined_partition.json')
part_attractors = get_attractors_list(data=data, partition_data=part)
print(f'clusters amount: {len(part.cluster_id_to_metabolites_mapping.keys())}')
print(f'smelling clusters amount: {len(part_attractors)}')
part_attractors.sort()
print(f'smelling clusters: {part_attractors}')

print('FloralClust-Louvain')
data_path = 'Data/combined_data.xls'
data = load_xcms.get_normilezed_xcms_df(data_path=data_path)
part = voc_louvain.get_partition(data_path, graph_threshold=0.9, addaptive_threshold=0.95)
part = get_partition_after_separation(data=data,data_path=data_path,part_data=part)
print(f'clusters amount: {len(part.cluster_id_to_metabolites_mapping.keys())}')
part_attractors = get_attractors_list(data=data, partition_data=part)
part_attractors.sort()
print(f'smelling clusters amount: {len(part_attractors)}')
print(f'smelling clusters: {part_attractors}')

print('FloralClust-Infomap')
data_path = 'Data/combined_data.xls'
data = load_xcms.get_normilezed_xcms_df(data_path=data_path)
part = voc_infomap.get_partition(data_path, graph_threshold=0.9, addaptive_threshold=0.95)
part = get_partition_after_separation(data=data,data_path=data_path,part_data=part)
print(f'clusters amount: {len(part.cluster_id_to_metabolites_mapping.keys())}')
part_attractors = get_attractors_list(data=data, partition_data=part)
part_attractors.sort()
print(f'smelling clusters amount: {len(part_attractors)}')
print(f'smelling clusters: {part_attractors}')



