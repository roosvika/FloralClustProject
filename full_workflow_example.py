from compound_hypothesises.RT_seperation import create_compounds_msps
from data_loading import load_xcms
import clustering.FloralClust.louvain as voc_louvain
import clustering.FloralClust.corr_infomap as voc_infomap
from floralclust_utils.smell import get_attractors_list
from summerize_results.barcharts_creation import create_flower_barcharts_from_partition

data_path = 'Data/combined_data.xls'
data = load_xcms.get_normilezed_xcms_df(data_path=data_path)
combined_part = voc_louvain.get_partition(data_path, graph_threshold=0.9, addaptive_threshold=0.95, dest_dir='Results/CombinedDataResults', file_name='combined_FloralClust_T_09_AT_0.95')
# combined_part = partition.create_from_path('Results/CombinedDataResults/combined_FloralClust_T_09_AT_0.95_partition.json')
data = load_xcms.get_normilezed_xcms_df(data_path=data_path)
smelling_clusters = get_attractors_list(data=data, partition_data=combined_part)
create_compounds_msps(data_path=data_path,
                      clusters_list=smelling_clusters,
                      part_data=combined_part,
                      dest='Results/CombinedDataResults/SmellingMSPS')


create_flower_barcharts_from_partition(data=data, dest='Results/CombinedDataResults/Barchars',
                                           part_data=combined_part)


# data_path = 'Data/combined_data.xls'
# data = load_xcms.get_normilezed_xcms_df(data_path=data_path)
#  combined_part = voc_infomap.get_partition(data_path, graph_threshold=0.9, addaptive_threshold=0.95)
# # combined_part = voc_louvain.get_partition(data_path, graph_threshold=0.9, addaptive_threshold=0.95, dest_dir='Results/CombinedDataResults', file_name='combined_FloralClust_T_09_AT_0.95')
# # combined_part = partition.create_from_path('Results/CombinedDataResults/combined_FloralClust_T_09_AT_0.95_partition.json')
# data = load_xcms.get_normilezed_xcms_df(data_path=data_path)
# smelling_clusters = get_attractors_list(data=data, partition_data=combined_part)
# create_compounds_msps(data_path=data_path,
#                       clusters_list=smelling_clusters,
#                       part_data=combined_part,
#                       dest='Results/CombinedDataResults/InfomapBased/SmellingMSPS')
#

# create_flower_barcharts_from_partition(data=data, dest='Results/CombinedDataResults/Barchars',
#                                            part_data=combined_part)
