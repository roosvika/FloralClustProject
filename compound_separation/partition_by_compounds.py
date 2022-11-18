from clustering.partition import PartitionData
import pandas as pd
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler


def get_partition_after_separation(data, data_path, part_data, eps=0.5):
    # clusters_list = get_attractors_list(data, part_data)
    clusters_list = list(part_data.cluster_id_to_metabolites_mapping.keys())
    id_to_cluster_mapping={}
    cluster_to_metabolites_mapping = {}
    old_ids_to_new = {}
    cluster_ids = list(part_data.cluster_id_to_metabolites_mapping.keys())
    final_cluster_to_metabolites = {}
    for new_id in range(1, len(cluster_ids)+1):
        cluster_to_metabolites_mapping[new_id] = part_data.cluster_id_to_metabolites_mapping[cluster_ids[new_id-1]]
        old_ids_to_new[cluster_ids[new_id-1]] = new_id
    max_cluter_id = max(cluster_to_metabolites_mapping.keys())
    free_cluster_id = 1
    smelling_rt = []
    compounds_num = 0
    for cluster_id in clusters_list:
        cluster_metabolites_list = part_data.cluster_id_to_metabolites_mapping[cluster_id]
        metabolites_volume_df = data.loc[cluster_metabolites_list]
        flower_headers = [header for header in metabolites_volume_df.columns if not ('Blank' in header)]
        flower_volumes = metabolites_volume_df[flower_headers]
        flower_volumes['rtmed'] = pd.read_excel(data_path).loc[flower_volumes.index]['rtmed']
        flower_volumes['mz'] = pd.read_excel(data_path).loc[flower_volumes.index]['mzmed'].round(2)
        std_slc = StandardScaler()
        X_std = std_slc.fit_transform([[x] for x in flower_volumes['rtmed']])
        clt = DBSCAN(eps=eps, min_samples=2)
        model = clt.fit(X_std)
        clusters = pd.DataFrame(model.fit_predict(X_std))
        flower_volumes['metabolite_id'] = list(clusters[0])
        metabolites_in_cluster = set(list(clusters[0]))
        for metab_id in metabolites_in_cluster:
            matb_hyposesys = flower_volumes[flower_volumes['metabolite_id'] == metab_id]
            if len(matb_hyposesys) > 1 and metab_id>-1:
                final_cluster_to_metabolites[free_cluster_id] = list(matb_hyposesys.index)
                free_cluster_id+=1
    for cluster_id in final_cluster_to_metabolites.keys():
        for metab in final_cluster_to_metabolites[cluster_id]:
            id_to_cluster_mapping[metab] = cluster_id
    part = PartitionData(cluster_id_to_metabolites_mapping=final_cluster_to_metabolites,
                         id_to_cluster_mapping=id_to_cluster_mapping)
    return part



# data_path = '../Data/combined_data.xls'
# data = load_xcms.get_normilezed_xcms_df(data_path=data_path)
# # combined_part = clustering.louvain.get_partition(data_path, threshold=0.9)
# combined_part = voc_louvain.get_partition(data_path, graph_threshold=0.9, addaptive_threshold=0.95)
# data = load_xcms.get_normilezed_xcms_df(data_path=data_path)
# part_after_separation = get_partition_after_separation(data =data, data_path=data_path,
#                       part_data=combined_part)
# smelling_clusters = get_attractors_list(data=data, partition_data=combined_part)
# smelling_clusters = get_attractors_list(data=data, partition_data=part_after_separation)
# print(smelling_clusters)
#
# smelling_RTs = []
# data = pd.read_excel(data_path)
# data = data.sort_values(by=['featureidx'])
# for cluster_id in smelling_clusters:
#     metabs = np.mean((data.loc[part_after_separation.cluster_id_to_metabolites_mapping[cluster_id]])['rtmed']).round(2)
#     smelling_RTs+=[metabs]
# print(smelling_RTs)