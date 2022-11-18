import data_loading.load_xcms as load_xcms
from clustering.partition import get_sub_part_data
from compound_separation.partition_by_compounds import get_partition_after_separation
import clustering.FloralClust.corr_infomap as voc_infomap
from clustering_similarity import similarity
import pandas as pd
import numpy as np
from floralclust_utils.smell import get_attractors_list

res_df = pd.DataFrame(
    columns=['Method', "graph_threshold", "addaptive_threshold", 'Seven Adaptive Clusters', 'Eight Adaptive Clusters',
             'Seven Adaptive Smelling Clusters', 'Eight Adaptive Smelling Clusters',
             'Seven RT Clusters', 'Eight RT Clusters',
             'Seven RT Smelling Clusters', 'Eight RT Smelling Clusters',
             'RI', 'ARI', 'NMI', 'Purity', 'Smelling RI', 'Smelling ARI', 'Smelling NMI', 'Smelling Purity'])

thresholds_options  = [round(x, 2) for x in np.arange(0.5, 1, 0.05)]
thresholds_pairs = [(x,y) for x in thresholds_options for y in thresholds_options]


for (graph_threshold, addaptive_threshold) in thresholds_pairs:
    res = {}
    print(f'graph threshold: {graph_threshold}, addaptive threshold: {addaptive_threshold}')
    data_path = 'Data/Erucaria_7-19.xls'
    data = load_xcms.get_normilezed_xcms_df(data_path=data_path)
    seven_part = voc_infomap.get_partition(data_path, graph_threshold=graph_threshold, addaptive_threshold=addaptive_threshold)
    res['Seven Adaptive Clusters'] = len(seven_part.cluster_id_to_metabolites_mapping.keys())
    res['Seven Adaptive Smelling Clusters'] = len(get_attractors_list(data=data, partition_data=seven_part))
    seven_part = get_partition_after_separation(data=data,data_path=data_path,part_data=seven_part)
    res['Seven RT Clusters'] = len(seven_part.cluster_id_to_metabolites_mapping.keys())
    res['Seven RT Smelling Clusters'] = len(get_attractors_list(data=data, partition_data=seven_part))
    seven_part_smelling = get_sub_part_data(partition=seven_part, clusters_list=get_attractors_list(data=data,
                                                                                                    partition_data=seven_part))

    data_path = 'Data/Erucaria_8-20.xls'
    data = load_xcms.get_normilezed_xcms_df(data_path=data_path)
    eight_part = voc_infomap.get_partition(data_path, graph_threshold=graph_threshold, addaptive_threshold=addaptive_threshold)
    res['Eight Adaptive Clusters'] = len(eight_part.cluster_id_to_metabolites_mapping.keys())
    res['Eight Adaptive Smelling Clusters'] = len(get_attractors_list(data=data, partition_data=eight_part))
    eight_part = get_partition_after_separation(data=data,data_path=data_path,part_data=eight_part)
    res['Eight RT Clusters'] = len(eight_part.cluster_id_to_metabolites_mapping.keys())
    res['Eight RT Smelling Clusters'] = len(get_attractors_list(data=data, partition_data=eight_part))
    # print(f'eight smelling: {get_attractors_list(data=data, partition_data=eight_part)}')
    eight_part_smelling = get_sub_part_data(partition=eight_part, clusters_list=get_attractors_list(data=data,
                                                                                                    partition_data=eight_part))

    ri = similarity.culculate_RI(seven_part, eight_part)
    # # print(f'RI: {ri}')
    ari = similarity.calculate_ARI(seven_part, eight_part)
    # # print(f'ARI: {ari}')
    nmi = similarity.calculate_nmi(seven_part, eight_part)
    # # print(f'NMI: {nmi}')
    purity = similarity.culculate_purity(seven_part, eight_part)
    # print(f'Purity: {purity}')
    # print('Smelling')
    smelling_ri = similarity.culculate_RI(seven_part_smelling, eight_part_smelling)
    print(f'RI: {smelling_ri}')
    smelling_ari = similarity.calculate_ARI(seven_part_smelling, eight_part_smelling)
    print(f'ARI: {smelling_ari}')
    smelling_nmi = similarity.calculate_nmi(seven_part_smelling, eight_part_smelling)
    print(f'NMI: {smelling_nmi}')
    smelling_purity = similarity.culculate_purity(seven_part_smelling, eight_part_smelling)
    print(f'Purity: {smelling_purity}')
    print('*' * 50)
    res['Method'] = 'FloralClust-Infomap'
    res['graph_threshold'] = graph_threshold
    res['addaptive_threshold'] = addaptive_threshold
    res['RI'] = ri
    res['ARI'] = ari
    res['NMI'] = nmi
    res['Purity'] = purity
    res['Smelling RI'] = smelling_ri
    res['Smelling ARI'] = smelling_ari
    res['Smelling NMI'] = smelling_nmi
    res['Smelling Purity'] = smelling_purity
    res_df = res_df.append(res, ignore_index=True)
res_df.to_excel('Results/similarity_different_parameters_infomap.xlsx')