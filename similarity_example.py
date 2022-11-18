import data_loading.load_xcms as load_xcms
from clustering.partition import get_sub_part_data, create_from_path
import clustering.FloralClust.louvain as voc_louvain
import clustering.FloralClust.corr_infomap as voc_infomap
import clustering.louvain as louvain
from clustering_similarity import similarity
import clustering.corr_infomap as corr_infomap
import pandas as pd

from floralclust_utils.smell import get_attractors_list

res_df = pd.DataFrame(columns=['Method','RI', 'ARI', 'NMI', 'Purity', 'SmellingRI', 'SmellingARI', 'SmellingNMI', 'SmellingPurity'])
# RAMClust
print('RAMClust')
data_path = 'Data/Erucaria_7-19.xls'
data = load_xcms.get_normilezed_xcms_df(data_path=data_path)
seven_part = create_from_path('ramclust_results/seven_partition.json')
seven_attractors = get_attractors_list(data=data, partition_data=seven_part)
print(f'seven smelling: {seven_attractors}')
seven_part_smelling = get_sub_part_data(partition=seven_part, clusters_list=seven_attractors)

data_path = 'Data/Erucaria_8-20.xls'
data = load_xcms.get_normilezed_xcms_df(data_path=data_path)
eight_part = create_from_path('ramclust_results/eight_partition.json')
eight_attractors = get_attractors_list(data=data, partition_data=eight_part)
print(f'eight smelling: {eight_attractors}')
eight_part_smelling = get_sub_part_data(partition=eight_part, clusters_list=eight_attractors)


ri=similarity.culculate_RI(seven_part, eight_part)
print(f'RI: {ri}')
ari=similarity.calculate_ARI(seven_part, eight_part)
print(f'ARI: {ari}')
nmi=similarity.calculate_nmi(seven_part, eight_part)
print(f'NMI: {nmi}')
purity = similarity.culculate_purity(seven_part, eight_part)
print(f'Purity: {purity}')
print('Smelling')
smelling_ri = similarity.culculate_RI(seven_part_smelling, eight_part_smelling)
print(f'RI: {smelling_ri}')
smelling_ari = similarity.calculate_ARI(seven_part_smelling, eight_part_smelling)
print(f'ARI: {smelling_ari}')
smelling_nmi = similarity.calculate_nmi(seven_part_smelling, eight_part_smelling)
print(f'NMI: {smelling_nmi}')
smelling_purity = similarity.culculate_purity(seven_part_smelling, eight_part_smelling)
print(f'Purity: {smelling_purity}')
print('*' * 50)
res={}
res['Method']='RAMClust'
res['RI']=ri
res['ARI']=ari
res['NMI']=nmi
res['Purity']=purity
res['SmellingRI']=smelling_ri
res['SmellingARI']=smelling_ari
res['SmellingNMI']=smelling_nmi
res['SmellingPurity']=smelling_purity
res_df = res_df.append(res, ignore_index=True)


# Louvain
print('Louvain')
data_path = 'Data/Erucaria_7-19.xls'
data = load_xcms.get_normilezed_xcms_df(data_path=data_path)
seven_part = louvain.get_partition(data_path, threshold=0.9)
print(f'seven smelling: {get_attractors_list(data=data, partition_data=seven_part)}')
seven_part_smelling = get_sub_part_data(partition=seven_part, clusters_list=get_attractors_list(data=data,
                                                                                                partition_data=seven_part))

data_path = 'Data/Erucaria_8-20.xls'
data = load_xcms.get_normilezed_xcms_df(data_path=data_path)
eight_part = louvain.get_partition(data_path, threshold=0.9)
print(f'eight smelling: {get_attractors_list(data=data, partition_data=eight_part)}')
eight_part_smelling = get_sub_part_data(partition=eight_part, clusters_list=get_attractors_list(data=data,
                                                                                                partition_data=eight_part))
ri=similarity.culculate_RI(seven_part, eight_part)
print(f'RI: {ri}')
ari=similarity.calculate_ARI(seven_part, eight_part)
print(f'ARI: {ari}')
nmi=similarity.calculate_nmi(seven_part, eight_part)
print(f'NMI: {nmi}')
purity = similarity.culculate_purity(seven_part, eight_part)
print(f'Purity: {purity}')
print('Smelling')
smelling_ri = similarity.culculate_RI(seven_part_smelling, eight_part_smelling)
print(f'RI: {smelling_ri}')
smelling_ari = similarity.calculate_ARI(seven_part_smelling, eight_part_smelling)
print(f'ARI: {smelling_ari}')
smelling_nmi = similarity.calculate_nmi(seven_part_smelling, eight_part_smelling)
print(f'NMI: {smelling_nmi}')
smelling_purity = similarity.culculate_purity(seven_part_smelling, eight_part_smelling)
print(f'Purity: {smelling_purity}')
print('*' * 50)
res={}
res['Method']='Louvain'
res['RI']=ri
res['ARI']=ari
res['NMI']=nmi
res['Purity']=purity
res['SmellingRI']=smelling_ri
res['SmellingARI']=smelling_ari
res['SmellingNMI']=smelling_nmi
res['SmellingPurity']=smelling_purity
res_df = res_df.append(res, ignore_index=True)

# Infomap
print('Infomap')
data_path = 'Data/Erucaria_7-19.xls'
data = load_xcms.get_normilezed_xcms_df(data_path=data_path)
seven_part = corr_infomap.get_partition(data_path, graph_threshold=0.9)
print(f'seven smelling: {get_attractors_list(data=data, partition_data=seven_part)}')
seven_part_smelling = get_sub_part_data(partition=seven_part, clusters_list=get_attractors_list(data=data,
                                                                                                partition_data=seven_part))

data_path = 'Data/Erucaria_8-20.xls'
data = load_xcms.get_normilezed_xcms_df(data_path=data_path)
eight_part = corr_infomap.get_partition(data_path, graph_threshold=0.9)
print(f'eight smelling: {get_attractors_list(data=data, partition_data=eight_part)}')
eight_part_smelling = get_sub_part_data(partition=eight_part, clusters_list=get_attractors_list(data=data,
                                                                                                partition_data=eight_part))

ri=similarity.culculate_RI(seven_part, eight_part)
print(f'RI: {ri}')
ari=similarity.calculate_ARI(seven_part, eight_part)
print(f'ARI: {ari}')
nmi=similarity.calculate_nmi(seven_part, eight_part)
print(f'NMI: {nmi}')
purity = similarity.culculate_purity(seven_part, eight_part)
print(f'Purity: {purity}')
print('Smelling')
smelling_ri = similarity.culculate_RI(seven_part_smelling, eight_part_smelling)
print(f'RI: {smelling_ri}')
smelling_ari = similarity.calculate_ARI(seven_part_smelling, eight_part_smelling)
print(f'ARI: {smelling_ari}')
smelling_nmi = similarity.calculate_nmi(seven_part_smelling, eight_part_smelling)
print(f'NMI: {smelling_nmi}')
smelling_purity = similarity.culculate_purity(seven_part_smelling, eight_part_smelling)
print(f'Purity: {smelling_purity}')
print('*' * 50)
res={}
res['Method']='Infomap'
res['RI']=ri
res['ARI']=ari
res['NMI']=nmi
res['Purity']=purity
res['SmellingRI']=smelling_ri
res['SmellingARI']=smelling_ari
res['SmellingNMI']=smelling_nmi
res['SmellingPurity']=smelling_purity
res_df = res_df.append(res, ignore_index=True)

# FloralClustLouvain
print('VOClust-Louvain')
data_path = 'Data/Erucaria_7-19.xls'
data = load_xcms.get_normilezed_xcms_df(data_path=data_path)
seven_part = voc_louvain.get_partition(data_path, graph_threshold=0.9, addaptive_threshold=0.95)
print(f'seven smelling: {get_attractors_list(data=data, partition_data=seven_part)}')
seven_part_smelling = get_sub_part_data(partition=seven_part, clusters_list=get_attractors_list(data=data,
                                                                                                partition_data=seven_part))

data_path = 'Data/Erucaria_8-20.xls'
data = load_xcms.get_normilezed_xcms_df(data_path=data_path)
eight_part = voc_louvain.get_partition(data_path, graph_threshold=0.9, addaptive_threshold=0.95)
print(f'eight smelling: {get_attractors_list(data=data, partition_data=eight_part)}')
eight_part_smelling = get_sub_part_data(partition=eight_part, clusters_list=get_attractors_list(data=data,
                                                                                                partition_data=eight_part))


ri=similarity.culculate_RI(seven_part, eight_part)
print(f'RI: {ri}')
ari=similarity.calculate_ARI(seven_part, eight_part)
print(f'ARI: {ari}')
nmi=similarity.calculate_nmi(seven_part, eight_part)
print(f'NMI: {nmi}')
purity = similarity.culculate_purity(seven_part, eight_part)
print(f'Purity: {purity}')
print('Smelling')
smelling_ri = similarity.culculate_RI(seven_part_smelling, eight_part_smelling)
print(f'RI: {smelling_ri}')
smelling_ari = similarity.calculate_ARI(seven_part_smelling, eight_part_smelling)
print(f'ARI: {smelling_ari}')
smelling_nmi = similarity.calculate_nmi(seven_part_smelling, eight_part_smelling)
print(f'NMI: {smelling_nmi}')
smelling_purity = similarity.culculate_purity(seven_part_smelling, eight_part_smelling)
print(f'Purity: {smelling_purity}')
print('*' * 50)
res={}
res['Method']='FloralClust-Louvain'
res['RI']=ri
res['ARI']=ari
res['NMI']=nmi
res['Purity']=purity
res['SmellingRI']=smelling_ri
res['SmellingARI']=smelling_ari
res['SmellingNMI']=smelling_nmi
res['SmellingPurity']=smelling_purity
res_df = res_df.append(res, ignore_index=True)


# FloralClustInfomap
print('FloralClust-Infomap')
data_path = 'Data/Erucaria_7-19.xls'
data = load_xcms.get_normilezed_xcms_df(data_path=data_path)
seven_part = voc_infomap.get_partition(data_path, graph_threshold=0.9, addaptive_threshold=0.95)
print(f'seven smelling: {get_attractors_list(data=data, partition_data=seven_part)}')
seven_part_smelling = get_sub_part_data(partition=seven_part, clusters_list=get_attractors_list(data=data,
                                                                                                partition_data=seven_part))

data_path = 'Data/Erucaria_8-20.xls'
data = load_xcms.get_normilezed_xcms_df(data_path=data_path)
eight_part = voc_infomap.get_partition(data_path, graph_threshold=0.9, addaptive_threshold=0.95)
print(f'eight smelling: {get_attractors_list(data=data, partition_data=eight_part)}')
eight_part_smelling = get_sub_part_data(partition=eight_part, clusters_list=get_attractors_list(data=data,
                                                                                                partition_data=eight_part))


ri=similarity.culculate_RI(seven_part, eight_part)
print(f'RI: {ri}')
ari=similarity.calculate_ARI(seven_part, eight_part)
print(f'ARI: {ari}')
nmi=similarity.calculate_nmi(seven_part, eight_part)
print(f'NMI: {nmi}')
purity = similarity.culculate_purity(seven_part, eight_part)
print(f'Purity: {purity}')
print('Smelling')
smelling_ri = similarity.culculate_RI(seven_part_smelling, eight_part_smelling)
print(f'RI: {smelling_ri}')
smelling_ari = similarity.calculate_ARI(seven_part_smelling, eight_part_smelling)
print(f'ARI: {smelling_ari}')
smelling_nmi = similarity.calculate_nmi(seven_part_smelling, eight_part_smelling)
print(f'NMI: {smelling_nmi}')
smelling_purity = similarity.culculate_purity(seven_part_smelling, eight_part_smelling)
print(f'Purity: {smelling_purity}')
print('*' * 50)
res={}
res['Method']='FloralClust-Infomap'
res['RI']=ri
res['ARI']=ari
res['NMI']=nmi
res['Purity']=purity
res['SmellingRI']=smelling_ri
res['SmellingARI']=smelling_ari
res['SmellingNMI']=smelling_nmi
res['SmellingPurity']=smelling_purity
res_df = res_df.append(res, ignore_index=True)



print(res_df)
res_df.to_csv('Results/clustering_similarity_modified.csv')