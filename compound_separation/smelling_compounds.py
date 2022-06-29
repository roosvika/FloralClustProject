from smell import get_attractors_list
import pandas as pd
from IPython.utils.path import ensure_dir_exists
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler
import numpy as np
import os


def create_smelling_compounds_hypothesises(data, data_path, part_data, dest, eps=0.5):
    clusters_list = get_attractors_list(data, part_data)
    smelling_rt = []
    compounds_num = 0
    ensure_dir_exists(dest)
    for cluster_id in clusters_list:
        cluster_metabolites_list = part_data.cluster_id_to_metabolites_mapping[cluster_id]
        metabolites_volume_df = data.loc[cluster_metabolites_list]
        flower_headers = [header for header in metabolites_volume_df.columns if not ('Blank' in header)]
        flower_volumes = metabolites_volume_df[flower_headers]
        flower_volumes['rtmed'] = pd.read_excel(data_path).loc[flower_volumes.index]['rtmed']
        flower_volumes['mz'] = pd.read_excel(data_path).loc[flower_volumes.index]['mzmed'].round(2)
        std_slc = StandardScaler()
        X_std = std_slc.fit_transform([[x] for x in flower_volumes['rtmed']])
        clt = DBSCAN(eps=eps, min_samples=1)
        model = clt.fit(X_std)
        clusters = pd.DataFrame(model.fit_predict(X_std))
        flower_volumes['metabolite_id'] = list(clusters[0])
        metabolites_in_cluster = set(list(clusters[0]))
        for metab_id in metabolites_in_cluster:
            matb_hyposesys = flower_volumes[flower_volumes['metabolite_id'] == metab_id]
            if len(matb_hyposesys) > 1:
                rt = np.mean(matb_hyposesys['rtmed']).round(2)
                matb_hyposesys.to_csv(os.path.join(dest, str(cluster_id) + '_' + str(rt) + '.csv'))
                smelling_rt.append(rt)
                compounds_num += 1
                for col in [col for col in flower_volumes.columns if 'sample' in col]:
                    msp_file_path = dest + '/' + str(rt) + '_' + str(cluster_id) + '_' + col.split('.')[0] + '.msp'
                    msp_file_path = msp_file_path.replace(':', '-')
                    with open(msp_file_path,
                              'w') as f:
                        f.write('NAME: ' + str(cluster_id) + '_' + str(rt) + '\n')
                        f.write('IONTYPE: positive\n\n')
                        f.write('Num Peaks: ' + str(len(matb_hyposesys)) + '\n')
                        for feature in matb_hyposesys.iterrows():
                            f.write(str(feature[1]['mz']) + '  ' + str(feature[1][col]) + '\n')
                        f.close()
    smelling_rt = list(set(smelling_rt))
    smelling_rt.sort()
    print(f'{compounds_num} smelling compounds detected')
    print(f'Retention times: {smelling_rt}')
    print(f'MSP and CSV files saved to {dest}')




