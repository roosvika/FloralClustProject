import os

import pandas as pd
import numpy as np
from IPython.utils.path import ensure_dir_exists
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler

from data_loading.load_xcms import get_xcms_df
from sklearn.neighbors import NearestNeighbors
from matplotlib import pyplot as plt


def create_compounds_msps(data_path, clusters_list, part_data , dest, eps=0.1):
    data = get_xcms_df(data_path)
    ensure_dir_exists(dest)
    for cluster_id in clusters_list:
        cluster_metabolites_list = part_data.cluster_id_to_metabolites_mapping[cluster_id]
        metabolites_volume_df = data.loc[cluster_metabolites_list]
        flower_headers = [header for header in metabolites_volume_df.columns if not ('blank' in header.lower())]
        flower_volumes = metabolites_volume_df[flower_headers]
        flower_volumes['rtmed'] = pd.read_excel(data_path).loc[flower_volumes.index]['rtmed']
        flower_volumes['mz'] = pd.read_excel(data_path).loc[flower_volumes.index]['mzmed'].round(2)
        std_slc = StandardScaler()
        X_std = std_slc.fit_transform([[x] for x in flower_volumes['rtmed']])
        clt = DBSCAN(eps=0.5, min_samples=1)
        model = clt.fit(X_std)
        clusters = pd.DataFrame(model.fit_predict(X_std))
        flower_volumes['metabolite_id'] = list(clusters[0])
        metabolites_in_cluster = set(list(clusters[0]))
        for metab_id in metabolites_in_cluster:
            matb_hyposesys = flower_volumes[flower_volumes['metabolite_id'] == metab_id]
            if len(matb_hyposesys) > 1 and metab_id>-1:
                rt = np.mean(matb_hyposesys['rtmed']).round(2)
                matb_hyposesys.to_csv(os.path.join(dest, str(cluster_id) + '_' + str(rt) + '.csv'))
                for col in [col for col in flower_volumes.columns if 'sample' in col]:
                    with open(dest+'/'+ str(rt)+'_'+str(cluster_id) + '_' + col.split('.')[0].replace(':','_') + '.msp', 'w') as f:
                        f.write('NAME: ' + str(cluster_id) + '_' + str(rt) + '\n')
                        f.write('IONTYPE: positive\n\n')
                        f.write('Num Peaks: ' + str(len(matb_hyposesys)) + '\n')
                        for feature in matb_hyposesys.iterrows():
                            f.write(str(feature[1]['mz']) + '  ' + str(feature[1][col]) + '\n')
                        f.close()


