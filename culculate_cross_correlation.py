import numpy as np
from scipy.interpolate import interp1d
import statsmodels.api as sm

flower_timeline = [8, 11, 14, 17, 20]
# dict: keys - cluster ids, values - volumes during the day
flower_clusters_volumes = {}
flower_clusters_volumes['2'] = [372.67, 322.58, 219.36, 122.24, 122.24]
flower_clusters_volumes['3'] = [48.98, 63.16, 49.68, 21.34, 21.34]
flower_clusters_volumes['4'] = [45.40, 45.02, 20.78, 8.80, 8.80]
flower_clusters_volumes['7'] = [11.87, 5.29, 4.44, 5.58, 5.58]
flower_clusters_volumes['11'] = [17.19, 28.09, 10.47, 4.02, 4.02]
flower_clusters_volumes['13'] = [0.41, 0.64, 0.38, 0.17, 0.17]
flower_clusters_volumes['31'] = [0.21, 0.3, 0.63, 0.49, 0.49]

bugs_amount = {}
bugs_amount["Hymenopterans"] = [0, 3, 10.5, 9, 1.5, 0]
bugs_amount["Diptera"] = [0, 3.5, 7.5, 1, 1, 0]
bugs_amount["Coleoptera"] = [0, 1, 5.5, 4.5, 1, 0]
bugs_amount["Lepidoptera"] = [0, 0.5, 1, 1.5, 0.5, 0]

for bug_name, bugs in bugs_amount.items():
    best_corr = 0
    best_lag = 0
    best_cluster = 0
    for cluster_id, volumes in flower_clusters_volumes.items():
        x = np.array(flower_timeline)
        y = np.array(volumes)
        f = interp1d(x, y)
        flower = [f(8), f(10), f(12), f(14), f(16), f(18)]
        cross_corr = sm.tsa.stattools.ccf(bugs, flower, adjusted=False)
        max_corr = max(cross_corr)
        if max_corr > best_corr and np.where(cross_corr == max_corr)[0][0]<1:
            best_corr = max_corr
            best_cluster = cluster_id
            best_lag = np.where(cross_corr == max_corr)[0][0]
    print(f'bug: {bug_name}')
    print(f'cluster: {best_cluster}')
    print(f'cross correlation: {best_corr}')
    print(f'lag: {best_lag}')
    print('-' * 20)
