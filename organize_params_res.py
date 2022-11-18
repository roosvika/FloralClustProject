import pandas as pd

res_df = pd.DataFrame(
    columns=['Method', "graph_threshold", "addaptive_threshold", 'RI', 'ARI', 'NMI', 'Purity'])

with open('Results/params_results.txt') as f:
    for line in f:
        if 'graph threshold: ' in line:
            graph_threshold = float(line.split('graph threshold: ')[1].split(',')[0])
            addaptive_threshold = float(line.split('addaptive threshold: ')[1])
        if 'RI: ' in line and not ('ARI:' in line):
            ri = float(line.split('RI: ')[1])
        if 'ARI: ' in line:
            ari = float(line.split('ARI: ')[1])
        if 'NMI: ' in line:
            nmi = float(line.split('NMI: ')[1])
        if 'Purity:' in line:
            purity = float(line.split('Purity: ')[1])
            res = {}
            res['Method'] = 'FloralClust-Louvain'
            res['graph_threshold'] = graph_threshold
            res['addaptive_threshold'] = addaptive_threshold
            res['RI'] = ri
            res['ARI'] = ari
            res['NMI'] = nmi
            res['Purity'] = purity
            res_df = res_df.append(res, ignore_index=True)
res_df.to_excel('Results/similarity_different_parameters_louvain.xlsx')
