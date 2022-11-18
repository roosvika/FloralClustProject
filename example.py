import data_loading.load_xcms as load_xcms
import clustering
from compound_separation.smelling_compounds import create_smelling_compounds_hypothesises
import clustering.FloralClust.louvain as floral_louvain
from clustering_similarity import similarity
import clustering.corr_infomap as corr_infomap

data_path = 'Data/combined_data.xls'
data = load_xcms.get_normilezed_xcms_df(data_path=data_path)
# louvain_part = floral_louvain.get_partition(data_path=data_path, graph_threshold=0.9, addaptive_threshold=0.95)
# create_smelling_compounds_hypothesises(data, data_path, louvain_part, 'Results/ErucariaCombined/Smelling/MSP')




infomap_part = corr_infomap.get_partition(data_path=data_path,graph_threshold=0.8)
create_smelling_compounds_hypothesises(data, data_path, infomap_part, 'Results/ErucariaCombined/Infomap')
