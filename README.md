# FloralClustProject
The identification of floral volatiles that ensure pollination services is
an important step in creating and implementing ecologically-friendly
solutions to agricultural challenges. Typically, volatile organic compounds (VOCs) are identified based on Gas Chromatography-Mass
Spectrometry (GC-MS) targeted analysis, a time-consuming and subjective process relying on prior knowledge and manual curation. To
improve this, we propose FloralClust, a novel pseudo-targeted approach
for VOC classification based on an adaptive clustering algorithm that
considers flowers’ circadian rhythms as flower-insect communication
agents.

### Environment
This software was developed and supports Python 3.6 and Windows operating system.
### Installation
Installation should take up to 5 minutes.
1. Using .whl:

    1.1. Download .whl file

     malunally from: https://github.com/roosvika/FloralClustProject/tree/master/dist

     using wget:
     
      ```console
        foo@bar:~$ wget https://github.com/roosvika/FloralClustProject/blob/master/dist/FloralClust-1.0-py2.py3-none-any.whl
      ```    
   
    1.2. Install using pip
      ```console
        foo@bar:~$ pip install FloralClust-1.0-py2.py3-none-any.whl
      ```  
 2. Using git clone:
    ```console
        foo@bar:~$ git clone https://github.com/roosvika/FloralClustProject.git
        foo@bar:~$ cd FloralClustProject
        foo@bar:~$ pip install -r requirements.txt
    ```

### Input data file
The input file should be an agrregated XCMS (suffix .xls) file which can be generated using XCMS Online: https://xcmsonline.scripps.edu/landing_page.php?pgcontent=mainPage.
For each feature should appear the following metadata columns: featureidx,	name,	mzmed, mzmin,	mzmax,	rtmed,	rtmin,	rtmax,	npeaks.
Additionally - for each sample taken at a certain time interval should appear a column of the following format `'sample_<start time in 24 hour format>-<end time in 24 hour format>_<sample number>'`.

For instance sample_8:00-11:00_2. The values in those columns should be the generated intensity by XCMS online (flower samples). The number of replicants for time interval and the numer of time intervals 
are not limited. 
Additionaly for each time interval add one blank sample column. 
The header of those columns should be of the following format `'sample_<start time in 24 hour format>-<end time in 24 hour format>_blank'`

Example file:
![image](https://user-images.githubusercontent.com/62721219/202784744-3134b92a-7f1d-412f-8857-c87054dc6cfd.png)

The data files used in our proof of concept are located in https://github.com/roosvika/FloralClustProject/tree/master/Data.

### Run FloralClust
Relevent imports:
```python
from compound_hypothesises.RT_seperation import create_compounds_msps
from data_loading import load_xcms
import clustering.FloralClust.louvain as voc_louvain
import clustering.FloralClust.corr_infomap as voc_infomap
from floralclust_utils.smell import get_attractors_list
from summerize_results.barcharts_creation import create_flower_barcharts_from_partition
```

Load Data:
```python
data_path = 'Data/combined_data.xls'
data = load_xcms.get_normilezed_xcms_df(data_path=data_path)
```

Create partiotion using FloralClust-Louvain
(This might take several minutes)

```python
combined_part = voc_louvain.get_partition(data_path, graph_threshold=0.9, addaptive_threshold=0.95, dest_dir='Results/CombinedDataResults', file_name='combined_FloralClust_T_09_AT_0.95')
```

Create partiotion using FloralClust-Infomap
(This might take several minutes)
```python
combined_part = voc_infomap.get_partition(data_path, graph_threshold=0.9, addaptive_threshold=0.95, dest_dir='Results/CombinedDataResults', file_name='combined_FloralClust_T_09_AT_0.95')
```


Get only flower-insect communication related clusters:
```python
data = load_xcms.get_normilezed_xcms_df(data_path=data_path)
smelling_clusters = get_attractors_list(data=data, partition_data=combined_part)
```
Create files representing the flower-insect communication
```python
create_compounds_msps(data_path=data_path,
                      clusters_list=smelling_clusters,
                      part_data=combined_part,
                      dest='Results/CombinedDataResults/SmellingMSPS')
```
For each compound 2 files are generated: .csv and .msp.
The file names represent the cluster from stage 3 of FloralClust an the avarage RT.
The .csv files contain all the m/z intensity measure of the compounds which can be compared to different libraries and the .msp files can be an input for Nist program:
![image](https://user-images.githubusercontent.com/62721219/202797003-8c200655-90c1-4880-953a-1375edf4e35e.png)



Create barcharchs describing the dayly pattern of the different clusters emission pattern.
```python
create_flower_barcharts_from_partition(data=data, dest='Results/CombinedDataResults/Barchars',
                                           part_data=combined_part)
```

Results for the running example for the data file https://github.com/roosvika/FloralClustProject/blob/master/Data/combined_data.xls can be found in https://github.com/roosvika/FloralClustProject/tree/master/Results/CombinedDataResults.
The total run of the example might take some time (up to 10-15 minutes depending on the speed of your computer.)


