from setuptools import setup

setup(
    name='FloralClust',
    version='1.0',
    packages=['clustering', 'clustering.FloralClust', 'data_loading', 'floralclust_utils', 'summerize_results',
              'compound_separation', 'clustering_similarity', 'compound_hypothesises'],
    url='',
    license='',
    author='Viktoria Roos',
    author_email='golobvi@post.bgu.ac.il',
    description='Detection of flower-insect related volatiles in GC-MS samples',
    install_requires=[
        'infomap==1.0.1',
        'ipython==7.16.3',
        'matplotlib==3.3.4',
        'networkx==2.5.1',
        'numpy==1.19.5',
        'pandas==1.1.5',
        'python_louvain==0.16',
        'scikit_learn==0.24.2',
        'scipy==1.5.4',
        'setuptools==40.8.0',
        'statsmodels==0.12.2',
        'xlrd==2.0.1'
    ]
)
