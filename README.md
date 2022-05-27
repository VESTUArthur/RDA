# Rhythmic Data Analysis

Include different functions that allow the user to perform rhythmic data analysis (LS, ARS, JTK, Cosinor, RAIN).

## Installation

 Can be installed with ``pip`` :
```python
pip install -i https://test.pypi.org/simple/ rda-VESTUArthur
```
Then you can ``import`` the package with :
```python
from rda_package import rda
```
The setup.cfg file install automatically the dependencies necessary for the proper functioning of the package.

Otherwise, if the dependencies do not install automatically, you must manually install them :
```python
pip install numpy
pip install pandas
pip install plotly
pip install rpy2
pip install CosinorPy
pip install matplotlib
pip install matplotlib_venn
pip install seaborn
pip install statsmodels
pip install sklearn
```
To finish, don't forget to install R if it's not already installed : https://cran.r-project.org/bin/windows/base/

## Functions

`cycMouseLiverRNA(filename)`:
    Save MetaCycle dataset cycMouseLiverRNA
    
``cycMouseLiverProtein(filename):``
    Save MetaCycle dataset cycMouseLiverProtein
    
``menetRNASeqMouseLiver(filename):``
    Save RAIN dataset menetRNASeqMouseLiver
    
``meta2d_format(filename,sep=','):``
    Try format a given file to be usable with meta2d

``meta_JTK(filename,filestyle='csv',timepoints='line1')``
    Perform JTK analysis thanks to metacycle and store the result in the metaout folder

``meta_ARS(filename,filestyle='csv',timepoints='line1')``
    Perform ARS analysis thanks to metacycle and store the result in the metaout folder

``meta_LS(filename,filestyle='csv',timepoints='line1')``
    Perform LS analysis thanks to metacycle and store the result in the metaout folder
    
``meta2d(filename,filestyle='csv',timepoints='line1'):``
    Perform meta2d analysis (JTK,LS and if no replicates ARS) and store the result in the metaout folder

``cosinorpy(filename,sep=',', n_components = 2, period = 24, names = "",folder=None, **kwargs):``
    Perform Cosinor analysis and store the result in the cosinorpyout folder

``cosinor_pop(filename,sep,period):``
    Perform Cosinor analysis for population data and store the result in the cosinorpyout folder

``rain(filename,sample_rate=1,n_replicate=1,period=24):``
    Perform RAIN analysis and store the result in the rainout folder
    
``periodogram(df):``
    Plot the periodogram of a given dataset in cosinor format

``plot_meta2d(filename,pvalue_plot=False,amplitude_plot=False,period_plot=False,qvalue_plot=False,phase_plot=False):``
    Plot meta2d result in downloadable graphic table

``pv_to_file(filename):``
    Find all models results and save them in one file

``pv_dist(filename):``
    Plot and save in images folder pvalues distributions
   
``venn(filename):``
    Plot and save in images folder venn diagram

``synt_rhythmic_data(filename,half_rnd=False,n_test=1,n_components=1,noise=0.5,replicates=1):``
    Create test data  (rhythmic)

``synt_random_data(filename,n_test=1,replicates=1):``
    Create random test data  (non-rhythmic)

``make_metrics(filename,y=None,half_rnd=False):``
    Make metrics of an analyzed file

``plot_metrics(filename):``
    Plot metrics comparaison of ARS,JTK,LS,Meta2d,Cosinor,Rain

``file_rda(filename,metrics=False,half_rnd=True,n_components=1,replicates=1,sample_rate=2,y=None):``
    Perform meta2d,ARS,JTK,LS,Rain,Cosinor, make pv distribution, venn diagram and can plot metrics