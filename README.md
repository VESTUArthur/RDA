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
...and MetaCycle : https://github.com/gangwug/MetaCycle
```R
# install 'devtools' in R(>3.0.2)
install.packages("devtools")
# install MetaCycle
devtools::install_github('gangwug/MetaCycle')
```

...and Rain too : https://www.bioconductor.org/packages/release/bioc/html/rain.html
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("rain")
```

## Functions

`cycMouseLiverRNA(filename)`:
    Save MetaCycle dataset cycMouseLiverRNA.
    
``cycMouseLiverProtein(filename):``
    Save MetaCycle dataset cycMouseLiverProtein.
    
``menetRNASeqMouseLiver(filename):``
    Save RAIN dataset menetRNASeqMouseLiver.
    
``meta2d_format(filename,sep=','):``
    Try format a given file to be usable with meta2d.
    
``meta2d(filename,filestyle='csv',timepoints='line1',models=["ARS", "JTK", "LS"]):``
    Perform meta2d analysis (JTK, LS and/or, if no replicates, ARS) and store the result in the metaout folder.

``cosinorpy(filename,sep=',', n_components = [1,2,3], period = 24, folder=None, **kwargs):``
    Perform Cosinor analysis and store the result in the cosinorpyout folder.

``cosinor1py(filename,sep=',', period = 24,folder=None):``
    Perform 1 component Cosinor analysis and store the result in the cosinorpyout folder.

``cosinor_pop(filename,sep,period):``
    Perform Cosinor analysis for population data and store the result in the cosinorpyout folder.

``rain(filename,sample_rate=1,n_replicate=1,period=24):``
    Perform RAIN analysis and store the result in the rainout folder.
    
``periodogram(df):``
    Plot the periodogram of a given dataset in cosinor format.

``plot_meta2d(filename,pvalue_plot=False,amplitude_plot=False,period_plot=False,qvalue_plot=False,phase_plot=False):``
    Plot meta2d result in downloadable graphic table.

``pv_load(filename):``
    Find all models p-values and save them in one file.

``pv_dist(filename):``
    Plot and save in images folder p-values distributions.

``pv_venn(filename):``
    Plot and save in images folder venn diagram using p-values.

``qv_load(filename):``
    Find all models q-values and save them in one file.

``qv_dist(filename):``
    Plot and save in images folder q-values distributions.
   
``qv_venn(filename):``
    Plot and save in images folder venn diagram using q-values.

``synt_rhythmic_data(filename,half_rnd=False,n_test=1,n_components=1,noise=0.5,replicates=1):``
    Create test data.  (rhythmic)

``synt_random_data(filename,n_test=1,replicates=1):``
    Create random test data.  (non-rhythmic)

``make_metrics(filename,y=None,half_rnd=False,conf_matrix=True,pvalue=False,qvalue=True):``
    Make metrics of an analyzed file.

``plot_metrics(filename,qvalue=True,pvalue=False):``
    Plot metrics comparaison of ARS,JTK,LS,Meta2d,Cosinor,Rain.

``file_rda(filename,filestyle='csv',metrics=False,half_rnd=True,n_components=3,replicates=1,sample_rate=2,period=24,y=None,pvalue=False,qvalue=True):``
    Perform meta2d,ARS,JTK,LS,Rain,Cosinor, make pv distribution, venn diagram and can plot metrics.

``cosinor_read(filename,sep='\t'):``
    Read file in cosinor format, xlsx, csv or txt.

``export_csv(df,filename):``
    Write a dataframe in a csv in cosinor format.

``plot_data(df,filename=None):``
    Plot file or dataframe data.

``cosinor_peaks(df,filename):``
    Plot cosinor peaks of an analysed file.

``analysis(df,filename,lines='all',dt=None,time_unit_label='hours',T_cutoff = None):``
    pyBOAT signal analysis.

``plot_detrend(x,y,deg=[1,2,3,5]):``
    Plot detrend ploynomial curve.

``detrend(x,y,deg):`` 
    Detrend data.
