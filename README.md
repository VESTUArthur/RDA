# Rhythmic Data Analysis

Include different functions that allow the user to perform rhythmic data analysis (LS, ARS, JTK, Cosinor, RAIN)

## Installation

 Can be intalled with ``pip`` :
```python
pip install -i https://test.pypi.org/simple/ rda-VESTUArthur
```
then you can ``import`` the package with :
```python
from rda_package import rda
```
## Functions

`cycMouseLiverRNA(filename)`:
    Save MetaCycle dataset cycMouseLiverRNA
    
``cycMouseLiverProtein(filename):``
    Save MetaCycle dataset cycMouseLiverProtein
    
``menetRNASeqMouseLiver(filename):``
    Save RAIN dataset menetRNASeqMouseLiver
    
``meta2dFormat(filename,sep=','):``
    Try format a given file to be usable with meta2d
    
``meta2d(filename,filestyle='csv',timepoints='line1'):``
    Perform meta2d analysis (JTK,ARS,LS) and store the result in the metaout folder

``cosinorpy(filename,sep=',', n_components = 2, period = 24, names = "",folder=None, **kwargs):``
    Perform Cosinor analysis and store the result in the cosinorpyout folder

``cosinorPop(df,filename,period):``
    Perform Cosinor analysis for population data and store the result in the cosinorpyout folder

``rain(filename,sampleRate=1,nbReplicate=1,period=24):``
    Perform RAIN analysis and store the result in the rainout folder
    
``periodogram(df):``
    Plot the periodogram of a given dataset in cosinor format

``plotMeta2d(filename,pvalue_plot=False,amplitude_plot=False,period_plot=False,qvalue_plot=False,phase_plot=False):``
    Plot meta2d result in downloadable graphic table

``pValues(filename):``
    Plot and save in images folder pvalues distributions
   
``venn(filename):``
    Plot and save in images folder venn diagram

``cosinorTestData(df=file_parser.generate_test_data(phase = 0, n_components = 1, name="test1", noise=0.5, replicates = 3),phase = np.pi,n_components = 1,name="test2",noise=0.5, replicates = 3):``
    Load test data in cosinor format


