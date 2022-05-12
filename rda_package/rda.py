import numpy as np
import pandas as pd
import plotly.figure_factory as ff
import os
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from CosinorPy import file_parser, cosinor, cosinor1
import random
import re
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import statsmodels.stats.multitest as multi
from matplotlib_venn import venn2


def cycMouseLiverRNA(filename) :
        """  
        Save MetaCycle dataset cycMouseLiverRNA
        ...

        Parameters
        ----------
        filename : str
            name of the saved file
        """
        r = robjects.r
        r['library']('MetaCycle')
        r(f"""annot<-cycMouseLiverRNA[1]
            data<-cycMouseLiverRNA
            data$geneName<-NULL
            colnames(data) <- sub("..", "", colnames(data))
            head(data)
            write.csv(data,'{filename}')""")
    
def cycMouseLiverProtein(filename) :
        """  
        Save MetaCycle dataset cycMouseLiverProtein
        ...

        Parameters
        ----------
        filename : str
            name of the saved file
        """
        r = robjects.r
        r['library']('MetaCycle')
        r(f"""annot<-cycMouseLiverProtein[1]
            data<-cycMouseLiverRNA
            data$geneName<-NULL
            head(data)
            colnames(data)<-sub('..', '', colnames(data))
            write.csv(data,'{filename}')""")
            
def menetRNASeqMouseLiver(filename) :
        """  
        Save RAIN dataset menetRNASeqMouseLiver
        ...

        Parameters
        ----------
        filename : str
            name of the saved file
        """
        r = robjects.r
        r(f"""
            citation("Biobase")
            library(rain)
            data(menetRNASeqMouseLiver)
            data <- menetRNASeqMouseLiver
            rownames(data) <- NULL
            write.csv(data, "{filename}")""")

def meta2dFormat(filename,sep=','):
        """  
        Try format a given file to be usable with meta2d
        ...

        Parameters
        ----------
        filename : str
            name of the saved file
        sep : str
            separator of the file
        """
        df1 = pd.read_csv(filename, sep=sep)
        header = np.array(list(map(lambda s: re.sub('[^0-9_]','', s), df1.columns[1:])))
        print(df1.columns)
        print(header)
        reps = int(header[-1].split('_')[2])
        t1 = int(header[0].split("_")[1])
        t2 = int(header[reps].split("_")[1])
        max_time = int(header[-1].split("_")[1])
        dt = t2-t1
        measurements = max_time//dt
        shuffle = np.array([(np.arange(0,measurements*reps+1,reps))+i for i in range(reps)]).flatten()
        print(shuffle)
        x= ['id']
        print(list(map(lambda t: int(t.split('_')[1]), header[shuffle])))
        newnames=list(map(lambda t: int(t.split('_')[1]), header[shuffle]))
        x.extend(newnames)
        df1.columns=x
        df1.to_csv(f"{filename[:-4]}.csv", index=False)

def meta2d(filename,filestyle='csv',timepoints='line1'):
        """  
        Perform meta2d analysis (JTK,ARS,LS) and store the result in the metaout folder
        ...

        Parameters
        ----------
        filename : str
            name of the input file
        filestyle : str
            The data format of input file, must be "txt", or "csv", or a character vector containing field separator character(sep),
            quoting character (quote), and the character used for decimal points(dec, for details see read.table).
        timepoints : str, [] 
            a numeric vector corresponding to sampling time points of input time-series data; 
            if sampling time points are in the first line of input file, it could be set as a character sting-"Line1" or "line1".
        """
        r = robjects.r
        r['library']('MetaCycle')
        r(f"""meta2d('{filename}',timepoints='{timepoints}',filestyle = '{filestyle}')""")
        print('Meta2d Done :)')

   

def meta2dOut(filename):
        dff = pd.read_csv(f'metaout/meta2d_{filename}')
        return dff




def cosinorpy(filename,sep=',', n_components = 2, period = 24,folder=None, **kwargs):
        """  
        Perform Cosinor analysis and store the result in the cosinorpyout folder
        ...

        Parameters
        ----------
        filename : str
            name of the input file
        sep : str
            separator of the file
        n_components : int
            number of components in each models
        period : int
            period of the analysed data
        folder : str
            folder to store Plot if wanted
        """
        df=file_parser.read_csv(filename,sep)
        os.makedirs('cosinorpyout', exist_ok=True)
        df['test'] = df['test'].astype(str)
        df_results = pd.DataFrame(columns = ['test', 'period', 'n_components', 'p', 'q', 'p_reject', 'q_reject', 'RSS', 'R2', 'R2_adj', 'log-likelihood', 'amplitude', 'acrophase', 'mesor', 'peaks', 'heights', 'troughs', 'heights2'], dtype=float)
        if type(period) == int:
            period = [period]  
        if type(n_components) == int:
            n_components = [n_components]
        names = np.unique(df.test) 
        for test in names:
            for n_comps in n_components:
                for per in period:            
                    if n_comps == 0:
                        per = 100000
                    X, Y = np.array(df[df.test == test].x), np.array(df[df.test == test].y)    
                    if folder:                    
                        save_to = os.path.join(folder,test+'_compnts='+str(n_comps) +'_per=' + str(per))
                    else:
                        save_to = None              
                    results, statistics, rhythm_param, _, _ = cosinor.fit_me(X, Y, n_components = n_comps, period = per, name = test, save_to=save_to,plot=False,**kwargs)
                    try:
                        R2, R2_adj = results.rsquared,results.rsquared_adj
                    except:
                        R2, R2_adj = np.nan, np.nan
                    #TODO append change
                    df_results = df_results.append({'test': test, 
                                            'period': per,
                                            'n_components': n_comps,
                                            'p': statistics['p'], 
                                            'p_reject': statistics['p_reject'],
                                            'RSS': statistics['RSS'],
                                            'R2': R2, 
                                            'R2_adj': R2_adj,
                                            'ME': statistics['ME'],
                                            'resid_SE': statistics['resid_SE'],
                                            'log-likelihood': results.llf,        
                                            'amplitude': rhythm_param['amplitude'],
                                            'acrophase': rhythm_param['acrophase'],
                                            'mesor': rhythm_param['mesor'],
                                            'peaks': rhythm_param['peaks'],
                                            'heights': rhythm_param['heights'],
                                            'troughs': rhythm_param['troughs'],
                                            'heights2': rhythm_param['heights2']
                                            
                                            }, ignore_index=True)
                    if n_comps == 0:
                        break        
        df_results.q = multi.multipletests(df_results.p, method = 'fdr_bh')[1]
        df_results.q_reject = multi.multipletests(df_results.p_reject, method = 'fdr_bh')[1]  
        df_best_fits = cosinor.get_best_fits(df_results,criterium='RSS', reverse = False)
        df_best_fits.to_csv(f"cosinorpyout/COSINORresult_{filename}", index=False)   
        print('Cosinor Done :)') 
        return df_results
    
def cosinorpyPop(filename,sep,period):
        """  
        Perform Cosinor analysis on population data and store the result in the cosinorpyout folder
        ...

        Parameters
        ----------
        filename : str
            name of the input file
        sep : str
            separator of the file
        period : int
            period of the analysed data
        """
        df=file_parser.read_csv(filename,sep)
        os.makedirs('cosinorpyout', exist_ok=True)
        df_results = cosinor.population_fit_group(df,period=period,folder='cosinorpyout')
        df_best_fits = cosinor.get_best_fits(df_results,criterium='RSS', reverse = False)
        df_best_fits.to_csv(f"cosinorpyout/COSINORPopresult_{filename}", index=False)
    
def rain(filename,sampleRate=1,nbReplicate=1,period=24):
        """  
        Perform RAIN analysis and store the result in the rainout folder
        ...

        Parameters
        ----------
        filename : str
            name of the input file
        sampleRate : int
            the rate of the sample collection, interval between two sample
        nbReplicate : int
            number of replicate in the data set
        period : int
            period of the analysed data
        """
        r = robjects.r
        r['library']('rain')
        r(f"""data <- read.csv("{filename}", row.names = 1)
            head(data)
            sampleRate <- {sampleRate}
            nbReplicate <- {nbReplicate}
            period <- {period}
            res <- rain(t(data), deltat = sampleRate, 
                nr.series = nbReplicate, period = period, verbose = TRUE)
            dir.create("rainout", showWarnings = FALSE)
            write.csv(res, "rainout/RAINresult_{filename}")""")
        print('RAIN Done :)')

def periodogram(df):
        """  
        Plot the periodogram of a given dataset in cosinor format
        ...

        Parameters
        ----------
        df : DataFrame
            input dataframe in cosinor format
        """
        cosinor.periodogram_df(df)

def plotMeta2d(filename,pvalue_plot=False,amplitude_plot=False,period_plot=False,qvalue_plot=False,phase_plot=False):
        """  
        Plot the periodogram of a given dataset in cosinor format
        ...

        Parameters
        ----------
        filename : str
            name of the input file
        pvalue_plot : bool
            if True plot the Table
        amplitude_plot : bool
            if True plot the Table
        period_plot : bool
            if True plot the Table
        qvalue_plot : bool
            if True plot the Table
        phase_plot : bool
            if True plot the Table
        """
        dff = pd.read_csv(f'metaout/meta2d_{filename}')
        colnames=dff.columns.to_list()
        pvalue=['CycID']
        amplitude=['CycID']
        period=['CycID']
        qvalue=['CycID']
        phase=['CycID']
        for i in range(len(colnames)):
            if colnames[i][-2]=='u':
                pvalue.append(colnames[i])
            if colnames[i][-2]=='d':
                amplitude.append(colnames[i])
            if colnames[i][-2]=='o':
                period.append(colnames[i])
            if colnames[i][-2]=='.':
                qvalue.append(colnames[i])
            if colnames[i][-2]=='s':
                phase.append(colnames[i]) 
        if pvalue_plot == True:
            fig1 = ff.create_table(dff[pvalue])
            fig1.show()
            for i in range(1,len(pvalue)):
                hfig1= px.histogram(dff[pvalue[i]])
                hfig1.show()
        if amplitude_plot==True:
            fig2 = ff.create_table(dff[amplitude])
            fig2.show()
        if period_plot==True:
            fig3 = ff.create_table(dff[period])
            fig3.show()
        if qvalue_plot==True:
            fig4 = ff.create_table(dff[qvalue])
            fig4.show()
        if phase_plot==True:
            fig5 = ff.create_table(dff[phase])
            fig5.show()

def pValues(filename):
        """  
        Plot and save in images folder pvalues distributions
        ...

        Parameters
        ----------
        filename : str
            name of the input file
        """
        os.makedirs(f'images/{filename[:-4]}/dist', exist_ok=True)
        try:
            mout = pd.read_csv(f'metaout/meta2d_{filename}')
            colnames=mout.columns.to_list()
            pvalue=['CycID']
            for i in range(len(colnames)):
                if colnames[i][-2]=='u':
                    pvalue.append(colnames[i])
            pv=mout[pvalue]
        except:
            print('no meta2out of this file')
        try:
            cout=pd.read_csv(f"cosinorpyout/COSINORresult_{filename}")
            pv['Cosinor_pvalue']=cout['p']
        except:
            print('no Cosinorout of this file')
        try:
            rout=pd.read_csv(f"rainout/RAINresult_{filename}")
            pv["Rain_pvalue"]=rout["pVal"]
        except:
            print('no Rainout of this file')
        fig1 = ff.create_table(pv)
        fig1.show()
        for i in pv.columns[1:]:
            plt.hist(pv[i])
            plt.ylabel('count')
            plt.title(i)
            plt.savefig(f'images/{filename[:-4]}/dist/dist_{i}_{filename[:-4]}.png',facecolor='white')
            plt.show()
            print(i,': ok')
    
def venn(filename):
        """  
        Plot and save in images folder venn diagram
        ...

        Parameters
        ----------
        filename : str
            name of the input file
        """
        os.makedirs(f'images/{filename[:-4]}/venn', exist_ok=True)
        try:
            mout = pd.read_csv(f'metaout/meta2d_{filename}')
            colnames=mout.columns.to_list()
            pvalue=['CycID']
            for i in range(len(colnames)):
                if colnames[i][-2]=='u':
                    pvalue.append(colnames[i])
            pv=mout[pvalue]
        except:
            print('no meta2out of this file')
        try:
            cout=pd.read_csv(f"cosinorpyout/COSINORresult_{filename}")
            pv['Cosinor_pvalue']=cout['p']
        except:
            print('no Cosinorout of this file')
        try:
            rout=pd.read_csv(f"rainout/RAINresult_{filename}")
            pv["Rain_pvalue"]=rout["pVal"]
        except:
            print('no Rainout of this file')
        for col in range(1,len(pv.columns)):
            for coll in range(col,len(pv.columns)):
                print(range(1,len(pv)),'-----col',col,'coll',coll)
                l1=0
                for val in range(len(pv.values)):
                    if pv.iloc[[val], [col]].to_numpy()[0]<0.05 and pv.iloc[[val], [coll]].to_numpy()[0]<0.05:
                        l1+=1
                l2=pv[pv.iloc[:, [col]]<0.05].count()[col]
                print('cal2:',pv[pv.iloc[:, [col]]<0.05].count()[col])
                print('sel2:',pv[pv.iloc[:, [col]]<0.05].count())
                l3=pv[pv.iloc[:, [coll]]<0.05].count()[coll]
                print('cal3:',pv[pv.iloc[:, [coll]]<0.05].count()[coll])
                print('sel3:',pv[pv.iloc[:, [coll]]<0.05].count())
                print('col:',col,pv.columns[col],'coll:',coll,pv.columns[coll],'l2:',l2,'l3:',l3,'l1:',l1)
                if ((l1>0) or (l2>0) or (l3>0)) and (col!=coll):
                    venn2(subsets = (l2-l1,l3-l1,l1), set_labels = (pv.columns[col], pv.columns[coll]))
                    plt.savefig(f'images/{filename[:-4]}/venn/venn_{pv.columns[col]}_{pv.columns[coll]}_{filename[:-4]}.png',facecolor='white')
                    plt.show()
                else: 
                    print('no venn')
        return pv

def syntRhythmicData(filename):
        """  
        Load test data in cosinor format (rhythmic)
        ...

        Parameters
        ----------
        ...

        """
        replicates=3
        df = file_parser.generate_test_data(phase = 0, n_components = 1, name="test1", noise=0.5, replicates = 3)
        df2 = file_parser.generate_test_data(phase = np.pi,n_components = 1,name="test2",noise=0.5, replicates = 3)
        df=pd.concat([df,df2],ignore_index=True)
        dff=pd.DataFrame()
        for i in range(1,replicates+1):
            dfx= 'ZT_'+ df["x"].astype(int).astype(str) + f'_{i}' 
            dff=pd.concat([dff,pd.Series(dfx.unique())],ignore_index=True)
        df_new =pd.DataFrame(index=df.test.unique(),columns=dff.to_numpy().flatten()) 
        var = 0
        for i in range(len(df_new)):
            for col in df_new:
                df_new[col].iloc[i]= df.iloc[var].y
                var+=1  
        dfout=pd.DataFrame()
        for i in range(1,replicates+1):
            dfx=df["x"].astype(int).astype(str)
            print(dfx)
            dfout=pd.concat([dfout,pd.Series(dfx.unique())],ignore_index=True)
        print(dff.to_numpy().flatten())
        df_new.columns=dfout.to_numpy().flatten()
        df_new.to_csv(f"{filename[:-4]}.csv")
        return df_new

def syntRandomData(filename):
        """  
        Load test data in cosinor format (non-rhythmic)
        ...

        Parameters
        ----------
        ...

        """
        replicates=3
        df = file_parser.generate_test_data_group_random(phase = 0, n_components = 1, name="test1", noise=0.5, replicates = 3)
        df2 = file_parser.generate_test_data_group_random(phase = np.pi,n_components = 1,name="test2",noise=0.5, replicates = 3)
        df=pd.concat([df,df2],ignore_index=True)
        dff=pd.DataFrame()
        for i in range(1,replicates+1):
            dfx= 'ZT_'+ df["x"].astype(int).astype(str) + f'_{i}' 
            dff=pd.concat([dff,pd.Series(dfx.unique())],ignore_index=True)
        df_new =pd.DataFrame(index=df.test.unique(),columns=dff.to_numpy().flatten()) 
        var = 0
        for i in range(len(df_new)):
            for col in df_new:
                df_new[col].iloc[i]= df.iloc[var].y
                var+=1  
        dfout=pd.DataFrame()
        for i in range(1,replicates+1):
            dfx=df["x"].astype(int).astype(str)
            print(dfx)
            dfout=pd.concat([dfout,pd.Series(dfx.unique())],ignore_index=True)
        print(dff.to_numpy().flatten())
        df_new.columns=dfout.to_numpy().flatten()
        df_new.to_csv(f"{filename[:-4]}.csv")
        return df_new