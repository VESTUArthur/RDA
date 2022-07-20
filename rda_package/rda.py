
import numpy as np
import pandas as pd
import plotly.figure_factory as ff
import os
import rpy2.robjects as robjects
from CosinorPy import file_parser, cosinor, cosinor1
import random
import re
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import statsmodels.stats.multitest as multi
from matplotlib_venn import venn2
from sklearn.preprocessing import Binarizer
from sklearn.metrics import accuracy_score
from sklearn.metrics import f1_score
from sklearn.metrics import recall_score
from sklearn.metrics import roc_auc_score
from sklearn.metrics import precision_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import matthews_corrcoef
from pyboat import WAnalyzer

def cycMouseLiverRNA(filename) :
        """  
        Save MetaCycle dataset cycMouseLiverRNA
        ...

        Parameters
        ----------
        filename : str
            name of the saved file
        """
        os.makedirs(f'Out/{filename[:-4]}', exist_ok=True)
        r = robjects.r
        r['library']('MetaCycle')
        r(f"""annot<-cycMouseLiverRNA[1]
            data<-cycMouseLiverRNA
            data$geneName<-NULL
            colnames(data) <- sub("..", "", colnames(data))
            head(data)
            write.csv(data,'Out/{filename[:-4]}/{filename}')""")
    
def cycMouseLiverProtein(filename) :
        """  
        Save MetaCycle dataset cycMouseLiverProtein
        ...

        Parameters
        ----------
        filename : str
            name of the saved file
        """
        os.makedirs(f'Out/{filename[:-4]}', exist_ok=True)
        r = robjects.r
        r['library']('MetaCycle')
        r(f"""annot<-cycMouseLiverProtein[1]
            data<-cycMouseLiverProtein
            data$geneName<-NULL
            head(data)
            colnames(data)<-sub('..', '', colnames(data))
            write.csv(data,'Out/{filename[:-4]}/{filename}')""")
            
def menetRNASeqMouseLiver(filename) :
        """  
        Save RAIN dataset menetRNASeqMouseLiver
        ...

        Parameters
        ----------
        filename : str
            name of the saved file
        """
        os.makedirs(f'Out/{filename[:-4]}', exist_ok=True)
        r = robjects.r
        r(f"""
            citation("Biobase")
            library(rain)
            data(menetRNASeqMouseLiver)
            data <- menetRNASeqMouseLiver
            rownames(data) <- NULL
            write.csv(data, "Out/{filename[:-4]}/{filename}")""")

def meta2d_format(filename,sep=','):
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
        os.makedirs(f'Out/{filename[:-4]}', exist_ok=True)
        df1 = pd.read_csv(filename, sep=sep)
        header = np.array(list(map(lambda s: re.sub('[^0-9_]','', s), df1.columns[1:])))
        print(df1.columns)
        print(header)
        try:
            reps = int(header[-1].split('_')[2])
            t1 = int(header[0].split("_")[1])
            t2 = int(header[reps].split("_")[1])
            max_time = int(header[-1].split("_")[1])
        except:
            try:
                reps = int(header[-1].split('_')[1])
                t1 = int(header[0].split("_")[0])
                t2 = int(header[reps].split("_")[0])
                max_time = int(header[-1].split("_")[0])
            except: raise
        dt = t2-t1
        measurements = max_time//dt
        shuffle = np.array([(np.arange(0,measurements*reps+1,reps))+i for i in range(reps)]).flatten()
        x= ['id']
        #print(list(map(lambda t: int(t.split('_')[1]), header[shuffle])))
        newnames=list(map(lambda t: int(t.split('_')[1]), header[shuffle]))
        x.extend(newnames)
        df1.columns=x
        df1.to_csv(f"Out/{filename[:-4]}/{filename[:-4]}.csv", index=False)
        return df1


def meta2d(filename,filestyle='csv',timepoints='line1',models=["ARS", "JTK", "LS"]):
        """  
        Perform meta2d analysis (JTK,ARS,LS) and store the result in the metaout folder
        ...

        Parameters
        ----------
        filename : str
            name of the input file (with path)
        filestyle : str
            The data format of input file, must be "txt", or "csv", or a character vector containing field separator character(sep),
            quoting character (quote), and the character used for decimal points(dec, for details see read.table).
        timepoints : str, [] 
            a numeric vector corresponding to sampling time points of input time-series data; 
            if sampling time points are in the first line of input file, it could be set as a character sting-"Line1" or "line1".
        """
        os.makedirs(f'Out/{filename.split("/")[-1][:-4]}/metaout', exist_ok=True)
        print(filename.split('/')[-1][:-4])
        r = robjects.r
        r['library']('MetaCycle')
        r(f"""meta2d('{filename}',timepoints='{timepoints}',filestyle = '{filestyle}',cycMethod = c({str(models)[1:-1]}),outdir='Out/{filename.split('/')[-1][:-4]}/metaout')""")
        print('Meta2d Done :)')

def cosinorpy(filename,sep=',',n_components = [1,2,3], period = 24,folder=None, **kwargs):
        """ 
        Perform Cosinor analysis and store the result in the cosinorpyout folder
        ...

        Parameters
        ----------
        filename : str
            name of the input file (with path)
        sep : str
            separator of the file
        n_components : int
            number of components in each models
        period : int
            period of the analysed data
        folder : str
            folder to store Plot if wanted
        """
        if(filename[-4:]=='xlsx'):
            df=file_parser.read_excel(filename)
            print(df)
        else:
            df=file_parser.read_csv(filename,sep)
        filename=filename.split("/")[-1][:-4]
        os.makedirs(f'Out/{filename}/cosinorpyout', exist_ok=True)
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
                    df_dict=pd.DataFrame({'test': [test],
                                            'period': [per],
                                            'n_components': [n_comps],
                                            'p': [statistics['p']], 
                                            'p_reject': [statistics['p_reject']],
                                            'RSS': [statistics['RSS']],
                                            'R2': [R2], 
                                            'R2_adj': [R2_adj],
                                            'ME': [statistics['ME']],
                                            'resid_SE': [statistics['resid_SE']],
                                            'log-likelihood': [results.llf],        
                                            'amplitude':[ rhythm_param['amplitude']],
                                            'acrophase': [rhythm_param['acrophase']],
                                            'mesor': [rhythm_param['mesor']],
                                            'peaks': [rhythm_param['peaks']],
                                            'heights':[ rhythm_param['heights']],
                                            'troughs': [rhythm_param['troughs']],
                                            'heights2':[ rhythm_param['heights2']]
                                            })
                    df_results= pd.concat([df_results,df_dict],ignore_index=True)
                    if n_comps == 0:
                        break        
        df_results.q = multi.multipletests(df_results.p, method = 'fdr_bh')[1]
        df_results.q_reject = multi.multipletests(df_results.p_reject, method = 'fdr_bh')[1]  
        df_best_models = cosinor.get_best_fits(df_results,criterium='RSS', reverse = False)
        df_results_extended = cosinor.analyse_best_models(df, df_best_models, analysis="bootstrap")
        df_results_extended['peaks']=df_results['peaks']
        df_results_extended['heights']=df_results['heights']
        df_results_extended.to_csv(f"Out/{filename}/cosinorpyout/COSINORresult_{filename}.csv", index=False)   
        print('Cosinor Done :)') 
        return df_results_extended

    
def cosinor1py(filename,sep=',', period = 24,folder=None):
        """ 
        Perform Cosinor analysis and store the result in the cosinorpyout folder
        ...

        Parameters
        ----------
        filename : str
            name of the input file (with path)
        sep : str
            separator of the file
        n_components : int
            number of components in each models
        period : int
            period of the analysed data
        folder : str
            folder to store Plot if wanted
        """
        if(filename[-4:]=='xlsx'):
            df=file_parser.read_excel(filename)
            print(df)
        else:
            df=file_parser.read_csv(filename,sep)
        filename=filename.split("/")[-1][:-4]
        os.makedirs(f'Out/{filename}/cosinorpyout', exist_ok=True)
        df['test'] = df['test'].astype(str)
        df_results = cosinor1.fit_group(df,period=period, plot_on=False)
        df_results.to_csv(f"Out/{filename}/cosinorpyout/COSINOR1result_{filename}.csv", index=False)   
        print('Cosinor1 Done :)') 
        return df_results


def cosinorpy_pop(filename,sep,period):
        """  
        Perform Cosinor analysis on population data and store the result in the cosinorpyout folder
        ...

        Parameters
        ----------
        filename : str
            name of the input file (with path)
        sep : str
            separator of the file
        period : int
            period of the analysed data
        """
        os.makedirs(f'Out/{filename.split("/")[-1][:-4]}/cosinorpyout', exist_ok=True)
        df=file_parser.read_csv(filename,sep)
        df_results = cosinor.population_fit_group(df,period=period,folder='cosinorpyout')
        df_best_fits = cosinor.get_best_fits(df_results,criterium='RSS', reverse = False)
        df_best_fits.to_csv(f"Out/{filename.split('/')[-1][:-4]}/cosinorpyout/COSINORPopresult_{filename.split('/')[-1][:-4]}", index=False)
    
def rain(filename,sample_rate=1,n_replicate=1,period=24):
        """  
        Perform RAIN analysis and store the result in the rainout folder
        ...

        Parameters
        ----------
        filename : str
            name of the input file (with path)
        sample_rate : int
            the rate of the sample collection, interval between two sample
        n_replicate : int
            number of replicate in the data set
        period : int
            period of the analysed data
        """
        os.makedirs(f'Out/{filename.split("/")[-1][:-4]}/rainout', exist_ok=True)
        r = robjects.r
        r['library']('rain')
        r(f"""data <- read.csv("{filename}", row.names = 1)
            sampleRate <- {sample_rate}
            nbReplicate <- {n_replicate}
            period <- {period}
            res <- rain(t(data), deltat = sampleRate, 
                nr.series = nbReplicate, period = period, verbose = TRUE)
            dir.create("Out/{filename.split('/')[-1][:-4]}/rainout", showWarnings = FALSE)
            write.csv(res, "Out/{filename.split('/')[-1][:-4]}/rainout/RAINresult_{filename.split('/')[-1][:-4]}.csv")""")
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

def plot_meta2d(filename,pvalue_plot=False,amplitude_plot=False,period_plot=False,qvalue_plot=False,phase_plot=False):
        """  
        Plot meta2d result in downloadable graphic table
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
        dff = pd.read_csv(f'Out/{filename[:-4]}/metaout/meta2d_{filename}')
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

def pv_load(filename):
        """  
        Find all models results and save them in one file
        ...

        Parameters
        ----------
        filename : str
            name of the analyzed file 
        """
        filename=filename.split("/")[-1]
        pv=pd.DataFrame()
        try:
            mout = pd.read_csv(f'Out/{filename[:-4]}/metaout/meta2d_{filename}')
            colnames=mout.columns.to_list()
            pvalue=['CycID']
            for i in range(len(colnames)):
                if colnames[i][-2]=='u':
                    pvalue.append(colnames[i])
            pv=mout[pvalue]
        except:
            print('no metaout of this file')
        try:
            cout=pd.read_csv(f"Out/{filename[:-4]}/cosinorpyout/COSINORresult_{filename}")
            pv['Cosinor_pvalue']=cout['p']
        except:
            print('no Cosinorout of this file')
        try:
            cout=pd.read_csv(f"Out/{filename[:-4]}/cosinorpyout/COSINOR1result_{filename}")
            pv['Cosinor1_pvalue']=cout['p(amplitude)']
            pv['Cosinor1(p)_pvalue']=cout['p']
        except:
            print('no Cosinor1out of this file')
        try:
            rout=pd.read_csv(f"Out/{filename[:-4]}/rainout/RAINresult_{filename}")
            pv["Rain_pvalue"]=rout['pVal']
        except:
            print('no Rainout of this file')
        pv.to_csv(f"Out/{filename[:-4]}/pv_{filename[:-4]}.csv",index=False)
        return pv

def pv_dist(filename):
        """  
        Plot and save in images folder pvalues distributions
        ...

        Parameters
        ----------
        filename : str
            name of the analyzed file

        """
        pv=pd.DataFrame()
        filename=filename.split("/")[-1]
        os.makedirs(f'Out/{filename[:-4]}/dist', exist_ok=True)
        pv=pd.read_csv(f"Out/{filename[:-4]}/pv_{filename[:-4]}.csv")
        for i in pv.columns[1:]:
            plt.hist(pv[i])
            plt.ylabel('count')
            plt.title(i)
            plt.savefig(f'Out/{filename[:-4]}/dist/pv_dist_{i}_{filename[:-4]}.png',facecolor='white')
            #plt.show()
            plt.clf()
            print(i,': ok')
        print('pValue Done :)') 
        return pv

def pv_venn(filename):
        """  
        Plot and save in images folder venn diagram
        ...

        Parameters
        ----------
        filename : str
            name of the analyzed file
        """
        pv=pd.DataFrame()
        filename=filename.split("/")[-1]
        os.makedirs(f'Out/{filename[:-4]}/venn', exist_ok=True)
        pv=pd.read_csv(f"Out/{filename[:-4]}/pv_{filename[:-4]}.csv")
        for col in range(1,len(pv.columns)):
            for coll in range(col,len(pv.columns)):
                l1=0
                for val in range(len(pv.values)):
                    if pv.iloc[[val], [col]].to_numpy()[0]<0.05 and pv.iloc[[val], [coll]].to_numpy()[0]<0.05:
                        l1+=1
                l2=pv[pv.iloc[:, [col]]<0.05].count()[col]
                l3=pv[pv.iloc[:, [coll]]<0.05].count()[coll]
                if ((l1>0) or (l2>0) or (l3>0)) and (col!=coll):
                    venn2(subsets = (l2-l1,l3-l1,l1), set_labels = (pv.columns[col], pv.columns[coll]))
                    plt.savefig(f'Out/{filename[:-4]}/venn/venn_{pv.columns[col]}_{pv.columns[coll]}_{filename[:-4]}.png',facecolor='white')
                    plt.clf()
                else: 
                    print('no venn')
        print('venn Done :)') 
        return pv

def qv_load(filename):
        """  
        Find all models results and save them in one file
        ...

        Parameters
        ----------
        filename : str
            name of the analyzed file 
        """
        filename=filename.split("/")[-1]
        qv=pd.DataFrame()
        try:
            mout = pd.read_csv(f'Out/{filename[:-4]}/metaout/meta2d_{filename}')
            colnames=mout.columns.to_list()
            qv['CycID']=mout['CycID']
            qvalue=[]
            for i in range(len(colnames)):
                if colnames[i][-2]=='u':
                    qvalue.append(colnames[i])
            for col in qvalue:
                Qs = multi.multipletests(mout[qvalue][col], method = 'fdr_bh')[1]
                qv[f"{col.split('_')[0]}_qvalue"]=Qs
        except:
            print('no metaout of this file')
        try:
            cout=pd.read_csv(f"Out/{filename[:-4]}/cosinorpyout/COSINORresult_{filename}")
            qv['Cosinor_qvalue']=cout['q']
        except:
            print('no Cosinorout of this file')
        try:
            cout=pd.read_csv(f"Out/{filename[:-4]}/cosinorpyout/COSINOR1result_{filename}")
            qv['Cosinor1_qvalue']=cout['q(amplitude)']
            qv['Cosinor1(q)_qvalue']=cout['q']
        except:
            print('no Cosinor1out of this file')
        try:
            rout=pd.read_csv(f"Out/{filename[:-4]}/rainout/RAINresult_{filename}")
            Qs = multi.multipletests(rout['pVal'], method = 'fdr_bh')[1]
            qv["Rain_qvalue"]=Qs
        except:
            print('no Rainout of this file')
        qv.to_csv(f"Out/{filename[:-4]}/qv_{filename[:-4]}.csv",index=False)
        return qv



def qv_dist(filename):
        """  
        Plot and save in images folder pvalues distributions
        ...

        Parameters
        ----------
        filename : str
            name of the analyzed file

        """
        qv=pd.DataFrame()
        filename=filename.split("/")[-1]
        os.makedirs(f'Out/{filename[:-4]}/dist', exist_ok=True)
        qv=pd.read_csv(f"Out/{filename[:-4]}/qv_{filename[:-4]}.csv")
        for i in qv.columns[1:]:
            plt.hist(qv[i])
            plt.ylabel('count')
            plt.title(i)
            plt.savefig(f'Out/{filename[:-4]}/dist/qv_dist_{i}_{filename[:-4]}.png',facecolor='white')
            #plt.show()
            plt.clf()
            print(i,': ok')
        print('pValue Done :)') 
        return qv


def qv_venn(filename):
        """  
        Plot and save in images folder venn diagram
        ...

        Parameters
        ----------
        filename : str
            name of the analyzed file
        """
        qv=pd.DataFrame()
        filename=filename.split("/")[-1]
        os.makedirs(f'Out/{filename[:-4]}/venn', exist_ok=True)
        qv=pd.read_csv(f"Out/{filename[:-4]}/qv_{filename[:-4]}.csv", )
        for col in range(1,len(qv.columns)):
            for coll in range(col,len(qv.columns)):
                l1=0
                for val in range(len(qv.values)):
                    if qv.iloc[[val], [col]].to_numpy()[0]<0.05 and qv.iloc[[val], [coll]].to_numpy()[0]<0.05:
                        l1+=1
                l2=qv[qv.iloc[:, [col]]<0.05].count()[col]
                l3=qv[qv.iloc[:, [coll]]<0.05].count()[coll]
                if ((l1>0) or (l2>0) or (l3>0)) and (col!=coll):
                    venn2(subsets = (l2-l1,l3-l1,l1), set_labels = (qv.columns[col], qv.columns[coll]))
                    plt.savefig(f'Out/{filename[:-4]}/venn/venn_{qv.columns[col]}_{qv.columns[coll]}_{filename[:-4]}.png',facecolor='white')
                    plt.clf()
                else: 
                    print('no venn')
        print('venn Done :)') 
        return qv

def synt_rhythmic_data(filename,half_rnd=False,n_test=1,n_components=1,noise=0.5,replicates=1):
        """  
        Create test data  (rhythmic)
        ...

        Parameters
        ----------
        filename : str
            name of the output file
        half_rnd : bool
            make half of the data rhythmic and the other half non-rhythmic
        n_test : int
            number of line in the dataset
        n_components : int
            number of components in the cosinor data generator
        noise : int
            % of gaussian noise added
        replicates : int
            number of replicate in the dataset
        """
        os.makedirs(f'Out/{filename[:-4]}', exist_ok=True)
        if(n_components==1):
            df_rhd=file_parser.generate_test_data_group(N=n_test,n_components = n_components, noise=noise, replicates = replicates)
        else:
           df_rhd=file_parser.generate_test_data_group(N=n_test,n_components = n_components, noise=noise, replicates = replicates, 
                amplitudes = [1, np.random.random(), np.random.random()], phase = [2*np.pi*np.random.random(), 2*np.pi*np.random.random(), 2*np.pi*np.random.random()])
        df_str_col=pd.DataFrame()
        for i in range(1,replicates+1):
            str_col= 'ZT_'+ df_rhd["x"].astype(int).astype(str) + f'_{i}' 
            df_str_col=pd.concat([df_str_col,pd.Series(str_col.unique())],ignore_index=True)
        df_rhd_res =pd.DataFrame(index=df_rhd.test.unique(),columns=df_str_col.to_numpy().flatten()) 
        var = 0
        for i in range(len(df_rhd_res)):
            for col in df_rhd_res:
                df_rhd_res[col].iloc[i]= df_rhd.iloc[var].y
                var+=1  
        df_int_col=pd.DataFrame()
        for i in range(1,replicates+1):
            int_col=df_rhd["x"].astype(int).astype(str)
            df_int_col=pd.concat([df_int_col,pd.Series(int_col.unique())],ignore_index=True)
        df_rhd_res.columns=df_int_col.to_numpy().astype(int).flatten()
        df_results = df_rhd_res
        if(half_rnd==True):
            half=n_test//2
            df_rnd_res=df_rhd_res.iloc[:half,]
            print(df_rnd_res)
            df_rhd_res=df_rhd_res.iloc[half:,]
            print('---------ok-----------')
            print(df_rhd_res)
            #df_rnd_res.columns=random.sample(list(df_int_col.to_numpy().astype(int).flatten()), len(df_int_col.to_numpy().flatten()))
            df_rnd_res.columns=random.sample(list(df_int_col.to_numpy().astype(int).flatten()), len(df_int_col.to_numpy().flatten()))   
            df_rnd_res = df_rnd_res.sort_index(axis=1)
            print(df_rnd_res.columns) 
            df_results=pd.concat([df_rhd_res,df_rnd_res],ignore_index=True)
        print(df_results.columns)          
        df_results.to_csv(f"Out/{filename[:-4]}/{filename[:-4]}.csv")
        return df_results

def synt_random_data(filename,n_test=1,replicates=1):
        """  
        Create random test data  (non-rhythmic)
        ...

        Parameters
        ----------
        filename : str
            name of the output file
        n_test : int
            number of line in the dataset
        replicates : int
            number of replicate in the dataset
        """
        os.makedirs(f'Out/{filename[:-4]}', exist_ok=True)
        df=file_parser.generate_test_data_group_random(N=n_test, replicates = replicates)
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
            #print(dfx)
            dfout=pd.concat([dfout,pd.Series(dfx.unique())],ignore_index=True)
        #print(dff.to_numpy().flatten())
        #df_new.columns=dfout.to_numpy().flatten()
        print(df_new.columns)
        df_new.columns=random.sample(list(dfout.to_numpy().astype(int).flatten()), len(dfout.to_numpy().flatten()))  
        print(df_new.columns)
        df_new = df_new.sort_index(axis=1)
        df_new.to_csv(f"Out/{filename[:-4]}/{filename[:-4]}.csv")
        #file_parser.export_csv(df,f"Out/{filename[:-4]}/{filename[:-4]}.csv")
        return df_new

def make_metrics(filename,y=None,half_rnd=False,conf_matrix=True,pvalue=False,qvalue=True):
        """  
        Make metrics of an analyzed file
        ...

        Parameters
        ----------
        filename : str
            name of the analyzed file
        y : list
            target
        half_rnd : bool
            set y first half to 1 and the other to 0
        conf_matrix : bool
            plot confusion matrix and save them in Out/filename/confusion_matrix/
        pvalue : bool
            make metrics using pvalues and save them in Out/filename/
        qvalue : bool
            make metrics using qvalues and save them in Out/filename/
        """
        filename=filename.split("/")[-1]
        vals=[]
        if(qvalue==True):
            vals.append('qv')
        if(pvalue==True):    
            vals.append('pv')
        pv=pd.DataFrame()
        for val in vals:
            if(val=='qv'):
                pv = qv_load(filename)
            elif(val=='pv'):    
                pv = pv_load(filename)
            pv = pv.dropna()
            #print(pv.iloc[:,1:])
            transformer = Binarizer(threshold=0.05).fit_transform(pv.iloc[:,1:])
            label = pv.iloc[:,[0]]
            x =pd.DataFrame(transformer, columns=pv.columns[1:]).applymap(lambda x: 1 if x==0 else 0)

            if(half_rnd==True):
                y = pd.DataFrame([1] * (len(x)//2), columns=['y'])
                y = pd.concat([y,pd.DataFrame([0] * (len(x)//2), columns=['y'])],ignore_index=True)
            pv['y'] = y
            df_results = pd.DataFrame(columns = ['auc','precision','recall','f1','accuracy','model','mcc'], dtype=float)
            for col in x :
                try:
                    if conf_matrix==True:
                        os.makedirs(f'Out/{filename[:-4]}/confusion_matrix', exist_ok=True)
                        cm = confusion_matrix(y,x[col])
                        group_names = ['True Neg','False Pos','False Neg','True Pos']
                        group_counts = ['{0:0.0f}'.format(value) for value in cm.flatten()]
                        group_percentages = ['{0:.2%}'.format(value) for value in cm.flatten()/np.sum(cm)]
                        labels = [f'{v1}\n{v2}\n{v3}' for v1, v2, v3 in zip(group_names,group_counts,group_percentages)]
                        labels = np.asarray(labels).reshape(2,2)
                        #plt.figure(figsize=(10, 10))
                        sns.heatmap(cm, annot=labels, fmt='',cmap='Blues')
                        plt.suptitle(f'{val} Metrics {col.split("_")[0]}  {filename[:-4]}',fontsize =20)
                        plt.savefig(f"Out/{filename[:-4]}/confusion_matrix/{val}_{filename[:-4]}_{col.split('_')[0]}_confusion_matrix.png", bbox_inches="tight", facecolor='white')
                        plt.show()
                    #print(col, ' accuracy:', accuracy_score(y,x[col]), ' precision:', precision_score(y,x[col]), ' recall:', recall_score(y,x[col]), ' f1:',f1_score(y,x[col]), ' auc:', roc_auc_score(y,x[col]),'mcc:',matthews_corrcoef(y,x[col]))
                    df_results= pd.concat([df_results,pd.DataFrame({'auc': [roc_auc_score(y,x[col])], 
                                                    'precision': [precision_score(y,x[col])],
                                                    'recall': [recall_score(y,x[col])],
                                                    'f1': [f1_score(y,x[col])], 
                                                    'accuracy': [accuracy_score(y,x[col])],
                                                    'model': [col.split('_')[0]],
                                                    'mcc': [matthews_corrcoef(y,x[col])]
                                                    })],ignore_index=True)
                except:
                    df_results= pd.concat([df_results,pd.DataFrame({'auc': [0], 
                                                    'precision': [precision_score(y,x[col])],
                                                    'recall': [recall_score(y,x[col])],
                                                    'f1': [f1_score(y,x[col])], 
                                                    'accuracy': [accuracy_score(y,x[col])],
                                                    'model': [col.split('_')[0]],
                                                    'mcc': [matthews_corrcoef(y,x[col])]
                                                    })],ignore_index=True)
            if(val=='qv'):
                df_results.to_csv(f"Out/{filename[:-4]}/qv_metrics_{filename[:-4]}.csv",index=False)
            elif(val=='pv'):    
                df_results.to_csv(f"Out/{filename[:-4]}/pv_metrics_{filename[:-4]}.csv",index=False)
        print('metrics Done :)')
        return pv



def plot_metrics(filename,qvalue=True,pvalue=False):
        """  
        Plot metrics comparaison of ARS,JTK,LS,Meta2d,Cosinor,Rain and save them in Out/filename/
        ...

        Parameters
        ----------
        filename : str
            name of the analyzed file
        qvalue : bool
            plot metrics using qvalues and save them in Out/filename/
        pvalue : bool
            plot metrics using qvalues and save them in Out/filename/
        """
        filename=filename.split("/")[-1]
        vals=[]
        if(qvalue==True):
            vals.append('qv')
        if(pvalue==True):    
            vals.append('pv')
        for val in vals:
            if(val=='qv'):
                df_metrics = pd.read_csv(f"Out/{filename[:-4]}/qv_metrics_{filename[:-4]}.csv")
            elif(val=='pv'):    
                df_metrics = pd.read_csv(f"Out/{filename[:-4]}/pv_metrics_{filename[:-4]}.csv")
            ncols = 2
            nrows = 3
            sns.set_style("white")
            flatui = ['#d9d9d9','#bdbdbd','#969696','#737373','#525252','#252525']
            fig, axes = plt.subplots(ncols = ncols, nrows = nrows, sharey=False)
            axes = axes.flatten()         
            fig.set_size_inches(10, 10)
            metrics = ["precision", "f1", "recall", "accuracy", "auc","mcc"]
            for ax, metric in zip(axes, metrics):
                #sns.barplot(data=df_after_ind, x='model', y=metric, ax=ax, ci=95) # ci=95 --> 95% confidence interval
                sns.barplot(data=df_metrics, x='model', y=metric, ax=ax, ci=68,palette=flatui) # ci=68 --> standard error!
                ax.set_ylabel(metric)
            plt.suptitle(f'{val} Metrics  {filename[:-4]}')
            fig.subplots_adjust(top=0.95)
            if(val=='qv'):
                plt.savefig(f"Out/{filename[:-4]}/{filename[:-4]}_qv_metrics.png", bbox_inches="tight", facecolor='white')
                plt.show()
            elif(val=='pv'):    
                plt.savefig(f"Out/{filename[:-4]}/{filename[:-4]}_pv_metrics.png", bbox_inches="tight", facecolor='white')
                plt.show()

def file_rda(filename,filestyle='csv',metrics=False,half_rnd=True,n_components=3,replicates=1,sample_rate=2,period=24,y=None,pvalue=False,qvalue=True):
        """  
        Perform meta2d,ARS,JTK,LS,Rain,Cosinor, make pv distribution, venn diagram and can plot metrics and save them in Out/filename/
        ...

        Parameters
        ----------
        filename : str
            name (with path) of the analyzed file
        filestyle : str
            type of file, csv or txt
        metrics : bool
            if true make metric
        y : list
            target
        half_rnd : bool
            create target data for half of the data rhythmic and the other half non-rhythmic
        n_components : int
            number of components in the cosinor data generator
        replicates : int
            number of replicate in the dataset
        sample_rate : int
            the rate of the sample collection, interval between two sample
        period : int
            period of the signal
        qvalue : bool
            make and plot metrics using qvalues and save them in Out/filename/
        pvalue : bool
            make and plot metrics using qvalues and save them in Out/filename/

        """
        print(filename)
        meta2d(filename,filestyle)
        rain(filename,sample_rate=sample_rate,n_replicate=replicates,period=period)
        cosinorpy(filename,n_components=n_components)
        cosinor1py(filename)
        if(pvalue==True):
            pv_load(filename)
            pv_dist(filename)
            pv_venn(filename)
        if(qvalue==True):
            qv_load(filename)
            qv_dist(filename)
            qv_venn(filename)
        if(metrics==True):
            make_metrics(filename,y,half_rnd,pvalue=pvalue,qvalue=qvalue)
            plot_metrics(filename,pvalue=pvalue,qvalue=qvalue)

def cosinor_read(filename,sep='\t'):
        """  
        Read file in cosinor format, xlsx, csv or txt.
        ...

        Parameters
        ----------
        filename : str
            name of the input file
        sep : str
            separator of the file
        """
        filestyle = filename[-4:]
        if(filestyle=='xlsx'):
            df = file_parser.read_excel(filename)
        else:
            df = file_parser.read_csv(filename,sep)
        return df

def export_csv(df,filename):
        """  
        Write a dataframe in a csv in cosinor format.
        ...

        Parameters
        ----------
        filename : str
            saved file name
        df : DataFrame
            dataframe to save in cosinor format
        """
        filestyle = filename[-4:]
        if(filestyle=='xlsx'):
            file_parser.export_csv(df,f'{filename[:-4]}csv')
        else:
            file_parser.export_csv(df,filename)

def plot_data(df,filename=None):
        """  
        Plot file or dataframe data.
        ...

        Parameters
        ----------
        filename : str
            input file name
        df : DataFrame
            input dataframe
        """
        names = df.test.unique()
        for name in names:
            dff=df[df['test']==name]
            sns.lineplot(data=dff,x='x',y='y',hue='test')
            if(filename!=None):
                plt.savefig(f"{filename}_{name}.png",facecolor='white')
            plt.show()

def cosinor_peaks(df,filename):
        """  
        Plot cosinor peaks of an analysed file.
        ...

        Parameters
        ----------
        filename : str
            input file name
        df : DataFrame
            intput dataframe
        """
        names = df.test.unique()
        print(names)
        path=f"Out/{filename[:-4]}/cosinorpyout/COSINORresult_{filename[:-4]}.csv"
        df_peaks=pd.read_csv(path)
        for name in names:
            dff=df[df['test']==name]
            df_peak=df_peaks[df_peaks['test']==name]
            x_peak= str(df_peak['peaks'].iloc[0]).strip('][').split(' ')
            y_peak=str(df_peak['heights'].iloc[0]).strip('][').split(' ')
            while("" in x_peak) :
                x_peak.remove("")
            while("" in y_peak) :
                y_peak.remove("")
            print(x_peak)
            x_peak= [float(x) for x in x_peak]
            print(y_peak)
            y_peak= [float(y) for y in y_peak]
            fig, ax = plt.subplots()
            sns.lineplot(data=dff,x='x',y='y',hue='test',ax=ax)
            ax.plot(x_peak,y_peak,'or')
            print(x_peak,y_peak)
            if(filename!=None):
                plt.savefig(f"{filename}_{name}.png",facecolor='white')
            plt.show()

def analysis(df,filename,lines='all',dt=None,time_unit_label='hours',T_cutoff = None):
        """  
        pyBOAT signal analysis.
        ...

        Parameters
        ----------
        filename : str
            input file name
        df : DataFrame
            input dataframe
        lines : str or list or int
            number of lines analysed
        dt : int
            delta time
        time_unit_label : str
            time unit in minutes or hours
        T_cutoff : int
            filter cut off
        """
        filename=filename.split('/')[-1][:-4]
        print(filename)
        ridge=pd.DataFrame()
        if lines == 'all':
            lines = range(len(df))
        if type(lines) == int:
            lines = [lines]  
        #print("lines:",lines)
        if dt==None:
            dt = int(df.columns[1]) - int(df.columns[0])
        #print("dt",dt)
        if T_cutoff == None:
            T_cutoff = int(2*df.columns[-1])
        #print('T_cutoff:',int(2*df.columns[-1]))
        for x in lines:
                #print(f'line: {x}')
                signal = df.iloc[x][1:].interpolate(method ='linear', limit_direction ='forward').to_list()
                t = df.columns[1:].astype(int).to_list()
                periods = t
                #print('periods:',periods)
                wAn= WAnalyzer(periods=periods,dt=dt, time_unit_label=time_unit_label)
                plt.ion()
                plt.suptitle(f'pyBOAT line{x} {filename}')
                modulus, ransform = wAn.compute_spectrum(signal,T_c=T_cutoff)
                ridge_tmp = wAn.get_maxRidge()
                ridge_tmp['line'] = x
                ridge = pd.concat([ridge,ridge_tmp.groupby(['line']).mean().drop('time',axis=1)])
                os.makedirs(f'Out/{filename}/analysis', exist_ok=True)
                plt.savefig(f"Out/{filename}/analysis/plt_line{x}_{filename}.png",facecolor='white')
        ridge.to_csv(f"Out/{filename}/analysis/ridge_{filename}.csv")
        return ridge

def plot_detrend(x,y,deg=[1,2,3,5]):
        """  
        Plot detrend ploynomial curve.
        ...

        Parameters
        ----------
        x: list or series
            time points
        y : list or series
            values
        deg : list or int
            ploynomial curve degree
        """
        if type(deg)!=list:
            deg=[deg]
        for deg_val in deg:
            model= np.polyfit(x,y,deg_val)
            predicted = np.polyval(model, x)
            yp=y - predicted
            plt.plot(x, y,'r')
            plt.plot(x, predicted,'b')
            plt.plot(x, yp,'g')
            plt.title(f'Detrended Residual deg={deg_val}')
            plt.show()

def detrend(x,y,deg):  
        """  
        Detrend data.
        ...

        Parameters
        ----------
        x: list or series
            time points
        y : list or series
            values
        deg : int
            ploynomial curve degree
        """
        model= np.polyfit(x,y,deg)
        predicted = np.polyval(model, x)
        yp=y - predicted
        print(f'Detrended Residual deg={deg}')
        return yp