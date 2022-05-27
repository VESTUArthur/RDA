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
from sklearn.preprocessing import Binarizer
from sklearn.metrics import accuracy_score
from sklearn.metrics import f1_score
from sklearn.metrics import recall_score
from sklearn.metrics import roc_auc_score
from sklearn.metrics import precision_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import matthews_corrcoef

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


def meta2d(filename,filestyle='csv',timepoints='line1',models=("ARS", "JTK", "LS")):
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
        r(f"""meta2d('{filename}',timepoints='{timepoints}',filestyle = '{filestyle}',cycMethod = c("ARS", "JTK", "LS"),outdir='Out/{filename.split('/')[-1][:-4]}/metaout')""")
        print('Meta2d Done :)')

def cosinorpy(filename,sep=',',folder_in='', n_components = 2, period = 24,folder=None, **kwargs):
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
        os.makedirs(f'Out/{filename[:-4]}/cosinorpyout', exist_ok=True)
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

                    """df_results= pd.concat([df_results,pd.DataFrame.from_dict({'test': test, 
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
                                            })])"""


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
        df_best_fits.to_csv(f"Out/{filename.split('/')[-1][:-4]}/cosinorpyout/COSINORresult_{filename.split('/')[-1][:-4]}.csv", index=False)   
        print('Cosinor Done :)') 
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
        os.makedirs(f'Out/{filename[:-4]}/cosinorpyout', exist_ok=True)
        df=file_parser.read_csv(filename,sep)
        os.makedirs('cosinorpyout', exist_ok=True)
        df_results = cosinor.population_fit_group(df,period=period,folder='cosinorpyout')
        df_best_fits = cosinor.get_best_fits(df_results,criterium='RSS', reverse = False)
        df_best_fits.to_csv(f"Out/{filename[:-4]}/cosinorpyout/COSINORPopresult_{filename}", index=False)
    
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
        os.makedirs(f'Out/{filename[:-4]}/rainout', exist_ok=True)
        r = robjects.r
        r['library']('rain')
        r(f"""data <- read.csv("{filename}", row.names = 1)
            sampleRate <- {sample_rate}
            nbReplicate <- {n_replicate}
            period <- {period}
            res <- rain(t(data), deltat = sampleRate, 
                nr.series = nbReplicate, period = period, verbose = TRUE)
            dir.create("Out/{filename[:-4]}/rainout", showWarnings = FALSE)
            write.csv(res, "Out/{filename[:-4]}/rainout/RAINresult_{filename}")""")
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

def pv_to_file(filename):
        """  
        Find all models results and save them in one file
        ...

        Parameters
        ----------
        filename : str
            name of the analyzed file 
        """
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
            print('no meta2out of this file')
        try:
            cout=pd.read_csv(f"Out/{filename[:-4]}/cosinorpyout/COSINORresult_{filename}")
            pv['Cosinor_pvalue']=cout['p']
        except:
            print('no Cosinorout of this file')
        try:
            rout=pd.read_csv(f"Out/{filename[:-4]}/rainout/RAINresult_{filename}")
            pv["Rain_pvalue"]=rout["pVal"]
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
        try:
            mout = pd.read_csv(f'Out/{filename[:-4]}/metaout/meta2d_{filename}')
            colnames=mout.columns.to_list()
            pvalue=['CycID']
            for i in range(len(colnames)):
                if colnames[i][-2]=='u':
                    pvalue.append(colnames[i])
            pv=mout[pvalue]
        except:
            print('no meta2out of this file')
        try:
            cout=pd.read_csv(f"Out/{filename[:-4]}/cosinorpyout/COSINORresult_{filename}")
            pv['Cosinor_pvalue']=cout['p']
        except:
            print('no Cosinorout of this file')
        try:
            rout=pd.read_csv(f"Out/{filename[:-4]}/rainout/RAINresult_{filename}")
            pv["Rain_pvalue"]=rout["pVal"]
        except:
            print('no Rainout of this file')
        pv.to_csv(f"Out/{filename[:-4]}/pv_{filename[:-4]}.csv",index=False)
        #fig1 = ff.create_table(pv)
        #fig1.show()
        for i in pv.columns[1:]:
            plt.hist(pv[i])
            plt.ylabel('count')
            plt.title(i)
            plt.savefig(f'Out/{filename[:-4]}/dist/dist_{i}_{filename[:-4]}.png',facecolor='white')
            #plt.show()
            plt.clf()
            print(i,': ok')
        print('pValue Done :)') 
        return pv
    
def venn(filename):
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
        try:
            mout = pd.read_csv(f'Out/{filename[:-4]}/metaout/meta2d_{filename}')
            colnames=mout.columns.to_list()
            pvalue=['CycID']
            for i in range(len(colnames)):
                if colnames[i][-2]=='u':
                    pvalue.append(colnames[i])
            pv=mout[pvalue]
        except:
            print('no meta2out of this file')
        try:
            cout=pd.read_csv(f"Out/{filename[:-4]}/cosinorpyout/COSINORresult_{filename}")
            pv['Cosinor_pvalue']=cout['p']
        except:
            print('no Cosinorout of this file')
        try:
            rout=pd.read_csv(f"Out/{filename[:-4]}/rainout/RAINresult_{filename}")
            pv["Rain_pvalue"]=rout["pVal"]
        except:
            print('no Rainout of this file')
        pv.to_csv(f"Out/{filename[:-4]}/pv_{filename[:-4]}.csv", index=False)
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

def make_metrics(filename,y=None,half_rnd=False,conf_matrix=True):
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
        """
        pv=pd.DataFrame()
        try:
            mout = pd.read_csv(f'Out/{filename[:-4]}/metaout/meta2d_{filename}',sep=',')
            colnames=mout.columns.to_list()
            pvalue=['CycID']
            for i in range(len(colnames)):
                if colnames[i][-2]=='u':
                    pvalue.append(colnames[i])
            pv=mout[pvalue]
        except:
            print('no meta2out of this file')
        try:
            cout=pd.read_csv(f"Out/{filename[:-4]}/cosinorpyout/COSINORresult_{filename}")
            pv['Cosinor_pvalue']=cout['p']
        except:
            print('no Cosinorout of this file')
        try:
            rout=pd.read_csv(f"Out/{filename[:-4]}/rainout/RAINresult_{filename}")
            pv["Rain_pvalue"]=rout["pVal"]
        except:
            print('no Rainout of this file')
        pv = pv.dropna()
        print(pv.iloc[:,1:])
        transformer = Binarizer(threshold=0.05).fit_transform(pv.iloc[:,1:])
        label = pv.iloc[:,[0]]
        x =pd.DataFrame(transformer, columns=pv.columns[1:]).applymap(lambda x: 1 if x==0 else 0)
        if(half_rnd==True):
            y = pd.DataFrame([1] * (len(x)//2), columns=['y'])
            y = pd.concat([y,pd.DataFrame([0] * (len(x)//2), columns=['y'])],ignore_index=True)
        pv['y'] = y
        pv.to_csv(f"Out/{filename[:-4]}/pv_{filename[:-4]}.csv",index=False)
        df_results = pd.DataFrame(columns = ['auc','precision','recall','f1','accuracy','model'], dtype=float)
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
                    plt.suptitle(f'Metrics {col.split("_")[0]}  {filename[:-4]}',fontsize =20)
                    plt.savefig(f"Out/{filename[:-4]}/confusion_matrix/{filename[:-4]}_{col.split('_')[0]}_confusion_matrix.png", bbox_inches="tight", facecolor='white')
                    plt.show()
                print(col, ' accuracy:', accuracy_score(y,x[col]), ' precision:', precision_score(y,x[col]), ' recall:', recall_score(y,x[col]), ' f1:',f1_score(y,x[col]), ' auc:', roc_auc_score(y,x[col]),'mcc:',matthews_corrcoef(y,x[col]))
                df_results = df_results.append({'auc': roc_auc_score(y,x[col]), 
                                                'precision': precision_score(y,x[col]),
                                                'recall': recall_score(y,x[col]),
                                                'f1': f1_score(y,x[col]), 
                                                'accuracy': accuracy_score(y,x[col]),
                                                'model': col.split('_')[0],
                                                'mcc': matthews_corrcoef(y,x[col])
                                                }, ignore_index=True)
            except:
                print(col, ' accuracy:', accuracy_score(y,x[col]), ' precision:', precision_score(y,x[col]), ' recall:', recall_score(y,x[col]), ' f1:',f1_score(y,x[col]), ' auc:', 0)
                tn, fp, fn, tp = confusion_matrix(y,x[col]).ravel()
                df_results = df_results.append({'auc': 0, 
                                                'precision': precision_score(y,x[col]),
                                                'recall': recall_score(y,x[col]),
                                                'f1': f1_score(y,x[col]), 
                                                'accuracy': accuracy_score(y,x[col]),
                                                'model': col.split('_')[0],
                                                'mcc': matthews_corrcoef(y,x[col])
                                                }, ignore_index=True)
                
        print('metrics Done :)') 
        df_results.to_csv(f"Out/{filename[:-4]}/metrics_{filename[:-4]}.csv",index=False)
        return df_results

def plot_metrics(filename):
        """  
        Plot metrics comparaison of ARS,JTK,LS,Meta2d,Cosinor,Rain
        ...

        Parameters
        ----------
        filename : str
            name of the analyzed file
        """
        df_metrics = pd.read_csv(f"Out/{filename[:-4]}/metrics_{filename[:-4]}.csv")
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
        plt.suptitle(f'Metrics  {filename[:-4]}')
        fig.subplots_adjust(top=0.95)
        plt.savefig(f"Out/{filename[:-4]}/{filename[:-4]}_metrics.png", bbox_inches="tight", facecolor='white')
        plt.show()

def file_rda(filename,path='',filestyle='csv',metrics=False,half_rnd=True,n_components=1,replicates=1,sample_rate=2,period=24,y=None):
        """  
        Perform meta2d,ARS,JTK,LS,Rain,Cosinor, make pv distribution, venn diagram and can plot metrics
        ...

        Parameters
        ----------
        filename : str
            name (with path) of the analyzed file
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
        """
        if(path!=''):
            file=f"{path}/{filename}"
        else:
            file=filename
        meta2d(file,filestyle)
        rain(file,sample_rate=sample_rate,n_replicate=replicates,period=period)
        cosinorpy(file,n_components=n_components)
        pv_dist(filename)
        venn(filename)
        if(metrics==True):
            make_metrics(filename,y,half_rnd)
            plot_metrics(filename)

def cosinor_read(filename,sep='\t'):
    filestyle = filename[-4:]
    if(filestyle=='xlsx'):
        df = file_parser.read_excel(filename)
    else:
        df = file_parser.read_csv(filename,sep)
    return df

def export_csv(df,filename):
    filestyle = filename[-4:]
    if(filestyle=='xlsx'):
        file_parser.export_csv(df,f'{filename[:-4]}csv')
    else:
       file_parser.export_csv(df,filename)

def plot_data(df,filename=None):
    names = df.test.unique()
    for name in names:
        dff=df[df['test']==name]
        sns.lineplot(data=dff,x='x',y='y',hue='test')
        if(filename!=None):
            plt.savefig(f"{filename}_{name}.png",facecolor='white')
        plt.show()