# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 10:21:49 2024

@author: Jack
"""
import pandas as pd
import os

#%% Different split
inputloc = r"WorkSpace/SSIT_Data_V1.csv"

df_input = pd.read_csv(inputloc)

groups = df_input['Classification']
print(groups)
tps = [i[7:] for i in groups] 

df_input['timepoint'] = tps
    
    
exportpath = 'SSITData.csv'

df_input.to_csv(exportpath)
    






os.chdir(r'/home/formanj/scRNAseq_model/UpDownProject/Story')
#%% Split SSIT preprocessing into useable format
inputloc = r"SSITData.csv"

df_input = pd.read_csv(inputloc)

groups = set(df_input['Group'])
Groupnames = ['Group_'+str(i) for i in groups]

c = ['Time', 'UMAP_X1', 'UMAP_X2', 'PCA_X1', 'PCA_X2'] + Groupnames
df_output = pd.DataFrame(columns=c)

for i, r in df_input.iterrows():
    crow = [r['Time'], r['UMAP X1'], r['UMAP X2'], r['PCA X1'], r['PCA X2']] + [0] * len(groups)
    df_output.loc[len(df_output)] = crow
    df_output['Group_'+str(r['Group'])].iloc[i] = 1
    
    
exportpath = 'SSITData.csv'

df_output.to_csv(exportpath)
    
    
    
    
    
















































