#imports
from eqcorrscan import Tribe #reading tgz files with templates
import obspy
from obspy import Stream
import csv #for reading and saving data
import numpy as np
import pandas as pd
from glob import glob #looking through files
import scipy.signal as sp
import yaml #for config file
import matplotlib.pyplot as plt #for plotting
from tqdm import trange
from time import time #for time for code to run


#set parameters
with open('/home/smocz/expand_redpy/scripts/config.yaml') as file:
    config = yaml.load(file, Loader=yaml.FullLoader)
    
volc_list_names = config['volc_list_names']
vv = config['vv']
volc = volc_list_names[vv]
print(volc)
homedir = config['homedir']

# volc = "Rainier" #can change volc name to look at specific volcano without reloading config.yaml
# K = "7"

for volc in volc_list_names:
    for K in [5,6,7,8]:
        tdf = pd.read_csv(f'/home/smocz/expand_redpy_new_files/kmeans_K_{K}_Volcano_{volc}.csv') #read csv

        long_cl_list = [] #list of clusters from template names
        for t_name in tdf['template_name']: #for each template name
            long_cl_list.append(t_name[-3:]) #append the cluster ID

        cl_list = np.unique(long_cl_list) #sort and get rid of duplicates

        cdf = pd.DataFrame(data={'Cluster_ID':[],'percent_confidence':[],'Kmeans':[],'total_number_templates':[],'group_ties':[]}) #make dataframe to hold stats

        for cl in cl_list:
            print('cluster',cl)

            tdf_rows = [] #list of lists, each list within is a row of the dataframe
            for index, rows in tdf.iterrows(): #itterate each row of tdf
                row =[rows.template_name, rows.Kmeans] #record the template name and Kmeans group
                tdf_rows.append(row) #save to tdf_rows

            N_list = []
            for row in tdf_rows: #for each row of dataframe
                if row[0][-3:] == cl: #if the last 3 characters of the template name are the correct cluster ID
                    N_list.append(row[0]) #save the template name

            N = len(N_list) #number of templates for a cluster

            group_list = np.unique(tdf['Kmeans'].values.tolist())

            nt_list = [] #list of number of rows/templates for each group, same index as group_list
            for gr in group_list: #for each K number/number of groups
                nr = [] #list of rows that fit the critera
                for row in tdf_rows: #for each row
                    if row[0][-3:]==cl and row[1]==gr: #if the cluster is correct and the K number is correct
                        nr.append(row) #save to nr
                nt = len(nr) #number of rows/templates that fit the criteria
                nt_list.append(nt)


            n = max(nt_list) #number of templates in the most reliable group
            tie = ''
            if nt_list.count(max(nt_list)) > 1:
                print('nt_list has more than one of this max')
                tie = '*'
            n_idx = nt_list.index(n) #index can be applied to group_list

            print(f'{n/N*100} % confidence that cluster {cl} is part of kmeans group {group_list[n_idx]}. Based on {N}\
            total templates')
            print('------------------')

            cdf.loc[len(cdf)+1] = [cl,n/N*100,group_list[n_idx],N,tie] #add this result

        #     break
        cdf.to_csv(f'{homedir}confidence_K_{K}_Volcano_{volc}.csv')