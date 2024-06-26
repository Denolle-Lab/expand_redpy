#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import h5py
import yaml
import csv
from time import time
from glob import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import signal

#read config file for parameters
with open('/home/smocz/expand_redpy/scripts/config.yaml') as file:
    config = yaml.load(file, Loader=yaml.FullLoader)

smooth_length = config['smooth_length']
fs = config['fs']
tb = config['tb']
ta = config['ta']
fqmin = config['fqmin']
fqmax = config['fqmax']
chan = config['chan']
homedir = config['homedir']
readdir = config['readdir']
minsta = config['minsta']
grid_length = float(config['grid_length'])
grid_height = float(config['grid_height'])
step = config['step']
t_step = config['t_step']
vs_min = config['vs_min']
vs_max = config['vs_max']
vs_step = config['vs_step']
volc_lat_lon = config['volc_lat_lon']
volc_list_names = config['volc_list_names']
nlta = config['nlta']

vv = config['vv']
vv=3

iend = 3520 #end of possible pick times
istart = 176 #beginning of possible pick times (currently excluding taper)


print(volc_list_names[vv])

volc_md = pd.read_csv(readdir+'Volcano_Metadata.csv')
#associate network and station
volc_md['netsta'] = volc_md['Network'].astype(str)+'.'+volc_md['Station'].astype(str)


#get clusterid from template name
def getcl_id(t_name_str): #for normalized
    t_cl = int(t_name_str.split('_')[-1])
    return t_cl

#get station name from template name
def getnet_sta(template_str): #for normalized
    t_net = template_str.split('_')[0]
    t_sta = template_str.split('_')[1]
    netsta = t_net+'_'+t_sta
    return t_net, t_sta


# In[ ]:


#get some info
v = volc_md[volc_md['Volcano_Name'] == volc_list_names[vv]]['netsta'].values.tolist() #list of network and station per volc
zz = chan[-2:].lower() #the last two letters of channel names (essentially the letters in chan)
csv_name = f'{homedir}locations/{volc_list_names[vv]}_ELEP_normalized_picktimes.csv' #name of the csv for picktimes at this volcano
h5_name = f'{homedir}h5/{volc_list_names[vv]}_ELEP_smb_pred.h5' #name of h5 file for smb_pred


# In[ ]:


smbs = []
temps = []

with h5py.File(h5_name, "r") as f: #pull in fingerprints
    for ff in f.keys():
#         print(ff)
        if ff.startswith('smb_pred'):
            smbs.append(f[ff][()])
        if ff.startswith('template'):
            temps.append(f[ff][()])


# In[ ]:


with open(csv_name, 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Network','Station','Cluster_ID','Template_Name','SMB_peak']) #,'SMB_peak_MBF'
    file.close()


# In[ ]:


for smb, template in zip(smbs,temps): #go through each cluster
    
    time0 = time()
    
    print('---')
    cl = getcl_id(str(template[0])[2:-1])
    print(f'for cluster {cl}')
    
    dup_list = [] #a list to check if any smb has been saved twice
    
    temps_p = [] #list of template picks for the cluster
    rows = [] #list of csv rows for this cluster

    for ss,s in enumerate(smb):
        smb_pred = s[0][0]
        t_name = str(template[ss])[2:-1]
        
        if t_name in dup_list: #if we have already seen this template/smb
            continue #skip to next one
            
        else:
        
            dup_list.append(t_name) #record seeing this template/smb

#             #plot
            plt.plot(smb_pred);
            print(f'before mbf: {smb_pred}')
            plt.title(t_name)
            plt.show()
            plt.close()
            
            sfs = fs
            istart = int(istart) #t_before*sfs - t_around*sfs
            iend = int(iend) #np.min((t_before*sfs + t_around*sfs,smb_pred.shape[1]))
            print(f'istart {istart} or {istart/sfs} seconds; iend {iend} or {iend/sfs} seconds')

                
            print(f'after mbf: {smb_pred}')

            if np.isnan(smb_pred[istart:iend]).any(): #if smb_pred is just nans, skip
                continue
            else:

                ### from mbf_elep_fun.py ###

                peaks = signal.find_peaks(smb_pred[istart:iend],distance=5*sfs, height=0.03)

                if len(peaks[0]) == 0:
                    peaks = signal.find_peaks(smb_pred[istart:iend],distance=5*sfs)
                    
                print('peaks:')
                print(peaks[0]) ## PEAKS, NOT PICKS. they tell you how long after istart ###
                
                
                ### END from mbf_elep_fun.py ###


                picks = [i+istart for i in peaks[0]] # ADDED ISTART
                temps_p.append(picks)
                csv_picks=' '.join([str(i/fs)for i in picks]) #formatted for saving in csv
                
                net, sta = getnet_sta(t_name)

                #append picks to csv
                row = [net, sta, cl, t_name, csv_picks] ### APPEND PICKS NOT PEAKS
                rows.append(row)

    ### FILTERING PEAKS ###

    one_peak = [] #list of n and y if the template only has one peak/picktime
    for p in temps_p: #for saved picks
        if len(p)==1:
            one_peak.append('y')
        else:
            one_peak.append('n')
            
            
    if one_peak.count('y')/len(one_peak) ==1: #if there is only one peak per template
        #write info to csv as is
        for row in rows:
            with open(csv_name, 'a', newline='') as file:
                writer = csv.writer(file)
                writer.writerow(row)
                file.close()

    if one_peak.count('y')/len(one_peak) >=0.5 and one_peak.count('y')/len(one_peak) <1: #if 75% of templates have one peak

        one_p_value = [] #list of peak values when there is only one
        for peaks in temps_p:
            if len(peaks)==1:
                one_p_value.append(peaks[0])
            else:
                continue

        for en,peaks in enumerate(temps_p):
            if len(peaks)>1:
                sums = []#list of sum of differences between possible peaks and the "confirmed" peaks
                for peak in peaks:
                    sums.append(sum([abs(peak-v) for v in one_p_value]))
                closest_peak = peaks[sums.index(min(sums))] #find which peak has the smallest summed distance

            else:
                continue

            rows[en][-1] = str(closest_peak/fs) #update the corresponding row's pick for csv

        #write updated info to csv
        for row in rows:
            with open(csv_name, 'a', newline='') as file:
                writer = csv.writer(file)
                writer.writerow(row)
                file.close()


    if one_peak.count('y')/len(one_peak) <0.5:

        for en,peaks in enumerate(temps_p):
            if len(peaks)>1:
                rows[en][-1] = 'UNCERTAIN'

        #write updated info to csv
        for row in rows:
            with open(csv_name, 'a', newline='') as file:
                writer = csv.writer(file)
                writer.writerow(row)
                file.close()
    time1 = time()
    
    print(f'{time1-time0} seconds for this cluster')
        
#         break
#     break


# In[ ]:




