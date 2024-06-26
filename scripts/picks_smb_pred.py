
import requests
import h5py
import yaml
import csv
import math
import eqcorrscan
from eqcorrscan import Tribe
from time import time
import obspy
from obspy import UTCDateTime, Trace
import pandas as pd
from glob import glob
import numpy as np
from obspy.signal.cross_correlation import *
import matplotlib.pyplot as plt
from geopy import distance
from tqdm import trange

import torch
plt.rcParams.update({'font.size': 10})
# from utils import *


import seisbench.models as sbm
device = torch.device("cpu")

from ELEP.elep.ensemble_statistics import ensemble_statistics
from ELEP.elep.ensemble_coherence import ensemble_semblance 
from ELEP.elep.ensemble_learners import ensemble_regressor_cnn
from mbf_elep_func import apply_mbf
from ELEP.elep import mbf, mbf_utils
from ELEP.elep import trigger_func
from ELEP.elep.trigger_func import picks_summary_simple

from ELEP.elep.mbf_utils import make_LogFq, make_LinFq, rec_filter_coeff, create_obspy_trace
from ELEP.elep.mbf import MB_filter as MBF


# In[2]:


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

tape = 0.05 #percent of waveform to taper on either side

# istart = config['nlta']*fs
# print(istart)


print(volc_list_names[vv])


# In[3]:


volc_md = pd.read_csv(readdir+'Volcano_Metadata.csv')
#associate network and station
volc_md['netsta'] = volc_md['Network'].astype(str)+'.'+volc_md['Station'].astype(str)


#get clusterid from template name
def getcl_id(t_name_str): #for normalized
    t_cl = int(t_name_str.split('_')[-1])
    return t_cl


# In[5]:


# ML picker parameters
paras_semblance = {'dt':0.025, 'semblance_order':4, 'window_flag':True, 
                   'semblance_win':0.5, 'weight_flag':'max'}
p_thrd, s_thrd = 0.01, 0.05

# download models
pretrain_list = ["pnw","ethz","instance","scedc","stead","geofon","neic"]
pn_pnw_model = sbm.EQTransformer.from_pretrained('pnw')
pn_ethz_model = sbm.EQTransformer.from_pretrained("ethz")
pn_instance_model = sbm.EQTransformer.from_pretrained("instance")
pn_scedc_model = sbm.EQTransformer.from_pretrained("scedc")
pn_stead_model = sbm.EQTransformer.from_pretrained("stead")
pn_geofon_model = sbm.EQTransformer.from_pretrained("geofon")
pn_neic_model = sbm.EQTransformer.from_pretrained("neic")

#list of models to run through
list_models = [pn_pnw_model,pn_ethz_model,pn_scedc_model,pn_neic_model,pn_geofon_model,pn_stead_model,pn_instance_model]

pn_pnw_model.to(device); #imodel 0
pn_ethz_model.to(device); #imodel 1
pn_scedc_model.to(device); #imodel 2
pn_neic_model.to(device); #imodel 3
pn_geofon_model.to(device); #imodel 4
pn_stead_model.to(device); #imodel 5
pn_instance_model.to(device); #imodel 6


# Find picktimes

# In[6]:


#pull in h5 for volcano

all_temps = []
all_waves = []

for filepath in glob(f'/home/smocz/expand_redpy_new_files/h5/normalized_{volc_list_names[vv].lower()}_templates_*.h5'):
    net = filepath.split('_')[-2]
    with h5py.File(filepath, "r") as f: #pull in fingerprints
        template_name = f["template_name"][()]
        waveforms = f["waveforms"][()]
#         print(f.keys()) #print what data is in this file
    [all_temps.append(i) for i in template_name]
    [all_waves.append(i) for i in waveforms]
    
all_waves = np.array(all_waves)
all_temps = [str(i)[2:-1] for i in all_temps]


# In[7]:


#get some info
v = volc_md[volc_md['Volcano_Name'] == volc_list_names[vv]]['netsta'].values.tolist() #list of network and station per volc
clid = np.unique([getcl_id(i) for i in all_temps]) #list of cluster ids
cllen = len(str(max(clid))) #length of the largest cluster ID, used for zfill
zz = chan[-2:].lower() #the last two letters of channel names (essentially the letters in chan)
csv_name = f'{homedir}locations/{volc_list_names[vv]}_ELEP_normalized_picktimes.csv' #name of the csv for picktimes at this volcano
h5_name = f'{homedir}h5/{volc_list_names[vv]}_ELEP_smb_pred.h5' #name of h5 file for smb_pred

# In[8]:


# create csv for volcano
with open(csv_name, 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Network','Station','Cluster_ID','Template_Name','SMB_peak']) #,'SMB_peak_MBF'
    file.close()




### PULL IN TEMPLATES ###
# cl_trange = trange(max(clid), desc="Finding picktimes for each cluster", leave=True)
cl_trange=range(0,max(clid)+1) #max(clid)+1
for cl in cl_trange:
#     print('------') #print a divider
#     print("cluster:",str(cl).zfill(cllen)) #print the cluster ID
    
    temps_w = [] #list of templates waveforms for this cluster
    temps_n = [] #list of template names
    temps_p = [] #list of template picks (list of arrays)
    preds = [] #pred data for looking
    
    stopwatch0=time() #note the time
    
    s_trange = trange(len(v), desc=f"Finding picktimes for each station at cluster {cl}", leave=True) #has a progress bar
#     s_trange = range(0,len(v)) #no progress bar

    rows = [] #list of rows to append to csv
    for s in s_trange: #loop through stations that have a template for this cluster
        net, sta =  v[s].split('.') #add specific network per station
        
        for tt,t in enumerate(all_temps): #go through each template
            if t.split('_')[0]==net and t.split('_')[1]==sta and t.split('_')[-1]==str(cl): #if the template is for this station and cluster
                
                ### PREPARE DATA ###
#                 Trace(all_waves[tt]).plot(); #show before preparation
                wave = all_waves[tt].copy() #create copy
                
                t_tapered = Trace(wave).taper(tape) #make waveform a trace to taper it
                padded_wave = np.hstack((t_tapered.data[:],np.zeros(2480))) #take data from tapered trace and pad end 
                #with zeros so that len(t_trace)=6000 and will fit in the nueral network
                t_trace = Trace(padded_wave,{'sampling_rate':fs}) #make back into a Trace, and set sampling rate
#                 print(len(t_trace))
                
                temps_w.append(t_trace) #append trace
                temps_n.append(t) #append name
                
#                 t_trace.plot(); #plot trace after preparation

    
                ### FIND PICKS! ###

                #picking params
                evt_data = Stream(traces=[t_trace])
                sta_available = [sta]
                list_models = list_models
                twin = len(t_trace)-1


                dt = 1/fs; fs = fs
                nfqs = 5
                nt = 6000; nc = 3
                iend = len(wave)-1 #end of possible pick times
                istart = math.ceil(tape*(len(wave)-1)) #beginning of possible pick times (currently excluding taper)
                
                fq_list = make_LogFq(fqmin, fqmax, dt, nfqs)
                coeff_HP, coeff_LP = rec_filter_coeff(fq_list, dt)
                MBF_paras = {'f_min':fqmin, 'f_max':fqmax, 'nfqs':nfqs, 'frequencies':fq_list, 'CN_HP':coeff_HP, 'CN_LP':coeff_LP,                     'dt':dt, 'fs':fs, 'nt':nt, 'nc':nc, 'npoles': 2}

                paras_semblance = {'dt':dt, 'semblance_order':2, 'window_flag':True, 
                                   'semblance_win':0.5, 'weight_flag':'max'}
                
                #find picktimes!
                peaks,smb_pred = apply_mbf(evt_data, sta_available, list_models, MBF_paras, paras_semblance, istart, iend) #smb_peak,smb_peak_mbf,
                
                
                print(peaks[0]) # THESE PICKS ARE OFFSET BY ISTART
                
                picks = [i+istart for i in peaks[0]] # ADDED ISTART
                temps_p.append(picks)
                csv_picks=' '.join([str(i/fs)for i in picks]) #formatted for saving in csv (i/fs = pick time in seconds)
                
                preds.append([smb_pred])
        
                #append picks to csv
                row = [net, sta, cl, t, csv_picks] ### APPEND PICKS NOT PEAKS
                rows.append(row)

    if len(temps_n)==0:
        print(f'no templates for cluster {cl}')
        continue
                
    ### FILTERING PEAKS ###

    one_peak = [] #list of n and y if the template only has one peak/picktime
    for p in temps_p: #for each template's picks
        if len(p)==1: #if only one possible pick
            one_peak.append('y')
        else: 
            one_peak.append('n')
            
    if one_peak.count('y')/len(one_peak) ==1: #if there is only one peak per template for all templates in this cluster
        #write info to csv as is
        for row in rows:
            with open(csv_name, 'a', newline='') as file:
                writer = csv.writer(file)
                writer.writerow(row)
                file.close()
            
    if one_peak.count('y')/len(one_peak) >=0.5 and one_peak.count('y')/len(one_peak) <1: #if 50-75% of templates have one peak
        
        one_p_value = [] #list of peak values for templates with only one possible peak AKA "confirmed" peaks
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
                
                
    if one_peak.count('y')/len(one_peak) <0.5: #if less than 50% of templates in a cluster have "confirmed" peaks
        
        for en,peaks in enumerate(temps_p): #for each list of peaks in the cluster
            if len(peaks)>1: #if more than one possible peak
                rows[en][-1] = 'UNCERTAIN' #label the pick for that template waveform as "uncertain"
                
        #write updated info to csv
        for row in rows:
            with open(csv_name, 'a', newline='') as file:
                writer = csv.writer(file)
                writer.writerow(row)
                file.close()
                
    #save smb_pred to h5
    try:
        open(h5_name)
        with h5py.File(h5_name, "a") as f: #append info
            f.create_dataset(f"smb_pred_cl_{cl}", data=np.array(preds))
            f.create_dataset(f"template_names_cl_{cl}", data=temps_n)
    except:
        with h5py.File(h5_name, "w") as f: #create file
            f.create_dataset(f"smb_pred_cl_{cl}", data=np.array(preds))
            f.create_dataset(f"template_names_cl_{cl}", data=temps_n)
    
#     break

