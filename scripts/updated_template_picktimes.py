#!/usr/bin/env python
# coding: utf-8

# Created August 17, 2023
# 
# Purpose is to update location grid search in the following ways:
# 1. use templates whose constituent waveforms were normalized before stacking (fixes most false features)
# 2. use elep to find picktimes instead of envelope cross correlation
# 3. base velocity model of grid search on p and s picktimes when possible
# 4. compare locations to ComCat potential matches from the <a href="https://assets.pnsn.org/red/">redpy website</a> when possible

# In[1]:


import requests
import h5py
import yaml
import csv
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

vv = config['vv']


print(volc_list_names[vv])


# In[3]:


volc_md = pd.read_csv(readdir+'Volcano_Metadata.csv')
#associate network and station
volc_md['netsta'] = volc_md['Network'].astype(str)+'.'+volc_md['Station'].astype(str)


# In[4]:


#color map for plotting
def get_cmap(n, name='viridis'): #hsv
#     Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
#     RGB color; the keyword argument name must be a standard mpl colormap name.
    return plt.cm.get_cmap(name, n)

# define function to predict synthetic arrival times
def travel_time(t0, x, y, vs, sta_x, sta_y):
    dist = np.sqrt((sta_x - x)**2 + (sta_y - y)**2)
    tt = t0 + dist/vs
    return tt

# define function to compute residual sum of squares
def error(synth_arrivals,arrivals):
    res = arrivals - synth_arrivals   #make sure arrivals are in the right order, maybe iterate through keys
    res_sqr = res**2
    rss = np.sum(res_sqr)
    return rss

# define function to iterate through grid and calculate travel time residuals
def gridsearch(t0,x_vect,y_vect,sta_x,sta_y,vs,arrivals):
    rss_mat = np.zeros((len(t0),len(x_vect),len(y_vect)))
    rss_mat[:,:,:] = np.nan
    for i in range(len(t0)): 
        for j in range(len(x_vect)):
            for k in range(len(y_vect)):
                for m in range(len(vs)): #parameterize velocity
                    synth_arrivals = []
                    for h in range(len(sta_x)):
                        tt = travel_time(t0[i],x_vect[j],y_vect[k],vs[m],sta_x[h],sta_y[h]) 
                    #add vs in nested loop, vector 1000-5000, per cluster to account for p and s waves
                        synth_arrivals.append(tt)
                    rss = error(np.array(synth_arrivals),np.array(arrivals))
                    rss_mat[i,j,k] = rss
    return rss_mat

# define function to convert the location index into latitude and longitude
def location(x_dist, y_dist, start_lat, start_lon):
    bearing = 90-np.rad2deg(np.arctan(y_dist/x_dist))
    dist = np.sqrt((x_dist)**2 + (y_dist)**2)
    d = distance.geodesic(meters = dist)
    loc_lat = d.destination(point=[start_lat,start_lon], bearing=bearing)[0]
    loc_lon = d.destination(point=[start_lat,start_lon], bearing=bearing)[1]
    return loc_lat, loc_lon, d

# define function to find diameter in meters of the error on the location
def error_diameter(new_array):
    min_idx = np.min(new_array[:,1])
    max_idx = np.max(new_array[:,1])
    difference = max_idx-min_idx
    diameter_m = difference*1000
    return diameter_m 

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


# In[8]:


#create csv for volcano
with open(csv_name, 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Network','Station','Cluster_ID','Template_Name','SMB_peak','SMB_peak_MBF'])
    file.close()


# In[17]:


### PULL IN TEMPLATES ###
cl_trange = trange(max(clid), desc="Finding picktimes for each cluster", leave=True)
for cl in cl_trange:
#     print('------') #print a divider
#     print("cluster:",str(cl).zfill(cllen)) #print the cluster ID
    
    temps_w = [] #list of templates waveforms for this cluster
    temps_n = [] #list of template names
    temps_p = [] #list of template picks (list of arrays)
    
    stopwatch0=time() #note the time
    
    s_trange = trange(len(v), desc=f"Finding picktimes for each station at cluster {cl}", leave=True)
    for s in s_trange: #loop through stations that have a template for this cluster
        net, sta =  v[s].split('.') #add specific network per station
        
        for tt,t in enumerate(all_temps): #go through each template
            if t.split('_')[0]==net and t.split('_')[1]==sta and t.split('_')[-1]==str(cl): #if the template is for this station and cluster
                
                ### PREPARE DATA ###
#                 Trace(all_waves[tt]).plot(); #show before preparation
                wave = all_waves[tt].copy() #create copy
                
                t_tapered = Trace(wave).taper(0.05) #make waveform a trace to taper it
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
                fq_list = make_LogFq(fqmin, fqmax, dt, nfqs)
                coeff_HP, coeff_LP = rec_filter_coeff(fq_list, dt)
                MBF_paras = {'f_min':fqmin, 'f_max':fqmax, 'nfqs':nfqs, 'frequencies':fq_list, 'CN_HP':coeff_HP, 'CN_LP':coeff_LP,                     'dt':dt, 'fs':fs, 'nt':nt, 'nc':nc, 'npoles': 2}

                paras_semblance = {'dt':dt, 'semblance_order':2, 'window_flag':True, 
                                   'semblance_win':0.5, 'weight_flag':'max'}
                
                #find picktimes!
                smb_peak,smb_peak_mbf = apply_mbf(evt_data, sta_available,                 list_models, MBF_paras, paras_semblance)
        
                #append picks to csv
                row = [net, sta, cl, t, smb_peak[0], smb_peak_mbf[0]]
                with open(csv_name, 'a', newline='') as file:
                    writer = csv.writer(file)
                    writer.writerow(row)
                    file.close()
                
    
#         break
#     break