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
import matplotlib.pyplot as plt

import tsfel
import random
from datetime import timedelta
import calendar
from tsfel import time_series_features_extractor

#%config InlineBackend.figure_format = "png"

#from Feature_Extraction import compute_hibert

import warnings

# Ignore all warnings
warnings.filterwarnings("ignore")

# displaying all columns from pandas dataframe
# Set display options to show all columns
pd.set_option('display.max_columns', None)

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

# vv = config['vv']
vv=3


print(volc_list_names[vv])


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

print('first template:',all_temps[0])

#color map for plotting
def get_cmap(n, name='viridis'): #hsv
#     Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
#     RGB color; the keyword argument name must be a standard mpl colormap name.
    return plt.cm.get_cmap(name, n)

#get clusterid from template name
def getcl_id(t_name_str): #for normalized
    t_cl = int(t_name_str.split('_')[-1])
    return t_cl

#get station name from template name
def getnet_sta(t_name_str): #for normalized
    t_net = t_name_str.split('_')[0]
    t_sta = t_name_str.split('_')[1]
    return t_net,t_sta

eq_Z = all_waves

cfg_file = tsfel.get_features_by_domain()

# Extract features for earthquakes
features_eqz = pd.DataFrame([])
for i in range(len(all_temps)+1):
    try:
    
        df = time_series_features_extractor(cfg_file, eq_Z[i], fs=40,)
        df['template'] = all_temps[i]
        features_eqz = pd.concat([features_eqz,df])
        
    except:
        pass
    
features_eqz.to_csv(f'{homedir}{volc_list_names[vv]}_tsfel_features.csv',index=False)