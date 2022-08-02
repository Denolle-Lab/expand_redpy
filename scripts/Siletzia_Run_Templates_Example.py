#!/usr/bin/env python
# coding: utf-8

# In[6]:


# import everything

import numpy as np
import obspy
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
import matplotlib.pyplot as plt
import h5py
import sys
sys.path.append("/data/wsd01/pnwstore/")
import eqcorrscan
from eqcorrscan.core.match_filter import match_filter
from eqcorrscan.core.match_filter.tribe import Tribe
from time import time
import csv
import pandas as pd
from glob import glob

from obspy.core.utcdatetime import UTCDateTime


# In[7]:


volc_md = pd.read_csv('Volcano_Metadata.csv') # read metadata file to create dataframe of labels

# create lists of stations used at each volcano/for each data file

Baker_sta = volc_md[volc_md['Volcano_Name'] == 'Mt_Baker']['Station'].values.tolist()
Hood_sta = volc_md[volc_md['Volcano_Name'] == 'Mt_Hood']['Station'].values.tolist() # missing from Volcano_Metadata.csv
St_Helens_sta = volc_md[volc_md['Volcano_Name'] == 'Mt_St_Helens']['Station'].values.tolist()
Newberry_sta = volc_md[volc_md['Volcano_Name'] == 'Newberry']['Station'].values.tolist() # missing from Volcano_Metadata.csv
Rainier_sta = volc_md[volc_md['Volcano_Name'] == 'Mt_Rainier']['Station'].values.tolist()

#create list of volcanoes

volc_list_names = ['Baker','Hood','St_Helens','Newberry','Rainier'] # list of names of each volcano
volc_sta = [Baker_sta,Hood_sta,St_Helens_sta,Newberry_sta,Rainier_sta] # lists of stations connected to respective volcanoes


# In[8]:


#Filter based on distance and then stack all earthquakes in a group to create a template

# Template metadata
fqmin = 1
fqmax = 10
fs=40
prepick_len = 0.3
trig_int = 6

homedir = '/home/smocz/redpy_expand_new_files/' #home directory or directory to save new files to
datatdir = '/home/smocz/2019_data/' #directory to get data from

year = 2019

# tribe = eqcorrscan.core.match_filter.tribe.Tribe(templates = stack_templates)


# ### Run the templates over a year - loop through multiple tribes - integrate saving as a csv - Jul 14, 2022

# Save as .py files to run on terminal - divide up into each volcano

# In[ ]:

# ### Running one test template

# Station HOOD on Hood for cluster 17

# In[5]:



# In[ ]:


# loops through julian days
# run over a year

# for vv,v in enumerate(volc_sta):
v = Hood_sta
vv = 1
for s in range(0,7): #usually 0,len(Hood_sta)
    station = v[s]
    print(station)
    print(volc_list_names[vv])
    try:
        T = Tribe().read(*glob(homedir+'/templates/Volcano_'+volc_list_names[vv]+'_Station_'+station+'_Channel_*.tgz'))
        print(T)
    except:
        print('No Tribe for Station')
        continue
    try:
        obspy.read(glob(datatdir+str(year)+'/*/'+station+'.*.'+str(year)+'.*')[0])
    except:
        print('No Data for Station')
        continue
    with open(homedir+'/detections/'+volc_list_names[vv]+'_'+v[s]+'_'+str(year)+'_detections.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["ID", "Template_Name", "Detection_Time"])
        file.close()
    for i in range(1,366): # normally 1,366
        parties = []
        t0=time()
        st = obspy.read(*glob(datatdir+str(i).zfill(3)+'/'+station+'.*.'+str(year)+'.*'))
        st.detrend(type='demean')
        st.resample(fs)
        st.filter(type='bandpass',freqmin=fqmin,freqmax=fqmax)
        st.merge(fill_value=0)
        t1=time()
        print("it tooks %2f s to download data" %(t1-t0))
        print(st)
        print(str(year)+str(i).zfill(3))
        for ii in range(0,len(T.templates)):
            T.templates[ii].prepick = prepick_len 
        if len(st)==0: continue
        try:
            party = T.detect(stream=st,starttime=st[0].stats.starttime,endtime=st[0].stats.endtime,threshold=0.6, threshold_type='absolute',xcorr_func = 'fmf',trig_int=trig_int, plot=False, return_stream=False, ignore_bad_data=True,overlap='calculate')
        except:
            continue
        party.decluster(metric='avg_cor',trig_int=trig_int) #had to add trig_int, it is minimum detection separation in seconds
        t2=time()
        print("it tooks %2f s to launch the party" %(t2-t1))
        print(party)
        print(len(party))
        if len(party) > 0: 
            print(party[0])
            print(party.families)
            parties.append(party)
            for ii in range(0,len(parties[0].families)):
                for iii in range(0,len(parties[0].families[ii].detections)):
                    row = [str(parties[0].families[ii].detections[iii].id),
                           str(parties[0].families[ii].detections[iii].template_name),
                           str(parties[0].families[ii].detections[iii].detect_time)]
                    with open(homedir+'/detections/'+volc_list_names[vv]+'_'+v[s]+'_'+str(year)+'_detections.csv', 'a',
                              newline='') as file:
                        writer = csv.writer(file)
                        writer.writerow(row)
                        file.close()
#         break
#    break


# In[ ]: