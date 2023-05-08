#!/usr/bin/env python
# coding: utf-8

# # currently set up for siletzia, change data pathway for cascadia
# 
# Created: July 14, 2022
# Updated: October 25, 2022

# In[2]:


# import everything
import yaml
import numpy as np
import obspy
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
import matplotlib.pyplot as plt
import h5py
# %matplotlib inline
# import sys
# sys.path.append("/data/wsd01/pnwstore/")
import eqcorrscan
from eqcorrscan.core.match_filter import match_filter
from eqcorrscan.core.match_filter.tribe import Tribe
from time import time
import csv
import pandas as pd
from glob import glob
from obspy.core.utcdatetime import UTCDateTime
client = Client('IRIS')


# In[3]:


with open('/home/smocz/expand_redpy/scripts/config_Helens.yaml') as file:
    config = yaml.load(file, Loader=yaml.FullLoader)
# Template metadata
fqmin = config['fqmin']
fqmax = config['fqmax']
fs = config['fs']
prepick_len = config['prepick_len']
trig_int = config['trig_int']
thr_t = config['thr_t']
thr = config['thr']
xfunc = config['xfunc']
plot = config['plot']
r_st = config['r_st']
i_b_d = config['i_b_d']
overlap = config['overlap']
metric = config['metric']
cores = config['cores']

homedir = config['homedir']
datadir = config['datadir']
readdir = config['readdir']

# year = config['year']
years = config['years']
vv = config['vv']
chan = config['chan']

# tribe = eqcorrscan.core.match_filter.tribe.Tribe(templates = stack_templates)


# In[4]:


volc_md = pd.read_csv(readdir+'Volcano_Metadata.csv') # read metadata file to create dataframe of labels
volc_md['netsta'] = volc_md['Network'].astype(str)+'.'+volc_md['Station'].astype(str)
# create lists of stations used at each volcano/for each data file

Baker_sta = volc_md[volc_md['Volcano_Name'] == 'Baker']['netsta'].values.tolist()
Hood_sta = volc_md[volc_md['Volcano_Name'] == 'Hood']['netsta'].values.tolist() # missing from Volcano_Metadata.csv
St_Helens_sta = volc_md[volc_md['Volcano_Name'] == 'St_Helens']['netsta'].values.tolist()
Newberry_sta = volc_md[volc_md['Volcano_Name'] == 'Newberry']['netsta'].values.tolist() # missing from Volcano_Metadata.csv
Rainier_sta = volc_md[volc_md['Volcano_Name'] == 'Rainier']['netsta'].values.tolist()

#create list of volcanoes

volc_list_names = ['Baker','Hood','Newberry','Rainier','St_Helens'] # list of names of each volcano
volc_sta = [Baker_sta,Hood_sta,Newberry_sta,Rainier_sta,St_Helens_sta] # lists of stations connected to respective volcanoes
# print(volc_sta)


# ### Run the templates over a year - loop through multiple tribes - integrate saving as a csv - Oct 25, 2022

# Save as .py files to run on terminal - divide up into each volcano

# In[ ]:
t00 = time()

# loops through julian days
# run over several years
for year in years:
    print(f"YEAR: {year}")
    # for vv,v in enumerate(volc_sta):
    v = volc_sta[vv]
    for s in range(0,len(v)): #usually 0,len(Hood_sta)
        net, sta =  v[s].split('.') #add specific network per station
        print(net,sta)
        print(volc_list_names[vv])
        try:
            T = Tribe().read(*glob(homedir+'templates/Volcano_'+volc_list_names[vv]+'_Network_'+net+'_Station_'+sta+'_Channel_*.tgz'))
            print(T)
        except:
            print('No tgz for Station')
            continue
#         try:
#             obspy.read(glob(f'{datadir}{year}/{net}/{sta}/{sta}.{net}.{year}.*')[0]).select(channel=chan)
#         except:
#             print('No Data for Station')
#             continue
        if UTCDateTime(volc_md[volc_md['netsta']==v[s]]['Starttime'].values.tolist()[0]).year <= year: #if the year of the startime of this station is less than or equal to the year we are running
            if UTCDateTime(volc_md[volc_md['netsta']==v[s]]['Endtime'].values.tolist()[0]).year >= year: #if the year of the endtime is greater than or equal to the year we are running
                with open(homedir+'detections/'+volc_list_names[vv]+'_'+v[s]+'_'+str(year)+'_detections.csv', 'w', newline='') as file: # if both are true, create csv
                    writer = csv.writer(file)
                    writer.writerow(["ID", "Template_Name", "Detection_Time"])
                    file.close()
            else: #if endtime statement is NOT true
                continue #skip this station for this year
        else: #if starttime statement is NOT true
            continue #skip this station for this year
        for i in range(1,366): # normally 1,366
            print('------')
            parties = []
            t0=time()
#             st = obspy.read(*glob(f'{datadir}{year}/*/{sta}/{sta}.{net}.{year}.{str(i).zfill(3)}')).select(channel=chan)
            sst = UTCDateTime(year=year,julday=i,hour=0,minute=0,second=0)
            try:
                st = client.get_waveforms(network=net, station=sta, channel=chan, location='*',starttime=sst, endtime=sst+86400)
            except:
                print(f'No data available for request. {net}.{sta} day {i}')
                continue            
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
                party = T.detect(stream=st,starttime=st[0].stats.starttime,endtime=st[-1].stats.endtime,threshold=thr, threshold_type=thr_t,xcorr_func = xfunc,trig_int=trig_int, plot=plot, return_stream=r_st, ignore_bad_data=i_b_d,overlap=overlap,cores=cores)
            except:
                print('can\'t detect')
                continue
            party.decluster(metric=metric,trig_int=trig_int) #had to add trig_int, it is minimum detection separation in seconds
            t2=time()
            print(f"it tooks {t2-t1} s to launch the party")
            print(party)
            print('detections:',len(party))
            if len(party) > 0: 
                print(party[0])
                print(party.families)
                parties.append(party)
                for ii in range(0,len(parties[0].families)):
                    for iii in range(0,len(parties[0].families[ii].detections)):
                        row = [str(parties[0].families[ii].detections[iii].id),
                               str(parties[0].families[ii].detections[iii].template_name),
                               str(parties[0].families[ii].detections[iii].detect_time)]
                        with open(homedir+'detections/'+volc_list_names[vv]+'_'+v[s]+'_'+str(year)+'_detections.csv', 'a',
                                  newline='') as file:
                            writer = csv.writer(file)
                            writer.writerow(row)
                            file.close()
    #         break
    #    break
t01 = time()
print(f'number of days to run detections at this volcano: {(t01-t00)/86400}')


# In[ ]:




