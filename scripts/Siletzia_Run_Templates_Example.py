# import everything

import numpy as np
import obspy
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
import matplotlib.pyplot as plt
import h5py
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


volc_md = pd.read_csv('Volcano_Metadata.csv') # read metadata file to create dataframe of labels

# create lists of stations used at each volcano/for each data file

Baker_sta = volc_md[volc_md['Volcano_Name'] == 'Mt_Baker']['Station'].values.tolist()
Hood_sta = volc_md[volc_md['Volcano_Name'] == 'Mt_Hood']['Station'].values.tolist() # missing from Volcano_Metadata.csv
St_Helens_sta = volc_md[volc_md['Volcano_Name'] == 'Mt_St_Helens']['Station'].values.tolist()
Newberry_sta = volc_md[volc_md['Volcano_Name'] == 'Newberry']['Station'].values.tolist() # missing from Volcano_Metadata.csv
Rainier_sta = volc_md[volc_md['Volcano_Name'] == 'Mt_Rainier']['Station'].values.tolist()

#create list of volcanoes

volc_list_names = ['Baker','Hood','Newberry','Rainier','St_Helens'] # list of names of each volcano
volc_sta = [Baker_sta,Hood_sta,Newberry_sta,Rainier_sta,St_Helens_sta] # lists of stations connected to respective volcanoes



# to replace the volcano loop,
vv = #read directory
# vv = 3 #if running multiple volcano scripts at the same times in the same directory, just set for each script
v = volc_sta[vv]

# loops through julian days
# run over a year

# for vv,v in enumerate(volc_sta): #can be replaced with v and vv above

for s in range(11,len(v)): #usually 0,len(v)
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
        obspy.read(glob(f'{datadir}{year}/*/{station}/{station}.*.{year}.*')[0]).select(channel=chan)
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
        st = obspy.read(*glob(f'{datadir}{year}/*/{sta}/{sta}.*.{year}.{str(i).zfill(3)}')).select(channel=chan)
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
            party = T.detect(stream=st,starttime=st[0].stats.starttime,endtime=st[0].stats.endtime,threshold=thr, threshold_type=thr_t,xcorr_func = xfunc,trig_int=trig_int, plot=plot, return_stream=r_st, ignore_bad_data=i_b_d,overlap=overlap)
        except:
            continue
        party.decluster(metric=metric,trig_int=trig_int) #had to add trig_int, it is minimum detection separation in seconds
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




