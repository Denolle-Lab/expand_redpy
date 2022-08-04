#!/usr/bin/env python
# coding: utf-8

# # Make Templates from REDpy Clusters using EQCorrScan
# Created: Jul 14, 2022
# 
# Updated: Jul 14, 2022

# ### Set Up

# Import everything

# In[1]:


import numpy as np
import obspy
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
import matplotlib.pyplot as plt
import h5py
get_ipython().run_line_magic('matplotlib', 'inline')
import sys
sys.path.append("/data/wsd01/pnwstore/")
import eqcorrscan
from eqcorrscan.core.match_filter import match_filter
from eqcorrscan.core.match_filter.tribe import Tribe
import pandas as pd

from time import time

from pnwstore.mseed import WaveformClient
from obspy.core.utcdatetime import UTCDateTime

client = WaveformClient()
client2 = Client('IRIS')


# Parameters

# In[2]:


fs = 40 #sampling rate in hertz
tb = 3 #time in seconds before origin time
ta = 10 #time in seconds after origin time
fqmin = 1 #minimum frequency for bandpass filter
fqmax = 10 #maximum frequency for bandpass filter
nbucket = 1 # breaking up matrix, if having keyerror for nbucket, just set it to 1
chan = '*HZ' #channel to get waveforms from


# Read REDpy Catalogs and Volcano Metadata File

# In[3]:


Baker = pd.read_csv('Baker_catalog.csv')
Hood = pd.read_csv('Hood_catalog.csv')


St_Helens = pd.read_csv('MountStHelens_catalog.csv')

# Combining borehole and local catalogs with St_Helens

Helens_Borehole = pd.read_csv('MSHborehole_catalog.csv')
Helens_Borehole['Clustered'] += 2000 
# Cluster 0 in Helens_Borehole is now Cluster 2000 in St_Helens
Helens_Local = pd.read_csv('MSHlocal_catalog.csv')
Helens_Local['Clustered'] += 3000
# Cluster 0 in Helens_Local is now Cluster 3000 in St_Helens

# Use St_Helens to access all three St Helens catalogs
St_Helens = pd.concat([St_Helens,Helens_Borehole,Helens_Local])

Newberry = pd.read_csv('Newberry_catalog.csv')
Rainier = pd.read_csv('Rainier_catalog.csv')

volc_md = pd.read_csv('Volcano_Metadata.csv')
# read metadata file to create dataframe of labels


# Create Lists of Stations for Each Volcano Using volc_md

# In[4]:


Baker_sta = volc_md[volc_md['Volcano_Name'] == 'Mt_Baker']['Station'].values.tolist()
Hood_sta = volc_md[volc_md['Volcano_Name'] == 'Mt_Hood']['Station'].values.tolist() 
St_Helens_sta = volc_md[volc_md['Volcano_Name'] == 'Mt_St_Helens']['Station'].values.tolist()
Newberry_sta = volc_md[volc_md['Volcano_Name'] == 'Newberry']['Station'].values.tolist() 
Rainier_sta = volc_md[volc_md['Volcano_Name'] == 'Mt_Rainier']['Station'].values.tolist()


# Create Lists of Volcano Information

# In[5]:


volc_list = [Baker,Hood,St_Helens,Newberry,Rainier] # list of dataframes for each volcano
volc_list_names = ['Baker','Hood','St_Helens','Newberry','Rainier'] # list of names of each volcano
volc_sta = [Baker_sta,Hood_sta,St_Helens_sta,Newberry_sta,Rainier_sta] # lists of stations connected to respective volcanoes


# ### Making the Templates

# Define update_data

# In[ ]:


def update_data(data, streamdata, ibucket):
    streamdata = np.expand_dims(streamdata, axis = 0)

    if ibucket in data:
        data[ibucket] = np.concatenate((data[ibucket], streamdata), axis = 0)
    else:
        data[ibucket] = streamdata
    return data


# Make and Save Templates

# In[ ]:


# QUESTION: how often should data and meta reset? Every volcano? 
# currently the seisbench is named by volcano and station - so I think we could put it between t0 and cid
# changed most (all?) .stats.station to v[s]

for vv,v in enumerate(volc_sta): #vv is the number in the list, v is the station list for current volcano
    clid = volc_list[vv]['Clustered'].values.tolist() #find the largest cluster ID for a volcano to set range
    for s in range(0,len(v)): #loop through stations
        t0 = time() #record time
        #Create Dictionary data and DataFrame meta
        data = {}
        meta = pd.DataFrame(columns = [
            "source_id", "source_origin_time", "source_latitude_deg", "source_longitude_deg", "source_type",
            "source_depth_km", "split", "source_magnitude", "station_network_code", "trace_channel", 
            "station_code", "station_location_code", "station_latitude_deg",  "station_longitude_deg",
            "station_elevation_m", "trace_name", "trace_sampling_rate_hz", "trace_start_time",
            "trace_S_arrival_sample", "trace_P_arrival_sample", "CODE"])
        
        cid = [] #cid will become a list of cluster IDs that have templates
        st3=obspy.Stream() #st3 will become a stream of traces, each trace being the template/stack for a cluster
        for cl in range(0,(clid[-1]+1)): #cl is cluster ID
        # clid[-1] is the highest cluster ID, and add 1 so that the last cluster id is not skipped
            t2=time()
            sst=obspy.Stream() #sst will have every trace for each cluster
            print('------------') #divider for clarity
            print('Cluster ID: '+str(cl)+' Volcano: '+volc_list_names[vv]+' Station: '+v[s]) #keeping track of what is currently running
            for ii,i in enumerate(volc_list[vv][volc_list[vv]['Clustered'] == cl]['datetime']): #i is each datetime from cl at this volcano
                stt=UTCDateTime(i)-tb #starttime
                et=UTCDateTime(i)+ta #endtime
                utct=UTCDateTime(i) #datetime from REDpy catalog
                try:
                    st = client.get_waveforms(network='*',station=v[s],location='*',channel=chan, year=utct.year, doy=utct.julday)
                    st = st.detrend(type = 'demean')
                    st.filter(type='bandpass',freqmin=fqmin,freqmax=fqmax)
                    st.resample(fs) #get same sampling rate for all events
                    st.trim(starttime=stt,endtime=et)
                    st.merge(fill_value = 0)
                    if len(st[0].data) >= round((ta+tb)*fs)-2: #if there is enough data to contribute to making a stack
                        sst.append(st[0])

                except:
                    pass
        #         break
            print(sst) #see how many traces were found
            st2=sst.copy() #copy of sst after appending is finished/all waveforms for a cluster are gathered, reference for shifting
            st4=st2.copy() #st4 will become the aligned version of st2
            for i in range(0, len(st2)):
                print('Working on shifting') # Will not print if len(st2)==0
                # Also serves as a divider between shift and cc values for different waveforms
                xcor = obspy.signal.cross_correlation.correlate(st2[i].data[:],sst[0].data[:],200)
                index = np.argmax(xcor)
                cc = round(xcor[index],9)
                shift = 200-index
                print(shift,cc) #the shift and cross correlation values before shifting/alignment. Perfect shift is 0, perfect cross correlation is 1
                if shift<0: st4[i].data[:shift]=st2[i].data[-shift:]
                if shift>0: st4[i].data[shift:]=st2[i].data[:-shift]
                xcor = obspy.signal.cross_correlation.correlate(st4[i].data[:],sst[0].data[:],200)
                index = np.argmax(xcor)
                cc = round(xcor[index],9)
                shift = 200-index
                print(shift,cc) #shift and cross correlation values after alignment. shift should be 0
            print('shift complete')
            print( )
            sst2=st4.copy() #sst2 is a copy of st4 once alignment is finished
            sst2.stack() #stack aligned waveforms
            if len(sst)<2: #some clusters will have 1 or 0 traces due to unavailable data
                print('sst not long enough')
                t4=time()
                print(str(t4-t2)+' seconds to attempt to find waveforms')
                continue

            print('length of sst: '+str(len(sst))) #should be 2 or higher
            st3.append(sst2[0])

            cid.append('rp'+volc_list_names[vv][:2].lower()+str(cl).zfill(len(str(clid[-1])))) #append the clusterID to cid. rp for REDpy, volc_list_names[vv][:2].lower() for the volcano name
            # writing seisbench h5py

            lat = volc_md[volc_md['Station'] == v[s]]['Latitude'].values[0]
            lon = volc_md[volc_md['Station'] == v[s]]['Longitude'].values[0]
           # fill in metadata 
            meta = meta.append({"source_id": cid[-1], "source_origin_time": '', 
                "source_latitude_deg": "%.3f" % 0, "source_longitude_deg": "%.3f" % 0, 
                "source_type": 'unknown',
                "source_depth_km": "%.3f" % 0, "source_magnitude": 0,
                "station_network_code": sst[0].stats.network, "trace_channel": sst[0].stats.channel, 
                "station_code": v[s], "station_location_code": sst[0].stats.location,
                "station_latitude_deg": lat,  "station_longitude_deg": lon,
                "station_elevation_m": 0,
                "trace_p_arrival_sample": 0, "CODE": v[s].lower()+'bhz'+cid[-1]}, ignore_index = True)

            # fill in data
            ibucket = np.random.choice(list(np.arange(nbucket) + 1))
            data = update_data(data, st3[-1], ibucket)
            print('ibucket: ',ibucket,data[ibucket])
            t3=time()
            print(str(t3-t2)+' seconds to make the stack for this cluster')
# BREAK FOR THE CLUSTER LOOP
    #         break
        t1=time()
        print(str(t1-t0)+'seconds to make stacks for one station for all clusters with enough data')
        # integrate tribe
        temp_list = [] #will be a list of streams,each stream has one stack as a trace
        for i in range(0,len(st3)):
            temp = obspy.Stream(st3[i]) #give each trace in st3 its own stream
            temp_list.append(temp) #append each stream to temp_list
            print('appending temp list')
        print(temp_list)

        print('got temp list!')


        template_stream = [] #will be a list of streams that have been filtered to become templates
        for i in range(0,len(temp_list)):
            template_stream.append(eqcorrscan.core.match_filter.template.Template(
                name=v[s].lower()+'bhz'+cid[i],st=temp_list[i], \
                    lowcut=fqmin, highcut=fqmax, samp_rate=fs, filt_order=4, process_length=200))
            #filter each stream from temp_list and append to template_stream
            print('appending template stream')
        print('made template stream!')

        tribe = eqcorrscan.core.match_filter.tribe.Tribe(templates = template_stream)
        #make a tribe from template_stream, each template will be made from a stream in template_stream
        print('making tribe...')

        if len(st3) > 0:
            print('Writing tribe in progress...')
            tribe.write('Volcano_' + volc_list_names[vv] + '_Station_' + v[s] + '_Channel_' + st3[0].stats.channel)
            print('Wrote tribe sucessfully!')
            # add to file
            meta.to_csv("/data/whd02/Data_rp/metadata_"+volc_list_names[vv]+"_"+v[s]+".csv",sep = ',', index=False)
            f = h5py.File("/data/whd02/Data_rp/waveforms_"+volc_list_names[vv]+"_"+v[s]+".hdf5",'a') #appending mode
                #If the file does not exist, it creates a new file for writing.
            # need to define f in order to close it in order to open it in mode w
            if f: f.close()
            f = h5py.File("/data/whd02/Data_rp/waveforms_"+volc_list_names[vv]+"_"+st3[0].stats.station+".hdf5",'w') #writing mode
            f['/data_format/component_order'] ='ZNE'
            print(range(nbucket))
            for b in range(nbucket):
                f['/data/bucket%d' % (b + 1)] = data[b + 1]
            f.close()
            
            #put saving cid here
            
            print('Saved!')
# BREAK FOR THE STATION LOOP
    #     break
# BREAK FOR THE VOLCANO LOOP
#     break

