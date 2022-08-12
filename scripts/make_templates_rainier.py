#!/usr/bin/env python
# coding: utf-8

# # Make Templates from REDpy Clusters using EQCorrScan
# Created: Jul 14, 2022
# 
# Updated: Aug 8, 2022

# ### Set Up

# Import everything

# In[ ]:


import yaml
import numpy as np
import obspy
from obspy import UTCDateTime
from obspy import Trace
from obspy.clients.fdsn import Client
from obspy.signal.trigger import classic_sta_lta, plot_trigger, trigger_onset
import matplotlib.pyplot as plt
import h5py
# %matplotlib inline
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

# In[ ]:


with open('/home/smocz/expand_redpy/scripts/config_rainier.yaml') as file:
    config = yaml.load(file, Loader=yaml.FullLoader)

fs = config['fs']
tb = config['tb']
ta = config['ta']
fqmin = config['fqmin']
fqmax = config['fqmax']
nbucket = config['nbucket']
chan = config['chan']
homedir = config['homedir']
datadir = config['datadir']
readdir = config['readdir']
f_o = config['f_o']
p_len = config['p_len']

vv = config['vv']

snr_t = config['snr_t']
thr_on = config['thr_on']
thr_off = config['thr_off']
nsta = config['nsta']
nlta = config['nlta']
pr = config['pr']

print(vv)


# Read REDpy Catalogs and Volcano Metadata File

# In[ ]:


Baker = pd.read_csv(readdir+'Baker_catalog.csv')
Hood = pd.read_csv(readdir+'Hood_catalog.csv')


St_Helens = pd.read_csv(readdir+'MountStHelens_catalog.csv')

# Combining borehole and local catalogs with St_Helens

Helens_Borehole = pd.read_csv(readdir+'MSHborehole_catalog.csv')
Helens_Borehole['Clustered'] += 2000 
# Cluster 0 in Helens_Borehole is now Cluster 2000 in St_Helens
Helens_Local = pd.read_csv(readdir+'MSHlocal_catalog.csv')
Helens_Local['Clustered'] += 3000
# Cluster 0 in Helens_Local is now Cluster 3000 in St_Helens

# Use St_Helens to access all three St Helens catalogs
St_Helens = pd.concat([St_Helens,Helens_Borehole,Helens_Local])

Newberry = pd.read_csv(readdir+'Newberry_catalog.csv')
Rainier = pd.read_csv(readdir+'Rainier_catalog.csv')

volc_md = pd.read_csv(readdir+'Volcano_Metadata.csv')
# read metadata file to create dataframe of labels


# Associate networks and stations

# In[ ]:


volc_md['netsta'] = volc_md['Network'].astype(str)+'.'+volc_md['Station'].astype(str)


# Create Lists of Stations for Each Volcano Using volc_md

# In[ ]:


Baker_sta = volc_md[volc_md['Volcano_Name'] == 'Baker']['netsta'].values.tolist()
Hood_sta = volc_md[volc_md['Volcano_Name'] == 'Hood']['netsta'].values.tolist() 
St_Helens_sta = volc_md[volc_md['Volcano_Name'] == 'St_Helens']['netsta'].values.tolist()
Newberry_sta = volc_md[volc_md['Volcano_Name'] == 'Newberry']['netsta'].values.tolist() 
Rainier_sta = volc_md[volc_md['Volcano_Name'] == 'Rainier']['netsta'].values.tolist()


# Create Lists of Volcano Information

# In[ ]:


volc_list = [Baker,Hood,Newberry,Rainier,St_Helens] # list of dataframes for each volcano
volc_list_names = ['Baker','Hood','Newberry','Rainier','St_Helens'] # list of names of each volcano
volc_sta = [Baker_sta,Hood_sta,Newberry_sta,Rainier_sta,St_Helens_sta] # lists of stations connected to respective volcanoes


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

# changed most (all?) .stats.station to sta

# for vv,v in enumerate(volc_sta): #vv is the number in the list, v is the station list for current volcano
v = volc_sta[vv]
clid = volc_list[vv]['Clustered'].values.tolist() #find the largest cluster ID for a volcano to set range
print(len(v))
for s in range(0,len(v)): #loop through stations
    net, sta =  v[s].split('.') #add specific network per station
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
        print('Cluster ID: '+str(cl)+' Volcano: '+volc_list_names[vv]+' Station: '+sta+' Network: '+net) #keeping track of what is currently running
        for ii,i in enumerate(volc_list[vv][volc_list[vv]['Clustered'] == cl]['datetime']): #i is each datetime from cl at this volcano
            print(i)
            stt=UTCDateTime(i)-tb-nlta #starttime
            et=UTCDateTime(i)+ta #endtime
            utct=UTCDateTime(i) #datetime from REDpy catalog
            try:
                st = client.get_waveforms(network=net,station=sta,location='*',
                                          channel=chan, year=utct.year, doy=utct.julday)
                st = st.detrend(type = 'demean')
                st.filter(type='bandpass',freqmin=fqmin,freqmax=fqmax)
                st.resample(fs) #get same sampling rate for all events
                st.trim(starttime=stt,endtime=et)
                st.merge(fill_value = 0)
            except:
                print('Could not find waveform')
                continue #if no waveform can be found, move onto next time
            if len(st)==0:
                print('no data for this time')
                continue
            print(len(st[0].data))
            if len(st[0].data) == round((ta+tb+nlta)*fs)+1: #if there is enough data to contribute to making a stack
                sst.append(st[0])
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

        lat = volc_md[volc_md['Station'] == sta]['Latitude'].values[0]
        lon = volc_md[volc_md['Station'] == sta]['Longitude'].values[0]
       # fill in metadata 
        meta = meta.append({"source_id": cid[-1], "source_origin_time": '', 
            "source_latitude_deg": "%.3f" % 0, "source_longitude_deg": "%.3f" % 0, 
            "source_type": 'unknown',
            "source_depth_km": "%.3f" % 0, "source_magnitude": 0,
            "station_network_code": net, "trace_channel": sst[0].stats.channel, 
            "station_code": sta, "station_location_code": sst[0].stats.location,
            "station_latitude_deg": lat,  "station_longitude_deg": lon,
            "station_elevation_m": 0,
            "trace_p_arrival_sample": 0, "CODE": sta.lower()+'bhz'+cid[-1]}, ignore_index = True)

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
    stack_list = [] #will be a list of streams,each stream has one stack as a trace
    for i in range(0,len(st3)):
        temp = obspy.Stream(st3[i]) #give each trace in st3 its own stream
        stack_list.append(temp) #append each stream to stack_list
        print('appending stack list')
    print(stack_list)

    print('got stack list!')

    #filter for snr with temp list, remember to keep cid in order
    cidu = [] #updated cid
    temp_list = []
    for ll,l in enumerate(stack_list):
        cft = classic_sta_lta(l[0].data, int(nsta * fs), int(nlta * fs))
        print('-------------')
        print(l[0].stats.starttime+nlta)
        plot_trigger(l[0], cft, thr_on, thr_off)
        print('C ID:', cid[ll])
        on_off = np.array(trigger_onset(cft, thr_on, thr_off))
        if not np.any(on_off):
            print('NO SIGNAL FOUND')
            continue
        index = []
        for ii,i in enumerate(on_off):
        #     print(i)
            if i[1]-i[0] == 0:
                index.append(ii)
        on_off = np.delete(on_off,index, axis=0)
        print('on_off')
        print(on_off)
        amps = []
        for ii,i in enumerate(on_off): #go through each possible signal window
            start = i[0]
            end = i[1]
            amps.append(Trace(data=l[0].data[start:end]).max()) #find max amplitude in that signal window
        # print(amps)
        amp = max(amps, key=abs) #find max amplitude of stream within sta/lta windows
        # print(amp)
        for ii,i in enumerate(on_off): #go through each possible signal window
            start = i[0]
            end = i[1]
            ooamp = Trace(data=l[0].data[start:end]).max() #find max amplitude in that signal window
            print('amp',amp,'ooamp',ooamp,'start',start,'end',end)
            if amp == ooamp: #if max amplitude of this window is the same as the overall max amplitude
                print('wow')
                trig_on = round(float(start / fs),4) #took out -nlta
                trig_off = round(float(end / fs),4) #get start and end times
                break #move on
        # show trigger on and off times, rounded to 4 decimal places
        print('Trigger on',trig_on-nlta,'seconds after tb time') #trig_on-nlta is picktime...?
        print('Trigger off',trig_off-nlta,'seconds after tb time')

        signal_window = l[0].copy()
        noise_window = l[0].copy()

        signal_window.trim(starttime=l[0].stats.starttime+trig_on,endtime=l[0].stats.starttime+trig_off)
        #i+trig_on-0.5 to include lead up to the signal
        noise_window.trim(starttime=l[0].stats.starttime,endtime=l[0].stats.starttime+nlta)

        snr = 20 * np.log(np.percentile(np.abs(signal_window.data),pr) 
                          / np.percentile(np.abs(noise_window.data),pr))/np.log(10)
        print('snr:',snr)
        if snr>=snr_t: #if SNR is greater than or equal to minimum, save it
            l.trim(starttime=l[0].stats.starttime+nlta,endtime=l[0].stats.endtime)
            cidu.append(cid[ll])
            temp_list.append(l)
    

    template_stream = [] #will be a list of streams that have been filtered to become templates
    for i in range(0,len(temp_list)):
        template_stream.append(eqcorrscan.core.match_filter.template.Template(
            name=sta.lower()+temp_list[i][0].stats.channel.lower()+cidu[i],st=temp_list[i], \
                lowcut=fqmin, highcut=fqmax, samp_rate=fs, filt_order=f_o, process_length=p_len))
        #filter each stream from temp_list and append to template_stream
        print('appending template stream')
    print('made template stream!')

    tribe = eqcorrscan.core.match_filter.tribe.Tribe(templates = template_stream)
    #make a tribe from template_stream, each template will be made from a stream in template_stream
#     print('making tribe...')

    if len(st3) > 0:
        print('Writing tribe in progress...')
        tribe.write(homedir+'/templates/Volcano_' + volc_list_names[vv] + '_Network_' +                     net + '_Station_' + sta + '_Channel_' + st3[0].stats.channel)
#         print('Wrote tribe sucessfully!')
        # add to file
        meta.to_csv("/data/whd02/Data_rp/metadata_"+volc_list_names[vv]+"_"+net+"_"+sta+".csv",sep = ',', index=False)
        f = h5py.File("/data/whd02/Data_rp/waveforms_"+volc_list_names[vv]+"_"+net+"_"+sta+".hdf5",'a') #appending mode
            #If the file does not exist, it creates a new file for writing.
        # need to define f in order to close it in order to open it in mode w
        if f: f.close()
        f = h5py.File("/data/whd02/Data_rp/waveforms_"+volc_list_names[vv]+"_"+net+"_"+sta+".hdf5",'w') #writing mode
#         f['/data_format/component_order'] ='ZNE' #dunno what this does
        print(range(nbucket))
        for b in range(nbucket):
            f['/data/bucket%d' % (b + 1)] = data[b + 1]
        f.close()

        print('Saved!')
# BREAK FOR THE STATION LOOP
#     break
# BREAK FOR THE VOLCANO LOOP
#     break


# In[ ]:




