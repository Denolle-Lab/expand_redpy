import yaml
import numpy as np
import math
import obspy
from obspy import UTCDateTime
from obspy import Trace
from obspy.clients.fdsn import Client
from obspy.signal.trigger import classic_sta_lta, plot_trigger, trigger_onset
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec #for organizing subplots
from tqdm import trange
                        # t = trange(len(waveforms), desc="Trimming Waveforms ", leave=True)
                        #    for i in t:
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

tt0=time()

with open('/home/smocz/expand_redpy/scripts/config.yaml') as file:
    config = yaml.load(file, Loader=yaml.FullLoader)

fs = config['fs']
tb = config['tb']
ta = config['ta']
fqmin = config['fqmin']
fqmax = config['fqmax']
nbucket = config['nbucket']
chan = config['chan']
homedir = config['homedir']
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

volc_md['netsta'] = volc_md['Network'].astype(str)+'.'+volc_md['Station'].astype(str)
Baker_sta = volc_md[volc_md['Volcano_Name'] == 'Baker']['netsta'].values.tolist()
Hood_sta = volc_md[volc_md['Volcano_Name'] == 'Hood']['netsta'].values.tolist() 
St_Helens_sta = volc_md[volc_md['Volcano_Name'] == 'St_Helens']['netsta'].values.tolist()
Newberry_sta = volc_md[volc_md['Volcano_Name'] == 'Newberry']['netsta'].values.tolist() 
Rainier_sta = volc_md[volc_md['Volcano_Name'] == 'Rainier']['netsta'].values.tolist()

volc_list = [Baker,Hood,Newberry,Rainier,St_Helens] # list of dataframes for each volcano
volc_list_names = ['Baker','Hood','Newberry','Rainier','St_Helens'] # list of names of each volcano
volc_sta = [Baker_sta,Hood_sta,Newberry_sta,Rainier_sta,St_Helens_sta] # lists of stations connected to respective volcanoes

# for vv,v in enumerate(volc_sta): #vv is the number in the list, v is the station list for current volcano
v = volc_sta[vv]
clid = volc_list[vv]['Clustered'].values.tolist() #find the largest cluster ID for a volcano to set range
for s in range(0,len(v)): #loop through stations
    norm_stack_list = [] #list of normalized stacks for every cluster
    norm_stack_names = [] #list of names for the normalized stacks
    
    net, sta =  v[s].split('.') #add specific network per station
    t0 = time() #record time
    cid = [] #cid will become a list of cluster IDs that have templates
    st3=obspy.Stream() #st3 will become a stream of traces, each trace being the template/stack for a cluster

    for cl in range(0,(clid[-1]+1)): #cl is cluster ID
    # clid[-1] is the highest cluster ID, and add 1 so that the last cluster id is not skipped
        t2=time()
        sst=obspy.Stream() #sst will have every trace for each cluster
        print('------------') #divider for clarity
        print('Cluster ID: '+str(cl)+' Volcano: '+volc_list_names[vv]+' Station: '+sta+' Network: '+net) #keeping track of what is currently running
        date_list = volc_list[vv][volc_list[vv]['Clustered'] == cl]['datetime'].values.tolist()
        for ii,i in enumerate(date_list): #i is each datetime from cl at this volcano
            stt=UTCDateTime(i)-tb-nlta #starttime
            et=UTCDateTime(i)+ta #endtime
            utct=UTCDateTime(i) #datetime from REDpy catalog
            
            sta_st = UTCDateTime(volc_md[volc_md['netsta']==v[s]]['Starttime'].values[0])
            sta_nd = UTCDateTime(volc_md[volc_md['netsta']==v[s]]['Endtime'].values[0])
            
            if UTCDateTime(i)<sta_st or UTCDateTime(i)>sta_nd:
#                 print(f'{i} before sta start {sta_st} or after sta end {sta_nd}')
                continue
            
#             print(i)
            
            ### PULL WAVEFORM ###
            
            st = client.get_waveforms(network=net,station=sta,location='*',
                                      channel=chan, year=utct.year, doy=utct.julday)
            try:
                st = st.detrend(type = 'demean')
                st.filter(type='bandpass',freqmin=fqmin,freqmax=fqmax)
                st.resample(fs) #get same sampling rate for all events
                st.trim(starttime=stt,endtime=et)
                st.merge(fill_value = 0)
            except:
                print('could not process')
                continue

            ### CHECK DATA ###
            
            if len(st)==0:
                print('no data for this time')
                continue
#             print(len(st[0].data))
            if len(st[0].data) == round((ta+tb+nlta)*fs)+1: #if there is enough data to contribute to making a stack
                sst.append(st[0])
    #         break
    
        ### SHIFT TRACES IN CL ###
    
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
            print('before alignment:',shift,cc) #the shift and cross correlation values before shifting/alignment. Perfect shift is 0, perfect cross correlation is 1
            if shift<0: st4[i].data[:shift]=st2[i].data[-shift:]
            if shift>0: st4[i].data[shift:]=st2[i].data[:-shift]
            xcor = obspy.signal.cross_correlation.correlate(st4[i].data[:],sst[0].data[:],200)
            index = np.argmax(xcor)
            cc = round(xcor[index],9)
            shift = 200-index
            print('after alignment',shift,cc) #shift and cross correlation values after alignment. shift should be 0
        print('shift complete')
        print( )
        
        
        sst2=st4.copy() #sst2 is a copy of st4 once alignment is finished
        sst3=st4.copy() #another copy for normalization
        sst2.stack() #stack aligned waveforms
        
        sst3 = sst3.normalize()
        sst3.stack()
        
        if len(sst)<2: #some clusters will have 1 or 0 traces due to unavailable data
            print('sst not long enough')
            t4=time()
            print(f'{str(t4-t2)} seconds to attempt to find {len(date_list)} waveforms')
            continue

        print('length of sst: '+str(len(sst))) #should be 2 or higher
        st3.append(sst2[0])

        cid.append('rp'+volc_list_names[vv][:2].lower()+str(cl).zfill(len(str(clid[-1])))) #append the clusterID to cid. rp for REDpy, volc_list_names[vv][:2].lower() for the volcano name
        # writing seisbench h5py

        lat = volc_md[volc_md['netsta']==v[s]]['Latitude'].values[0]
        lon = volc_md[volc_md['netsta']==v[s]]['Longitude'].values[0]
        
        ### PLOTTING ###
        
#         if len(sst) <=4:
#             height = len(sst)
#         if len(sst)>4 and len(sst)<100:
#             height = math.ceil(len(sst)/2)
#         if len(sst) >=100:
#             height = math.ceil(len(sst)/5)
        
#         fig, ax0 = plt.subplots(figsize=(6,height))
#         gs = GridSpec(height, 2, figure=fig) #make GridSpec for formatting subplots, based on number of waveforms
#         ax1 = fig.add_subplot(gs[0:height,0:1])
#         ax2 = fig.add_subplot(gs[0:1,1:2])
#         ax3 = fig.add_subplot(gs[1:2,1:2])
        
#         yscale = 2 #how far to space waveforms from eachother
#         wavecolor = 'black'
#         for ww, wave in enumerate(sst):
#             ax1.plot(wave.data[:]/np.max(np.abs(wave.data))+yscale+(yscale*ww),color=wavecolor,linewidth=.5)
#         ax1.tick_params(axis='both', which='both', bottom=False,left=False,labelbottom=False,labelleft=False)

            
#         #plot un-normalized stack
#         ax2.plot(sst2[0].data[:],color='red',linewidth=.5)
#         ax2.set_title('Stack (no normalization)')
#         ax2.tick_params(axis='x', which='both', bottom=False,labelbottom=False)

        
#         #plot normalized stack
#         ax3.plot(sst3[0].data[:],color='orange',linewidth=.5)
#         ax3.set_title('Normalized Stack')
#         ax3.tick_params(axis='x', which='both', bottom=False,labelbottom=False)
        
#         fig.suptitle(f'Cluster {cl} at {v[s]} on {volc_list_names[vv]}')
#         fig.set_tight_layout(True)
#         fig.delaxes(ax0) #remove unused ax
#         fig.savefig(f'{homedir}stack_plots/Volcano_{volc_list_names[vv]}_netsta_{v[s]}_cl_{str(cl).zfill(len(str(clid[-1])))}.svg')
        
        ### SNR on the normalized stack (sst3) ###

        cft = classic_sta_lta(sst3[0].data, int(nsta * fs), int(nlta * fs)) #basis for stalta
#         print('-------------')
#         print(sst3[0].stats.starttime+nlta)
#         plot_trigger(sst3[0], cft, thr_on, thr_off) #print stalta plots
#         print('C ID:', cid[ll])
        on_off = np.array(trigger_onset(cft, thr_on, thr_off)) #x value of vlines in the plot (plot_trigger is NOT needed for this to work)
        if not np.any(on_off):
            print('NO SIGNAL FOUND')
            continue
        index = []
        for ii,i in enumerate(on_off): #for each signal detection in on_off
        #     print(i)
            if i[1]-i[0] == 0: #if the start and end times are same, mark
                index.append(ii)
        on_off = np.delete(on_off,index, axis=0) #remove marked signal windows by index
        print('on_off')
        print(on_off)
        amps = []
        for ii,i in enumerate(on_off): #go through each possible signal window
            start = i[0]
            end = i[1]
            amps.append(Trace(data=sst3[0].data[start:end]).max()) #find max amplitude in that signal window
        # print(amps)
        amp = max(amps, key=abs) #find max amplitude of stream within sta/lta windows
        # print(amp)
        for ii,i in enumerate(on_off): #go through each possible signal window
            start = i[0]
            end = i[1]
            ooamp = Trace(data=sst3[0].data[start:end]).max() #find max amplitude in that signal window
            print('amp',amp,'ooamp',ooamp,'start',start,'end',end)
            if amp == ooamp: #if max amplitude of this window is the same as the overall max amplitude
#                 print('wow')
                trig_on = round(float(start / fs),4) #took out -nlta
                trig_off = round(float(end / fs),4) #get start and end times
                break #move on
        # show trigger on and off times, rounded to 4 decimal places
        print('Trigger on',trig_on-nlta,'seconds after tb time') #trig_on-nlta is picktime...?
        print('Trigger off',trig_off-nlta,'seconds after tb time')

        signal_window = sst3[0].copy()
        noise_window = sst3[0].copy()
        

        signal_window.trim(starttime=sst3[0].stats.starttime+trig_on,endtime=sst3[0].stats.starttime+trig_off)
        #i+trig_on-0.5 to include lead up to the signal?
        noise_window.trim(starttime=sst3[0].stats.starttime,endtime=sst3[0].stats.starttime+nlta) 
        #the noise window is the part that can't be considered in the sta/lta anyway

        snr = 20 * np.log(np.percentile(np.abs(signal_window.data),pr) 
                          / np.percentile(np.abs(noise_window.data),pr))/np.log(10)
        print('snr:',snr)
        if snr<snr_t:
            print('SNR too low')
            continue
            
        ### APPENDING ###
        
        #append normalized waveform data to a list
        norm_stack_list.append(sst3[0].data[:])
        norm_stack_names.append(f'{net}_{sta}_rp{volc_list_names[vv][:2].lower()}_{cl}')
                
#         break
    ### Save Station info to h5 ###
    with h5py.File(f"{homedir}h5/normalized_{volc_list_names[vv].lower()}_templates_{net}_{sta}.h5", "w") as f:
        f.create_dataset("waveforms", data=np.array(norm_stack_list))
        f.create_dataset("template_name", data=norm_stack_names)
        
#     break

tt1 = time()
print(f'{(tt1-tt0)/60} minutes to get normalized stacks and plots for {volc_list_names[vv]}')
