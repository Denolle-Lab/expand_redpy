#imports
from eqcorrscan import Tribe
import obspy
from obspy import Stream
import csv
import numpy as np
import pandas as pd
from glob import glob
import scipy.signal as sp
import yaml
import matplotlib.pyplot as plt

#set parameters
with open('/home/smocz/expand_redpy/scripts/config.yaml') as file:
        config = yaml.load(file, Loader=yaml.FullLoader)
            
volc_list_names = config['volc_list_names']
vv = config['vv']
volc = volc_list_names[vv]
print(volc)
filename = f'Volcano_{volc}_Network_*_Station_*_Channel_*.tgz' #name of .tgz files to glob
homedir = config['homedir']
path = homedir+'templates/'

winlen = 40

fs = config['fs']
len_data = 10000
fqmin = config['fqmin']
fqmax = config['fqmax']

# spectrogram parameters
sgramMode='magnitude'
sgramScaling='spectrum'

# frequency/time resolution
nperseg = 64 #affects time resolution
noverlap = nperseg/4
nfft = 1024

def maxwindow(st,winlen): #provide stream object (one trace each) and window length in s
    fs = st[0].stats.sampling_rate
    stt = st[0].stats.starttime

    numwindows = round(len(st[0])/fs)-winlen+1 #how many windows will be cycled (+1 is for range)
    print('number of windows:',numwindows-1)

    sum_list = []
    for i in range(0,numwindows):
        tr = st[0].copy() #create copy to avoid trimming st[0]
        tt = tr.trim(stt+i,stt+i+winlen) #trim to a window the length of winlen

        amps_list = []
        for ii in tt: #for each amplitude/point in the trace
            amps_list.append(abs(ii)) #append the amplitude 
        sum_list.append(sum(amps_list)) #sum and save the amplitude for this possible window

    max_inx = sum_list.index(max(sum_list)) #find the index of the window with maximum sum(abs(amplitudes))

    print('maximum amp window index:',max_inx) #

    final_tr = st[0].copy()
    tr_trim = final_tr.trim(stt+max_inx,stt+max_inx+winlen) #use index to find the window
    print('Max Window:')
    tr_trim.plot();

    return(Stream(tr_trim))
