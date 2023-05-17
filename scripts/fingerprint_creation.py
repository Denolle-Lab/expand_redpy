#imports
from eqcorrscan import Tribe #reading tgz files with templates
import obspy
from obspy import Stream, Trace
import csv #for reading and saving data
import numpy as np
import pandas as pd
from glob import glob #looking through files
import scipy.signal as sp
import yaml #for config file
import matplotlib.pyplot as plt #for plotting
from tqdm import trange
from time import time #for time for code to run
import h5py

from specufex import BayesianNonparametricNMF, BayesianHMM #nmf and hmm functions
from sklearn.cluster import KMeans #kmeans clustering function

from matplotlib.gridspec import GridSpec #for plotting clusters

#set parameters
with open('/home/smocz/expand_redpy/scripts/config.yaml') as file:
    config = yaml.load(file, Loader=yaml.FullLoader)
    
volc_list_names = config['volc_list_names']
vv = config['vv']
volc = volc_list_names[vv]
print(volc)
filename = f'Volcano_{volc}_Network_*_Station_*_Channel_*.tgz' #name of .tgz files to glob
homedir = config['homedir']
print(homedir)
path = homedir+'templates/'

winlen = config['winlen']

fs = config['fs']
fqmin = config['fqmin']
fqmax = config['fqmax']

# spectrogram parameters
sgramMode='magnitude'
sgramScaling='spectrum'

# frequency/time resolution
nperseg = config['nperseg'] 
noverlap = nperseg/config['dnoverlap']
nfft = config['nfft']

#Kmeans
K=config['K'] # number of clusters to fit

#NMF
batches_nmf = config['batches_nmf']
batch_size = config['batch_size']

#HMM
num_states = config['num_states']
batches_hmm = config['batches_hmm']

savedir = '/home/smocz/expand_redpy_new_files/h5/' #directory to save h5 file in
h5name = 'curated_set.h5'

##########################################################################################

# create window with max possible sum(abs(amplitudes))
# goal is to align stream windows and provide shorter windows without as much noise
def maxwindow(st,winlen,fs): #provide stream object (one trace each), window length in s, and sampling rate
    fs_ = st[0].stats.sampling_rate #if sampling rate is in metadata, this can read it
    if fs_ > 1.: fs=fs_
    else: st[0].stats.sampling_rate = fs
    stt = st[0].stats.starttime #if starttime is in metadata, this can read it

    print(len(st[0]))
    print(fs)
    numwindows = round(len(st[0])/fs)-winlen+1 #how many windows will be cycled (+1 is for range)
    print('number of windows:',numwindows-1)

    sum_list = [] #list of sum of abs value of amplitudes for each window
    for i in range(0,numwindows): #i is index of windows
        tw = st[0].copy() #create copy to avoid trimming st[0]
        tt = tw.trim(stt+i,stt+i+winlen) #trim to a window the length of winlen
#         print(stt+i)

        amps_list = [] #list of absolute value amplitudes
        for ii in tt: #for each amplitude/point (ii) in the trace (tt)
            amps_list.append(abs(ii)) #append the amplitude 
        sum_list.append(sum(amps_list)) #sum and save the amplitude for this possible window
    
    max_sec = sum_list.index(max(sum_list)) #find the index of the window with maximum sum(abs(amplitudes))
    
    print('maximum amp window index:',max_sec) #print max_sec or the # second after beginning of template that 
    #the maximum window occurs
    
    final_tw = st[0].copy() #make a copy for trimming
    
    #adding 3s before to account for pwave/beginning of waveform
    if max_sec-3 <= 0: #if 3s before max_sec is 0 or negative
        sec = 0 #sec is 0
    else: #if 3s before max_sec is >0
        sec = max_sec-3 #sec is 3s before max_sec
    
    tw_trim = final_tw.trim(stt+max_sec,stt+max_sec+winlen) #use sec to trim the window corretly
    print('Max Window:') #show the maximum window
    tw_trim.plot();
    
    return(Stream(tw_trim)) #return the stream of the maximum window

##########################################################################################

with h5py.File(f"{savedir}{h5name}", "r") as f: #read file
    waveforms = f["waveforms"][()] #pull in waveforms
    template_name = f["template_name"].asstr()[()] #pull in waveform id
    group_id = f["group_id"][()] #pull in group id

print(f"{len(waveforms)} waveforms in file")

#shorten window length around max amplitude

waveforms_n = []
for wave in waveforms:
    wave_ = Trace(wave)
    wave_.plot();
    new_wave = maxwindow(Stream(traces=[wave_]),winlen,fs)
    waveforms_n.append(np.array(new_wave[0].data))
#     break


tdf = pd.DataFrame({"template_name":template_name,"waveform":list(waveforms_n)})
tdf.head()

##########################################################################################

#calculate raw spectrograms with scipy
fSTFT, tSTFT, STFT_raw = sp.spectrogram(
    x=np.stack(tdf["waveform"].values),
    fs=fs,
    nperseg=nperseg,
    noverlap=noverlap,
    nfft=nfft,
    scaling=sgramScaling,
    axis=-1,
    mode=sgramMode
)

#quality check for NaN
np.isnan(STFT_raw).any()

#bandpass

freq_slice = np.where((fSTFT >= fqmin) & (fSTFT <= fqmax))
fSTFT   = fSTFT[freq_slice]
STFT_0 = STFT_raw[:,freq_slice,:].squeeze()

# STFT_0 = STFT_raw

#normalize spectrogram and make values nonnegative for nmf

normConstant = np.median(STFT_0, axis=(1,2))
STFT_norm = STFT_0 / normConstant[:,np.newaxis,np.newaxis]  # norm by median
del STFT_0
STFT_dB = 20*np.log10(STFT_norm, where=STFT_norm != 0) # convert to dB (database?)
del STFT_norm
STFT = np.maximum(0, STFT_dB) # make sure nonnegative
del STFT_dB

tdf["stft"] = list(STFT)
tdf.head()

#quality check
bad_idx = tdf["stft"][tdf["stft"].apply(lambda x: np.isnan(x).any())].index
print(f"Bad spectrograms: \n{tdf.loc[bad_idx].template_name}")
tdf = tdf.drop(bad_idx).sort_values("template_name")

nmf = BayesianNonparametricNMF(np.stack(tdf["stft"].values).shape)
print(np.stack(tdf["stft"].values).shape)

t = trange(batches_nmf, desc="NMF fit progress ", leave=True)
for i in t:
    idx = np.random.randint(len(tdf["stft"].values), size=batch_size)
    nmf.fit(tdf["stft"].iloc[idx].values)
    t.set_postfix_str(f"Patterns: {nmf.num_pat}")
    break
    
# #get activation matrix (Vs) from nmf
# Vs = nmf.transform(tdf["stft"].values)

# hmm = BayesianHMM(nmf.num_pat, nmf.gain, num_state=num_states, Neff=50000)

# t = trange(batches_hmm, desc="HMM fit progress ", leave=True)
# for i in t:
#     idx = np.random.randint(Vs.shape[0], size=1)
#     hmm.fit(Vs[idx])
    
# fingerprints, As, gams = hmm.transform(Vs) #create fingerprints

# ##########################################################################################

# with h5py.File(f"{savedir}new_{h5name}", "w") as f:
#     f.create_dataset("waveforms", data=waveforms_n)
#     f.create_dataset("group_id", data=group_id)
#     f.create_dataset("template_name", data=template_name)
#     f.create_dataset("fingerprints", data=fingerprints)