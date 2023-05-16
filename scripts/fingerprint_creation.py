#imports
from eqcorrscan import Tribe #reading tgz files with templates
import obspy
from obspy import Stream
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

with h5py.File(f"{savedir}{h5name}", "r") as f:
    waveforms = f["waveforms"][()]
    template_name = f["template_name"].asstr()[()]
    group_id = f["group_id"][()]

print(f"{len(waveforms)} waveforms in file")

tdf = pd.DataFrame({"template_name":template_name,"waveform":list(waveforms)})
tdf.head()

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

t = trange(batches_nmf, desc="NMF fit progress ", leave=True)
for i in t:
    idx = np.random.randint(len(tdf["stft"].values), size=batch_size)
    nmf.fit(tdf["stft"].iloc[idx].values)
    t.set_postfix_str(f"Patterns: {nmf.num_pat}")
    
#get activation matrix (Vs) from nmf
Vs = nmf.transform(tdf["stft"].values)

hmm = BayesianHMM(nmf.num_pat, nmf.gain, num_state=num_states, Neff=50000)

t = trange(batches_hmm, desc="HMM fit progress ", leave=True)
for i in t:
    idx = np.random.randint(Vs.shape[0], size=1)
    hmm.fit(Vs[idx])
    
fingerprints, As, gams = hmm.transform(Vs) #create fingerprints


with h5py.File(f"{savedir}{h5name}", "w") as f:
    f.create_dataset("waveforms", data=waveforms)
    f.create_dataset("group_id", data=group_id_list)
    f.create_dataset("template_name", data=name_list)
    f.create_dataset("fingerprints", data=fingerprints)