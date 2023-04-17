#imports and parameters
import h5py #to save and read h5 files
import numpy as np #for arrays and math
import pandas as pd #for dataframes and csv
import scipy.signal as sp #for calculating spectrograms
import matplotlib.pyplot as plt #for plots
from sklearn.cluster import KMeans #for grouping

from specufex import BayesianNonparametricNMF, BayesianHMM #for NMF and HMM
from tqdm import trange #for a progress meter

homedir = '/home/smocz/expand_redpy_new_files/templates/' #directory to with .tgz file
filename = 'Volcano_Rainier_Network_UW_Station_RCM_Channel_HHZ.tgz' #name of .tgz file
savedir = '/home/smocz/expand_redpy_new_files/h5/' #directory to save h5 file in
h5name = '' #name of h5 file

#stream filtering
fs = 40 #sampling rate
fqmin = 1 #minimum frequency for bandpass
fqmax = 10 #maximum frequency for bandpass

# spectrogram parameters
sgramMode='magnitude'
sgramScaling='spectrum'

# frequency/time resolution
nperseg = 64
noverlap = nperseg/4
nfft = 1024

#NMF
batches_nmf = 100000
batch_size = 1

#HMM
num_states = 8
batches_hmm = 5000

data_files = ["specufex_data_full_tgz.h5"] #list of data file names "specufex_data_selected.h5","specufex_data.h5"

############################

for filename in data_files: #for each h5 file
    #load waveforms from h5 file
    with h5py.File(f"/home/smocz/expand_redpy/notebooks/tutorials/{filename}", "r") as f:
        waveforms = f["waveforms"][()]
        template_name = f["template_name"].asstr()[()]

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
    
    #frequency bandpass
    freq_slice = np.where((fSTFT >= fqmin) & (fSTFT <= fqmax))
    fSTFT   = fSTFT[freq_slice]
    STFT_0 = STFT_raw[:,freq_slice,:].squeeze()
    
    #normalize
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
    

    #nmf
    nmf = BayesianNonparametricNMF(np.stack(tdf["stft"].values).shape)

    t = trange(batches_nmf, desc="NMF fit progress ", leave=True)
    for i in t:
        idx = np.random.randint(len(tdf["stft"].values), size=batch_size)
        nmf.fit(tdf["stft"].iloc[idx].values)
        t.set_postfix_str(f"Patterns: {nmf.num_pat}")

    #show activation matrix
    Vs = nmf.transform(tdf["stft"].values)
    
    hmm = BayesianHMM(nmf.num_pat, nmf.gain, num_state=num_states, Neff=50000)

    t = trange(batches_hmm, desc="HMM fit progress ", leave=True)
    for i in t:
        idx = np.random.randint(Vs.shape[0], size=1)
        hmm.fit(Vs[idx])
        
    fingerprints, As, gams = hmm.transform(Vs)


    # save data to an hdf5
    with h5py.File(f"/home/smocz/expand_redpy/notebooks/tutorials/{filename}", "w") as f:
        f["Vs"] = Vs
        f["waveforms"] = waveforms
        f["template_name"] = template_name
        f["fingerprints"] = fingerprints
        f["hmm_EB"] = hmm.EB