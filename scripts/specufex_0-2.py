#!/usr/bin/env python
# coding: utf-8

# Adapted from https://github.com/Specufex/specufex Specufex code
# 
# Created Feb 1, 2023
# Updated mAR 9, 2023

# In[ ]:


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

from specufex import BayesianNonparametricNMF, BayesianHMM #nmf and hmm functions
from sklearn.cluster import KMeans #kmeans clustering function


# In[ ]:


#set parameters
with open('/home/smocz/expand_redpy/scripts/config.yaml') as file:
    config = yaml.load(file, Loader=yaml.FullLoader)
    
volc_list_names = config['volc_list_names']

homedir = config['homedir']
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
K_list = [5,6,7,8]

#NMF
batches_nmf = config['batches_nmf']
batch_size = config['batch_size']

#HMM
num_states = config['num_states']
batches_hmm = config['batches_hmm']


# ## functions

# In[ ]:


# create window with max possible sum(abs(amplitudes))
# goal is to align stream windows and provide shorter windows without as much noise
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
    
    print('maximum amp window index:',max_inx) #print max_inx or the # second after beginning of template that 
    #the maximum window occurs
    
    final_tr = st[0].copy() #make a copy for trimming
    
    #adding 3s before to account for pwave/beginning of waveform
    if max_inx-3 <= 0: #if 3s before max_inx is 0 or negative
        inx = 0 #inx is 0
    else: #if 3s before max_inx is >0
        inx = max_inx-3 #inx is 3s before max_inx
    
    tr_trim = final_tr.trim(stt+max_inx,stt+max_inx+winlen) #use inx to trim the window corretly
    print('Max Window:') #show the maximum window
    tr_trim.plot();
    
    return(Stream(tr_trim)) #return the stream of the maximum window


# ### Create spectrogram from templates

# Create a dataframe of trimmed template waveforms

for vv in [0,1,2]:
    volc = volc_list_names[vv]
    print(volc)
    filename = f'Volcano_{volc}_Network_*_Station_*_Channel_*.tgz' #name of .tgz files to glob


    tdf = pd.DataFrame(data={'template_name':[],'waveform':[]}) #make a two column dataframe with labels

    for f_name in glob(f'{path}{filename}'): #loop thru tgz files at a volcano
        T = Tribe().read(f_name)

        for t in T: #loop thru tgz file
            print(f'{t.name} template original stream:') #show original stream

            t.st.plot();

            max_stream = maxwindow(st=t.st,winlen=winlen) #maxwindow func

            net = f_name.split('_')[6].lower() #f_name.split('_')[6].lower() is the network name

            tdf.loc[len(tdf)+1] = [f'{net}.{t.name}',max_stream[0].data] #add this trimmed
            #template to the df at index len+1

    #         max_stream[0].spectrogram() #show spectrogram
            print('---') #barrier

    #         break #only one template
    #     break #when break is on, only first station will be added


    # In[ ]:


    tdf_spare = tdf.copy() #create a copy of the unaltered df for easier testing
    # tdf = tdf_spare #make tdf refer to this copy


    # In[ ]:


    #test tdf
    tdf


    # In[ ]:


    #to see if trim length is not consistent - had some problems with trim at one point
    l_list = []
    for i in range(1,len(tdf['waveform'])):
        l = len(tdf['waveform'][i]) #number of data points in a waveform
        if l != (winlen*fs)+1: #if it is not 1601 (40fs times 40s plus 1)
            print('inx',i,'l=',l) #show the index and number of data points
            raise Exception('Some waveforms have not been properly trimmed') #raise an error if the trim lengths are not the same
    # print(np.unique(l_list))


    # Calculate spectrograms from trimmed template df

    # In[ ]:


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


    # In[ ]:


    #quality check
    np.isnan(STFT_raw).any()


    # In[ ]:


    #bandpass filter, optional
    freq_slice = np.where((fSTFT >= fqmin) & (fSTFT <= fqmax))
    fSTFT   = fSTFT[freq_slice]
    STFT_0 = STFT_raw[:,freq_slice,:].squeeze()


    # In[ ]:


    normConstant = np.median(STFT_0, axis=(1,2))
    STFT_norm = STFT_0 / normConstant[:,np.newaxis,np.newaxis]  # norm by median
    del STFT_0
    STFT_dB = 20*np.log10(STFT_norm, where=STFT_norm != 0) # convert to dB
    del STFT_norm
    STFT = np.maximum(0, STFT_dB) # make sure nonnegative
    del STFT_dB

    tdf["stft"] = list(STFT)
    tdf.head()


    # In[ ]:


    #quality check
    bad_idx = tdf["stft"][tdf["stft"].apply(lambda x: np.isnan(x).any())].index
    print(f"Bad spectrograms: \n{tdf.loc[bad_idx].template_name}")
    tdf = tdf.drop(bad_idx).sort_values("template_name")

    #this can mess with index sometimes (?)


    # In[ ]:


    #plotting example spectrogram
    n_spectrogram = 0 # index of spectrogram to plot

    f, ax = plt.subplots(1,2, figsize=(10,5))

    ax[0].pcolormesh(tSTFT,fSTFT,STFT_raw[n_spectrogram,freq_slice,:].squeeze())
    ax[0].set_xlabel("Timestep")
    ax[0].set_ylabel("Frequency (Hz)")
    ax[0].set_title("Original spectrogram")

    ax[1].pcolormesh(tSTFT,fSTFT, STFT[n_spectrogram])
    ax[1].set_xlabel("Timestep")
    ax[1].set_title("Normalized spectrogram")


    # In[ ]:


    #save spectrogram df as a csv
    #tdf.to_csv(f'{homedir}{volc}_{winlen}s_window_spectrograms.csv',index=False)

    #right now this is kind of useless, it to_csv doesn't save the full waveform 
    #or stft arrays, so they would have to be modified to be saved or saved in some other way


    # ### Run Specufex

    # NMF

    # In[ ]:


    nmf = BayesianNonparametricNMF(np.stack(tdf["stft"].values).shape)


    # In[ ]:


    t = trange(batches_nmf, desc="NMF fit progress ", leave=True)
    for i in t:
        idx = np.random.randint(len(tdf["stft"].values), size=batch_size)
        nmf.fit(tdf["stft"].iloc[idx].values)
        t.set_postfix_str(f"Patterns: {nmf.num_pat}")


    # In[ ]:


    plt.pcolormesh(nmf.EW@np.diag(nmf.EA[0]))
    plt.xlabel("NMF pattern number")
    plt.xticks(range(0,nmf.num_pat,2), range(0,nmf.num_pat,2))
    plt.ylabel("Frequency (Hz)")
    plt.show()


    # Activation Matrices

    # In[ ]:


    Vs = nmf.transform(tdf["stft"].values)
    # save Vs to an hdf5
    # with h5py.File("data/geysers/Vs.h5", "w") as f:
    #     f["Vs"] = Vs


    # In[ ]:


    f, axes = plt.subplots(1,2,figsize=(15,5))

    axes[0].pcolormesh(Vs[0])
    axes[0].set_xlabel("Timestep")
    axes[0].set_ylabel("NMF pattern")
    # axes[0].set_yticks(range(0,nmf.num_pat,2), range(0,nmf.num_pat,2))
    # axes[0].set_title("Activation matrix")

    # axes[1].pcolormesh(tdf["stft"].iloc[0])
    # axes[1].set_xlabel("Timestep")
    # axes[1].set_ylabel("Frequency")
    # axes[1].set_title("Normalized spectrogram")

    # plt.show()


    # HMM

    # In[ ]:


    hmm = BayesianHMM(nmf.num_pat, nmf.gain, num_state=num_states, Neff=50000)


    # In[ ]:


    t = trange(batches_hmm, desc="HMM fit progress ", leave=True)
    for i in t:
        idx = np.random.randint(Vs.shape[0], size=1)
        hmm.fit(Vs[idx])


    # In[ ]:


    plt.figure(figsize=(10,5))
    plt.imshow(hmm.EB, origin="lower")
    plt.ylabel("HMM state")
    plt.xlabel("Frequency pattern")
    _=plt.yticks(range(0,num_states,5), range(0, num_states,5))


    # Fingerprints

    # In[ ]:


    fingerprints, As, gams = hmm.transform(Vs)


    # In[ ]:


    plt.imshow(fingerprints[0])


    # ### Cluster with Kmeans

    # In[ ]:


    # convert fingerprints from 2D array to 1D array
    fingerprints_ = fingerprints.reshape((fingerprints.shape[0], fingerprints.shape[1]**2))

    for K in K_list:
        #Predicted labels
        y_pred = KMeans(n_clusters=K, random_state=42).fit_predict(fingerprints_)


        # ### Saving the clustering results

        # In[ ]:


        print((y_pred))


        # In[ ]:


        for i in np.unique(list(y_pred)): #list how many templates are in each kmeans cluster
            print("Kmeans",i,"Count",list(y_pred).count(i))


        # In[ ]:
        save_tdf = tdf.copy() #create a copy of tdf for saving

        save_tdf['Kmeans'] = list(y_pred) #add Kmeans cluster to df
        test = save_tdf.drop(labels=["waveform","stft"],axis=1) #drop waveform and stft for saving
        test.reset_index(drop=True) #reset index (gets altered when editing df), drop=True removes old index

        #save to csv
        test.to_csv(f'/home/smocz/expand_redpy_new_files/kmeans_K_{K}_Volcano_{volc}.csv',index=False)