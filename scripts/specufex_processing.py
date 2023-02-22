#applies NMF and HMM to a csv containing a df of spectrograms as created in specufex_testing.ipynb

#imports, adapted from specufex geysers_tutorial.ipynb
import datetime

import h5py
from matplotlib import pyplot as plt
from matplotlib import cm
import numpy as np
import pandas as pd
import yaml
import scipy.signal as sp
from sklearn.cluster import KMeans
import seaborn as sns
from tqdm import trange

from specufex import BayesianNonparametricNMF, BayesianHMM

#read config.yaml
with open('/home/smocz/expand_redpy/scripts/config.yaml') as file:
        config = yaml.load(file, Loader=yaml.FullLoader)

volc_list_names = config['volc_list_names']
vv = config['vv']
volc = volc_list_names[vv]
homedir = config['homedir']


#run NMF

