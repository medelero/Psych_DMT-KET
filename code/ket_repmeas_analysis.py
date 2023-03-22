#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 23:34:54 2023

@author: medelero
"""


import glob
from scipy import io
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

path = '/home/medelero/Desktop/papers/SciRep_Psychodelics/data/ketamine_measures/EC'
paths = glob.glob(path+'/*')
paths.sort()

slopes_keta = []
slopes_awake = []

ssd_keta = []
ssd_awake = []
psd_ket = np.zeros((22, 62, 129))
psd_awake = np.zeros((22, 62, 129))

for index, i in enumerate(paths):

    mat = io.loadmat(i)

    if i[-20:][4:7] == '007':
        slopes_keta.append(mat['aperiodics'][:, 2])
        ssd_keta.append(mat['ssd_pattern'][0])
        psd_ket[index] = mat['psd']
    else: 
        slopes_awake.append(mat['aperiodics'][:, 2])
        ssd_awake.append(mat['ssd_pattern'][0])
        psd_awake[index] = mat['psd']

        

#%%
slopes_keta = np.array(slopes_keta)
slopes_awake = np.array(slopes_awake)

pow_ket = np.array(ssd_keta)
pow_awake = np.array(ssd_awake)

#%% Plot SSD Alpha Patterns
import mne

path = '/media/medelero/TOSHIBA EXT/EEG_ketamine/Farnes_et_al_PLOS_ONE_Dryad/spontaneous/EO/210_20161207_0002eyesOpen_afterICA.set'

raw_eo = mne.io.read_epochs_eeglab(path)
pattern = mne.EvokedArray(data=np.array([pow_awake[9], pow_ket[9]]).T, 
                          info=raw_eo.info)
pattern.plot_topomap(np.arange(0, 0.005, 0.003), units=dict(eeg='A.U.'),
                     time_format='')
#%%

df = pd.DataFrame({'Awake EC': np.mean(slopes_awake, 1),
                   'Ketamine EC':  np.mean(slopes_keta, 1)})


jitter_2 = 0.05
np.random.seed(3)
df_jitter_2 = pd.DataFrame(np.random.normal(loc=0, scale=jitter_2, size=df.values.shape), columns=df.columns
                           )
    
df_jitter_2 += np.arange(len(df.columns))
# Define pre-settings
w = 6
h = 6
title_size = 20
xlab_size = 15
ylab_size = 20
labels = ['Awake EC', 'Ketamine EC']

df_long = pd.melt(df,  value_vars=['Awake EC', 'Ketamine EC'])

# Create empty figure and plot the individual datapoints
fig, ax = plt.subplots(figsize=(5,5))

for col in df:
    ax.plot(df_jitter_2[col], df[col], 'o', alpha=0.7, zorder=2, ms=10, mew=1.5)

for idx in df.index:
    ax.plot(df_jitter_2.loc[idx,['Awake EC', 'Ketamine EC']], df.loc[idx,['Awake EC', 'Ketamine EC']], color = 'black', linewidth = 2, linestyle = '-',alpha = .3)    
    sns.pointplot(x='variable', y='value', ci=95, data=df_long, join=False, scale=0.01, color = 'black', capsize = .03)    
    sns.violinplot(x='variable', y='value', data=df_long, hue = 'variable', split = True, inner = 'quartile', cut=1, alpha=.01)

    
#Additonal settings
    ax.set_xticks(range(len(df.columns)))
    ax.set_xticklabels((labels), size= xlab_size)
    ax.set_xlim(-1, len(df.columns))
    ax.set_ylabel('1/f Slope', size = ylab_size)
    ax.set_title('All Elecs', size = title_size)
    ax.legend_.remove()
    sns.despine()
    plt.setp(ax.collections, alpha=.01) 
   
plt.tight_layout()
sns.despine()

#%%

df = pd.DataFrame({'Awake EC': np.mean(slopes_awake, 1),
                   'Ketamine EC':  np.mean(slopes_keta, 1)})


jitter_2 = 0.05
np.random.seed(3)
df_jitter_2 = pd.DataFrame(np.random.normal(loc=0, scale=jitter_2, size=df.values.shape), columns=df.columns
                           )
    
df_jitter_2 += np.arange(len(df.columns))
# Define pre-settings
w = 6
h = 6
title_size = 20
xlab_size = 15
ylab_size = 20
labels = ['Awake EC', 'Ketamine EC']

df_long = pd.melt(df,  value_vars=['Awake EC', 'Ketamine EC'])

# Create empty figure and plot the individual datapoints
fig, ax = plt.subplots(figsize=(5,5))

for col in df:
    ax.plot(df_jitter_2[col], df[col], 'o', alpha=0.7, zorder=2, ms=10, mew=1.5)

for idx in df.index:
    ax.plot(df_jitter_2.loc[idx,['Awake EC', 'Ketamine EC']], df.loc[idx,['Awake EC', 'Ketamine EC']], color = 'black', linewidth = 2, linestyle = '-',alpha = .3)    
    sns.pointplot(x='variable', y='value', ci=95, data=df_long, join=False, scale=0.01, color = 'black', capsize = .03)    
    sns.violinplot(x='variable', y='value', data=df_long, hue = 'variable', split = True, inner = 'quartile', cut=1, alpha=.01)

    
#Additonal settings
    ax.set_xticks(range(len(df.columns)))
    ax.set_xticklabels((labels), size= xlab_size)
    ax.set_xlim(-1, len(df.columns))
    ax.set_ylabel('1/f Slope', size = ylab_size)
    ax.set_title('All Elecs', size = title_size)
    ax.legend_.remove()
    sns.despine()
    plt.setp(ax.collections, alpha=.01) 
   
plt.tight_layout()
sns.despine()

#%%

alpha_anter_awake = pow_awake[:, -7] / pow_awake[:, 2]
alpha_anter_ket = pow_ket[:, -7] / pow_ket[:, 2]

df = pd.DataFrame({'Awake EC': alpha_anter_awake,
                   'Ketamine EC':  alpha_anter_ket})


jitter_2 = 0.05
np.random.seed(3)
df_jitter_2 = pd.DataFrame(np.random.normal(loc=0, scale=jitter_2, size=df.values.shape), columns=df.columns
                           )
    
df_jitter_2 += np.arange(len(df.columns))
# Define pre-settings
w = 6
h = 6
title_size = 20
xlab_size = 15
ylab_size = 20
labels = ['Awake EC', 'Ketamine EC']

df_long = pd.melt(df,  value_vars=['Awake EC', 'Ketamine EC'])

# Create empty figure and plot the individual datapoints
fig, ax = plt.subplots(figsize=(5,5))

for col in df:
    ax.plot(df_jitter_2[col], df[col], 'o', alpha=0.7, zorder=2, ms=10, mew=1.5)

for idx in df.index:
    ax.plot(df_jitter_2.loc[idx,['Awake EC', 'Ketamine EC']], df.loc[idx,['Awake EC', 'Ketamine EC']], color = 'black', linewidth = 2, linestyle = '-',alpha = .3)    
    sns.pointplot(x='variable', y='value', ci=95, data=df_long, join=False, scale=0.01, color = 'black', capsize = .03)    
    sns.violinplot(x='variable', y='value', data=df_long, hue = 'variable', split = True, inner = 'quartile', cut=1, alpha=.01)

    
#Additonal settings
    ax.set_xticks(range(len(df.columns)))
    ax.set_xticklabels((labels), size= xlab_size)
    ax.set_xlim(-1, len(df.columns))
    ax.set_ylabel('1/f Slope', size = ylab_size)
    ax.set_title('All Elecs', size = title_size)
    ax.legend_.remove()
    sns.despine()
    plt.setp(ax.collections, alpha=.01) 
   
plt.tight_layout()
sns.despine()