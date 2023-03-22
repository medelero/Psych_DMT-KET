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

path = '/home/medelero/Desktop/papers/SciRep_Psychodelics/data/dmt_measures'
paths = glob.glob(path+'/*')
paths.sort()

slopes_dmt = []
slopes_ec = []
slopes_eo = []

ssd_dmt = []
ssd_ec = []
ssd_eo = []
for i in paths:
    mat = io.loadmat(i)
    slopes_dmt.append(mat['dmt_aperiodics'][:,2])
    slopes_ec.append(mat['ec_aperiodics'][:,2])
    slopes_eo.append(mat['eo_aperiodics'][:,2])

    ssd_dmt.append(mat['dmt_ssd_pattern'][0])
    ssd_ec.append(mat['ec_ssd_pattern'][0])
    ssd_eo.append(mat['eo_ssd_pattern'][0])
    
slopes_dmt = np.array(slopes_dmt)
slopes_ec = np.array(slopes_ec)
slopes_eo = np.array(slopes_eo)

pow_dmt = np.mean(np.array(ssd_dmt), 0)
pow_eo = np.mean(np.array(ssd_eo), 0)
pow_ec = np.mean(np.array(ssd_ec), 0)

#%% Plot SSD Alpha Patterns
import mne

path = '/media/medelero/TOSHIBA EXT/DMT_cleaneeg/CleanEEG/'
subj = 'S01'

raw_eo = mne.io.read_epochs_eeglab(path + subj + '-EO_ICA_pruned.set')
pattern = mne.EvokedArray(data=np.array([pow_ec, pow_dmt]).T, 
                          info=raw_eo.info)
pattern.plot_topomap(np.arange(0, 0.004, 0.002), units=dict(eeg='A.U.'),
                     time_format='')

#%%

df = pd.DataFrame({'Eyes Open': np.mean(slopes_eo, 1),
                   'DMT':  np.mean(slopes_dmt, 1)})


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
labels = ['Eyes Open', 'DMT']

df_long = pd.melt(df,  value_vars=['Eyes Open', 'DMT'])

# Create empty figure and plot the individual datapoints
fig, ax = plt.subplots(figsize=(5,5))

for col in df:
    ax.plot(df_jitter_2[col], df[col], 'o', alpha=0.7, zorder=2, ms=10, mew=1.5)

for idx in df.index:
    ax.plot(df_jitter_2.loc[idx,['Eyes Open', 'DMT']], df.loc[idx,['Eyes Open', 'DMT']], color = 'black', linewidth = 2, linestyle = '-',alpha = .3)    
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

df = pd.DataFrame({'Eyes Open': np.mean(slopes_eo, 0),
                   'DMT':  np.mean(slopes_dmt, 0)})


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
labels = ['Eyes Open', 'DMT']

df_long = pd.melt(df,  value_vars=['Eyes Open', 'DMT'])

# Create empty figure and plot the individual datapoints
fig, ax = plt.subplots(figsize=(5,5))

for col in df:
    ax.plot(df_jitter_2[col], df[col], 'o', alpha=0.7, zorder=2, ms=10, mew=1.5)

for idx in df.index:
    ax.plot(df_jitter_2.loc[idx,['Eyes Open', 'DMT']], df.loc[idx,['Eyes Open', 'DMT']], color = 'black', linewidth = 2, linestyle = '-',alpha = .3)    
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

df = pd.DataFrame({'Eyes Closed': np.mean(slopes_ec, 1),
                   'DMT':  np.mean(slopes_dmt, 1)})


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
labels = ['Eyes Closed', 'DMT']

df_long = pd.melt(df,  value_vars=['Eyes Closed', 'DMT'])

# Create empty figure and plot the individual datapoints
fig, ax = plt.subplots(figsize=(5,5))

for col in df:
    ax.plot(df_jitter_2[col], df[col], 'o', alpha=0.7, zorder=2, ms=10, mew=1.5)

for idx in df.index:
    ax.plot(df_jitter_2.loc[idx,['Eyes Closed', 'DMT']], df.loc[idx,['Eyes Closed', 'DMT']], color = 'black', linewidth = 2, linestyle = '-',alpha = .3)    
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

df = pd.DataFrame({'Eyes Closed': np.mean(slopes_ec, 0),
                   'DMT':  np.mean(slopes_dmt, 0)})


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
labels = ['Eyes Closed', 'DMT']

df_long = pd.melt(df,  value_vars=['Eyes Closed', 'DMT'])

# Create empty figure and plot the individual datapoints
fig, ax = plt.subplots(figsize=(5,5))

for col in df:
    ax.plot(df_jitter_2[col], df[col], 'o', alpha=0.7, zorder=2, ms=10, mew=1.5)

for idx in df.index:
    ax.plot(df_jitter_2.loc[idx,['Eyes Closed', 'DMT']], df.loc[idx,['Eyes Closed', 'DMT']], color = 'black', linewidth = 2, linestyle = '-',alpha = .3)    
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

