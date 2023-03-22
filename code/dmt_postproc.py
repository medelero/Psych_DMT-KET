#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  5 14:50:52 2023

@author: medelero
"""

import mne
import os
import numpy as np
import fooof
from scipy import io
from neurodsp import spectral
from mne.decoding import SSD

#%% Load Path

path = '/media/medelero/TOSHIBA EXT/DMT_cleaneeg/CleanEEG/'

subj_id = []

for i in np.arange(1, 36, 1):
    if i < 10:
        subj_id.append('S0' + str(i))
    else: subj_id.append('S' + str(i))
    
subj_id = np.delete(subj_id, [31, 23, 16, 8, 5, 2]) #idx of missing subjects

#%%

for subj in subj_id:

    raw_eo = mne.io.read_epochs_eeglab(path + subj + '-EO_ICA_pruned.set')
    raw_ec = mne.io.read_epochs_eeglab(path + subj + '-EC_ICA_pruned.set')
    raw_dmt = mne.io.read_epochs_eeglab(path + subj + '-DMT_ICA_pruned.set')
    
    fs = int(raw_eo.info['sfreq'])

# Alpha Freq params

    freqs_sig = 8, 13
    freqs_noise = 6, 15

# Calculate Spectro-Spatial Decomposition (SSD) spatial filters

# EO
    ssd_eo = SSD(info=raw_eo.info, reg='oas', 
              filt_params_signal=dict(l_freq=freqs_sig[0], h_freq=freqs_sig[1],
                                      l_trans_bandwidth=1, h_trans_bandwidth=1),
              filt_params_noise=dict(l_freq=freqs_noise[0], h_freq=freqs_noise[1],
                                     l_trans_bandwidth=1, h_trans_bandwidth=1))
    ssd_eo.fit(X=raw_eo.get_data())

    ssd_eo_sources = ssd_eo.transform(X=raw_eo.get_data())
    
    pattern_eo = ssd_eo.patterns_[:5]

    # Get psd of SSD-filtered signals.
    eo_ssd_psd, freqs = mne.time_frequency.psd_array_welch(ssd_eo_sources, 
                                                    sfreq=fs, 
                                                    n_per_seg=500,
                                                    n_overlap=100)

# EC

    ssd_ec = SSD(info=raw_ec.info, reg='oas', 
              filt_params_signal=dict(l_freq=freqs_sig[0], h_freq=freqs_sig[1],
                                      l_trans_bandwidth=1, h_trans_bandwidth=1),
              filt_params_noise=dict(l_freq=freqs_noise[0], h_freq=freqs_noise[1],
                                     l_trans_bandwidth=1, h_trans_bandwidth=1))
    ssd_ec.fit(X=raw_ec.get_data())

    ssd_ec_sources = ssd_ec.transform(X=raw_ec.get_data())
    
    pattern_ec = ssd_ec.patterns_[:5]

    # Get psd of SSD-filtered signals.
    ec_ssd_psd, freqs = mne.time_frequency.psd_array_welch(ssd_ec_sources, 
                                                    sfreq=fs, 
                                                    n_per_seg=500,
                                                    n_overlap=100)

# DMT

    ssd_dmt = SSD(info=raw_dmt.info, reg='oas', 
              filt_params_signal=dict(l_freq=freqs_sig[0], h_freq=freqs_sig[1],
                                      l_trans_bandwidth=1, h_trans_bandwidth=1),
              filt_params_noise=dict(l_freq=freqs_noise[0], h_freq=freqs_noise[1],
                                     l_trans_bandwidth=1, h_trans_bandwidth=1))
    ssd_dmt.fit(X=raw_dmt.get_data())
    ssd_dmt_sources = ssd_dmt.transform(X=raw_dmt.get_data())
    pattern_dmt = ssd_dmt.patterns_[:5]

    # Get psd of SSD-filtered signals.
    dmt_ssd_psd, freqs = mne.time_frequency.psd_array_welch(ssd_dmt_sources, 
                                                    sfreq=fs, 
                                                    n_per_seg=500,
                                                    n_overlap=100)
    
# Slope

    psd_eo, f = mne.time_frequency.psd_array_welch(np.concatenate(raw_eo.get_data(),1), 
                                                   sfreq=fs)
    fg = fooof.FOOOFGroup(aperiodic_mode='knee')
    fg.fit(f, psd_eo, freq_range=[0.5, 40])
    eo_aperiodics = fg.get_params('aperiodic_params')
    eo_peaks = fg.get_params('peak_params')
    eo_peaks = eo_peaks[eo_peaks[:, 0] < 15]

    psd_ec, f = mne.time_frequency.psd_array_welch(np.concatenate(raw_ec.get_data(),1), 
                                                   sfreq=fs)    
    fg = fooof.FOOOFGroup(aperiodic_mode='knee')
    fg.fit(f, psd_ec, freq_range=[0.5, 40])   
    ec_aperiodics = fg.get_params('aperiodic_params')
    ec_peaks = fg.get_params('peak_params')
    ec_peaks = ec_peaks[ec_peaks[:, 0] < 15]

    psd_dmt, f = mne.time_frequency.psd_array_welch(np.concatenate(raw_dmt.get_data(),1), 
                                                    sfreq=fs)    
    fg = fooof.FOOOFGroup(aperiodic_mode='knee')
    fg.fit(f, psd_dmt, freq_range=[0.5, 40])    
    dmt_aperiodics = fg.get_params('aperiodic_params')
    dmt_peaks = fg.get_params('peak_params')
    dmt_peaks = dmt_peaks[dmt_peaks[:, 0] < 15]
    
    data_dic = {}
    data_dic['dmt_psd'] = psd_dmt
    data_dic['dmt_peak_params'] = dmt_peaks
    data_dic['dmt_aperiodics'] = dmt_aperiodics
    data_dic['dmt_ssd_pattern'] = pattern_dmt
    data_dic['dmt_ssd_sources'] = ssd_dmt_sources

    data_dic['ec_psd'] = psd_ec    
    data_dic['ec_peak_params'] = ec_peaks
    data_dic['ec_aperiodics'] = ec_aperiodics
    data_dic['ec_ssd_pattern'] = pattern_ec
    data_dic['ec_ssd_sources'] = ssd_ec_sources
    
    data_dic['eo_psd'] = psd_eo
    data_dic['eo_peak_params'] = eo_peaks
    data_dic['eo_aperiodics'] = eo_aperiodics
    data_dic['eo_ssd_pattern'] = pattern_eo
    data_dic['eo_ssd_sources'] = ssd_eo_sources
    
    io.savemat(subj+'_measures.mat', data_dic)
    print('subject saved correctly')
    
    