#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 20:15:27 2023

@author: medelero
"""

import mne
import os
import numpy as np
import fooof
from scipy import io
from mne.decoding import SSD

#%% Eyes Close

path = '/media/medelero/TOSHIBA EXT/EEG_ketamine/Farnes_et_al_PLOS_ONE_Dryad/spontaneous/EC'

eeg_files = np.sort(os.listdir(path))

for f in eeg_files:
    subj_id_1 = f[:3]
    subj_id_2 = f[14:17]
    
    subj_id = subj_id_1 + '_' + subj_id_2
    
    raw = mne.io.read_epochs_eeglab(path + '/' + f[:-4] + '.set')
    fs = raw.info['sfreq']

# Alpha Freq params

    freqs_sig = 8, 13
    freqs_noise = 6, 15

# Calculate Spectro-Spatial Decomposition (SSD) spatial filters

# EO
    ssd = SSD(info=raw.info, reg='oas', 
              filt_params_signal=dict(l_freq=freqs_sig[0], h_freq=freqs_sig[1],
                                      l_trans_bandwidth=1, h_trans_bandwidth=1),
              filt_params_noise=dict(l_freq=freqs_noise[0], h_freq=freqs_noise[1],
                                     l_trans_bandwidth=1, h_trans_bandwidth=1))
    ssd.fit(X=raw.get_data())

    ssd_sources = ssd.transform(X=raw.get_data())
    
    pattern = ssd.patterns_[:5]

    # Get psd of SSD-filtered signals.
    ssd_psd, freqs = mne.time_frequency.psd_array_welch(ssd_sources, 
                                                    sfreq=fs, 
                                                    n_per_seg=500,
                                                    n_overlap=100)

    psd, f = mne.time_frequency.psd_array_welch(np.concatenate(raw.get_data(),1), 
                                                   sfreq=fs)
    fg = fooof.FOOOFGroup(aperiodic_mode='knee')
    fg.fit(f, psd, freq_range=[0.5, 40])
    aperiodics = fg.get_params('aperiodic_params')
    peaks = fg.get_params('peak_params')
    peaks = peaks[peaks[:, 0] < 15]

    data_dic = {}
    data_dic['psd'] = psd
    data_dic['peak_params'] = peaks
    data_dic['aperiodics'] = aperiodics
    data_dic['ssd_pattern'] = pattern
    data_dic['ssd_sources'] = ssd_sources
    
    io.savemat(subj_id+'_measures.mat', data_dic)
    print('subject saved correctly')
    
#%% Eyes Open

path = '/media/medelero/TOSHIBA EXT/EEG_ketamine/Farnes_et_al_PLOS_ONE_Dryad/spontaneous/EO'

eeg_files = np.sort(os.listdir(path))

for f in eeg_files:
    subj_id_1 = f[:3]
    subj_id_2 = f[14:17]
    
    subj_id = subj_id_1 + '_' + subj_id_2
    
    raw = mne.io.read_epochs_eeglab(path + '/' + f[:-4] + '.set')
    fs = raw.info['sfreq']

# Alpha Freq params

    freqs_sig = 8, 13
    freqs_noise = 6, 15

# Calculate Spectro-Spatial Decomposition (SSD) spatial filters

# EO
    ssd = SSD(info=raw.info, reg='oas', 
              filt_params_signal=dict(l_freq=freqs_sig[0], h_freq=freqs_sig[1],
                                      l_trans_bandwidth=1, h_trans_bandwidth=1),
              filt_params_noise=dict(l_freq=freqs_noise[0], h_freq=freqs_noise[1],
                                     l_trans_bandwidth=1, h_trans_bandwidth=1))
    ssd.fit(X=raw.get_data())

    ssd_sources = ssd.transform(X=raw.get_data())
    
    pattern = ssd.patterns_[:5]

    # Get psd of SSD-filtered signals.
    ssd_psd, freqs = mne.time_frequency.psd_array_welch(ssd_sources, 
                                                    sfreq=fs, 
                                                    n_per_seg=500,
                                                    n_overlap=100)

    psd, f = mne.time_frequency.psd_array_welch(np.concatenate(raw.get_data(),1), 
                                                   sfreq=fs)
    fg = fooof.FOOOFGroup(aperiodic_mode='knee')
    fg.fit(f, psd, freq_range=[0.5, 40])
    aperiodics = fg.get_params('aperiodic_params')
    peaks = fg.get_params('peak_params')
    peaks = peaks[peaks[:, 0] < 15]

    data_dic = {}
    data_dic['psd'] = psd
    data_dic['peak_params'] = peaks
    data_dic['aperiodics'] = aperiodics
    data_dic['ssd_pattern'] = pattern
    data_dic['ssd_sources'] = ssd_sources
    
    io.savemat(subj_id+'_measures.mat', data_dic)
    print('subject saved correctly')