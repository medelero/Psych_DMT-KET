#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 11:23:10 2023

@author: medelero
"""

import pandas as pd

# Load data
scales = pd.read_csv('/home/medelero/Desktop/papers/SciRep_Psychodelics/data/dmt_scales_results_replaced.csv')

# Drop some noisy subjects
scales = scales[scales['Unnamed: 0'].str.contains('S03|S06|S09|S17|S24|S32')==False]
scales = scales.drop(columns='Unnamed: 0')

# Def function
def is_what_percent_of(num_a, num_b):
    return (num_a / num_b) * 100

# Calculate % change
perc_change = np.zeros((slopes_dmt.shape))

for i in range(slopes_dmt.shape[0]):
    for j in range(slopes_dmt.shape[1]):
        perc_change[i, j] = is_what_percent_of(slopes_dmt[i,j], slopes_ec[i, j])

