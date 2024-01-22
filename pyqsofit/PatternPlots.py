#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 11:03:04 2024

@author: emilytemple
"""

"""
This code is used to analyze the output csv files created in the
ReadPyqso file.
The purpose of this is to search for any trends within the outlier sample in 
regards to FWHM, velocity shifts, kurtosis, asymmetries, etc...
"""

import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd


#Here we read in the data from the csv files
#the 'df' is to indicate it is a dataframe

OutlierLineProfile_df = pd.read_csv('OutliersLineProfiles_output_data2.0.csv')
OutlierLineShifts_df = pd.read_csv('OutliersLineShifts_output_data2.0.csv')

ControlLineProfile_df = pd.read_csv('ControlLineProfiles_output_data2.0.csv')
ControlLineShifts_df = pd.read_csv('ControlLineShifts_output_data2.0.csv')

plt.figure()
#plt.semilogy()
plt.scatter((OutlierLineShifts_df['VShift Center BHA']), OutlierLineShifts_df['Kurtosis BHA'], label='Outliers')
plt.scatter((ControlLineShifts_df['VShift Center BHA']), ControlLineShifts_df['Kurtosis BHA'], label= 'Control')
plt.legend()
plt.xlabel('BHB Center VShift')
plt.ylabel('BHB Kurtosis')
