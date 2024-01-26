#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 09:48:29 2024

@author: emilytemple

This script will be used to quantity an appropriate SNR for the sample.
This means you can use this SNR as criteria to weed out bad datasets. 

The authors in the Liu et al. 2019 paper have SNR listed, 
however I will calculate my own based on the data. 
I also hope to cross-check my calculation with the authors'.

"""
import statistics as stat
import numpy as np
import matplotlib.pyplot as plt
import os, re
from astropy.io import fits 


#first we need to extrapolate the information from the files
#this includes the flux, wavelengths, and error
#we can use the error to help compute SNR I think... 

#load in the .npy data files saved by pyqso 
#need to define the paths to where everything is
#I want to loops this eventually, try single source first

# sourcename='2031-53848-0060'
path1 = '/Users/emilytemple/documents/rsmbh-agn-fit2.0/pyqsofit/' 

# data = fits.open(os.path.join(path1+'Data/Outliers/spec-'+sourcename+'.fits'))
# lam = 10**data[1].data['loglam']                           # OBS wavelength (A)
# flux = data[1].data['flux']                           # OBS flux (erg/s/cm^2/A)
# err = 1./np.sqrt(data[1].data['ivar'])                          # 1 sigma error
# z = data[2].data['z'][0]                                             # Redshift
    

# # okay now we simply use the above data to calculate SNR
# # we can average it over the wavelength after, but let's see the spread first 
# # divide flux/error, look into error specifics 

# SNR_calc = flux/err

# plt.figure()
# plt.plot(lam,flux, label='spec')
# plt.plot(lam,SNR_calc, color='orange', label='snr')
# plt.legend()

'''
Okay so this is the SNR calculation from the direct SDSS data and error.
I wonder if we can do this on a "cleaner" spectrum...

I also want to see the mean SNR and do some stat things
'''







#alright now we are going to go through ALL the outliers
#and cut out "bad" SNR
Data_path = '/Users/emilytemple/documents/rsmbh-agn-fit2.0/pyqsofit/Data/Outliers'
Data_list_names = os.listdir(Data_path) 

#Pattern to extract the desired number sequence
pattern = r'spec-(\d{4}-\d{5}-\d{4})\.fits'

# Initialize an empty list to store the extracted sequences
Data_list = []

# Iterate through the file names and extract the sequences
for name in Data_list_names:
    match = re.search(pattern, name)
    if match:
        number_sequence = match.group(1)
        Data_list.append(number_sequence)

def SNR_cut(sourcename):
    
    data = fits.open(os.path.join(path1+'Data/Outliers/spec-'+sourcename+'.fits'))
    lam = 10**data[1].data['loglam']                           # OBS wavelength (A)
    flux = data[1].data['flux']                           # OBS flux (erg/s/cm^2/A)
    err = 1./np.sqrt(data[1].data['ivar'])                          # 1 sigma error
    
    SNR_calc = flux/err #SNR calc
    
    #saving each SNR calculation by sourcename
    path_snr_save = '/Users/emilytemple/documents/rsmbh-agn-fit2.0/pyqsofit/Fit Results/Line Properties/OutlierSNR'
    np.save(path_snr_save+'/'+sourcename, SNR_calc)
    
    plt.figure()
    plt.plot(lam,flux, label='spec')
    plt.plot(lam,SNR_calc, color='orange', label='snr')
    plt.title(f'{sourcename}')
    plt.xlabel('wavelength (Ang)')
    plt.ylabel('flux')
    plt.legend()
    return
    
loop = True
if loop: 
    for source in Data_list:
        SNR_cut(source)
else:
    SNR_cut('2031-53848-0060') 







