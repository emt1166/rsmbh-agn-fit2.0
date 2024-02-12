#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 17:08:15 2024

@author: emilytemple

This script will be used to assess the effect of SNR on the BLR measurements. 
It "messes up" a given spectrum to degrade its SNR. 

This is done using the already given SDSS wavelength bins & gaussian variables. 

The plan is as follows:
Want to generate spectra with an increased sigma of the gaussian used to determine 
the flux of a given wavelength bin.
    1) Investigate the SDSS binning process (gaussian values for each flux in a bin)
    2) Code a way to pick out a random number asssociated with a gaussian paramater
    3) Use this random number to "mess up" the spectra and create noise
    4) See how this added noise (reduced SNR) effects the measurements
    
The goal of this is to determine what may be affecting our final measurements. 
A lot of these lower luminosity quasar spectra (Liu et al. 2019) are lower SNR by nature. 

----------------------------------------------------------------------------------------------------------
Turns out, there is a built-in function called random.gauss(). 
Inputs of random.gauss() are mu, sigma.
mu is the mean
sigma is the standard deviation

This function can be used to choose a random number in a gaussian distribution, and therefore
can be used to manipulate SDSS spectra. 

Here, our mu, the mean, is going to be the binned flux measurement. 
The sigma, is going to be the standard deviation, error, or uncertainty defined by SDSS. 
In this case, it is a "one sigma" uncertainty (or error). 

The plan:
    1) Use random.gauss() to pick out a value to replace the recorded flux value with
    2) Loop this process through every wavelength bin, effectively replacing all flux values
    with new ones, and introducing an increased S/N
    3) Compare the line profile properties to see what is happening when S/N is changed 


"""
#code
import statistics as stat
import numpy as np
import matplotlib.pyplot as plt
import os, re
from astropy.io import fits
from astropy.table import Table
import pandas as pd
import random

#first we need to extrapolate the information from the files
#this includes the flux, wavelengths, and error
#we can use the error to help compute SNR 
path1 = '/Users/emilytemple/documents/rsmbh-agn-fit2.0/pyqsofit/' 

#the path to where the data is stored
Data_path = '/Users/emilytemple/documents/rsmbh-agn-fit2.0/pyqsofit/Data/Outliers'
Data_list_names = os.listdir(Data_path) 


#Will input a LOOP here eventually -----------------------------------------------------------------------------------

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

#sourcename='1624-53386-0032'
#this source is one of the highest SNR

def Input_SN(sourcename):

    data = fits.open(os.path.join(path1+'Data/Outliers/spec-'+sourcename+'.fits'))
    lam = 10**data[1].data['loglam']                           # OBS wavelength (A)
    flux = data[1].data['flux']                           # OBS flux (erg/s/cm^2/A)
    err = 1./np.sqrt(data[1].data['ivar'])                          # 1 sigma error
    
    #SNR_calc = flux/err #SNR calc if needed
    
    new_flux = []
    i = 0
    while i < len(flux):
        #coeff_err = np.random.randint(1,11)#sets how noisy it is, will change to np.random prob
        f = random.gauss(flux[i],1*err[i])
        #print(i, f)
        new_flux.append(f)
        if i > len(flux):
            break
        i +=1 
    
    '''
    Now we need to save these as a usable spectrum to do the rest of the pyqso
    fitting process, just like with the organic spectrum from SDSS.
    
    I think all you have to do is input the new_flux array and let it run... 
    I am going to save the new_flux array somewhere and have it pull from that. 
    Much easier than manipulating everything at once, just change the paths. 
    
    
    I will make a new folder to save all of the new_flux results in, and make sure to
    identify them using the sourcename. Outliers New Flux is the folder, in the Data folder. 
    '''
    
    #saving the new_flux arrays to be used
    
    #setting the path to the new flux arrays (introducing random excess noise)
    new_flux_path = '/Users/emilytemple/documents/rsmbh-agn-fit2.0/pyqsofit/Data/Outliers New Flux'
    #specific placement for the degree of noise introduced (defined by coeff_err)
    snr_degree1 = new_flux_path+'/1'
    snr_degree5 = new_flux_path+'/1-5'
    snr_degree10 = new_flux_path+'/1-10'
    snr_degree20 = new_flux_path+'/10-20'
    TEST = new_flux_path+'/TEST'
    #saving the file to the appropriate folder, with ID name
    np.savetxt(TEST+'/'+sourcename+'.new_flux.txt', new_flux)
    
    plt.figure()
    plt.plot(lam,flux, label='Orig. Spec')
    plt.plot(lam, new_flux, label='Manip. Spec')
    plt.xlim(6200,7000)
    plt.legend()
    
    
    
    return 

# Add choice of loop ------------------------------------------------------------------------------------------------
loop = False
if loop: 
    for source in Data_list:
        Input_SN(source)
else:
    Input_SN('0332-52367-0639')







