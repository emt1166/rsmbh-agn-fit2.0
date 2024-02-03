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

"""

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

sourcename='1624-53386-0032'
#this source is one of the highest SNR

data = fits.open(os.path.join(path1+'Data/Outliers/spec-'+sourcename+'.fits'))
lam = 10**data[1].data['loglam']                           # OBS wavelength (A)
flux = data[1].data['flux']                           # OBS flux (erg/s/cm^2/A)
err = 1./np.sqrt(data[1].data['ivar'])                          # 1 sigma error

SNR_calc = flux/err #SNR calc

'''
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
'''


#print(random.gauss(flux,err))

# plt.figure()
# plt.plot(lam,flux, label='OG Spectrum')
# plt.plot(lam, random.gauss(flux,err), label='Manipulated Spectrum')
# plt.xlim(6000,7000)
# plt.legend()

#need to make it so that the random.gauss() is NOT the same for each flux value
#this will require a simple for loop


new_flux = []
i = 0
while i < len(flux):
    coeff_err = np.random.randint(1,11)#sets how noisy it is
    f = random.gauss(flux[i],coeff_err*err[i])
    new_flux.append(f)
    if i > len(flux):
        break
    i +=1 

#np.savetxt(path1+'new_flux.txt', new_flux)
plt.figure()
plt.plot(lam,flux, label='Orig. Spec')
plt.plot(lam, new_flux, label='Manip. Spec')
plt.xlim(6300,7000)
plt.legend()
#okay now we have it so the SN is changed by different amounts across the spec
#what if I do different multiples of the error?? (10*error gives a BAD SNR...)
'''








