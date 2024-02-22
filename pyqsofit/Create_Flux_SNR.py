#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 12:48:25 2024

@author: emilytemple

This script is used to create the different flux arrays for each individual object.
It uses the sigma given by SDSS and a gaussian distribution to choose a random
"new" value of the flux.

This was originally in the SNR_CompareFits.py, but it didn't make any sense to put
all of that in one file. 
Run this file first to create all the "new" flux values per sourcename.
Right now, wecreate 5 iterations of new flux values. This will eventually be increased
to 500 for a much better sample size. 


"""

import os,timeit
import numpy as np
from PyQSOFit import QSOFit
from astropy.io import fits
import matplotlib.pyplot as plt
import warnings
import re
import pandas as pd
import scipy.integrate as integrate
import random

# Use custom matplotlib style to make Yue happy
QSOFit.set_mpl_style()
# Ignore warnings?
warnings.filterwarnings("ignore")

# -----------------------------------------------------------------------------------------
#choose which object to look at MOVED TO BOTTOM W/ LOOP CHOICE
#sourcename = '0813-52354-0561'
# ----------------------------------------------------------------------------------------


# Setting the file paths that the code uses
# The path of the source code file and qsopar.fits
path1 = '/Users/emilytemple/documents/rsmbh-agn-fit2.0/pyqsofit/'
# save fits path, unused
path3 = path1+'/Fit Results/QA Other'
# The path of the SNR New Flux Comparison Data
pathF = '/Users/emilytemple/documents/rsmbh-agn-fit2.0/pyqsofit/Data/Outliers New Flux/LOOP1/'

# loop through all the files, create their directories, and save results
Data_path = '/Users/emilytemple/documents/rsmbh-agn-fit2.0/pyqsofit/Data/Outliers'
Data_list_names = os.listdir(Data_path)

#Pattern to extract the desired number sequence, plate-mjd-fiber
pattern = r'spec-(\d{4}-\d{5}-\d{4})\.fits'

# Initialize an empty list to store the extracted sequences
Data_list = []

# Iterate through the file names and extract the sequences
for name in Data_list_names:
    match = re.search(pattern, name)
    if match:
        number_sequence = match.group(1)
        Data_list.append(number_sequence)
        

# ------------------------------------------------------------------------------------------
def Input_SN(sourcename): #gives 5 different SN files per sourcename

    #Creating directories automatically from the sourcename
    #if not done so already
        
    if not os.path.exists(pathF+sourcename):
        os.makedirs(pathF+sourcename)
        
    if not os.path.exists(pathF+sourcename+'/new_flux'):
        os.makedirs(pathF+sourcename+'/new_flux')

    #path to save ALL outputs, everything goes here
    path2 = pathF+sourcename+'/'
    
    data = fits.open(os.path.join(path1+'Data/Outliers/spec-'+sourcename+'.fits'))
    lam = 10**data[1].data['loglam']                           # OBS wavelength (A)
    flux = data[1].data['flux']                           # OBS flux (erg/s/cm^2/A)
    err = 1./np.sqrt(data[1].data['ivar'])                          # 1 sigma error
    
    #SNR_calc = flux/err #SNR calc if needed
    for n in range(5): #specify range for how many new_flux files we want to use
        #n +=1 #do i need this? don't think so?
        new_flux = []
        i = 0
        while i < len(flux):
            f = random.gauss(flux[i],1*err[i]) #pulling new flux value from gaussian distn 
            new_flux.append(f)
            if i > len(flux): #making new flux values across the entire length of the flux array
                break
            i +=1 
        #save new fluxes to this path, will be used in the SNR fit compare script    
        np.savetxt(pathF+sourcename+'/new_flux/'+sourcename+'.'+f'{n}'+'new_flux.txt', new_flux)
        
# ------------------------------------------------------------------------------------------
# Call the loop to generate the files, can loop it through the directory as well
# Loop over the values using the function
# you can CHOOSE to loop or to specify a spectrum 
loop = True

if loop:
    for source in Data_list:
        Input_SN(source)
else:
    Input_SN('0813-52354-0561')  
# generates the different flux arrays and saves them




