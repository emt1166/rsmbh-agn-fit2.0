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

NOTE: there is an option to turn OFF the loop for the SNR calculation
"""
import statistics as stat
import numpy as np
import matplotlib.pyplot as plt
import os, re
from astropy.io import fits 
import pandas as pd

#first we need to extrapolate the information from the files
#this includes the flux, wavelengths, and error
#we can use the error to help compute SNR 
path1 = '/Users/emilytemple/documents/rsmbh-agn-fit2.0/pyqsofit/' 

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
    new_flux_path = '/Users/emilytemple/documents/rsmbh-agn-fit2.0/pyqsofit/Data/Outliers New Flux/10-20/'+sourcename+'.new_flux.txt'
    synthetic_flux = np.loadtxt(new_flux_path)
    
    
    coeff_err = np.random.randint(1,2) #sets how noisy it is, might need to coordinate this
    SNR_calc = flux/err #SNR calc
    #print(len(lam))
       
    #saving each SNR calculation by sourcename
    path_snr_save = '/Users/emilytemple/documents/rsmbh-agn-fit2.0/pyqsofit/Fit Results/Line Properties/OutlierSNR'
    np.save(path_snr_save+'/'+sourcename, SNR_calc)
    
    
    #hbeta_range = lam[1970:2580] #for reference
    #halpha_range = lam[3200:3780] #for reference
    
    SNR_file_path = path_snr_save+'/'+sourcename
    SNR_file = np.load(SNR_file_path+'.npy')
    
    #calculating mean SNR for h alpha and h beta regions 
    #print(r'The H beta mean SNR is', stat.mean(SNR_file[1970:2580]))
    #print(r'The H alpha mean SNR is', stat.mean(SNR_file[3260:3770]))
    print(r'The total SNR is', stat.mean(SNR_file))
    
    
    # #need a way to sort the SNR of all outliers
    # #save the calculations into their own file
    # #organized by sourcename, and all three SNR calculations
    # h_beta_mean_snr = stat.mean(SNR_file[1970:2580])
    # h_alpha_mean_snr = stat.mean(SNR_file[3260:3770])
    # total_snr = stat.mean(SNR_file)
    # sourcename = f'{sourcename}'
    
    # # Creating a DataFrame
    # new_data = {
    #     'Sourcename': f'{sourcename}',
    #     'H beta mean SNR': [h_beta_mean_snr],
    #     'H alpha mean SNR': [h_alpha_mean_snr],
    #     'Total SNR': [total_snr]
    # }
    
    # # Check if the output CSV file exists
    # output_csv_path = 'SNR_output.csv'
    # if os.path.exists(output_csv_path):
    #     # If it exists, read the existing DataFrame
    #     existing_df = pd.read_csv(output_csv_path)
    #     # Append the new data to the existing DataFrame
    #     df = pd.concat([existing_df, pd.DataFrame(new_data)], ignore_index=True)
    # else:
    #     # If it doesn't exist, create a new DataFrame with the new data
    #     df = pd.DataFrame(new_data)

    # # Writing to CSV file
    # df.to_csv(output_csv_path, index=False)

   
    plt.figure() #plots spec and SNR calc to see comparison 
    plt.plot(lam,flux, label='spec')
    plt.plot(lam,SNR_calc, color='orange', label='snr')
    plt.title(f'{sourcename}')
    plt.xlabel('wavelength (Ang)')
    plt.ylabel('flux')
    plt.legend()
    return

#allowing CHOICE of loop, say False if you only want to look at a particular source    
loop = True
if loop: 
    for source in Data_list:
        SNR_cut(source)
else:
    SNR_cut('2974-54592-0229') 

#----------------------------------------------------------------------------------------------------------------------
# Now we need to take a look at our SNR and decide an appropriate cut
# this will also involve some trial and error and is kind of subjective... 
# I will try to use stats to do it (maybe quartile ranges?)




'''
Now that the initial calculation is done, we have some options:
    1) Focus on Halpha and Hbeta regions of interest and do SNR cut
    2) Average the SNR over the entire spectrum and use that

    Option 1 seems to be harder, but would work for our purposes perhaps. 
    Taking a look at the graphs, there are some objects that have a significant range in SNR. 
    This means averaging it would not be the best idea.... 
I will attempt Option #1. 
'''

'''
Okay because the length for wavelengths and error are the same, I am sure I can isolate Halpha and Hbeta 
I only want the range of SNR calculations from the following regions (VALUES of lam)
Halpha:
    8000 to 9000
Hbeta:
    6000 to 7000
These ranges are NOT redshift corrected, but that doesn't matter here
I am using the raw spec to calculate a preliminary SNR first. 

#by trial and error we get the following approx ranges:
hbeta = lam[1970:2600]
halpha = lam[3200:3780]
# these should be the same approx ranges in the SNR_calc array. 
# I will average over those

I INCORPORATED THIS INTO THE ABOVE LOOP
'''









