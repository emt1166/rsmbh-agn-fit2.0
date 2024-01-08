#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 13:06:43 2024

@author: emilytemple
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import pandas as pd
import os
import re

# IMPORTANT
# Remember to run the 'runscript' first to create all directories
# Also run the Fitting Script!
# Results from Fitting-Script.py are saved as plateID-MJD-fiberID_"".npy files


# Allowing a loop through directories to make this process faster
# Once you download all the spectral files from SDSS, this should be able to 
# loop through all the files, create their directories, and save results

Data_path = '/Users/emilytemple/documents/rsmbh-agn-fit2.0/pyqsofit/Data'
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
# -----------------------------------------------------------------------------
def lineprofilecalc_func(sourcename):
    
    # sourcename is just the plateID-MJD-fiberID.
    # sourcename = '0646-52523-0075'
    
    # Path of stored PyQSOFit fit results
    path = 'Fit Results/Line Properties/'+sourcename+'/Fit Data/'
    
    # Path for saving results from this code
    path2 = 'Fit Results/Line Properties/'+sourcename+'/'+'Line Profile Plots/'
    
    # -----------------------------------------------------------------------------
    
    # Obtaining line profile data result components from Fitting-Script.py
    bl = np.load(path+sourcename+'_BLData.npy')
    bl_profile = np.load(path+sourcename+'_BLSpectrum.npy')
    nl = np.load(path+sourcename+'_NLData.npy')
    data = np.load(path+sourcename+'_Data.npy')             # Flux from SDSS .fits file
    data_contFeII_sub = np.load(path+sourcename+'_DataCFeII.npy')  
    wavelength = np.load(path+sourcename+'_Wavelength.npy')
    
    # Converting .npy files into 2D dataframes to make velocity shift calculations
    # MUCH easier
    # BL Data
    bl_matrix = np.vstack((wavelength, bl)).T
    bl_data = pd.DataFrame(bl_matrix, columns=['Wavelength','Flux'])
    # NL Data
    nl_matrix = np.vstack((wavelength, nl)).T
    nl_data = pd.DataFrame(nl_matrix, columns=['Wavelength','Flux'])
    
    
    # -----------------------------------------------------------------------------
    
    # H_alpha Line Complex
    
    # If H_alpha is not present in your spectrum (z > 0.40), 
    # the velocity shift calculations can be skipped by answering 'no' to "bha_vs"
    
    bha_vs = 'yes'
    if(bha_vs == 'yes'):
        # Section 1: finds the peak velocity shift by obtaining the corresponding
        # wavelength of the max flux value of the BL profile
        # ---
        # BLs - change according to length of wings
        bha_data = bl_data.loc[(bl_data['Wavelength'] >= 6170) & \
                               (bl_data['Wavelength'] <= 6900)]
        bha_wave = bha_data['Wavelength']
        bha_flux = bha_data['Flux']
        # Finding peak 
        peak_bha_flux = bha_flux.max()
        peak_bha_wave = bha_data.loc[bha_data['Flux'] == peak_bha_flux]
        peak_bha = int(peak_bha_wave['Wavelength'])
        # ---
        # NLs 
        nha_data = nl_data.loc[(nl_data['Wavelength'] >= 6400) & \
                               (nl_data['Wavelength'] <=6800)]
        nha_flux = nha_data['Flux']
        # Finding peak: should be ~ vacuum wavelength and will be used for 
        # calculating all three types of velocity shifts
        peak_nha_flux = nha_flux.max()
        peak_nha_wave = nha_data.loc[nha_data['Flux'] == peak_nha_flux]
        peak_nha = int(peak_nha_wave['Wavelength'])
        # ---
        # Peak velocity shift 
        bha_pvs = ((peak_bha**2 - peak_nha**2) / (peak_bha**2 + peak_nha**2)) * 299792
        
        # Section 2: finds the centroid velocity shift by normalizing the BL 
        # profile and finding the cumulative distribution function (CDF) to obtain 
        # the corresponding wavelength at 50% of the CDF
        # ---
        # Normalizing BL profile and finding CDF
        bha_flux_rescaled = (bha_flux - bha_flux.min()) / (bha_flux.max() - bha_flux.min())
        bha_flux_area = integrate.trapezoid(bha_flux_rescaled, bha_wave)
        bha_flux_new = bha_flux_rescaled/bha_flux_area   
        bha_flux_norm = bha_flux_new / np.sum(bha_flux_new)
        bha_cdf = np.cumsum(bha_flux_norm)
        # ---
        # Finding centroid
        bha_ctr = np.interp(0.5, bha_cdf, bha_wave)
        # Ploting CDF 
        plot_bha_cdf = 'no'
        if(plot_bha_cdf == 'yes'):
            fig1 = plt.figure(figsize=(16,10))
            plt.plot(bha_wave, bha_cdf, linewidth=2, c='#DB1D1A')
            plt.title(sourcename+': CDF', fontsize=30)
            plt.ylabel('Probability', fontsize=20)
            plt.xlabel(r'Wavelength ($\rm \AA$)', fontsize=20)
            plt.xticks(fontsize=20)
            plt.yticks(fontsize=20)
            plt.tick_params(which='major', length=16, width=2)
            plt.tick_params(which='minor', length=8, width=1)
            plt.savefig(path2+sourcename+'_BHA_CDF.pdf')
            plt.axvline(bha_ctr, c='grey', linestyle='--', label='Centroid')
            plt.text(bha_ctr+5, 0.2, r'H$\alpha$ BL Centroid = {:.2f} $\AA$'.format(bha_ctr))
            plt.legend(fontsize=14)
        # ---
        # Centroid velocity shift
        bha_cvs = ((bha_ctr**2 - peak_nha**2) / (bha_ctr**2 + peak_nha**2)) * 299792 
    
        # Section 3: finding C80 shift by removing 10% of the area at each end of
        # the BL profile and obtaining the corresponding wavelength at the 
        # center.
        # ---
        bha_10 = np.interp(0.05, bha_cdf, bha_wave)
        bha_90 = np.interp(0.95, bha_cdf, bha_wave)
        bha_C80_data = bl_data.loc[(bl_data['Wavelength'] >= bha_10) & \
                                   (bl_data['Wavelength'] <= bha_90)]
        bha_C80_wave = bha_C80_data['Wavelength']
        bha_C80_flux = bha_C80_data['Flux']
        # Normalizing BL profile and finding CDF
        bha_C80_flux_rescaled = (bha_C80_flux - (bha_C80_flux.min())) / \
            (bha_C80_flux.max() - bha_C80_flux.min())
        bha_C80_area = integrate.trapezoid(bha_C80_flux_rescaled, bha_C80_wave)
        bha_C80_flux_new = bha_C80_flux_rescaled/bha_C80_area
        bha_C80_flux_norm = bha_C80_flux_new / np.sum(bha_C80_flux_new)
        bha_cdf_C80 = np.cumsum(bha_C80_flux_norm)
        # ---
        # Finding line center
        bha_C80_ctr = np.interp(0.5, bha_cdf_C80, bha_C80_wave)
        # C80 velocity shift
        bha_C80_vs = ((bha_C80_ctr**2 - peak_nha**2) / (bha_C80_ctr**2 + peak_nha**2)) * 299792 
        
        # Section 4: finding the Whittle 1985 line profile parameters; i.e. inter-
        # percentile velocity width (IPV), asymmetry (A), shift (S), and kurtosis
        # (K)
        # ---
        # Calculating FWHM
        bha_FWHM_calc = peak_bha_flux / 2
        bha_FWHM_range = bha_data.loc[(bha_data['Flux'] >= bha_FWHM_calc)]
        bha_FWHM_wave_min = bha_FWHM_range['Wavelength'].iloc[0]
        bha_FWHM_wave_max = bha_FWHM_range['Wavelength'].iloc[-1]
        bha_FWHM_wave = bha_FWHM_wave_max - bha_FWHM_wave_min    
        bha_FWHM = bha_FWHM_wave / peak_nha * 299792
        # ---
        # Calculating line parameters    
        bha_a = bha_C80_ctr - bha_10
        bha_b = (bha_C80_ctr - bha_90) * (-1)
        bha_IPV = (bha_a + bha_b) / (peak_nha) * 299792
        bha_A = (bha_a - bha_b) / (bha_a + bha_b)                          # asymmetry     
        bha_S = (bha_a - bha_b) / (peak_nha) * 299792
        bha_K = (1.397*bha_FWHM) / ((bha_a + bha_b) / (bha_C80_ctr) * 299792) # kurtosis
        
        # Section 5: Plotting EVERYTHING onto one plot
        # ---
        fig2 = plt.figure(figsize=(16,10))
        plt.plot(wavelength, nl, c='#0B6E4F', linewidth=2.25, label='NL Profile')
        plt.plot(wavelength, bl, c='#CB2C2A', linewidth=2.25, label='BL Profile')
        plt.plot(wavelength, bl_profile, c='b', linewidth=2.25, label='Cleaned BL Profile')
        plt.plot(wavelength, data_contFeII_sub, c='k', linewidth=2.25, label='Data - Continuum - FeII')
        # Plotting corresponding NL peak wavelength and the three wavelengths of 
        # the calculated velocity shifts as vertical lines 
        plt.axvline(peak_nha, c='#146746', linestyle=':', linewidth=1.75, label='NL Peak')
        plt.axvline(peak_bha, c='#CB2C2A', linestyle=':', linewidth=1.75, label='BL Peak')
        plt.axvline(bha_ctr, c='#629FD0', linestyle='--', linewidth=1.75, label='BL Centroid')
        plt.axvline(bha_C80_ctr, c='#9B469B', linestyle='--', linewidth=1.75, label='BL C80')
        # Plot details
        plt.xlabel(r'$\rm Rest \, Wavelength$ ($\rm \AA$)', fontsize=16)
        plt.xlim(6170, 6900)                                  # match with bha_data
        plt.ylabel(r'$\rm f_{\lambda}$ ($\rm 10^{-17} erg\;s^{-1}\;cm^{-2}\;\AA^{-1}$)', fontsize=16)
        plt.title(sourcename+r': H$\alpha$ Line Complex', fontsize=25)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        plt.tick_params(which='major', length=12, width=1)
        plt.legend(fontsize=12)
        plt.savefig(path2+sourcename+'_BHA_ProfileCalc.pdf')
        
    
    # -----------------------------------------------------------------------------
    
    # H_beta Line Complex
    
    # If H_beta is not present in your spectrum (z > 0.89), 
    # the velocity shift calculations can be skipped by answering 'no' to "bhb_vs"
    
    bhb_vs = 'yes'
    if(bhb_vs == 'yes'):
        # Section 1: finds the peak velocity shift by obtaining the corresponding
        # wavelength of the max flux value of the BL profile
        # ---
        # BLs - change according to length of wings
        bhb_data = bl_data.loc[(bl_data['Wavelength'] >= 4640) & \
                               (bl_data['Wavelength'] <= 5150)]
        bhb_wave = bhb_data['Wavelength']
        bhb_flux = bhb_data['Flux']
        # Finding peak 
        peak_bhb_flux = bhb_flux.max()
        peak_bhb_wave = bhb_data.loc[bhb_data['Flux'] == peak_bhb_flux]
        peak_bhb = int(peak_bhb_wave['Wavelength'])
        # ---
        # NLs 
        nhb_data = nl_data.loc[(nl_data['Wavelength'] >= 4640) & \
                               (nl_data['Wavelength'] <= 4900)]
        nhb_flux = nhb_data['Flux']
        # Finding peak: should be ~ vacuum wavelength and will be used for 
        # calculating all three types of velocity shifts
        peak_nhb_flux = nhb_flux.max()
        peak_nhb_wave = nhb_data.loc[nhb_data['Flux'] == peak_nhb_flux]
        peak_nhb = int(peak_nhb_wave['Wavelength'])
        # ---
        # Peak velocity shift 
        bhb_pvs = ((peak_bhb**2 - peak_nhb**2) / (peak_bhb**2 + peak_nhb**2)) * 299792
        
        # Section 2: finds the centroid velocity shift by normalizing the BL 
        # profile and finding the cumulative distribution function (CDF) to obtain 
        # the corresponding wavelength at 50% of the CDF
        # ---
        # Normalizing BL profile and finding CDF
        bhb_flux_rescaled = (bhb_flux - bhb_flux.min()) / (bhb_flux.max() - bhb_flux.min())
        bhb_flux_area = integrate.trapezoid(bhb_flux_rescaled, bhb_wave)
        bhb_flux_new = bhb_flux_rescaled/bhb_flux_area   
        bhb_flux_norm = bhb_flux_new / np.sum(bhb_flux_new)
        bhb_cdf = np.cumsum(bhb_flux_norm)
        # ---
        # Finding centroid
        bhb_ctr = np.interp(0.5, bhb_cdf, bhb_wave)
        # Ploting CDF 
        plot_bhb_cdf = 'no'
        if(plot_bhb_cdf == 'yes'):
            fig4 = plt.figure(figsize=(16,10))
            plt.plot(bhb_wave, bhb_cdf, linewidth=2, c='#DB1D1A')
            plt.title(sourcename+': CDF', fontsize=30)
            plt.ylabel('Probability', fontsize=20)
            plt.xlabel(r'Wavelength ($\rm \AA$)', fontsize=20)
            plt.xticks(fontsize=20)
            plt.yticks(fontsize=20)
            plt.tick_params(which='major', length=16, width=2)
            plt.tick_params(which='minor', length=8, width=1)
            plt.savefig(path2+sourcename+'_BHA_CDF.pdf')
            plt.axvline(bhb_ctr, c='grey', linestyle='--', label='Centroid')
            plt.text(bhb_ctr+5, 0.2, r'H$\alpha$ BL Centroid = {:.2f} $\AA$'.format(bhb_ctr))
            plt.legend(fontsize=14)
        # ---
        # Centroid velocity shift
        bhb_cvs = ((bhb_ctr**2 - peak_nhb**2) / (bhb_ctr**2 + peak_nhb**2)) * 299792 
        
        # Section 3: finding C80 shift by removing 10% of the area at each end of
        # the BL profile and obtaining the corresponding wavelength at the 
        # center.
        # ---
        bhb_10 = np.interp(0.05, bhb_cdf, bhb_wave)
        bhb_90 = np.interp(0.95, bhb_cdf, bhb_wave)
        bhb_C80_data = bl_data.loc[(bl_data['Wavelength'] >= bhb_10) & \
                                   (bl_data['Wavelength'] <= bhb_90)]
        bhb_C80_wave = bhb_C80_data['Wavelength']
        bhb_C80_flux = bhb_C80_data['Flux']
        # Normalizing BL profile and finding CDF
        bhb_C80_flux_rescaled = (bhb_C80_flux - (bhb_C80_flux.min())) / \
            (bhb_C80_flux.max() - bhb_C80_flux.min())
        bhb_C80_area = integrate.trapezoid(bhb_C80_flux_rescaled, bhb_C80_wave)
        bhb_C80_flux_new = bhb_C80_flux_rescaled/bhb_C80_area
        bhb_C80_flux_norm = bhb_C80_flux_new / np.sum(bhb_C80_flux_new)
        bhb_cdf_C80 = np.cumsum(bhb_C80_flux_norm)
        # ---
        # Finding line center
        bhb_C80_ctr = np.interp(0.5, bhb_cdf_C80, bhb_C80_wave)
        # C80 velocity shift
        bhb_C80_vs = ((bhb_C80_ctr**2 - peak_nhb**2) / (bhb_C80_ctr**2 + peak_nhb**2)) * 299792 
        
        # Section 4: finding the Whittle 1985 line profile parameters; i.e. inter-
        # percentile velocity width (IPV), asymmetry (A), shift (S), and kurtosis
        # (K)
        # ---
        # Calculating FWHM
        bhb_FWHM_calc = peak_bhb_flux / 2
        bhb_FWHM_range = bhb_data.loc[(bhb_data['Flux'] >= bhb_FWHM_calc)]
        bhb_FWHM_wave_min = bhb_FWHM_range['Wavelength'].iloc[0]
        bhb_FWHM_wave_max = bhb_FWHM_range['Wavelength'].iloc[-1]
        bhb_FWHM_wave = bhb_FWHM_wave_max - bhb_FWHM_wave_min    
        bhb_FWHM = bhb_FWHM_wave / peak_nhb * 299792
        # ---
        # Calculating line parameters    
        bhb_a = bhb_C80_ctr - bhb_10
        bhb_b = (bhb_C80_ctr - bhb_90) * (-1)
        bhb_IPV = (bhb_a + bhb_b) / (peak_nhb) * 299792
        bhb_A = (bhb_a - bhb_b) / (bhb_a + bhb_b)                          # asymmetry     
        bhb_S = (bhb_a - bhb_b) / (peak_nhb) * 299792
        bhb_K = (1.397*bhb_FWHM) / ((bhb_a + bhb_b) / (bhb_C80_ctr) * 299792) # kurtosis
        
        # Section 5: Plotting EVERYTHING onto one plot
        # ---
        fig5 = plt.figure(figsize=(16,10))
        plt.plot(wavelength, nl, c='#0B6E4F', linewidth=2.25, label='NL Profile')
        plt.plot(wavelength, bl, c='#CB2C2A', linewidth=2.25, label='BL Profile')
        plt.plot(wavelength, bl_profile, c='b', linewidth=2.25, label='Cleaned BL Profile')
        plt.plot(wavelength, data_contFeII_sub, c='k', linewidth=2.25, label='Data - Continuum - FeII')
        # Plotting corresponding NL peak wavelength and the three wavelengths of 
        # the calculated velocity shifts as vertical lines 
        plt.axvline(peak_nhb, c='#146746', linestyle=':', linewidth=1.75, label='NL Peak')
        plt.axvline(peak_bhb, c='#CB2C2A', linestyle=':', linewidth=1.75, label='BL Peak')
        plt.axvline(bhb_ctr, c='#629FD0', linestyle='--', linewidth=1.75, label='BL Centroid')
        plt.axvline(bhb_C80_ctr, c='#9B469B', linestyle='--', linewidth=1.75, label='BL C80')
        # Plot details
        plt.xlabel(r'$\rm Rest \, Wavelength$ ($\rm \AA$)', fontsize=16)
        plt.xlim(4640, 5150)                                  # match with bhb_data
        plt.ylabel(r'$\rm f_{\lambda}$ ($\rm 10^{-17} erg\;s^{-1}\;cm^{-2}\;\AA^{-1}$)', fontsize=16)
        plt.title(sourcename+r': H$\beta$ Line Complex', fontsize=25)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        plt.tick_params(which='major', length=12, width=1)
        plt.legend(fontsize=12)
        plt.savefig(path2+sourcename+'_BHB_ProfileCalc.pdf')
        
    
    
    # -----------------------------------------------------------------------------
    
    # Velocity Plots
    
    # This part of the code will overplot the normalized BL profiles on a velocity
    # scale for direct comparison
    
    plot_vel = 'yes'
    if(plot_vel == 'yes'):
        # Converting wavelengths to velocities
        bha_vel = (bha_wave - peak_nha) / (peak_nha) * 299792
        bhb_vel = (bhb_wave - peak_nhb) / (peak_nhb) * 299792
        fig3 = plt.figure(figsize=(16,10))
        plt.plot(bha_vel, bha_flux_norm, c='purple', linewidth=2, label=r'H$\alpha$ BL Profile')
        plt.plot(bhb_vel, bhb_flux_norm, c='darkgoldenrod', linewidth=2, label=r'H$\beta$ BL Profile')
        plt.axvline(0, c='k', linestyle=(0, (5, 5)), linewidth=1)
        # Plot details
        plt.xlabel(r'Velocity (km s$^{-1}$)', fontsize=16)                  
        plt.ylabel('Normalized Flux', fontsize=16)
        plt.title(sourcename, fontsize=25)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        plt.tick_params(which='major', length=12, width=1)
        plt.legend(fontsize=12)
        plt.savefig(path2+sourcename+'_BL_VProfile.pdf')
    
    
    # -----------------------------------------------------------------------------
    
    # Saving calculated velocity shifts and Whittle 1985 profile parameters into a 
    # separate .txt file
    
    with open(path2+sourcename+"_BLV-Shifts.txt","w") as f:
        print('Calculated BL Velocity Shifts and Line Profile Parameters:', file=f)
        print(" ", file=f)
        print("   Peak H_alpha Velocity Shift =", '%.2f' % bha_pvs, 'km s^-1', file=f)
        print("   Center H_alpha Velocity Shift =", '%.2f' % bha_cvs, 'km s^-1', file=f)
        print("   C80 H_alpha Velocity Shift =", '%.2f' % bha_C80_vs, 'km s^-1', file=f)
        print("   H_alpha IPV width =", '%.2f' % bha_IPV, file=f)
        print("   H_alpha Asymmetry =", '%.5f' % bha_A, file=f)
        print("   H_alpha Shift =", '%.2f' % bha_S, file=f)
        print("   H_alpha Kurtosis =", '%.5f' % bha_K, file=f)
        print(" ", file=f)
        print("   Peak H_beta Velocity Shift =", '%.2f' % bhb_pvs, 'km s^-1', file=f)
        print("   Center H_beta Velocity Shift =", '%.2f' % bhb_cvs, 'km s^-1', file=f)
        print("   C80 H_beta Velocity Shift =", '%.2f' % bhb_C80_vs, 'km s^-1', file=f)
        print("   H_beta IPV width =", '%.2f' % bhb_IPV, file=f)
        print("   H_beta Asymmetry =", '%.5f' % bhb_A, file=f)
        print("   H_beta Shift =", '%.2f' % bhb_S, file=f)
        print("   H_beta Kurtosis =", '%.5f' % bhb_K, file=f)
        
# ----------------------------------------------------------------------------------
# Deciding to loop through all values or not
loop = False
if loop:
    for source in Data_list:
        lineprofilecalc_func(source)
else:
    lineprofilecalc_func('0332-52367-0639')
    
    
    