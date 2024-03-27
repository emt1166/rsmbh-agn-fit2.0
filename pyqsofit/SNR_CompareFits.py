#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 12:52:10 2024

@author: emilytemple

This script will be used to compare the line profile calculations between successive iterations. 
Each iteration will have a changing signal-to-noise using the error in the fits files from SDSS. 

The first way I will implement this, is by using random.gauss() to pick a random value within the
one sigma error to "replace" the flux values with, making the spectrum noisier. 


This script is meant to be used with PyQSOFit, but has modifications for my purposes. 
I will likely be plotting the iterations on top of each other.
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
import time

# Use custom matplotlib style to make Yue happy
QSOFit.set_mpl_style()
# Ignore warnings?
warnings.filterwarnings("ignore")

# -----------------------------------------------------------------------------------------
#choose which object to look at
#sourcename = '0813-52354-0561' #Should be taken care of via loop now
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
def Fit_Function(sourcename):
    
    #Creating directories automatically from the sourcename
    #if not done so already, should be done already tbh
    
    if not os.path.exists(pathF+sourcename):
        os.makedirs(pathF+sourcename)
       
    #path to save ALL outputs, everything goes here, no more separate folders
    path2 = pathF+sourcename+'/'
    
    spec = 'spec-'+sourcename+'.fits'
    data = fits.open(os.path.join(path1+'Data/Outliers/'+spec))
    lam = 10**data[1].data['loglam']                           # OBS wavelength (A)
    
    #FLUX data will be changing due to input of noise via new_flux
    flux = data[1].data['flux']                           # OBS flux (erg/s/cm^2/A)
    err = 1./np.sqrt(data[1].data['ivar'])                          # 1 sigma error
    z = data[2].data['z'][0]                                             # Redshift
    
    # Optional information... 
    ra = data[0].header['plug_ra']                                             # RA 
    dec = data[0].header['plug_dec']                                          # DEC
    plateid = data[0].header['plateid']                             # SDSS plate ID
    mjd = data[0].header['mjd']                                          # SDSS MJD
    fiberid = data[0].header['fiberid']                             # SDSS fiber ID
    
    
    # ------------------------------------------------------------------------------------------
    # FITTING, NEED TO LOOP OVER THE SAME SOURCE (SNR COMPARE)
    # PyQSOFit - fitting portion
    
    # Preparing spectrum data, change FLUX as needed for SNR stuff 
    n = 3
    for n in range(5):
        new_flux = np.loadtxt(pathF+sourcename+'/new_flux/'+sourcename+'.'+f'{n}'+'new_flux.txt')
        
        q_mle = QSOFit(lam, new_flux, err, z, ra=ra, dec=dec, plateid=plateid, mjd=mjd, 
                        fiberid=fiberid, path=path1)
        start = timeit.default_timer()
        
        # Do the fitting. NOTE: Change any arguments accordingly
        # This has been edited to reflect version 2.0
        # the host is decomposed in this fit 
        
        q_mle.Fit(name=None, nsmooth=1, and_mask=False, or_mask=False, reject_badpix=True, deredden=True,
              wave_range=None,  wave_mask=None, decompose_host=True, host_prior=True, host_prior_scale=0.2,
              host_line_mask=True, decomp_na_mask=True, qso_type='CZBIN1',  
              npca_qso=10, host_type='PCA', npca_gal=5,
              Fe_uv_op=True, poly=True, BC=False, initial_guess=None,  
              rej_abs_conti=False, n_pix_min_conti=100, linefit=True,  rej_abs_line=False,
              MC=False, MCMC=False, nsamp=200, param_file_name='qsopar.fits', nburn=20, nthin=10, 
              epsilon_jitter=0., save_result=False, save_fits_name=None, save_fits_path=path3,  
              plot_fig=True, save_fig=False, plot_corner=False, verbose=True,  
              kwargs_plot={
                  'save_fig_path': '.',  
                  'broad_fwhm'   : 1200  # km/s, lower limit that code decide if a line component belongs to broad component
              },
              kwargs_conti_emcee={},
              kwargs_line_emcee={})
        
        end = timeit.default_timer()
        print('Fitting finished in : '+str(np.round(end-start))+'s')
        
        
        # Plotting the figures -----------------------------------------------------------------------
        # Path of line complex plots for SNR Comparison 
         
        # Plotting H_alpha coomplex
        plot_Ha = 'yes'  
        if(plot_Ha =='yes'):
            # Plotting broad H_alpha, NII, and SII line complex
            fig1 = plt.figure(figsize=(16,12))
            for p in range(len(q_mle.gauss_result)//3):
                if q_mle.CalFWHM(q_mle.gauss_result[3*p+2]) < 1000:          # FWHM max
                    color = 'g' # narrow
                else:
                    color = 'r' # broad
                plt.plot(q_mle.wave, q_mle.Onegauss(np.log(q_mle.wave), q_mle.gauss_result
                                                [p*3:(p+1)*3]), color=color)  
            # Plot total line model
            plt.plot(q_mle.wave, q_mle.Manygauss(np.log(q_mle.wave), q_mle.gauss_result), 
                              'b', lw=2)
            plt.plot(q_mle.wave, q_mle.line_flux,'k')
            plt.xlim(6300, 7000)
            plt.xlabel(r'$\rm Rest \, Wavelength$ ($\rm \AA$)', fontsize = 20)
            plt.ylabel(r'$\rm f_{\lambda}$ ($\rm 10^{-17} erg\;s^{-1}\;cm^{-2}\;\AA^{-1}$)'
                        , fontsize = 20)
            plt.title(r'Broad $H{\alpha}+[NII]+[SII]$', fontsize = 30)
            #plt.savefig(path2+sourcename+'_BroadHa_LineComplex.pdf')
        
        
        # Plotting broad H_beta and [OIII] line complex
        # EMT Note: Need caveat due to [OIII] wings frequently fit as BLR, not NLR... --------
        plot_Hb = 'yes'
        if(plot_Hb =='yes'):
            fig2 = plt.figure(figsize=(16,12))
            for p in range(len(q_mle.gauss_result)//3):
                if q_mle.CalFWHM(q_mle.gauss_result[3*p+2]) < 1000:  
                    color = 'g' # narrow
                else:
                    color = 'r' # broad
                plt.plot(q_mle.wave, q_mle.Onegauss(np.log(q_mle.wave), q_mle.gauss_result
                                                [p*3:(p+1)*3]), color=color)
            # Plot total line model
            plt.plot(q_mle.wave, q_mle.Manygauss(np.log(q_mle.wave), q_mle.gauss_result), 
                      'b', lw=2)
            plt.plot(q_mle.wave, q_mle.line_flux,'k')
            plt.xlim(4640, 5300)
            plt.xlabel(r'$\rm Rest \, Wavelength$ ($\rm \AA$)', fontsize = 20)
            plt.ylabel(r'$\rm f_{\lambda}$ ($\rm 10^{-17} erg\;s^{-1}\;cm^{-2}\;\AA^{-1}$)'
                        , fontsize = 20)
            plt.title(r'Broad $H{\beta}+[OIII]$', fontsize = 30)
            #plt.savefig(path2+sourcename+'_BroadHb_LineComplex.pdf')
            
            
        # -----------------------------------------------------------------------------------    
        # PyQSOFit calculates FWHM, Sigma, EW, Area, and SNR for each broad and narrow
        # component of emission lines. That information is obtained and then saved to a 
        # separate .txt file
        
        # Broad H_alpha
        fwhm_bha, sigma_bha, ew_bha, peak_bha, area_bha, snr_bha = q_mle.line_prop_from_name(
            'Ha_br','broad')
        # Narrow H_alpha
        fwhm_nha, sigma_nha, ew_nha, peak_nha, area_nha, snr_nha = q_mle.line_prop_from_name(
            'Ha_na', 'narrow')
        # Broad H_beta 
        fwhm_bhb, sigma_bhb, ew_bhb, peak_bhb, area_bhb, snr_bhb = q_mle.line_prop_from_name(
            'Hb_br', 'broad')
        # Narrow H_beta
        fwhm_nhb, sigma_nhb, ew_nhb, peak_nhb, area_nhb, snr_nhb = q_mle.line_prop_from_name(
            'Hb_na', 'narrow')
        # [OIII]5007 - core
        fwhm_oIII5, sigma_oIII5, ew_oIII5, peak_oIII5, area_oIII5, snr_oIII5 = q_mle.line_prop_from_name(
            'OIII5007c', 'narrow')
        # [OIII]4959 - core
        fwhm_oIII4, sigma_oIII4, ew_oIII4, peak_oIII4, area_oIII4, snr_oIII4 = q_mle.line_prop_from_name(
            'OIII4959c', 'narrow')
        
        #saving params identified by n 
        with open(path2+sourcename+'.'+f'{n}'+"_LineProperties.txt","w") as f:
            print('PyQSOFit Calculated Line Properties', file=f)
            print('', file=f)
            # Broad H_alpha
            print('Broad H_alpha:', file=f)
            print(' '  + "FWHM = " '%.2f' % fwhm_bha, "km s^-1", file=f)
            print(' '  + "sigma = " '%.2f' % sigma_bha, "km s^-1", file=f)
            print(' '  + "EW = ", '%.2f' % ew_bha, "Ang", file=f)
            print(' '  + "Peak = ", '%.2f' % peak_bha, "Ang", file=f)
            print(' '  + "Area = ", '%.2f' % area_bha, 
            "x 10^-17 erg s^-1 cm^-2", file=f)
            # Narrow H_alpha
            print('Narrow H_alpha:', file=f)
            print(' '  + "FWHM = " '%.2f' % fwhm_nha, "km s^-1", file=f)
            print(' '  + "sigma = " '%.2f' % sigma_nha, "km s^-1", file=f)
            print(' '  + "EW = ", '%.2f' % ew_nha, "Ang", file=f)
            print(' '  + "Peak = ", '%.2f' % peak_nha, "Ang", file=f)        
            print(' '  + "Area = ", '%.2f' % area_nha, 
            "x 10^-17 erg s^-1 cm^-2", file=f)
            # Broad H_beta
            print('Broad H_beta:', file=f)
            print(' '  + "FWHM = " '%.2f' % fwhm_bhb, "km s^-1", file=f)
            print(' '  + "sigma = " '%.2f' % sigma_bhb, "km s^-1", file=f)
            print(' '  + "EW = ", '%.2f' % ew_bhb, "Ang", file=f)
            print(' '  + "Peak = ", '%.2f' % peak_bhb, "Ang", file=f)
            print(' '  + "Area = ", '%.2f' % area_bhb, 
            "x 10^-17 erg s^-1 cm^-2", file=f)
            # Narrow H_beta
            print('Narrow H_beta:', file=f)
            print(' '  + "FWHM = " '%.2f' % fwhm_nhb, "km s^-1", file=f)
            print(' '  + "sigma = " '%.2f' % sigma_nhb, "km s^-1", file=f)
            print(' '  + "EW = ", '%.2f' % ew_nhb, "Ang", file=f)
            print(' '  + "Peak = ", '%.2f' % peak_nhb, "Ang", file=f)        
            print(' '  + "Area = ", '%.2f' % area_nhb, 
            "x 10^-17 erg s^-1 cm^-2", file=f)
            # Narrow [OIII]5007
            print('Narrow [OIII]5007:', file=f)
            print(' '  + "FWHM = " '%.2f' % fwhm_oIII5, "km s^-1", file=f)
            print(' '  + "sigma = " '%.2f' % sigma_oIII5, "km s^-1", file=f)
            print(' '  + "EW = ", '%.2f' % ew_oIII5, "Ang", file=f)
            print(' '  + "Peak = ", '%.2f' % peak_oIII5, "Ang", file=f)
            print(' '  + "Area = ", '%.2f' % area_oIII5, 
            "x 10^-17 erg s^-1 cm^-2", file=f)
            # Narrow [OIII]4959
            print('Narrow [OIII]4959:', file=f)
            print(' '  + "FWHM = " '%.2f' % fwhm_oIII4, "km s^-1", file=f)
            print(' '  + "sigma = " '%.2f' % sigma_oIII4, "km s^-1", file=f)
            print(' '  + "EW = ", '%.2f' % ew_oIII4, "Ang", file=f)
            print(' '  + "Peak = ", '%.2f' % peak_oIII4, "Ang", file=f)
            print(' '  + "Area = ", '%.2f' % area_oIII4, 
            "x 10^-17 erg s^-1 cm^-2", file=f)
    
              
        #--------------------------------------------------------------------------------------
        # Cleaing data through model subtraction to obtain a clean BL profile of the
        # spectrum. This part of the code is optional
        
        # Data subtraction and plotting 
        data_subtraction ='yes'                             
        if(data_subtraction == 'yes'):
            # Obtaining narrow lines from the fitted spectrum
            n_lines = np.zeros(len(q_mle.wave))
            for p in range(len(q_mle.gauss_result)//3):
                if q_mle.CalFWHM(q_mle.gauss_result[3*p+2]) < 1000: 
                    na = q_mle.Onegauss(np.log(q_mle.wave), q_mle.gauss_result[p*3:(p+1)*3])
                    n_lines = n_lines + na
                    
            # Obtaining broad lines from the fitted spectrum
            b_lines = np.zeros(len(q_mle.wave))
            for p in range(len(q_mle.gauss_result)//3):
                if q_mle.CalFWHM(q_mle.gauss_result[3*p+2]) > 1000: 
                    ba = q_mle.Onegauss(np.log(q_mle.wave), q_mle.gauss_result[p*3:(p+1)*3])
                    b_lines = b_lines + ba
            
            # Calling the separate models from the fit
            # NOTE: FLUX used here ----------------->
            data = q_mle.flux                               # Flux from SDSS .fits file
            continuum_FeII = q_mle.f_conti_model            # FeII template + continuum
            wavelength = q_mle.wave
            
            # Skip the error results before obtaining fitted line flux
            if q_mle.MCMC == True:
                gauss_result = q_mle.gauss_result[::2]
            else:
                gauss_result = q_mle.gauss_result
                
            line = q_mle.Manygauss(np.log(q_mle.wave), gauss_result) + q_mle.f_conti_model
        
            # Performing data subtraction
            data_contFeII_sub = data - continuum_FeII
            data_sub = data - continuum_FeII - n_lines
        
            # Plotting cleaned data
            fig6 = plt.figure(figsize=(15,5))
            plt.plot(wavelength, data_sub, c='k', label='BL Profile')
            plt.xlabel(r'$\rm Rest \, Wavelength$ ($\rm \AA$)', fontsize=20)
            plt.ylabel(r'$\rm f_{\lambda}$ ($\rm 10^{-17} erg\;s^{-1}\;cm^{-2}\;\AA^{-1}$)', fontsize=20)
            plt.title(f'ra,dec = ({np.round(ra, 4)},{np.round(dec, 4)})   {sourcename}   z = {np.round(float(z), 4)}',
                  fontsize=20)
            plt.legend()
            #plt.savefig(path2+sourcename+'_BLProfile.pdf')
            
            # Saving subtracted data into file to use for further calculations
            np.save(path2+sourcename+'_DataCFeII', data_contFeII_sub)
            np.save(path2+sourcename+'_Data', data)
            np.save(path2+sourcename+'_Wavelength', wavelength)
            np.save(path2+sourcename+'_BLSpectrum', data_sub)
            np.save(path2+sourcename+'_NLData', n_lines)
            np.save(path2+sourcename+'_BLData', b_lines)
        
        #---------------------------------------------------------------------------------------------
        # LINE PROFILE CALC
        path = path2
        # Obtaining line profile data result components from Fitting-Script.py
        # We will do this for the new_flux and then pull from the original spec folders to plot
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
            # EMT Note: need to figure out way to pick out NHa specifically... not NII....------------->
            nha_data = nl_data.loc[(nl_data['Wavelength'] >= 6500) & \
                                    (nl_data['Wavelength'] <=6600)]
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
                #plt.savefig(path2+sourcename+'_BHA_CDF.pdf')
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
            #plt.savefig(path2+sourcename+'_BHA_ProfileCalc.pdf')
            
        
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
                #plt.savefig(path2+sourcename+'_BHA_CDF.pdf')
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
            #plt.savefig(path2+sourcename+'_BHB_ProfileCalc.pdf')
            
        
        
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
            #plt.savefig(path2+sourcename+'_BL_VProfile.pdf')
        
        # Saving calculated velocity shifts and Whittle 1985 profile parameters into a 
        # separate .txt file, identified by n
        
        with open(path2+sourcename+'.'+f'{n}'+"_BLV-Shifts.txt","w") as f:
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
            
    # -----------------------------------------------------------------------------
    return 


# -----------------------------------------------------------------------------------------
# Allow choice of loop
loop = True
start_time = time.time()

if loop:
    for source in Data_list:
        Fit_Function(source)
else:
    Fit_Function('0719-52203-0092')
end_time=time.time()

elapsed_time = end_time - start_time
print("End time:", time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(end_time)))
print("Total loop runtime:", round(elapsed_time, 3),'seconds')



