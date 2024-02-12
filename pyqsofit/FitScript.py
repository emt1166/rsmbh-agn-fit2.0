#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 11:34:49 2024

@author: emilytemple
"""

import os,timeit
import numpy as np
from PyQSOFit import QSOFit
from astropy.io import fits
import matplotlib.pyplot as plt
import warnings
import re

# Use custom matplotlib style to make Yue happy
QSOFit.set_mpl_style()
# Ignore warnings?
warnings.filterwarnings("ignore")

# Setting the file paths that the code uses
# The path of the source code file and qsopar.fits
path1 = '/Users/emilytemple/documents/rsmbh-agn-fit2.0/pyqsofit/'
# The path of fitting results - can customize, I created a new directory
path2 = '/Users/emilytemple/documents/rsmbh-agn-fit2.0/pyqsofit/Fit Results/'     
# The path of fitted figure - same as above      
path3 = '/Users/emilytemple/documents/rsmbh-agn-fit2.0/pyqsofit/Fit Results/QA Other/'
# The path of dust reddening map
path4 = '/Users/emilytemple/documents/rsmbh-agn-fit2.0/pyqsofit/sfddata/'
# EMT: fitting control or outliers, for organization purposes  
pathC = path2+'/Line Properties/Control/'
pathO = path2+'/Line Properties/Outliers/'

# -----------------------------------------------------------------------------
# Allowing a loop through directories to make this process faster
# Once you download all the spectral files from SDSS, this should be able to 
# loop through all the files, create their directories, and save results

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


# -----------------------------------------------------------------------------
#making the following a function
def fittingscript_func(sourcename):
    
    # Opening spectrum to be fitted. NOTE: SDSS fits files are saved as
    # spec-plateID-MJD-fiberID.fits. Sourcename is just the plateID-MJD-fiberID
    #sourcename = '0646-52523-0075'
    spec = 'spec-'+sourcename+'.fits'
    data = fits.open(os.path.join(path1+'Data/Outliers/'+spec))
    lam = 10**data[1].data['loglam']                           # OBS wavelength (A)
    flux = data[1].data['flux']                           # OBS flux (erg/s/cm^2/A)
    err = 1./np.sqrt(data[1].data['ivar'])                          # 1 sigma error
    z = data[2].data['z'][0]                                             # Redshift
    
    # Optional information... 
    ra = data[0].header['plug_ra']                                             # RA 
    dec = data[0].header['plug_dec']                                          # DEC
    plateid = data[0].header['plateid']                             # SDSS plate ID
    mjd = data[0].header['mjd']                                          # SDSS MJD
    fiberid = data[0].header['fiberid']                             # SDSS fiber ID
    
    
    
    #-----------------------------------------------------------------------------
    #Creating directories automatically
    #NOTE: 'Run.py' should be run BEFORE this
    
    if os.path.exists(pathO+sourcename):
        if not os.path.exists(pathO+sourcename+'/Fit Data/'):
            os.makedirs(pathO+sourcename+'/Fit Data/')
        if not os.path.exists(pathO+sourcename+'/Line Profile Plots/'):
            os.makedirs(pathO+sourcename+'/Line Profile Plots/')
    
    # -----------------------------------------------------------------------------
    
    # PyQSOFit - fitting portion
    # Preparing spectrum data
    q_mle = QSOFit(lam, flux, err, z, ra=ra, dec=dec, plateid=plateid, mjd=mjd, 
                   fiberid=fiberid, path=path1)
    start = timeit.default_timer()
    
    # Do the fitting. NOTE: Change arguments accordingly
    # This has been edited to reflect version 2.0
    
    q_mle.Fit(name=None, nsmooth=1, and_mask=False, or_mask=False, reject_badpix=True, deredden=True,
          wave_range=None,  wave_mask=None, decompose_host=True, host_prior=True, host_prior_scale=0.2,
          host_line_mask=True, decomp_na_mask=True, qso_type='CZBIN1',  
          npca_qso=10, host_type='PCA', npca_gal=5,
          Fe_uv_op=True, poly=True, BC=False, initial_guess=None,  
          rej_abs_conti=False, n_pix_min_conti=100, linefit=True,  rej_abs_line=False,
          MC=False, MCMC=False, nsamp=200, param_file_name='qsopar.fits', nburn=20, nthin=10, 
          epsilon_jitter=0., save_result=False, save_fits_name=None, save_fits_path=path2,  
          plot_fig=True, save_fig=False, plot_corner=False, verbose=True,  
          kwargs_plot={
              'save_fig_path': '.',  
              'broad_fwhm'   : 1200  # km/s, lower limit that code decide if a line component belongs to broad component
          },
          kwargs_conti_emcee={},
          kwargs_line_emcee={})

    end = timeit.default_timer()
    print('Fitting finished in : '+str(np.round(end-start))+'s')
    
    
    # -----------------------------------------------------------------------------
    
    # Obtaining fit result files saved under QA Other subdirectory in Fit Results
    # source code directory
    
    # data = fits.open(path3+source+'.fits') 
    
    
    # -----------------------------------------------------------------------------
    
    # Saving each line complex individually because PyQSOFit has been modified to 
    # not plot these under the full fit figure; therefore, we plot them here 
    # separately. In addition, we only plot H_alpha, H_beta, and MgII line 
    # complexes since those will be used for finding velocity shifts
    
    # Path of line complex plots
    path5 = pathO+sourcename+'/'
                       
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
        plt.savefig(path5+sourcename+'_BroadHa_LineComplex.pdf')
    
    
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
        plt.savefig(path5+sourcename+'_BroadHb_LineComplex.pdf')
    
    # For high z, MgII lines appear and can be used for calculating velocity shifts
    plot_MgII = 'no'                  
    if(plot_MgII =='yes'):
        fig4 = plt.figure(figsize=(16,12))
        for p in range(len(q_mle.gauss_result)//3):
            if q_mle.CalFWHM(q_mle.gauss_result[3*p+2]) < 1200: 
                color = 'g' # narrow
            else:
                color = 'r' # broad
            plt.plot(q_mle.wave, q_mle.Onegauss(np.log(q_mle.wave), 
            q_mle.gauss_result[p*3:(p+1)*3]), color=color) 
        # Plot total line model
        plt.plot(q_mle.wave, q_mle.Manygauss(np.log(q_mle.wave), q_mle.gauss_result),
                 'b', lw=2)
        plt.plot(q_mle.wave, q_mle.line_flux,'k')
        plt.xlim(2700, 3000)
        plt.xlabel(r'$\rm Rest \, Wavelength$ ($\rm \AA$)', fontsize = 20)
        plt.ylabel(r'$\rm f_{\lambda}$ ($\rm 10^{-17} erg\;s^{-1}\;cm^{-2}\;\AA^{-1}$)'
                       , fontsize = 20)
        plt.title(r'Broad MgII', fontsize = 30)
        plt.savefig(path5+sourcename+'_BroadMgII_LineComplex.pdf')
    
    
    
    # -----------------------------------------------------------------------------    
    
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
    # Broad H_gamma
    fwhm_bhg, sigma_bhg, ew_bhg, peak_bhg, area_bhg, snr_bhg = q_mle.line_prop_from_name(
        'Hg_br', 'broad')
    # Broad H_delta
    fwhm_bhd, sigma_bhd, ew_bhd, peak_bhd, area_bhd, snr_bhd = q_mle.line_prop_from_name(
        'Hd_br', 'broad')
    # [OIII]5007 - core
    fwhm_oIII5, sigma_oIII5, ew_oIII5, peak_oIII5, area_oIII5, snr_oIII5 = q_mle.line_prop_from_name(
        'OIII5007c', 'narrow')
    # [OIII]4959 - core
    fwhm_oIII4, sigma_oIII4, ew_oIII4, peak_oIII4, area_oIII4, snr_oIII4 = q_mle.line_prop_from_name(
        'OIII4959c', 'narrow')
    # Broad MgII
    fwhm_bmgII, sigma_bmgII, ew_bmgII, peak_bmgII, area_bmgII, snr_bmgII = q_mle.line_prop_from_name(
        'MgII_br', 'broad')
    # Narrow MgII
    fwhm_nmgII, sigma_nmgII, ew_nmgII, peak_nmgII, area_nmgII, snr_nmgII = q_mle.line_prop_from_name(
        'MgII_na', 'narrow')
    
    with open(path5+sourcename+"_LineProperties.txt","w") as f:
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
        # Broad H_gamma
        print('Broad H_gamma:', file=f)
        print(' '  + "FWHM = " '%.2f' % fwhm_bhg, "km s^-1", file=f)
        print(' '  + "sigma = " '%.2f' % sigma_bhg, "km s^-1", file=f)
        print(' '  + "EW = ", '%.2f' % ew_bhg, "Ang", file=f)
        print(' '  + "Peak = ", '%.2f' % peak_bhg, "Ang", file=f)
        print(' '  + "Area = ", '%.2f' % area_bhg, 
        "x 10^-17 erg s^-1 cm^-2", file=f)
        # Broad H_delta
        print('Broad H_delta:', file=f)
        print(' '  + "FWHM = " '%.2f' % fwhm_bhd, "km s^-1", file=f)
        print(' '  + "sigma = " '%.2f' % sigma_bhd, "km s^-1", file=f)
        print(' '  + "EW = ", '%.2f' % ew_bhd, "Ang", file=f)
        print(' '  + "Peak = ", '%.2f' % peak_bhd, "Ang", file=f)
        print(' '  + "Area = ", '%.2f' % area_bhd, 
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
        # Broad MgII
        print('Broad MgII:', file=f)
        print(' '  + "FWHM = " '%.2f' % fwhm_bmgII, "km s^-1", file=f)
        print(' '  + "sigma = " '%.2f' % sigma_bmgII, "km s^-1", file=f)
        print(' '  + "EW = ", '%.2f' % ew_bmgII, "Ang", file=f)
        print(' '  + "Peak = ", '%.2f' % peak_bmgII, "Ang", file=f)
        print(' '  + "Area = ", '%.2f' % area_bmgII, 
        "x 10^-17 erg s^-1 cm^-2", file=f)
        # Narrow MgII
        print('Narrow MgII:', file=f)
        print(' '  + "FWHM = " '%.2f' % fwhm_nmgII, "km s^-1", file=f)
        print(' '  + "sigma = " '%.2f' % sigma_nmgII, "km s^-1", file=f)
        print(' '  + "EW = ", '%.2f' % ew_nmgII, "Ang", file=f)
        print(' '  + "Peak = ", '%.2f' % peak_nmgII, "Ang", file=f)
        print(' '  + "Area = ", '%.2f' % area_nmgII, 
        "x 10^-17 erg s^-1 cm^-2", file=f)
          
        
    # -----------------------------------------------------------------------------
    
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
        plt.savefig(path5+sourcename+'_BLProfile.pdf')
        
        # Saving subtracted data into file to use for further calculations
        np.save(path5+'/Fit Data/'+sourcename+'_DataCFeII', data_contFeII_sub)
        np.save(path5+'/Fit Data/'+sourcename+'_Data', data)
        np.save(path5+'/Fit Data/'+sourcename+'_Wavelength', wavelength)
        np.save(path5+'/Fit Data/'+sourcename+'_BLSpectrum', data_sub)
        np.save(path5+'/Fit Data/'+sourcename+'_NLData', n_lines)
        np.save(path5+'/Fit Data/'+sourcename+'_BLData', b_lines)
    
    return
# -----------------------------------------------------------------------------

#Now we can choose to loop through it
# if loop is false, then you can specify the spectrum you want

loop = False
if loop: 
    for source in Data_list:
        fittingscript_func(source)
else:
    fittingscript_func('0332-52367-0639')    
    
# -----------------------------------------------------------------------------
# Want a way to fit a spectrum multiple times for comparison purposes
# will probably write an entirely new script tbh 

    
    






