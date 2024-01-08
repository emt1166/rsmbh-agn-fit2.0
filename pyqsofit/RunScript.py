#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 10:58:37 2024

@author: emilytemple
"""

# This version will execute the fit as done in RunScript.py, but with the 
# addition of calculating the broad and narrow line properties. They can be
# changed accoriding to q.line_result_name

import re
import os,timeit
import numpy as np
from PyQSOFit import QSOFit
from astropy.io import fits
import matplotlib.pyplot as plt
import warnings
from ppxf.ppxf import ppxf

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

# -----------------------------------------------------------------------------
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


# Making a function for this script so we can loop over it
def runscript_func(sourcename):
        
    
    # -----------------------------------------------------------------------------
    #Creating directories automatically from the sourcename
    #if not done so already
    
    
    if not os.path.exists(path2+'Line Properties/'+sourcename):
        os.makedirs(path2+'Line Properties/'+sourcename)
       
    
    # -----------------------------------------------------------------------------
    
    # Opening spectrum to be fitted
    # NOTE: Remember to change data line, or else you will be fitting the previous
    # fits file used
    
    data = fits.open(os.path.join(path1+'Data/spec-'+sourcename+'.fits'))
    lam = 10**data[1].data['loglam']                           # OBS wavelength (A)
    flux = data[1].data['flux']                           # OBS flux (erg/s/cm^2/A)
    err = 1./np.sqrt(data[1].data['ivar'])                          # 1 sigma error
    z = data[2].data['z'][0]                                             # Redshift
    
    # Optional
    ra = data[0].header['plug_ra']                                             # RA 
    dec = data[0].header['plug_dec']                                          # DEC
    plateid = data[0].header['plateid']                             # SDSS plate ID
    mjd = data[0].header['mjd']                                          # SDSS MJD
    fiberid = data[0].header['fiberid']                             # SDSS fiber ID
    
    
    # -----------------------------------------------------------------------------
    
    # Fitting the spectrum of data - line 59 in PyQSOFit.py
    # Prepare data
    q_mle = QSOFit(lam, flux, err, z, ra=ra, dec=dec, plateid=plateid, mjd=mjd, fiberid=fiberid, path=path1)

    # Double check the installation path with the PCA / Fe template files
    # print('install path:', q_mle.install_path)

    # Change it if you installed them somewhere else
    #q_mle.install_path = '...'

    start = timeit.default_timer()
    # Do the fitting
    start = timeit.default_timer()
    # Do the fitting

    q_mle.Fit(name=None,  # customize the name of given targets. Default: plate-mjd-fiber
          # prepocessing parameters
          nsmooth=1,  # do n-pixel smoothing to the raw input flux and err spectra
          and_mask=False,  # delete the and masked pixels
          or_mask=False,  # delete the or masked pixels
          reject_badpix=True,  # reject 10 most possible outliers by the test of pointDistGESD
          deredden=True,  # correct the Galactic extinction
          wave_range=None,  # trim input wavelength
          wave_mask=None,  # 2-D array, mask the given range(s)

          # host decomposition parameters
          decompose_host=True,  # If True, the host galaxy-QSO decomposition will be applied
          host_prior=True, # If True, the code will adopt prior-informed method to assist decomposition. Currently, only 'CZBIN1' and 'DZBIN1' model for QSO PCA are available. And the model for galaxy must be PCA too.
          host_prior_scale=0.2, # scale of prior panelty. Usually, 0.2 works fine for SDSS spectra. Adjust it smaller if you find the prior affect the fitting results too much.

          host_line_mask=True, # If True, the line region of galaxy will be masked when subtracted from original spectra.
          decomp_na_mask=True, # If True, the narrow line region will be masked when perform decomposition
          qso_type='CZBIN1', # PCA template name for quasar
          npca_qso=10, # numebr of quasar templates
          host_type='PCA', # template name for galaxy
          npca_gal=5, # number of galaxy templates
          
          # continuum model fit parameters
          Fe_uv_op=True,  # If True, fit continuum with UV and optical FeII template
          poly=True,  # If True, fit continuum with the polynomial component to account for the dust reddening
          BC=False,  # If True, fit continuum with Balmer continua from 1000 to 3646A
          initial_guess=None,  # Initial parameters for continuum model, read the annotation of this function for detail
          rej_abs_conti=False,  # If True, it will iterately reject 3 sigma outlier absorption pixels in the continuum
          n_pix_min_conti=100,  # Minimum number of negative pixels for host continuuum fit to be rejected.

          # emission line fit parameters
          linefit=True,  # If True, the emission line will be fitted
          rej_abs_line=False,
          # If True, it will iterately reject 3 sigma outlier absorption pixels in the emission lines

          # fitting method selection
          MC=False,
          # If True, do Monte Carlo resampling of the spectrum based on the input error array to produce the MC error array
          MCMC=False,
          # If True, do Markov Chain Monte Carlo sampling of the posterior probability densities to produce the error array
          nsamp=200,
          # The number of trials of the MC process (if MC=True) or number samples to run MCMC chain (if MCMC=True)

          # advanced fitting parameters
          param_file_name='qsopar.fits',  # Name of the qso fitting parameter FITS file.
          nburn=20,  # The number of burn-in samples to run MCMC chain
          nthin=10,  # To set the MCMC chain returns every n samples
          epsilon_jitter=0.,
          # Initial jitter for every initial guass to avoid local minimum. (Under test, not recommanded to change)

          # customize the results
          save_result=False,  # If True, all the fitting results will be saved to a fits file
          save_fits_name=None,  # The output name of the result fits
          save_fits_path=path3,  # The output path of the result fits
          plot_fig=True,  # If True, the fitting results will be plotted
          save_fig=False,  # If True, the figure will be saved
          plot_corner=True,  # Whether or not to plot the corner plot results if MCMC=True

          # debugging mode
          verbose=True,  # turn on (True) or off (False) debugging output

          # sublevel parameters for figure plot and emcee
          kwargs_plot={
              'save_fig_path': '.',  # The output path of the figure
              'broad_fwhm'   : 1200  # km/s, lower limit that code decide if a line component belongs to broad component
          },
          kwargs_conti_emcee={},
          kwargs_line_emcee={})

    end = timeit.default_timer()

    print(f'Fitting finished in {np.round(end - start, 1)}s')

    
    
    # q_mc = QSOFit(lam, flux, err, z, ra=ra, dec=dec, plateid=plateid, mjd=mjd, fiberid=fiberid, path=path1)

    # start = timeit.default_timer()
    # # Do the fitting
    
    # q_mc.Fit(name=None, nsmooth=1, deredden=True, reject_badpix=False, wave_range=None, \
    #          wave_mask=None, decompose_host=True, BC03=False, npca_gal=5, npca_qso=10, \
    #              Fe_uv_op=True, poly=True, rej_abs_conti=False, rej_abs_line=True, MC=True, nsamp=200, linefit=True, \
    #                  save_result=True, kwargs_plot={'save_fig_path': '.'}, save_fits_name=None, verbose=True)

    # end = timeit.default_timer()

    # print(f'Fitting finished in {np.round(end - start, 1)}s')
    
    # -----------------------------------------------------------------------------
    
    
    # Obtain fit result files
    # Path of Line Properties Results for each system 
    
    path5 = '/Users/emilytemple/documents/rsmbh-agn-fit2.0/pyqsofit/Fit Results/Line Properties/'+sourcename+'/'
    
    # this line does not play nice, EMT edit 1/8/24
    # data=fits.open(path3+sourcename+'.fits') 
    
    # -----------------------------------------------------------------------------
    
    
    # Printing separate plots for each line
    # NOTE: If you want to include or exclude plots of other line/line complexes, 
    # adjust the following code
    
    # Plotting broad H_alpha, NII, and SII line complex
    fig1 = plt.figure(figsize=(16,12))
    for p in range(len(q_mle.gauss_result)//3):
        if q_mle.CalFWHM(q_mle.gauss_result[3*p+2]) < 1000:  # < 1000 km/s narrow
            color = 'g' # narrow
        else:
            color = 'r' # broad
        plt.plot(q_mle.wave, q_mle.Onegauss(np.log(q_mle.wave), q_mle.gauss_result
                                            [p*3:(p+1)*3]), color=color)  
    # Plot total line model
    plt.plot(q_mle.wave, q_mle.Manygauss(np.log(q_mle.wave), q_mle.gauss_result), 
             'b', lw=2)
    plt.plot(q_mle.wave, q_mle.line_flux,'k')
    plt.xlim(6400, 7000)
    plt.xlabel(r'$\rm Rest \, Wavelength$ ($\rm \AA$)', fontsize = 20)
    plt.ylabel(r'$\rm f_{\lambda}$ ($\rm 10^{-17} erg\;s^{-1}\;cm^{-2}\;\AA^{-1}$)'
               , fontsize = 20)
    plt.title(r'Broad $H_{\alpha}+[NII]+[SII]$', fontsize = 30)
    plt.savefig(path5+sourcename+'_BroadHa_LineComplex.pdf')
    
    # Plotting broad H_beta and [OIII] line complex
    fig2 = plt.figure(figsize=(16,12))
    for p in range(len(q_mle.gauss_result)//3):
        if q_mle.CalFWHM(q_mle.gauss_result[3*p+2]) < 1000:  # < 1000 km/s narrow
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
    plt.title(r'Broad $H_{\beta}+[OIII]$', fontsize = 30)
    plt.savefig(path5+sourcename+'_BroadHb_LineComplex.pdf')
     
    
       
    # -----------------------------------------------------------------------------    
    
    
    # Line fitting results
    
    # Use the following two commands to see what lines are
    # in the results file
    #print(q_mle.line_result_name)
    #print(q_mle.line_result)
    
    # NOTE: If you want to include or exclude properties of other lines adjust this 
    # section (and the text file section) of the code
    
    # H_alpha
    #fwhmha, sigmaha, ewha, peakha, areaha = q_mle.line_prop_from_name(
    #    'Ha_br', 'broad')
    
    # H_beta 
    #fwhmhb, sigmahb, ewhb, peakhb, areahb = q_mle.line_prop_from_name(
    #    'Hb_br', 'broad')
    
    # Narrow [OIII]5007:
    #fwhmo5, sigmao5, ewo5, peako5, areao5 = q_mle.line_prop_from_name(
    #    'OIII5007c', 'narrow')
    
    # Narrow [OIII]4959:
    #fwhmo4, sigmao4, ewo4, peako4, areao4 = q_mle.line_prop_from_name(
    #    'OIII4959c', 'narrow')
    
    # H_gamma
    #fwhmhg, sigmahg, ewhg, peakhg, areahg = q_mle.line_prop_from_name(
    #    'Hg_br', 'broad')
    
    
    # -----------------------------------------------------------------------------
    
    
    # Saving line properties into separate text file
    # NOTE: If lines were included or excluded above, do the same here too
    """
    with open(path5+sourcename+"_LineProperties.txt","w") as f:
        print('PyQSOFit Calculated Line Properties', file=f)
        print('', file=f)
        # Broad Ha - Component 1
        print('Broad Ha:', file=f)
        print('   '  + "FWHM = " '%.2f' % fwhmha, "km s\u207B\u00B9", file=f)
        print('   '  + "\u03C3 = " '%.2f' % sigmaha, "km s\u207B \u00B9", file=f)
        print('   '  + "EW = ", '%.2f' % ewha, "\u212B", file=f)
        print('   '  + "Peak = ", '%.2f' % peakha, "\u212B", file=f)
        print('   '  + "Area = ", '%.2f' % areaha, 
              "x 10\u207B\u00B9\u2077 erg s\u207B\u00B9 cm\u207B\u00B2", file=f)
        # Broad Hb
        print('Broad Hb:', file=f)
        print('   '  + "FWHM = " '%.2f' % fwhmhb, "km s\u207B\u00B9", file=f)
        print('   '  + "\u03C3 = " '%.2f' % sigmahb, "km s\u207B \u00B9", file=f)
        print('   '  + "EW = ", '%.2f' % ewhb, "\u212B", file=f)
        print('   '  + "Peak = ", '%.2f' % peakhb, "\u212B", file=f)
        print('   '  + "Area = ", '%.2f' % areahb, 
              "x 10\u207B\u00B9\u2077 erg s\u207B\u00B9 cm\u207B\u00B2", file=f)
        # Narrow [OIII]5007
    #    print('Narrow [OIII]5007:', file=f)
    #    print('   '  + "FWHM = " '%.2f' % fwhmo5, "km s\u207B\u00B9", file=f)
    #    print('   '  + "\u03C3 = " '%.2f' % sigmao5, "km s\u207B \u00B9", file=f)
    #    print('   '  + "EW = ", '%.2f' % ewo5, "\u212B", file=f)
    #    print('   '  + "Peak = ", '%.2f' % peako5, "\u212B", file=f)
    #    print('   '  + "Area = ", '%.2f' % areao5, 
    #          "x 10\u207B\u00B9\u2077 erg s\u207B\u00B9 cm\u207B\u00B2", file=f)
    #    # Narrow [OIII]4959
    #    print('Narrow [OIII]4959:', file=f)
    #    print('   '  + "FWHM = " '%.2f' % fwhmo4, "km s\u207B\u00B9", file=f)
    #    print('   '  + "\u03C3 = " '%.2f' % sigmao4, "km s\u207B \u00B9", file=f)
    #    print('   '  + "EW = ", '%.2f' % ewo4, "\u212B", file=f)
    #    print('   '  + "Peak = ", '%.2f' % peako4, "\u212B", file=f)
    #    print('   '  + "Area = ", '%.2f' % areao4, 
    #          "x 10\u207B\u00B9\u2077 erg s\u207B\u00B9 cm\u207B\u00B2", file=f)
         # Broad Hb
        print('Broad Hg:', file=f)
        print('   '  + "FWHM = " '%.2f' % fwhmhg, "km s\u207B\u00B9", file=f)
        print('   '  + "\u03C3 = " '%.2f' % sigmahg, "km s\u207B \u00B9", file=f)
        print('   '  + "EW = ", '%.2f' % ewhg, "\u212B", file=f)
        print('   '  + "Peak = ", '%.2f' % peakhg, "\u212B", file=f)
        print('   '  + "Area = ", '%.2f' % areahg, 
          "x 10\u207B\u00B9\u2077 erg s\u207B\u00B9 cm\u207B\u00B2", file=f)"""
    
    
    # -----------------------------------------------------------------------------
    
    
    # Extracting models for the whole spectrum
    # Plotting models
    fig4 = plt.figure(figsize=(15,5))
    # Plot the quasar rest frame spectrum after removed the host galaxy component
    plt.plot(q_mle.wave, q_mle.flux, 'grey',label='Data')
    plt.plot(q_mle.wave, q_mle.err, 'r',label='Error')
    
    # Skip the error results before plotting
    if q_mle.MCMC == True:
        gauss_result = q_mle.gauss_result[::2]
    else:
        gauss_result = q_mle.gauss_result
    
    
    # To plot the whole model, we use Manygauss to show the line fitting results 
    # saved in gauss_result  
    plt.plot(q_mle.wave, q_mle.Manygauss(np.log(q_mle.wave), gauss_result) + 
             q_mle.f_conti_model, 'b', label='Line', lw=2)
    plt.plot(q_mle.wave, q_mle.f_conti_model, 'c', lw=2,label='Continuum+FeII')
    plt.plot(q_mle.wave, q_mle.PL_poly_BC, 'orange', lw=2,label='Continuum')
    plt.plot(q_mle.wave, q_mle.host, 'm', lw=2,label='Host')
    plt.legend()
    plt.xlim(3500, 8000)
    plt.xlabel(r'$\rm Rest \, Wavelength$ ($\rm \AA$)', fontsize=20)
    plt.ylabel(r'$\rm f_{\lambda}$ ($\rm 10^{-17} erg\;s^{-1}\;cm^{-2}\;\AA^{-1}$)', fontsize=20)
    plt.savefig(path5+sourcename+'_SpectrumModels.pdf')
    
    
    #print('optical Fe flux (10^(-17) erg/s/cm^2): ' + 
    #q_mcmc.conti_result[q_mcmc.conti_result_name=='Fe_flux_4435_4685'][0])
    Fe_flux_result, Fe_flux_type, Fe_flux_name = q_mle.Get_Fe_flux(np.array([4400,4900]))
    print('Fe flux within a specific range: \n'+Fe_flux_name[0]+'= '+str(Fe_flux_result[0]))
    
    return

#----------------------------------------------------------------------------------
# Loop over the values using the function
# you can CHOOSE to loop or to specify a spectrum 
loop = False

if loop:
    for source in Data_list:
        runscript_func(source)
else:
    runscript_func('0332-52367-0639')




