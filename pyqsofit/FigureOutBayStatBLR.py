#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 15:35:03 2024

@author: emilytemple
"""


import re
import os,timeit
import numpy as np
from PyQSOFit import QSOFit
from astropy.io import fits
import matplotlib.pyplot as plt
import warnings
from ppxf.ppxf import ppxf
import sys
sys.path.append('../')
from astropy.table import Table
warnings.filterwarnings("ignore")

# Path to where this run code file sits
path_ex='/Users/emilytemple/documents/rsmbh-agn-fit2.0/pyqsofit'
#EMT edited the range values to better constrain them
newdata = np.rec.array([
(6564.61, r'H$\alpha$', 6350, 6800, 'Ha_br',     3,   0.1, 0.0, 1e10,   5e-3, 0.004,  0.05,     0.015, 0, 0, 0, 0.05,    1),
(6564.61, r'H$\alpha$', 6500, 6600, 'Ha_na',     1,   0.1, 0.0, 1e10,   1e-3, 5e-4,   0.001,   0.01,  1, 1, 0, 0.002,   1),
(6549.85, r'H$\alpha$', 6500, 6600, 'NII6549',   1,   0.1, 0.0, 1e10,   1e-3, 2.3e-4, 0.001,   5e-3,  1, 2, 1, 0.001,   1),
(6585.28, r'H$\alpha$', 6500, 6630, 'NII6585',   1,   0.1, 0.0, 1e10,   1e-3, 2.3e-4, 0.001,   5e-3,  1, 2, 1, 0.003,   1),
(6718.29, r'H$\alpha$', 6690, 6780, 'SII6718',   1,   0.1, 0.0, 1e10,   1e-3, 2.3e-4, 0.001,   5e-3,  1, 1, 0, 0.001,   1),
(6732.67, r'H$\alpha$', 6700, 6790, 'SII6732',   1,   0.1, 0.0, 1e10,   1e-3, 2.3e-4, 0.001,   5e-3,  1, 1, 0, 0.001,   1),

# Liza Matrecito - Added the following lines - September 22, 2022:
#(6302.05, 'OI', 6270, 6340, 'OI6302',           1,    0.1, 0.0, 1e10,   1e-3, 3.3e-4, 0.0017,   0.01, 1, 1, 0, 0.001,    1),
(5877.29, 'HeI', 5800, 5920, 'HeI5877_na',      1,    0.1, 0.0, 1e10,   1e-3, 3.3e-4, 0.001,   0.01, 0, 0, 0, 0.001,    1), # Added on November 29, 2022 
(5877.29, 'HeI', 5800, 5920, 'HeI5877_br',      1,    0.1, 0.0, 1e10,   5e-3, 0.004, 0.05,      0.015, 0, 0, 0, 0.001,   1), # Added on November 30, 2022 
# end of Liza Matrecito edits

(4862.68, r'H$\beta$', 4650, 5100, 'Hb_br',       3,   0.1, 0.0, 1e10,   5e-3, 0.004,  0.05,     0.01, 0, 0, 0, 0.01,    1),
(4862.68, r'H$\beta$', 4800, 4910, 'Hb_na',       1,   0.1, 0.0, 1e10,   1e-3, 2.3e-4, 0.001,   0.01, 1, 1, 0, 0.002,   1),
(4960.30, r'H$\beta$', 4920, 5100, 'OIII4959c',   1,   0.1, 0.0, 1e10,   1e-3, 2.3e-4, 0.001,   0.01, 1, 1, 1, 0.33,    1),
(5008.24, r'H$\beta$', 4960, 5100, 'OIII5007c',   1,   0.1, 0.0, 1e10,   1e-3, 2.3e-4, 0.001,   0.01, 1, 1, 1, 1,       1),
(4960.30, r'H$\beta$', 4910, 5000, 'OIII4959w',   1,   0.1, 0.0, 1e10,   3e-3, 2.3e-4, 0.003,    0.01, 2, 2, 2, 0.33,    1),
(5008.24, r'H$\beta$', 4960, 5080, 'OIII5007w',   1,   0.1, 0.0, 1e10,   3e-3, 2.3e-4, 0.003,    0.01, 2, 2, 2, 1,       1),
#(4687.02, r'H$\beta$', 4640, 5100, 'HeII4687_br', 1, 0.1, 0.0, 1e10, 5e-3, 0.004,  0.05,   0.005, 0, 0, 0, 0.001, 1),
#(4687.02, r'H$\beta$', 4640, 4700, 'HeII4687_na', 1,   0.1, 0.0, 1e10, 1e-3, 2.3e-4, 0.0017, 0.005, 1, 1, 0, 0.001,      1),

# Liza Matrecito - Added the following lines - September 22, 2022:
(4341.68, r'H$\gamma$', 4200, 4500, 'Hg_br',       2,   0.1, 0.0, 1e10,   5e-3, 0.004,  0.05,     0.01, 0, 0, 0, 0.01,    1),
(4341.68, r'H$\gamma$', 4200, 4500, 'Hg_na',       1,   0.1, 0.0, 1e10,   1e-3, 2.3e-4, 0.0017,   0.01, 1, 1, 0, 0.002,   1),
#(4364.44, r'H$\gamma$', 4200, 4500, 'OIII4364c',   1,   0.1, 0.0, 1e10,   1e-3, 2.3e-4, 0.0017,   0.01, 1, 1, 0, 0.002,   1),
#(4364.44, r'H$\gamma$', 4200, 4500, 'OIII4364w',  1,   0.1, 0.0, 1e10,   3e-3, 2.3e-4, 0.004,    0.01, 2, 2, 0, 0.001,   1),

(4102.98, r'H$\delta$', 4000, 4300, 'Hd_br',   1,   0.1, 0.0, 1e10,   5e-3, 0.004,  0.05,     0.01, 0, 0, 0, 0.01,    1),
(4102.98, r'H$\delta$', 4000, 4300, 'Hd_na',   1,   0.1, 0.0, 1e10,   1e-3, 2.3e-4, 0.0011,   0.01, 1, 1, 0, 0.002,   1),
# end of Liza Matrecito edits

#(3934.78, 'CaII', 3900, 3960, 'CaII3934', 2, 0.1, 0.0, 1e10, 1e-3, 3.333e-4, 0.0017, 0.01, 99, 0, 0, -0.001, 1),

# Liza Matrecito - Added the following lines - September 22, 2022:
(3869.85, 'NeIII', 3700, 3980, 'NeIII3869', 1, 0.1, 0.0, 1e10, 1e-3, 3.333e-4, 0.0017, 0.01, 0, 0, 0, 0.001, 1),
# end of Liza Matrecito edits

(3728.48, 'OII', 3650, 3800, 'OII3728', 1, 0.1, 0.0, 1e10, 1e-3, 3.333e-4, 0.0017, 0.01, 1, 1, 0, 0.001, 1),
    
#(3426.84, 'NeV', 3380, 3480, 'NeV3426',    1, 0.1, 0.0, 1e10, 1e-3, 3.333e-4, 0.0017, 0.01, 0, 0, 0, 0.001, 1),
#(3426.84, 'NeV', 3380, 3480, 'NeV3426_br', 1, 0.1, 0.0, 1e10, 5e-3, 0.0025,   0.02,   0.01, 0, 0, 0, 0.001, 1),

#(2798.75, 'MgII', 2700, 2900, 'MgII_br', 1, 0.1, 0.0, 1e10, 5e-3, 0.004, 0.05,   0.0017, 0, 0, 0, 0.05, 1),
#(2798.75, 'MgII', 2700, 2900, 'MgII_na', 1, 0.1, 0.0, 1e10, 1e-3, 5e-4,  0.0017, 0.01,   1, 1, 0, 0.002, 1),

#(1908.73, 'CIII', 1700, 1970, 'CIII_br',     2,   0.1, 0.0, 1e10,   5e-3, 0.004, 0.05,     0.015, 99, 0, 0, 0.01,    1),
#(1908.73, 'CIII', 1700, 1970, 'CIII_na',     1,   0.1, 0.0, 1e10,   1e-3, 5e-4,  0.0017,   0.01,  1, 1, 0, 0.002,    1),
#(1892.03, 'CIII', 1700, 1970, 'SiIII1892',   1,   0.1, 0.0, 1e10,   2e-3, 0.001, 0.015,    0.003, 1, 1, 0, 0.005,    1),
#(1857.40, 'CIII', 1700, 1970, 'AlIII1857',   1,   0.1, 0.0, 1e10,   2e-3, 0.001, 0.015,    0.003, 1, 1, 0, 0.005,    1),
#(1816.98, 'CIII', 1700, 1970, 'SiII1816',    1,   0.1, 0.0, 1e10,   2e-3, 0.001, 0.015,    0.01,  1, 1, 0, 0.0002,   1),
#(1786.7,  'CIII', 1700, 1970, 'FeII1787',    1,   0.1, 0.0, 1e10,   2e-3, 0.001, 0.015,    0.01,  1, 1, 0, 0.0002,   1),
#(1750.26, 'CIII', 1700, 1970, 'NIII1750',    1,   0.1, 0.0, 1e10,   2e-3, 0.001, 0.015,    0.01,  1, 1, 0, 0.001,    1),
#(1718.55, 'CIII', 1700, 1900, 'NIV1718',     1,   0.1, 0.0, 1e10,   2e-3, 0.001, 0.015,    0.01,  1, 1, 0, 0.001,    1),

#(1549.06, 'CIV', 1500, 1700, 'CIV_br',        1,   0.1, 0.0, 1e10,   5e-3, 0.004, 0.05,      0.015, 0, 0, 0, 0.05 ,   1),
#(1549.06, 'CIV', 1500, 1700, 'CIV_na',        1,   0.1, 0.0, 1e10,   1e-3, 5e-4,  0.0017,    0.01,  1, 1, 0, 0.002,   1),
#(1640.42, 'CIV', 1500, 1700, 'HeII1640',      1,   0.1, 0.0, 1e10,   1e-3, 5e-4, 0.0017,     0.008, 1, 1, 0, 0.002,   1),
#(1663.48, 'CIV', 1500, 1700, 'OIII1663',      1,   0.1, 0.0, 1e10,   1e-3, 5e-4,   0.0017,   0.008, 1, 1, 0, 0.002,   1),
#(1640.42, 'CIV', 1500, 1700, 'HeII1640_br',   1,   0.1, 0.0, 1e10,   5e-3, 0.0025, 0.02,     0.008, 1, 1, 0, 0.002,   1),
#(1663.48, 'CIV', 1500, 1700, 'OIII1663_br',   1,   0.1, 0.0, 1e10,   5e-3, 0.0025, 0.02,     0.008, 1, 1, 0, 0.002,   1),

#(1402.06, 'SiIV', 1290, 1450, 'SiIV_OIV1',   1,   0.1, 0.0, 1e10,   5e-3, 0.002, 0.05,    0.015, 1, 1, 0, 0.05,    1),
#(1396.76, 'SiIV', 1290, 1450, 'SiIV_OIV2',   1,   0.1, 0.0, 1e10,   5e-3, 0.002, 0.05,    0.015, 1, 1, 0, 0.05,    1),
#(1335.30, 'SiIV', 1290, 1450, 'CII1335',     1,   0.1, 0.0, 1e10,   2e-3, 0.001, 0.015,   0.01,  1, 1, 0, 0.001,   1),
#(1304.35, 'SiIV', 1290, 1450, 'OI1304',      1,   0.1, 0.0, 1e10,   2e-3, 0.001, 0.015,   0.01,  1, 1, 0, 0.001,   1),

(1215.67, 'Lya', 1150, 1290, 'Lya_br', 1, 0.1, 0.0, 1e10, 5e-3, 0.004, 0.05,   0.02, 0, 0, 0, 0.05 , 1),
(1215.67, 'Lya', 1150, 1290, 'Lya_na', 1, 0.1, 0.0, 1e10, 1e-3, 5e-4,  0.0017, 0.01, 0, 0, 0, 0.002, 1)],

formats = 'float32,      a20,  float32, float32,      a20,  int32, float32, float32, float32, float32, float32, float32, float32,   int32,  int32,  int32,   float32, int32',
names  =  ' lambda, compname,   minwav,  maxwav, linename, ngauss,  inisca,  minsca,  maxsca,  inisig,  minsig,  maxsig,  voff,     vindex, windex,  findex,  fvalue,  vary')

# Header
hdr = fits.Header()
hdr['lambda'] = 'Vacuum Wavelength in Ang'
hdr['minwav'] = 'Lower complex fitting wavelength range'
hdr['maxwav'] = 'Upper complex fitting wavelength range'
hdr['ngauss'] = 'Number of Gaussians for the line'

# Can be set to negative for absorption lines if you want
hdr['inisca'] = 'Initial guess of line scale [in ??]'
hdr['minsca'] = 'Lower range of line scale [??]'
hdr['maxsca'] = 'Upper range of line scale [??]'

hdr['inisig'] = 'Initial guess of linesigma [in lnlambda]'
hdr['minsig'] = 'Lower range of line sigma [lnlambda]'  
hdr['maxsig'] = 'Upper range of line sigma [lnlambda]'

hdr['voff  '] = 'Limits on velocity offset from the central wavelength [lnlambda]'
hdr['vindex'] = 'Entries w/ same NONZERO vindex constrained to have same velocity'
hdr['windex'] = 'Entries w/ same NONZERO windex constrained to have same width'
hdr['findex'] = 'Entries w/ same NONZERO findex have constrained flux ratios'
hdr['fvalue'] = 'Relative scale factor for entries w/ same findex'

hdr['vary'] = 'Whether or not to vary the line parameters (set to 0 to fix the line parameters to initial values)'

# Save line info
hdu = fits.BinTableHDU(data=newdata, header=hdr, name='data')
hdu.writeto(os.path.join(path_ex, 'qsopar.fits'), overwrite=True)

# Print table
Table(newdata)

# create a header
hdr0 = fits.Header()
hdr0['Author'] = 'ET'
primary_hdu = fits.PrimaryHDU(header=hdr0)

conti_windows = np.rec.array([
    (1150., 1170.), 
    (1275., 1290.),
    (1350., 1360.),
    (1445., 1465.),
    (1690., 1705.),
    (1770., 1810.),
    (1970., 2400.),
    (2480., 2675.),
    (2925., 3400.),
    (3775., 3832.),
    (4000., 4050.),
    (4200., 4230.),
    (4435., 4640.),
    (5100., 5535.),
    (6005., 6035.),
    (6110., 6250.),
    (6800., 7000.),
    (7160., 7180.),
    (7500., 7800.),
    (8050., 8150.), # Continuum fitting windows (to avoid emission line, etc.)  [AA]
    ], 
    formats = 'float32,  float32',
    names =    'min,     max')

hdu2 = fits.BinTableHDU(data=conti_windows, name='conti_windows')

conti_priors = np.rec.array([
    ('Fe_uv_norm',  0.0,   0.0,   1e10,  1), # Normalization of the MgII Fe template [flux]
    ('Fe_uv_FWHM',  3000,  1200,  18000, 1), # FWHM of the MgII Fe template [AA]
    ('Fe_uv_shift', 0.0,   -0.01, 0.01,  1), # Wavelength shift of the MgII Fe template [lnlambda]
    ('Fe_op_norm',  0.0,   0.0,   1e10,  1), # Normalization of the Hbeta/Halpha Fe template [flux]
    ('Fe_op_FWHM',  3000,  1200,  18000, 1), # FWHM of the Hbeta/Halpha Fe template [AA]
    ('Fe_op_shift', 0.0,   -0.01, 0.01,  1), # Wavelength shift of the Hbeta/Halpha Fe template [lnlambda]
    ('PL_norm',     1.0,   0.0,   1e10,  1), # Normalization of the power-law (PL) continuum f_lambda = (lambda/3000)^-alpha
    ('PL_slope',    -1.5,  -5.0,  3.0,   1), # Slope of the power-law (PL) continuum
    ('Blamer_norm', 0.0,   0.0,   1e10,  1), # Normalization of the Balmer continuum at < 3646 AA [flux] (Dietrich et al. 2002)
    ('Balmer_Te',   15000, 10000, 50000, 1), # Te of the Balmer continuum at < 3646 AA [K?]
    ('Balmer_Tau',  0.5,   0.1,   2.0,   1), # Tau of the Balmer continuum at < 3646 AA
    ('conti_a_0',   0.0,   None,  None,  1), # 1st coefficient of the polynomial continuum
    ('conti_a_1',   0.0,   None,  None,  1), # 2nd coefficient of the polynomial continuum
    ('conti_a_2',   0.0,   None,  None,  1), # 3rd coefficient of the polynomial continuum
    # Note: The min/max bounds on the conti_a_0 coefficients are ignored by the code,
    # so they can be determined automatically for numerical stability.
    ],

    formats = 'a20,  float32, float32, float32, int32',
    names = 'parname, initial,   min,     max,     vary')

hdr3 = fits.Header()
hdr3['ini'] = 'Initial guess of line scale [flux]'
hdr3['min'] = 'FWHM of the MgII Fe template'
hdr3['max'] = 'Wavelength shift of the MgII Fe template'

hdr3['vary'] = 'Whether or not to vary the parameter (set to 0 to fix the continuum parameter to initial values)'


hdu3 = fits.BinTableHDU(data=conti_priors, header=hdr3, name='conti_priors')

measure_info = Table(
    [
        [[1350, 1450, 3000, 4200, 5100]],
        [[
            # [2240, 2650], 
            [4435, 4685],
        ]]
    ],
    names=([
        'cont_loc',
        'Fe_flux_range'
    ]),
    dtype=([
        'float32',
        'float32'
    ])
)
hdr4 = fits.Header()
hdr4['cont_loc'] = 'The wavelength of continuum luminosity in results'
hdr4['Fe_flux_range'] = 'Fe emission wavelength range calculated in results'

hdu4 = fits.BinTableHDU(data=measure_info, header=hdr4, name='measure_info')

hdu_list = fits.HDUList([primary_hdu, hdu, hdu2, hdu3, hdu4])
hdu_list.writeto(os.path.join(path_ex, 'qsopar.fits'), overwrite=True)

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
# Getting the data file needed by specifying the path and sourcename

Data_path = '/Users/emilytemple/documents/rsmbh-agn-fit2.0/pyqsofit/Data/Outliers'
#Data_list_names = os.listdir(Data_path)
sourcename='0332-52367-0639'

data = fits.open(os.path.join(path1+'Data/Outliers/spec-'+sourcename+'.fits'))
lam = 10**data[1].data['loglam']                           # OBS wavelength (A)
flux = data[1].data['flux']                           # OBS flux (erg/s/cm^2/A)
err = 1./np.sqrt(data[1].data['ivar'])                          # 1 sigma error
z = data[2].data['z'][0]                                             # Redshift

# Optional
ra = data[0].header['plug_ra']  # RA
dec = data[0].header['plug_dec']  # DEC
plateid = data[0].header['plateid']  # SDSS plate ID
mjd = data[0].header['mjd']  # SDSS MJD
fiberid = data[0].header['fiberid']  # SDSS fiber ID

#Specifying the max num of gaussians to try to fit BLR
ngauss_max = 5  # stop at 5 components
bic_last = np.inf #start with an "infinite" BIC

# Number of Gaussians to try loop
# typically 5 gaussians is the max number ever needed, prob a bit overkill
for ngauss in range(1, ngauss_max):

    print(fr'Fitting broad H$\alpha$ with {ngauss} components.')

    # Change the number of Gaussians for the Ha line in the line parameter file

  #  hdu_new = hdu.copy() #making new copy of param table ?? for what
    if 'Ha_br' in hdu.data['linename']:
        hdu.data['ngauss'][hdu.data['linename'] == 'Ha_br'] = ngauss
    hdu_list = fits.HDUList([primary_hdu, hdu, hdu2, hdu3, hdu4])
    hdu_list.writeto(os.path.join(path_ex, 'qsopar.fits'), overwrite=True)

    # Do the fitting
    q = QSOFit(lam, flux, err, z, ra=ra, dec=dec, plateid=plateid, mjd=mjd, fiberid=fiberid, path=path1)

    q.Fit(name=None, nsmooth=1, deredden=True, reject_badpix=False, wave_range=None, \
          wave_mask=None, decompose_host=True, host_prior=True, decomp_na_mask=True, npca_gal=5, npca_qso=10, qso_type='CZBIN1',\
          Fe_uv_op=True, poly=True, BC=False, rej_abs_conti=False, rej_abs_line=False, \
          initial_guess=None, MCMC=True, nburn=20, nsamp=200, nthin=10, linefit=True, \
          save_result=False, plot_fig=False, verbose=False)

    mask_Ha_bic = q.line_result_name == '2_line_min_chi2'

    bic = float(q.line_result[mask_Ha_bic][0])
    print(f'BIC =', bic)
    print(f'BIC Last=', bic_last)

    print(f'Delta BIC = {np.round(bic_last - bic, 1)}')

    # Stop condition of Delta BIC = 10 is a good rule of thumb
    if bic_last - bic < 10:
        print(bic_last)
        print(f'{ngauss-1} components is prefered')
        break

    bic_last = bic

    # Plot the result
    if q.MCMC:
        gauss_result = q.gauss_result[::2]
    else:
        gauss_result = q.gauss_result

    fig, ax = plt.subplots(1, 1, figsize=(8, 5))

    # Plot individual line components
    for p in range(len(gauss_result) // 3):
        if q.CalFWHM(gauss_result[3 * p + 2]) < 1200:  # < 1200 km/s narrow
            color = 'g'  # narrow
        else:
            color = 'r'  # broad
        ax.plot(q.wave, q.Onegauss(np.log(q.wave), gauss_result[p * 3:(p + 1) * 3]), color=color)

    # Plot total line model
    ax.plot(q.wave, q.Manygauss(np.log(q.wave), gauss_result), 'b', lw=2)
    ax.plot(q.wave, q.line_flux, 'k')
    ax.set_xlim(6300, 6900)
    ax.set_xlabel(r'$\rm Rest \, Wavelength$ ($\rm \AA$)', fontsize=20)
    ax.set_ylabel(r'$\rm f_{\lambda}$ ($\rm 10^{-17} erg\;s^{-1}\;cm^{-2}\;\AA^{-1}$)', fontsize=20)

    c = 1  # Ha
    ax.text(0.02, 0.80, r'$\chi ^2_\nu=$' + str(np.round(float(q.comp_result[c * 7 + 4]), 2)),
            fontsize=16, transform=ax.transAxes)

    plt.show()