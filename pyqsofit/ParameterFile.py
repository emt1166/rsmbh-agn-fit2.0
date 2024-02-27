#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 14:41:15 2023

@author: emilytemple

This file contains the parameters for each of our lines.
The parameters here have been edited to reflect the best fitting results from 
trial and error with many different spectra. 
They can be edited further to fit particularily difficult/weird spectra. 

In order, the columns for the file are:
lambda, label, lower bound, upper boun,  n_gauss,
initial flux scale guess, lower flux, upper flux,
initial guess of line sigma (ln lambda), lower range, upper range,
v offset from centroid, vindex(tying), windex(tying), findex(tying), fvalue, var.
"""
import numpy as np
import sys, os
sys.path.append('../')
from astropy.io import fits
import warnings
from astropy.table import Table
warnings.filterwarnings("ignore")

# Path to where this run code file sits
path_ex='/Users/emilytemple/documents/rsmbh-agn-fit2.0/pyqsofit'
#EMT edited the range values to better constrain them
newdata = np.rec.array([
(6564.61, r'H$\alpha$', 6350, 6800, 'Ha_br',     3,   0.1, 0.0, 1e10,   5e-3, 0.004,  0.05,     0.015, 0, 0, 0, 0.05,    1),
(6564.61, r'H$\alpha$', 6500, 6600, 'Ha_na',     1,   0.1, 0.0, 1e10,   1e-3, 5e-4,   0.001414,   0.01,  1, 1, 0, 0.002,   1),
(6549.85, r'H$\alpha$', 6500, 6600, 'NII6549',   1,   0.1, 0.0, 1e10,   1e-3, 2.3e-4, 0.001414,   5e-3,  1, 2, 1, 0.001,   1),
(6585.28, r'H$\alpha$', 6500, 6630, 'NII6585',   1,   0.1, 0.0, 1e10,   1e-3, 2.3e-4, 0.001414,   5e-3,  1, 2, 1, 0.003,   1),
(6718.29, r'H$\alpha$', 6690, 6780, 'SII6718',   1,   0.1, 0.0, 1e10,   1e-3, 2.3e-4, 0.001414,   5e-3,  1, 1, 0, 0.001,   1),
(6732.67, r'H$\alpha$', 6700, 6790, 'SII6732',   1,   0.1, 0.0, 1e10,   1e-3, 2.3e-4, 0.001414,   5e-3,  1, 1, 0, 0.001,   1),

# Liza Matrecito - Added the following lines - September 22, 2022:
#(6302.05, 'OI', 6270, 6340, 'OI6302',           1,    0.1, 0.0, 1e10,   1e-3, 3.3e-4, 0.0017,   0.01, 1, 1, 0, 0.001,    1),
(5877.29, 'HeI', 5800, 5920, 'HeI5877_na',      1,    0.1, 0.0, 1e10,   1e-3, 3.3e-4, 0.001,   0.01, 0, 0, 0, 0.001,    1), # Added on November 29, 2022 
(5877.29, 'HeI', 5800, 5920, 'HeI5877_br',      1,    0.1, 0.0, 1e10,   5e-3, 0.004, 0.05,      0.015, 0, 0, 0, 0.001,   1), # Added on November 30, 2022 
# end of Liza Matrecito edits

(4862.68, r'H$\beta$', 4650, 5100, 'Hb_br',       3,   0.1, 0.0, 1e10,   5e-3, 0.004,  0.05,     0.01, 0, 0, 0, 0.01,    1),
(4862.68, r'H$\beta$', 4800, 4910, 'Hb_na',       1,   0.1, 0.0, 1e10,   1e-3, 2.3e-4, 0.001,   0.01, 1, 1, 0, 0.002,   1),
(4960.30, r'H$\beta$', 4920, 5100, 'OIII4959c',   1,   0.1, 0.0, 1e10,   0.00002, 0.00001, 0.001414,   0.01, 1, 1, 1, 0.33,    1),
(5008.24, r'H$\beta$', 4960, 5100, 'OIII5007c',   1,   0.1, 0.0, 1e10,   0.00002, 0.00001, 0.001414,   0.01, 1, 1, 1, 1,       1),
(4960.30, r'H$\beta$', 4910, 5000, 'OIII4959w',   1,   0.1, 0.0, 1e10,   0.00002, 0.00001, 0.001414,    0.01, 2, 2, 2, 0.33,    1),
(5008.24, r'H$\beta$', 4960, 5080, 'OIII5007w',   1,   0.1, 0.0, 1e10,   3e-3, 2.3e-4, 0.0035,    0.01, 2, 2, 2, 1,       1),
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
hdr['inisca'] = 'Initial guess of line scale [in flux??]'
hdr['minsca'] = 'Lower range of line scale [flux??]'
hdr['maxsca'] = 'Upper range of line scale [flux??]'

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