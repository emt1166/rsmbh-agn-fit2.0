#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 20:29:02 2024

@author: emilytemple

This script will be used to create gaussians for the rsmbh and the core left behind.
This will be used to obtain some simple theoretical values for our plots. 
Namely, we will be looking at the kurtosis and the asymmetry. 

These gaussians are going to be simple, but need some sort of physical meaning. 
Need to keep track of the "velocity" of each component. 
This means the FWHM of the gaussians needs to make sense in accordance with the 
kick velocity and the broadening due to velocity dispersion. 
"""

import numpy as np
import matplotlib.pyplot as plt
from numpy import random 
import seaborn as sns
import math
from astropy import units as u
from astropy import constants as const

x = random.normal(loc=0, scale=3, size=(2,4))

#sns.displot(random.normal(size=100), kde=True) #can input kde=False
#plt.show()

#okay that's cool, but let me just define a gaussian 
def gauss(mean, sigma, x):
    g = (1/(math.sqrt(2*math.pi*sigma**2))) * math.e**((-1*(x - mean)**2)/(2* sigma**2))
    return g

#Alright here is our core 
v_dis = np.linspace(-20000,20000, 25000)
core = 700*gauss(0,1500,0.8*v_dis)
plt.plot(v_dis,core, label='core')

#test case recoil
rsmbh = 470*gauss(-2800,2500,0.8*v_dis)
plt.plot(v_dis,rsmbh, label='recoil')

#showing the entire profile (add them)
tot = rsmbh+core
plt.plot(v_dis,tot, label='combined')
plt.ylabel('Arbitrary Flux Units')
plt.xlabel('Velocity (km/s)')
plt.legend()


'''
After the models are made, we need to measure the properties of them. 
These include the velocity shifts, asymmetries, and kurtosis. The FWHM
can also probably be done. 
'''

# ------------------------------------------------------------------------------
'''
SWITCHING GEARS
Going to try to make the modeling a bit more physical, so we need to use real
astrophytsics concepts.
I need to compute the relative flux ratios based on the SMBH mass and the BLR
radius --> use H recombo theory to get emission. 
'''
G = const.G
c = const.c

def schwarz(M):
    R_s = (2*G*M)/c**2
    return R_s

#print(schwarz((1e6*u.solMass).si)) #testing
#these are the different radii for each BH mass
SMBH1 = schwarz((1e6*u.solMass).si)
SMBH2 = schwarz((1e7*u.solMass).si)
SMBH3 = schwarz((1e8*u.solMass).si)
SMBH4 = schwarz((1e9*u.solMass).si)





