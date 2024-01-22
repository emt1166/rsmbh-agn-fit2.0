#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  1 20:38:59 2023

@author: emtem
"""

"""
This code allows you to pick out specific values from a .txt file
All you need to do is specify the path and which entry you want to read

I created this to be used with pyQSOfit outputs.
"""

import numpy as np
import os
import csv

#set path to specific object
#NOTE: these change depending on your set-up
Object = '1305-52757-0281'
Object_path = '/Users/emilytemple/documents/rsmbh-agn-fit2.0/pyqsofit/Fit Results/Line Properties/Outliers'

#set path to the output .txt file
#NOTE this changes based on your set-up
FWHM_path = Object_path+'/'+Object+'/'+Object+'_LineProperties.txt'
VShift_path = Object_path+'/'+Object+'/Line Profile Plots/'+Object+'_BLV-Shifts.txt'

#Now we want to access a specific value in each file
#This will give us the specific value we want 

#For FWHM text file-----------------------------------------------------------
Object_dict = {}

with open(FWHM_path, 'r') as fp:
    current_h_line = ''
    for i, line in enumerate(fp):
        if i < 2:
            continue
        line_text = line.strip()
        if line_text.endswith(':'):
            h_line = line_text.strip(':')
            Object_dict[h_line] = {}
            current_h_line = h_line
        else:
            line_values = line_text.split()
            agn_prop = line_values[0]
            agn_value = line_values[2]
            Object_dict[current_h_line][agn_prop] = agn_value
            
#print(Object_dict)
#print(Object_dict['Broad H_beta']['FWHM'])
#--------------------------------------------------------------------------

#For VShift text file ------------------------------------------------------
VObject_dict = {}

with open(VShift_path, 'r') as fp:
    for i, line in enumerate(fp):
        if i < 2:
            continue
        line_text = line.strip()
        line_values = line_text.split(' = ')
        if len(line_values) < 2:
            continue
        agn_prop = line_values[0]
        agn_value = line_values[1].split()[0]
        VObject_dict[agn_prop] = agn_value

#print(VObject_dict)
#print(VObject_dict['Peak H_alpha Velocity Shift'])


#----------------------------------------------------------------------------
"""
Now that we know this works for a single object, we want multiple objects
Let's put the ID's of the ones we want in a list format
Then, we can loop throughout that list
"""

#Defining a list of objects
#I think if I can automate this as much as possible, it'll be easier

Object_list = []
Object_names = os.listdir(Object_path)
#print(Object_names)


for object_name in Object_names:
    if object_name != '.DS_Store' and not object_name.endswith('.pdf'):
        Object_list.append(object_name)

#print(Object_list)

import csv

# Create a list to store the data
#these are for the two different .txt files we have right now
LineProfiles_output_data = []
LineShifts_output_data = []

for n in Object_list:
    FWHM_path_all = Object_path+'/'+n+'/'+n+'_LineProperties.txt'
    VShift_path_all = Object_path+'/'+n+'/Line Profile Plots/'+n+'_BLV-Shifts.txt'

    with open(FWHM_path_all, 'r') as fp:
        current_h_line = ''
        for i, line in enumerate(fp):
            if i < 2:
                continue
            line_text = line.strip()
            if line_text.endswith(':'):
                h_line = line_text.strip(':')
                Object_dict[h_line] = {}
                current_h_line = h_line
            else:
                line_values = line_text.split()
                agn_prop = line_values[0]
                agn_value = line_values[2]
                Object_dict[current_h_line][agn_prop] = agn_value


    with open(VShift_path_all, 'r') as fp:
        for i, line in enumerate(fp):
            if i < 2:
                continue
            line_text = line.strip()
            line_values = line_text.split(' = ')
            if len(line_values) < 2:
                continue
            agn_prop = line_values[0]
            agn_value = line_values[1].split()[0]
            VObject_dict[agn_prop] = agn_value
        
        fwhm_data = {'Object ID': n, 'FWHM BHA': Object_dict['Broad H_alpha']['FWHM'],
                     'FWHM BHB': Object_dict['Broad H_beta']['FWHM']}
        LineProfiles_output_data.append(fwhm_data)
        
        vshift_data = {'Object ID': n, 'VShift Peak BHA': VObject_dict['Peak H_alpha Velocity Shift'],
                       'Kurtosis BHA': VObject_dict['H_alpha Kurtosis'], 
                       'VShift Center BHA':VObject_dict['Center H_alpha Velocity Shift'],
                       'BHA Asymmetry': VObject_dict['H_alpha Asymmetry'],
                       'VShift Peak BHB': VObject_dict['Peak H_beta Velocity Shift'],
                       'VShift Center BHB':VObject_dict['Center H_beta Velocity Shift'],
                       'Kurtosis BHB': VObject_dict['H_beta Kurtosis'],
                       'BHB Asymmetry': VObject_dict['H_beta Asymmetry']}
        LineShifts_output_data.append(vshift_data)

# Specify the CSV file path
#there are two here because I wanted to separate some of the params for now
#Can change this later
csv_file_path1 = 'OutliersLineProfiles_output_data2.0.csv'
csv_file_path2 = 'OutliersLineShifts_output_data2.0.csv'

# Write the data to the CSV file
with open(csv_file_path1, 'w', newline='') as csv_file:
    fieldnames = ['Object ID', 'FWHM BHA','FWHM BHB']
    writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
    # Write the header
    writer.writeheader()
    # Write the data
    writer.writerows(LineProfiles_output_data)


with open(csv_file_path2, 'w', newline='') as csv_file:
    fieldnames = ['Object ID', 'VShift Peak BHA', 'Kurtosis BHA', 'VShift Center BHA', 'BHA Asymmetry',
                  'VShift Peak BHB', 'VShift Center BHB','Kurtosis BHB', 'BHB Asymmetry']
    writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
    # Write the header
    writer.writeheader()
    # Write the data
    writer.writerows(LineShifts_output_data)

#------------------------------------------------------------------------------
        
"""
The .csv files produced from this can then be used to analyze the results
from the two .txt output files produced by the code. 

You can compare the FWHM values, various velocity shift values,
asymmetries, and kurtosis values. 

One thing to note is to be aware of what files are in the directory.
This means if you have 100+ quasar folders/objects, be aware that this code
will go through and read all of the corresponding .txt files. 
If there are objects you do NOT wish you gather data from, you must specify that
on your own. 

"""
           



