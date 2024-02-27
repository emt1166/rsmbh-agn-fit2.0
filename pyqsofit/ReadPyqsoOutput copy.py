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
import re

#set path to the output files
#NOTE: these change depending on your set-up

# The path of the SNR New Flux Comparison Data
pathF = '/Users/emilytemple/documents/rsmbh-agn-fit2.0/pyqsofit/Data/Outliers New Flux/LOOP1/'
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

# -----------------------------------------------------------------------------------------------
def output_file(sourcename):
    Output_path = pathF+sourcename+'/' #path to the output files per sourcename
    
    #eventually will input another loop here to go over allll flux iterations using j instead of n
    j = 0 #this will change based on which iteration of flux we want to look at
    FWHM_path = Output_path+sourcename+'.'+f'{j}'+"_LineProperties.txt"
    VShift_path = Output_path+sourcename+'.'+f'{j}'+"_BLV-Shifts.txt"
    
    #----------------------------------------------------------------------------
    
   
    """
    Let's put the ID's of the ones we want in a list format
    Then, we can loop throughout that list
    """
    
    #Defining a list of the sourcenames in Object_list
    
    Object_list = []
    Object_names = os.listdir(pathF)
    
    
    for object_name in Object_names:
        if object_name != '.DS_Store' and not object_name.endswith('.pdf'):
            Object_list.append(object_name)
    
    # Create a list to store the data
    #these are for the two different .txt files we have right now
    LineProfiles_output_data = []
    LineShifts_output_data = []
    
    for source in Object_list: #for the sourcenames in the object list
        FWHM_path_all = pathF+source+'/'+source+'.'+f'{j}'+"_LineProperties.txt"
        VShift_path_all = pathF+source+'/'+source+'.'+f'{j}'+"_BLV-Shifts.txt"
        
        Object_dict = {}
        
        with open(FWHM_path_all, 'r') as fp:
            current_h_line = ''
            for i, line in enumerate(fp):
                line_text = line.strip()
                if i < 2 or line_text == '':
                    continue
                if line_text.endswith(':'):
                    h_line = line_text.strip(':')
                    Object_dict[h_line] = {}
                    current_h_line = h_line
                else:
                    line_values = line_text.split()
                    agn_prop = line_values[0]
                    agn_value = line_values[2]
                    Object_dict[current_h_line][agn_prop] = agn_value
        
        # Extract relevant properties for FWHM
        fwhm_data = {
            'Object ID': source,
            'FWHM BHA': Object_dict.get('Broad H_alpha', {}).get('FWHM', ''),
            'FWHM BHB': Object_dict.get('Broad H_beta', {}).get('FWHM', '')
        }
        
        # Append the extracted data to a list
        LineProfiles_output_data.append(fwhm_data)
            
        VObject_dict = {}
        
        # Open and read the file
        with open(VShift_path_all, 'r') as fp:
            # Read line by line
            for i, line in enumerate(fp):
                # Skip the first two lines
                if i < 2:
                    continue
                
                # Remove leading/trailing whitespace
                line_text = line.strip()
                
                # Split the line by ' = ' to extract property and value
                line_values = line_text.split(' = ')
                
                # Skip lines without property-value pairs
                if len(line_values) < 2:
                    continue
                
                # Extract property and value
                agn_prop = line_values[0].strip()
                agn_value = line_values[1].split()[0]
                
                # Store the property and value in the dictionary
                VObject_dict[agn_prop] = agn_value
        
        # Extract data and append to a list
        vshift_data = {
            'Object ID': source,
            'VShift Peak BHA': VObject_dict.get('Peak H_alpha Velocity Shift', ''),
            'Kurtosis BHA': VObject_dict.get('H_alpha Kurtosis', ''), 
            'VShift Center BHA': VObject_dict.get('Center H_alpha Velocity Shift', ''),
            'BHA Asymmetry': VObject_dict.get('H_alpha Asymmetry', ''),
            'VShift Peak BHB': VObject_dict.get('Peak H_beta Velocity Shift', ''),
            'VShift Center BHB': VObject_dict.get('Center H_beta Velocity Shift', ''),
            'Kurtosis BHB': VObject_dict.get('H_beta Kurtosis', ''),
            'BHB Asymmetry': VObject_dict.get('H_beta Asymmetry', '')
        }
        
        # Append the extracted data to a list
        LineShifts_output_data.append(vshift_data)
        
        # Specify the CSV file path
        #there are two here because I wanted to separate some of the params for now
        #Can change this later
        csv_file_path1 = 'SN_LineProfiles_outputTEST.csv'
        csv_file_path2 = 'SN_LineShifts_outputTEST.csv'
        
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
#ADDING CHOICE OF LOOP NOT WORKING RIGHT NOW
loop = True
if loop:
    for source in Data_list:
        output_file(source)
else:
    output_file('0813-52354-0561') 

      
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
           



