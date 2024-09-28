#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 17 12:51:50 2023

@author: ibarisic
"""

from astropy.io import fits
from time import perf_counter as t # Used for benchmarking

# It will be handy to have a stack
stack = []

# Import the NIRSpec Clean package
#import nsclean as nc
#from .nsclean import *
from msa3d.nsclean import NSClean1

opmode = 'MOS' # Set this to 'MOS' or 'IFU'
detector = 'NRS1'
#cleaner = nc.NSClean1(opmode, detector)
cleaner = NSClean1(opmode, detector)


    
def clean(files, output):
    
    for file in files:
        
        #with fits.open(file, mode='update') as hdul:
        with fits.open(file) as hdul:
            hdul[1].data = modified_data(file)
            #hdul.flush() 
            hdul.writeto(output + file.split('/')[-1], overwrite='True')
            
            
def modified_data(file):
    
        with fits.open(file) as hdul:
            D  = hdul[1].data
            H0 = hdul[0].header
            H1 = hdul[1].header
        
            stack.append(t()) # Setup to benchmark
            D = cleaner.clean(D) # This is the line that does all the work
            print('Elapsed time (s) = ', t()-stack.pop()) # Benchmark it!
            #print(D)
            return D
            
