#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 15:36:52 2023

@author: ibarisic
"""

"""

prepare the msa.fits file ; replace negative source-ids with a valid id for a given slitlet-id

""" 

import astropy.io.fits as pyfits
from astropy.io import fits
import numpy as np



def transform(file, save=False, path=None):
    data = pyfits.open(file)
    shutter_info = data['SHUTTER_INFO'].data
    
    all_targets = np.unique(shutter_info['source_id'])
    targets = all_targets[np.where(all_targets > 0)]
        
    for target in targets:
        index = np.where(shutter_info['source_id']==target)[0][0]
        
        match_slitlet = shutter_info['slitlet_id'][index]
        if target < 100000:
            shutter_info['source_id'][shutter_info['slitlet_id'] == match_slitlet] = target
        else:
            shutter_info['source_id'][shutter_info['slitlet_id'] == match_slitlet] = float(str(target)[:-3])
        
    if save:
        output_file = path
    
        data.writeto(output_file, overwrite=True)
        data.close()
        

    
