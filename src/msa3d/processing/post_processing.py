#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 13:00:14 2023

@author: ibarisic
@co-contributor: tnanayakkara
"""

from importlib import resources
import astroscrappy
import os
import matplotlib.pyplot as plt
from astropy.io import fits
from jwst import datamodels
import numpy as np
import grizli.utils
from msaexp import msa
from astropy.table import Table
import glob


def get_pathloss_removal(path, gal_id, sci_wave, n_of_slits = 3):
    
    """
    This code will compute the pathloss vector for a uniformly illuminated slit. 
    So when the data is reduced as point sources with no pathloss correction applied, 
    applying the uniform correction will remove the flux calibration effect due to a centered point source., 
    
    """
    
    #### technically all are considered as uniform sources in DR. 
    #meta_file = msa.MSAMetafile(f'../meta_data_files/jw02565{obsid}001_01_msa.fits')
    meta_file = msa.MSAMetafile(path)

    msa_metadata_slit_info = Table.read(path, hdu=2).to_pandas()
    msa_metadata_source_info = Table.read(path, hdu=3).to_pandas()
    
    ### only necessary if applying/removing PS corrections. 
    shutter_xs = meta_file.shutter_table[meta_file.shutter_table['source_id']==gal_id]['estimated_source_in_shutter_x']
    shutter_ys = meta_file.shutter_table[meta_file.shutter_table['source_id']==gal_id]['estimated_source_in_shutter_y']
    
        
    from jwst import pathloss
    from jwst.datamodels.pathloss import PathlossModel

    #ref_model = PathlossModel('/Users/ibarisic/Downloads/GO_2136_ALL_DATA/test_stage1/reduced/pathloss/jwst_nirspec_pathloss_0002.fits')
    #pathloss_ref =  fits.open('/Users/ibarisic/Downloads/GO_2136_ALL_DATA/test_stage1/reduced/pathloss/jwst_nirspec_pathloss_0002.fits')

    ref_model_file = resources.files('MSA3D') / 'processing' / 'jwst_nirspec_pathloss_0002.fits'
    ref_model = PathlossModel(ref_model_file)
    pathloss_ref =  fits.open(ref_model_file)

    # for software versions below 1.12 and cycle 1 data, instead of 'OPEN' argument, n_of_slits was used
    # assuming class definition / argument use was altered in the later pipeline versions
    pathloss_wcs = pathloss.pathloss.get_aperture_from_model(ref_model, 'OPEN').uniform_wcs
    pathloss_model = pathloss.pathloss.get_aperture_from_model(ref_model, 'OPEN').uniform_data
    uniform = True
    

    pathloss_wavelength, pathloss_vector, is_inside_slitlet = pathloss.pathloss.calculate_pathloss_vector(
        pathloss_model, pathloss_wcs, np.median(shutter_xs)-0.5,  np.median(shutter_ys)-0.5, calc_wave=True)
    

    pathloss_wavelength = pathloss_wavelength/1e-6

    
    pl_interp = np.interp(sci_wave, pathloss_wavelength, pathloss_vector )
    
    print(pl_interp, uniform)

    
    return pl_interp, uniform


def round_down_to_nearest(number, base):
    return (number // base) * base


def extract_wave(item):
    s2d = datamodels.open(item)
    s2dsci = s2d.data
    wcsobj = s2d.meta.wcs
    y, x = np.mgrid[:s2dsci.shape[0], :s2dsci.shape[1]]
    det2sky = wcsobj.get_transform('detector', 'world')
    ra, dec, s2dwave = det2sky(x, y)
    return s2dwave[0]


def pathloss_dict(spectra, path_to_msa):
    all_wave, all_pathloss, sources = [], [], []
    mockwave, mockpathloss = [], []
    
    for item in spectra:
        SOURCE = int(item.split('/')[-1].split('_')[1][1:])
        if SOURCE > 100:
            wave = extract_wave(item)
            
            print(item, SOURCE)
            pathloss, value = get_pathloss_removal(path_to_msa, SOURCE, wave)
            
            all_wave.append(wave)
            all_pathloss.append(pathloss)
            sources.append(SOURCE)               
        
    path_dict = {}    
    for key, value in zip(sources, all_pathloss):
        path_dict[key] = value

            
    return path_dict


def clean_astroscrappy(spectra):
    
    with fits.open(spectra) as hdul:
        fits_spectra  = hdul[1].data.transpose()

        spectra_size = np.shape(fits_spectra)
        spectra_clean = np.zeros(np.shape(fits_spectra))
        spectra_mask = np.zeros(np.shape(fits_spectra))
        spectra_clean[:,:] = fits_spectra
        print(spectra_size)
        
        
        w_min, w_max = 0, 100
        w_limit = round_down_to_nearest(spectra_size[0], 100)
        print(w_limit)
        while w_max < w_limit:
            mask, clean_data = astroscrappy.detect_cosmics(fits_spectra[w_min:w_max, :], sigclip=0.6, sigfrac=0.001, niter = 3)
            spectra_mask[w_min:w_max, :] = mask
            spectra_clean[w_min:w_max, :] = clean_data
            w_min += 100
            w_max += 100
            
        new_path = spectra.split('.')[0]
       
    
    #a_ratio = 30
    #plt.figure(figsize=(10,15))
    #plt.imshow(np.transpose(fits_spectra[:, :]), origin='lower', aspect=a_ratio, vmin=-0.2, vmax=1)
    #plt.figure(figsize=(10,15))
    #plt.imshow(np.transpose(spectra_mask[:, :]), origin='lower', aspect=a_ratio, vmin=-0.2, vmax=1)
    #plt.figure(figsize=(10,15))
    #plt.imshow(np.transpose(spectra_clean[:, :]), origin='lower', aspect=a_ratio, vmin=-0.2, vmax=1)
    
    return spectra_clean


def do_pathloss_astroscrappy(spectra, path_dict):
        
    for item in spectra:
        source_id = int(item.split('/')[-1].split('_')[1][1:])
        if source_id > 100:
            
            path_name = item[:-5] + '_pathcorr_astrocorr.fits'
            print(source_id, item, path_dict[source_id])
            
            if os.path.exists(path_name):
                print('already known')
                with fits.open(item, mode = 'readonly') as hdul_orig:
                    with fits.open(item[:-5] + '_pathcorr_astrocorr.fits', mode='update') as hdul:
                        for i, extension in enumerate(hdul):
                            if i == 1:
                                #plt.figure(figsize = (10,15))
                                #plt.imshow(hdul_orig[1].data, aspect = 20, vmin=-0.2, vmax=1, origin = 'lower')
                            
                                data = clean_astroscrappy(item).transpose()
                            
                                hdul[1].data = data/path_dict[source_id]
                                #plt.figure(figsize = (10,15))
                                #plt.imshow(hdul[1].data, aspect = 20, vmin=-0.2, vmax=1, origin = 'lower')
                        hdul.flush()
                
    
            else:
                print('not known')
                with fits.open(item, mode = 'readonly') as hdul_orig:
                    #print(item[:-5] + '_pathcorr.fits')
    
                    hdul_new = fits.HDUList()
    
                    #with fits.open(item[:-5] + '_pathcorr.fits', mode='append', memmap=True) as hdul_new:
                    for extension in hdul_orig:
                        hdul_new.append(extension.copy())
    
                    for i, extension in enumerate(hdul_new):
                        if i == 1:
                            #for row in range(np.shape(hdul_new[1].data)[0]):
                            #    hdul_new[1].data[row,:] = hdul_new[1].data[row,:] / path_dict[source_id]
    
                            #plt.figure(figsize = (10,15))
                            #plt.imshow(hdul_orig[1].data, aspect = 20, vmin=-0.2, vmax=1, origin = 'lower')
                            
                            data = clean_astroscrappy(item).transpose()
                            
                            hdul_new[1].data = data/path_dict[source_id]
                            #plt.figure(figsize = (10,15))
                            #plt.imshow(hdul_new[1].data, aspect = 20, vmin=-0.2, vmax=1, origin = 'lower')
    
                    hdul_new.writeto(path_name, overwrite=True)




# running it 
'''
#using CYCLE 1 data (reduced via 1-14)
msa_path = '/Users/ibarisic/Downloads/GO_2136_ALL_DATA/JWST/test_preprocess/jw02136001001_01_msa.fits'
all_spec = np.sort(glob.glob('/Users/ibarisic/Downloads/GO_2136_ALL_DATA/JWST/test_preprocess/process/exp_0311n_nobar/newoutput_*_s2d.fits'))

# in the implementation the line below should only be executed once - to find the dictionary for all galaxies
path_dict = pathloss_dict(all_spec, msa_path)
# based on which you can get pathloss+astroscrappy for all
do_pathloss_astroscrappy(all_spec, path_dict)
'''