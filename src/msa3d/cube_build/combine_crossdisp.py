#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 19:00:29 2023

@author: ibarisic
"""


import numpy as np
from .build_7_exposures import Build7exp
import glob
import matplotlib.pyplot as plt
from astropy.io import fits
from jwst import datamodels
from astropy.coordinates import SkyCoord
from astropy import units as u
import os


def separation_positionA(ra1, dec1, ra2, dec2):
    c1 = SkyCoord(ra1*u.deg, dec1*u.deg, frame='icrs')
    c2 = SkyCoord(ra2*u.deg, dec2*u.deg, frame='icrs')
    return c1.separation(c2).arcsecond


def extract_ra_dec(item):
    s2d_ = datamodels.open(item)
    s2dsci_ = s2d_.data
    wcsobj_ = s2d_.meta.wcs
    newy, newx = np.mgrid[:s2dsci_.shape[0], :s2dsci_.shape[1]]
    det2sky = wcsobj_.get_transform('detector', 'world')
    ra_, dec_, s2dwave_ = det2sky(newx, newy)
    return ra_[-2,0], dec_[-2,0], ra_[-1,0], dec_[-1,0]


def psize_native_dict(path1, path2):
    PSIZE_NATIVE_DICT = {}
    
    path1, path2 = np.sort(path1), np.sort(path2)
    for i in range(len(path1)):
        name = path1[i].split('/')[-1].split('_')[1]
        idx = int(name[1:])
        
        if idx > 100:
        
            ra1, dec1, ra2, dec2 = extract_ra_dec(path1[i])        
            dst = separation_positionA(ra1, dec1, ra2, dec2)
            
            PSIZE_NATIVE_DICT[idx] = round(dst, 5)
        
    return PSIZE_NATIVE_DICT


def generate_groups(whole_list, rows, columns):
    group_array = [x[:] for x in [[0] * rows] * columns]
    for i in range(columns):
        for j in range(rows):   
            if j%2 == 0:
                group_array[i][j] = whole_list[(j*9)+i]
            elif j%2 == 1:
                group_array[-(i+1)][j] = whole_list[(j*9)+i]
    return group_array


def make_cube_dir(folder_path, source_id):
    if os.path.exists(folder_path + 'cubes/cube_' + str(source_id)):
        print('exists')
        new_path = folder_path + 'cubes/cube_' + str(source_id) + '/'
    else:
        os.makedirs(folder_path + 'cubes/cube_' + str(source_id))
        new_path = folder_path + 'cubes/cube_' + str(source_id) + '/'
    return new_path


def parent(path, l=0):
    parent_path = os.path.dirname(path)
    if l == 0:
        return parent_path
    return parent(parent_path, l-1)


def build_cross_disp(spectra, folder_destination, psize_dict, source_id, groups):
    new_path = make_cube_dir(folder_destination, source_id)
    builder = Build7exp(path_destination=new_path, extension=1, psize_native=psize_dict[source_id])
    # this is done once for one galaxy
    builder.sample_wave(spectra)
    builder.common_grid(spectra[0])
    
    for i in range(len(groups)):
        N = i#8
        print(f"this is iteration : {N}")
        groups[N]
        median_grid = builder.make_7(groups[N], crossdisp_index = str(N+1), save_data=True)


def build_disp(dispersion_count, medians, source_id, save):
    x, y = np.shape(fits.getdata(medians[0]))
    #dispersion_count = M
    dispersion_grid = np.zeros([y, x, dispersion_count])
    print(np.shape(dispersion_grid), np.shape(dispersion_grid[:, : , 0]))
    
    for i, item in enumerate(medians):
        data = fits.getdata(item)
        print(item)
        print(len(data[:,0]), len(data[0,:]))
        print(data, np.shape(data), len(data[0]))
        data_ = np.transpose(data)
        print(data_, np.shape(data_), len(data_[0]))
        #print(data[0,:], data[:,0])
        dispersion_grid[:, :, i] = data_
    
    if save:
        fits.writeto(os.path.dirname(medians[0]) + '/cube_medians_' + str(source_id) + '.fits', dispersion_grid, overwrite = True)
        
        
def build_cube(directory, save_to, extension, psize_dict, rows, columns, save):
    for source_idx in psize_dict.keys():
        #if source_idx == 2465:
        #    break
        
        if source_idx < 10000:
            paths = np.sort(glob.glob(directory + extension + '/newoutput_s'+'0'+str(source_idx)+'_s2d_pathcorr_astrocorr.fits'))
        else:
            paths = np.sort(glob.glob(directory + extension + '/newoutput_s'+str(source_idx)+'_s2d_pathcorr_astrocorr.fits'))
        
        print(paths)
        
        full_group = generate_groups(paths, rows, columns)
        print(source_idx)
        print(full_group)
    
        # for a given source build the medians
        build_cross_disp(paths, save_to, psize_dict, source_idx, full_group)
        # read in the medians of said source
        all_medians = np.sort(glob.glob(save_to + 'cubes/cube_'+str(source_idx)+'/median*.fits'))
        print(all_medians)
        # make a cube
        build_disp(columns, all_medians, source_idx, save)
        
        
        
'''
ex_path1 = glob.glob('/Users/ibarisic/Downloads/go-2136_download_central_2dets/JWST/cleaned_mos_latest/exp_0311r_barshadow/newoutput_s*_s2d.fits')
ex_path2 = glob.glob('/Users/ibarisic/Downloads/go-2136_download_central_2dets/JWST/cleaned_mos_latest/exp_03119_barshadow/newoutput_s*_s2d.fits')

psize_dict = psize_native_dict(ex_path1, ex_path2)
print(psize_dict)

# building a do_round of build_medians for one galaxy
#all_paths = np.sort(glob.glob('/Users/ibarisic/Downloads/go-2136_download_central_2dets/JWST/cleaned_mos_latest/exp_*_barshadow/newoutput_s'+'0'+str(source_idx)+'_s2d_pathcorr_astrocorr.fits'))
directory = '/Users/ibarisic/Downloads/go-2136_download_central_2dets/JWST/cleaned_mos_latest/'

#!!!!!!!!!!!!!!!!!
destination = '/Users/ibarisic/Downloads/GO_2136_ALL_DATA/JWST/test_preprocess/' # !!!!!!!!!!!

build_cube(directory, destination, 'exp_*_barshadow', psize_dict, 7, 9, False)

ex_path1 = glob.glob('/Users/ibarisic/Downloads/go-2136_download_central_2dets/JWST/cleaned_mos_latest/exp_*_barshadow/newoutput_s08512_s2d.fits')
full_group = generate_groups(ex_path1, 7, 9)
print(np.shape(full_group))
print(full_group)
'''







    
    



