#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 15:56:54 2023

@author: ibarisic
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy
import scipy.ndimage as snd
from scipy.interpolate import RegularGridInterpolator
from astropy.io import fits
from jwst import datamodels

import glob


class Build7exp:
    def __init__(self, path_destination, extension, psize_native, n_fine = 4, nspec = 7, offset = 4.5):
        # offset : 0.075" moved 6x in dispersion direction = 4.5" offset
        #self.path = path
        self.psize_native = psize_native
        self.fine_gridN = n_fine
        self.psize_sample = psize_native / n_fine
        self.NSPEC = nspec
        self.OFFSET = offset
        self.DITHER = 0.075
        self.destination = path_destination
        self.ext = extension
        print("initiated")
    
    
    def get_shape(self, path):
        data = fits.getdata(path, self.ext)   #[2:,:]#[:, 1740:1800]
        return np.shape(data)
        

    def interp_grid(self, path):
        data = fits.getdata(path, self.ext)
        shape = np.shape(data)

        print(shape)
        
        x = np.linspace(0, shape[0]-1, shape[0])
        #y = np.linspace(0, shape[1]-1, shape[1])
        y = self.get_ra_dec_wave(path)
        
        print(x, y, len(x), len(y))
        
        #newx = np.linspace(0, shape[0]-1, (shape[0]-1) * self.fine_gridN + 1)  # +1  -- to achieve 0.25 step
        #print(newx, len(newx))
        
        #X, Y = np.meshgrid(newx, y, indexing='ij')
    
        #print(self.X, self.Y, np.shape(self.X), np.shape(self.Y))
        
        
        if self.ext != 2 :
            interp = RegularGridInterpolator((x,y), data, bounds_error=False, fill_value=np.nan)
            output = interp((self.X, self.Y))

        else:
            interp = RegularGridInterpolator((x,y), data**2, bounds_error=False, fill_value=np.nan)
            output = np.sqrt(interp((self.X, self.Y)))

        return output
    
    
    def common_grid(self, path):
        #data = fits.getdata(path, ext=1)
        #shape = np.shape(data)
        shape = self.get_shape(path)
        
        x = np.linspace(0, shape[0]-1, shape[0])
        newx = np.linspace(0, shape[0]-1, (shape[0]-1) * self.fine_gridN + 1)  # +1  -- to achieve 0.25 step  np.arange(0, 0.27, 0.025) / (0.1009/4) 
        newy = np.arange(self.min_lam, self.max_lam, step = self.step_lam)
        
        print(newx, newy, len(newx), len(newy))
        self.X, self.Y = np.meshgrid(newx, newy, indexing='ij')
        
        self.shape_common = np.shape(self.X)
        
        #return X, Y
        
        
    def get_ra_dec_wave(self, path):
        s2d = datamodels.open(path) 
        s2data, wcs = s2d.data, s2d.meta.wcs
        
        y, x = np.mgrid[:s2d.shape[0], :s2d.shape[1]]
        
        det2sky = wcs.get_transform('detector', 'world')
        ra, dec, wave = det2sky(x, y)
        
        return wave[0]
    
    
    def sample_wave(self, all_paths):
        min_wave, max_wave, wave_step = [], [], []

        for item in all_paths:
            wave = self.get_ra_dec_wave(item)
            min_wave.append(wave[0])
            max_wave.append(wave[-1])
            wave_step.append(np.min(np.diff(wave)))
            
        self.min_lam = np.min(min_wave)
        self.max_lam = np.max(max_wave)
        self.step_lam = np.min(wave_step)
                
        #return min_wave, max_wave, wave_step
    
    
    def project_common(self, path, count, spec_grid, dx_pix, nx_interp):
        new_data = self.interp_grid(path)        
        #spec_grid[int(dx_pix[count]):int(dx_pix[count])+nx_interp,:,count] = new_data
        spec_grid[count, int(dx_pix[count]):int(dx_pix[count])+nx_interp, :] = new_data
        
        return spec_grid
    
    
    def make_7(self, all_paths, crossdisp_index, save_data=False, save_err=False, save_wght=False):
        gal_id = all_paths[0].split('/')[-1].split('_')[1]
        spec_grid, dx_pix, nx_interp = self.pixel_scaling(all_paths[0])
                
        #fig, axes = plt.subplots(2, 1, figsize = (10,10))
        for i in range(self.NSPEC):
            print(all_paths[i])
            self.project_common(all_paths[i], i, spec_grid, dx_pix, nx_interp)
            #axes[0].imshow(spec_grid[0, :, 1700:1950], vmin=-0.5, vmax=1.5, origin='lower')
            #axes[1].imshow(spec_grid[6, :, 1700:1950], vmin=-0.5, vmax=1.5, origin='lower')
        
        print(spec_grid)
        median_grid = np.median(spec_grid, axis=0)
        #mean_grid = np.mean(spec_grid, axis = 0)
        print(np.shape(median_grid))
        
        if save_data:
            fits.writeto(self.destination + 'spec_lam' + crossdisp_index + '_' + gal_id + '_all.fits', spec_grid, overwrite = True)
            fits.writeto(self.destination + 'median_lam' + crossdisp_index + '_' + gal_id + '_all.fits', median_grid, overwrite = True)
        if save_err:
            fits.writeto(self.destination + 'error' + crossdisp_index + '_' + gal_id + '_all.fits', spec_grid, overwrite = True)
        if save_wght:
            fits.writeto(self.destination + 'weights' + crossdisp_index + '_' + gal_id + '_all.fits', spec_grid, overwrite = True)
        
        
        return median_grid
    
    
    def pixel_scaling(self, path):
        # Pixel scaling
        shape = self.get_shape(path)
        scale_factor = self.psize_native / self.psize_sample
        
        nx = int((shape[0] + self.OFFSET) * scale_factor)  # pixels in spatial dimension. +0.4 arcsec to account for offsets.
        #ny = shape[1]      # pixels in spectral dimension
        ny = self.shape_common[1]
        
        print(nx, ny)
        #spec_grid = np.zeros([nx, ny, self.NSPEC])
        spec_grid = np.zeros([self.NSPEC, nx, ny])
        
        print('scale factor:', scale_factor)
        print('grid size for interpolation:', np.shape(spec_grid))
        print('number of exposures:', self.NSPEC)
        
        # Interpolation offsets: 0.075 arcsec per dither step
        #dx_arcsec = np.arange(self.NSPEC - 1 , -1, -1) * self.DITHER  # dithers along slit direction, in arcseconds
        dx_arcsec = np.arange(self.NSPEC) * self.DITHER
        dx_pix = dx_arcsec / self.psize_sample  # in pixels 
        print(dx_arcsec)
        print(dx_pix)
        
        nx_interp = np.shape(self.interp_grid(path))[0]
        print(f'this: {np.shape(self.interp_grid(path)), nx_interp}')
        
        return spec_grid, dx_pix, nx_interp
           
