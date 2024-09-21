#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 14:02:50 2024

@author: ibarisic
"""

import os
import numpy as np
import glob
import shutil
from pre_processing import stsci_stage1, clean_bkg
from processing import post_processing as pp
import cube_build.combine_crossdisp as cc
from cube_build.build_7_exposures import *
from processing import stsci_stage2 as st2
from multiprocessing import Pool
from processing import transform_msa_sourceids as trans_msa
#from processing.nsclean import *
import argparse



class MSA3D:
    def __init__(self, data_path, msa_path, N_exposures = 63, rows = 7, columns = 9, run_preprocess = False, run_process = False, run_postprocess = False, run_cubebuild = False):
        self.run_preprocess = run_preprocess
        self.run_process = run_process
        self.run_postprocess = run_postprocess
        self.run_cubebuild = run_cubebuild
        
        self.base_path = os.path.dirname(data_path[0])
        print(self.base_path)
        
        print('running postprocessing')
        
        if self.run_postprocess:
            mid_exposure = int(np.ceil(N_exposures/2)) - 1
            print(mid_exposure)
            
            if self.run_process or self.run_preprocess:
                data_path_2d = np.sort(glob.glob(self.base_path + '/reduction/process/exp*nobar/newoutput*_s2d.fits'))
                folder_path = np.sort(glob.glob(self.base_path + '/reduction/process/exp*nobar/'))
                print(len(data_path_2d), len(folder_path), os.path.dirname(data_path_2d[0]))
                
                
            else:
                # base_path = /Users/ibarisic/Downloads/GO_2136_ALL_DATA/JWST/reduction/process/exp*nobar/newoutput*_s2d.fits
                data_path_2d = data_path
                folder_path = np.sort(glob.glob(self.parent(data_path_2d[0], 1) + '/exp*nobar/'))
            
            single_exposure = np.sort(glob.glob(os.path.dirname(folder_path[mid_exposure]) + '/newoutput*_s2d.fits'))
            print(len(single_exposure), single_exposure[0])
            path_dict = pp.pathloss_dict(single_exposure, msa_path)
            pp.do_pathloss_astroscrappy(data_path_2d, path_dict)
            
        
        if self.run_cubebuild:
            if self.run_process:
                exposure_path = np.sort(glob.glob(self.base_path + '/reduction/process/exp*nobar/'))
                directory = self.base_path + '/reduction/process/'
            else:
                exposure_path = np.sort(glob.glob(self.parent(data_path[0], 1) + '/exp*nobar/'))
                directory = self.parent(exposure_path[0], 1) + '/'
                
            print('this is the ', exposure_path, directory)
                
            groups = cc.generate_groups(exposure_path, rows, columns)
            mid_col, mid_row = int(np.ceil(columns/2)-1), int(np.ceil(rows/2)-1)
            ex_path1, ex_path2 = glob.glob(groups[mid_col][mid_row] + 'newoutput_s*_s2d.fits'), glob.glob(groups[mid_col][mid_row-1] + 'newoutput_s*_s2d.fits')
            
            psize_dict = cc.psize_native_dict(ex_path1, ex_path2)
            print(psize_dict)
            
            save_to = self.parent(directory, 1) + '/'
            cc.build_cube(directory, save_to, 'exp_*_nobar', psize_dict, rows, columns, True)
            
            
            
    def parent(self, path, l=0):
        parent_path = os.path.dirname(path)
        if l == 0:
            return parent_path
        return self.parent(parent_path, l-1)
    
    
    
    
def st2_multiprocessing(data, members):
    reshaped_entry = np.array(data).reshape(63, 2)
    st2_groups = st2.divide_into_groups(reshaped_entry, members)
    print(st2_groups)
        
    workers = len(st2_groups)
    with Pool(processes=workers) as pool:
        #result = pool.map(st2.run_stage2, st2_groups)        
        pool.map(st2.run_stage2, st2_groups)        


      
        

        
def run(data_path, msa_path, run_process=False, run_postprocess=True, run_cubebuild=True, N_exposures = 63, rows = 7, columns = 9, N_gmembers = 9, fname = 'reduction', transform_msa = True, run_stage1 = False):
    #data_path = np.sort(glob.glob(data_path_))
    #print(data_path)
    
    if transform_msa:
        trans_msa.transform(msa_path, True, msa_path)
    
    
    if run_stage1:
        stage1 = stsci_stage1.PreProcess(data_path, fname = fname, do_stage1=True, do_nsclean=True)
        #stage1.do_round(data_path)
        store_stage1 = stage1.output_dir
        print(store_stage1)

        
    if run_process:
        if run_stage1:
            shutil.copy(msa_path, store_stage1)
            data_path = np.sort(glob.glob(os.path.join(store_stage1, 'jw*rate.fits')))
            print(data_path)
        st2_multiprocessing(data_path, N_gmembers)
        
        
    if run_postprocess or run_cubebuild:
        MSA3D(data_path = data_path, msa_path = msa_path, N_exposures=N_exposures, rows = rows, columns = columns, run_preprocess=run_stage1, run_process = run_process, run_postprocess=run_postprocess, run_cubebuild = run_cubebuild)
        
        
    
        
        
        
data_entries = np.sort(glob.glob('/Users/ibarisic/Downloads/GO_2136_ALL_DATA/JWST/jw*rate.fits')) 
msa_path = '/Users/ibarisic/Downloads/GO_2136_ALL_DATA/JWST/jw02136001001_01_msa.fits'

#trans_msa.transform(msa_path, True, msa_path)
#print('transformed MSA file')

'''
if __name__ == '__main__':
    
    # if stage 1
    #data_entries = glob.glob('/Users/ibarisic/Downloads/GO_2136_ALL_DATA/JWST/jw02136001001_0311n_00001_*_uncal.fits')
    
    # if stage 2 and up
    #data_entries = np.sort(glob.glob('/Users/ibarisic/Downloads/GO_2136_ALL_DATA/JWST/jw*rate.fits')) 
    #data_entries = np.sort(glob.glob('/Users/ibarisic/Downloads/GO_2136_ALL_DATA/JWST/reduction/pre_processing/jw*rate.fits')) 
    #data_entries = np.sort(glob.glob('/Users/ibarisic/Downloads/GO_2136_ALL_DATA/JWST/reduction2/test_run/reduction/process/exp*nobar/newoutput*_s2d.fits')) 

    data_entries = np.sort(glob.glob('/Users/ibarisic/Downloads/GO_2136_ALL_DATA/JWST/reduction/process/exp*nobar/newoutput*_s2d.fits')) 
    msa_path = '/Users/ibarisic/Downloads/GO_2136_ALL_DATA/JWST/jw02136001001_01_msa.fits'
    run(data_entries, msa_path)



    #parser = argparse.ArgumentParser(description="msa3d reduction")

    #parser.add_argument('data_path', type=str, help='Input file path')
    #parser.add_argument('msa_path', type=str, help='Input msa file path')
    #parser.add_argument('stage2', type=str, help='Run Stage2/3')
    #parser.add_argument('postprocess', type=str, help='Run postprocessing')
    #parser.add_argument('cube_build', type=str, help='Run cube build')
    #parser.add_argument('--verbose', '-v', action='store_true', help='Enable verbose mode')   

    #args = parser.parse_args()
    
    #print(f"Input file: {args.data_path}")
    #print(f"Input msa file: {args.msa_path}")
    #print(f"Run Stage2/3: {args.stage2}")
    #print(f"Run postprocessing: {args.postprocess}")
    #print(f"Run cube build: {args.cube_build}")
    #if args.verbose:
    #    print("Verbose mode enabled")
    
    #run(str(args.data_path), str(args.msa_path), run_process = args.stage2, run_postprocess = args.postprocess, run_cubebuild = args.cube_build)
'''    
    

        
            
    