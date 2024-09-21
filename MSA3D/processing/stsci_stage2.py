#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 14:00:29 2023

@author: ibarisic

This JWST processing algorithm is based on 

(A) STScI JWST Stage 2 reduction pipeline 

where we use custom set of steps to reduce rate files
    

(A) Copyright (C) 2020 Association of Universities for Research in Astronomy (AURA)

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    1. Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    2. Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.

    3. The name of AURA and its representatives may not be used to
      endorse or promote products derived from this software without
      specific prior written permission.
      


"""

import numpy as np
import glob
import os
import time
#print(os.environ['CRDS_PATH'])
os.environ["CRDS_SERVER_URL"] = "https://jwst-crds.stsci.edu"

# testing pmap choice
#os.environ["CRDS_CONTEXT"] = "jwst_1075.pmap"


#import asdf
#import json
#from jwst import datamodels
from jwst.associations import asn_from_list
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base

#import astropy.io.fits as pyfits
#from astropy.io import fits
#from astropy.utils.data import get_pkg_data_filename
#import matplotlib.pyplot as plt

#import jwst
from jwst.pipeline import Spec2Pipeline, Spec3Pipeline
from multiprocessing import Pool




class CustomPipeline:
    def __init__(self, path, fname, do_stage2=True, do_stage3=True, barshadow=False, pathloss=False):
        #self.det = None
        self.do_stage2 = do_stage2
        self.do_stage3 = do_stage3
        self.make_dir(path, fname)
       
    
    def do_round(self, path):
        if self.do_stage2:
            for item in path:
                self.spec_2(item)
        if self.do_stage3:
            self.spec3(path)
        
    
    def spec_2(self, path):
        spec2 = Spec2Pipeline()
        
        spec2_param_reffile = os.path.join(self.output_dir, 'calwebb_spec2.asdf')
        spec2.export_config(spec2_param_reffile)
        
        spec2.save_results = True
        spec2.output_dir = self.output_dir
        
        self.spec2_steps = [spec2.bkg_subtract, spec2.master_background_mos, 
                            spec2.pathloss, spec2.barshadow, 
                            spec2.wavecorr, spec2.srctype, 
                            spec2.straylight, spec2.fringe, spec2.residual_fringe, 
                            spec2.cube_build, spec2.extract_1d]
        
        #self.add_step(spec2, barshadow, pathloss)        
        self.skip(self.spec2_steps)
        
        # NSClean addition from the STScI pipeline
        #spec2.nsclean.skip = False
        #spec2.nsclean.save_results = True
        
     
        result = spec2.run(path)


    def skip(self, step):
        for item in step:
            print(item)
            item.skip = True
            
            
    def add_step(self, pipe, barshadow, pathloss):
        if barshadow:
            self.spec2_steps.append(pipe.barshadow)
        if pathloss:
            self.spec2_steps.append(pipe.pathloss)
            
            
    def spec3(self, path):
        spec3 = Spec3Pipeline()
        
        spec3_param_reffile = os.path.join(self.output_dir, 'calwebb_spec3.asdf')
        spec3.export_config(spec3_param_reffile)
        
        self.combine_asn()
        
        spec3.save_results = True
        spec3.output_dir = self.output_dir

        self.spec3_steps = [spec3.assign_mtwcs, spec3.master_background, 
                            spec3.mrs_imatch, spec3.outlier_detection,
                            spec3.cube_build, spec3.extract_1d, 
                            spec3.combine_1d, spec3.photom]

        self.skip(self.spec3_steps)
        
        
        asn = self.output_dir + '/manual_calwebb3_asn.json'
        spec3.run(asn)
        
    
    def combine_asn(self):
        input_files = np.sort(glob.glob(self.output_dir + '/jw*_cal.fits'))
        out_asn = asn_from_list.asn_from_list(input_files,  rule=DMS_Level3_Base, product_name="newoutput")

        output_asn = os.path.join(self.output_dir, "manual_calwebb3_asn.json")

        with open(output_asn, "w") as outfile:
            name, serialized = out_asn.dump(format='json')
            outfile.write(serialized)
            
    
    def make_dir(self, path, fname):
        if len(path)==1:
            dir_name = os.path.dirname(path)
            if os.path.split(dir_name)[-1] == 'nsclean':
                newpath = str(self.parent(path, 1)) + '/process/' + fname
            else:
                if os.path.split(dir_name)[-1] == 'reduction':
                    newpath = str(self.parent(path)) + '/process/' + fname
                elif os.path.split(self.parent(dir_name))[-1] == 'reduction':
                    newpath = str(self.parent(path, 1)) + '/process/' + fname
                else:
                    newpath = str(self.parent(path)) + '/reduction/process/' + fname
                    #newpath = str(self.parent(path, 1)) + '/process/' + fname

        else:
            dir_name = os.path.dirname(path[0])
            if os.path.split(dir_name)[-1] == 'nsclean':
                newpath = str(self.parent(path[0], 1)) + '/process/' + fname
            else:
                if os.path.split(dir_name)[-1] == 'reduction':
                    newpath = str(self.parent(path[0])) + '/process/' + fname
                elif os.path.split(self.parent(dir_name))[-1] == 'reduction':
                    newpath = str(self.parent(path[0], 1)) + '/process/' + fname
                else:
                    newpath = str(self.parent(path[0])) + '/reduction/process/' + fname
                    #newpath = str(self.parent(path[0], 1)) + '/process/' + fname

            
        print(newpath)
            
        if not os.path.exists(newpath):
            os.makedirs(newpath)
            
        self.output_dir = newpath
        
        
    def parent(self, path, l=0):
        parent_path = os.path.dirname(path)
        if l == 0:
            return parent_path
        return self.parent(parent_path, l-1)
        
        
    


def divide_into_groups(data, group_size):
    return [data[i:i + group_size] for i in range(0, len(data), group_size)]


def run_stage2(group):
    for pair in group:
        arg = pair[0].split('/')[-1][14:19]
        print('this is', pair, arg)
        pipe = CustomPipeline(pair, 'exp_'+arg+'_nobar')
        pipe.do_round(pair)
    
    

#arg = '0311n'
#file_name = ['/Users/ibarisic/Downloads/GO_2136_ALL_DATA/JWST/test_preprocess/jw02136001001_' + arg + '_00001_nrs1_rate.fits',
#             '/Users/ibarisic/Downloads/GO_2136_ALL_DATA/JWST/test_preprocess/jw02136001001_' + arg + '_00001_nrs2_rate.fits']

#pipe = CustomPipeline(file_name, 'exp_'+arg+'_nobar')
#pipe.do_round(file_name)

'''
if __name__ == '__main__':
    
    data_entries = np.sort(glob.glob('/Users/ibarisic/Downloads/GO_2136_ALL_DATA/JWST/jw*rate.fits'))
    reshaped_entry = np.array(data_entries).reshape(63, 2)
        
    groups = divide_into_groups(reshaped_entry, 9)
    print(groups)
    
    workers = len(groups)
    with Pool(processes=workers) as pool:
        result = pool.map(run_stage2, groups)
'''



