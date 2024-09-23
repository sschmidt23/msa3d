#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 15:18:33 2023

@author: ibarisic


This JWST pre-processing algorithm is based on 

(A) STScI JWST Stage 1 reduction pipeline and
(B) NSClean algorithm (Benjamin Rauscher) 

where in (A) we use custom set of steps to reduce uncalibrated files, 
and (B) in its original form.
    

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

import glob
import os
#print(os.environ['CRDS_PATH'])
os.environ["CRDS_SERVER_URL"] = "https://jwst-crds.stsci.edu"

from jwst.pipeline import Detector1Pipeline
#from .clean_bkg import *
from . import clean_bkg


class PreProcess:
    def __init__(self, path, fname, do_stage1=False, do_nsclean=False):
        self.preprocess_dir(path, fname)
        self.do_stage1 = do_stage1
        self.do_nsclean = do_nsclean

    
    def do_round(self, path):
        
        if self.do_stage1:
            for uncal in path:
                self.run_stage1(uncal)
                
        if self.do_nsclean:
            
            path_to_clean = glob.glob(self.output_dir + '/jw*_rate.fits')
            output_path = os.path.join(self.output_dir + '/nsclean/')
            self.make_dir(output_path)
            print(path_to_clean, output_path)
            
            clean_bkg.clean(path_to_clean, output_path)
                
        
        elif self.do_nsclean and not self.do_stage1:
            path_to_clean = glob.glob(str(os.path.dirname(path[0])) + '/jw*_rate.fits')
            output_path = os.path.join(self.output_dir + '/nsclean/')
            print(path_to_clean, output_path)
            self.make_dir(output_path)
            
            clean_bkg.clean(path_to_clean, output_path)
            
            
    def run_stage1(self, filename):
        det1 = Detector1Pipeline()
    
        det_param_reffile = os.path.join(self.output_dir, '/manual_calwebb_det1.asdf')
    
        det1.save_results = True
        det1.output_dir = self.output_dir
    
        det1.jump.expand_large_events = True
        det1.jump.sat_required_snowball = False
        det1.jump.min_jump_area = 10.0
        
        det1.run(filename)
        
        
    def preprocess_dir(self, path, fname):
        if len(path)==1:
            newpath = str(os.path.dirname(path)) + '/' + fname + '/pre_process'
        else:
            newpath = str(os.path.dirname(path[0])) + '/' + fname + '/pre_process'
            
        self.make_dir(newpath)
        self.output_dir = newpath
        
    
    def make_dir(self, trace):
        if not os.path.exists(trace):
            os.makedirs(trace)

    

'''
#level0_data_files = glob.glob('/Users/ibarisic/Downloads/GO_2136_ALL_DATA/test_stage1/*uncal.fits')
level0_data_files = glob.glob('/Users/ibarisic/Downloads/GO_2136_ALL_DATA/JWST/jw02136001001_0311n_00001_*_uncal.fits')
level0_data_files.sort()
print(level0_data_files)



stage1 = PreProcess(level0_data_files, fname = 'reduction', do_stage1=True, do_nsclean=True)
stage1.do_round(level0_data_files)
''' 
    