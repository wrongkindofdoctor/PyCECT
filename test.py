#!/usr/bin/env python
import ConfigParser
import sys, getopt, os 
import numpy as np 
import Nio 
import time
import re
from asaptools.partition import EqualStride, Duplicate,EqualLength
import asaptools.simplecomm as simplecomm 
import pyEnsLib

# This routine creates a summary file from an ensemble of AM4 history files
# Put output time frequency options (tslice=0,1,2) in a dictionary
output_dict=dict()
output_dict[0] = ['1x12m0d']
output_dict[1] = ['1x1m0d','1x1m0d']   
output_dict[2] = ['1x0m2d','1x0m8d','2x0m1d']

#output files
# output directory tree is root/platform.compiler-debug/prod/repro(-openmp)/duration_pelayout/history/
input_dir='/net2/Jessica.Liptak/'
in_files=[]
#
# Form ensembles, each missing one member; compute RMSZs and global means
#    #for each variable, we also do max norm also (currently done in pyStats)
#    tslice = opts_dict['tslice']
tslice=1
time_out=output_dict[tslice]
#print time_out

if(os.path.exists(input_dir)):
  # Get the list of files
  in_dirs_temp = os.listdir(input_dir)
  #print in_dirs_temp
  file_list =[]
  # get directories with output for desired time average
  for d in in_dirs_temp:
      for root,subdirs,files in os.walk(os.path.join(input_dir,d)):
          if len(files)>0:
             for f in files:
                 file_path=os.path.join(root,f)
                 for t in time_out:
                     if t in file_path and '19790101' in file_path:
                        #print file_path
                        file_list.append(file_path)
print file_list

          #print sd
#      for f in files: 
#          print f
#          extension = '.'.join(f.rsplit('.')[1:])
#          #print extension
##          if extension == 'nc':
#             for p in platform_list:
#                 for t in target_list:
#                     if p in root and t in root:
                        #print os.path.join(root,f)
#                        in_files.append(os.path.join(root,f))
  
#in_files=sorted(in_files_temp)
#print 'in_files=',in_files

