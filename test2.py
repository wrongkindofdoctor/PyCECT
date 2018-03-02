#!/usr/bin/env python
import sys,getopt,os
import numpy as np
import Nio
import time
import pyEnsLib
import json
import random
import glob
import re
from datetime import datetime
from asaptools.partition import EqualStride, Duplicate
import asaptools.simplecomm as simplecomm 

# This routine creates a summary file from an ensemble of AM4 history files
def main(argv):
    # Get command line stuff and store in a dictionary
    # Get command line stuff and store in a dictionary
    s="""verbose sumfile= indir= input_globs= tslice= nPC= sigMul= 
         minPCFail= minRunFail= numRunFile= printVarTest popens
         jsonfile= mpi_enable nbin= minrange= maxrange= outfile= 
         casejson= npick= pepsi_gm test_failure pop_tol= 
         pop_threshold= prn_std_mean lev= eet= json_case=
         mach= histfolder= outputfreq= usetiles="""
    optkeys = s.split()
    try: 
        opts, args = getopt.getopt(argv, "h", optkeys)
    except getopt.GetoptError:
        pyEnsLib.EnsSum_usage()
        sys.exit(2)

    # Put command line options in a dictionary - also set defaults
    opts_dict={}
    
    # Defaults
    opts_dict['input_globs'] = ''
    opts_dict['indir'] = ''
    opts_dict['tslice'] = 0
    opts_dict['nPC'] = 20
    opts_dict['sigMul'] = 2
    opts_dict['verbose'] = False
    opts_dict['minPCFail'] = 3
    opts_dict['minRunFail'] = 2
    opts_dict['numRunFile'] = 3
    opts_dict['printVarTest'] = False
    opts_dict['popens'] = False
    opts_dict['jsonfile'] = ''
    opts_dict['mpi_enable'] = False
    opts_dict['nbin'] = 40
    opts_dict['minrange'] = 0.0
    opts_dict['maxrange'] = 4.0
    opts_dict['outfile'] = 'testcase.result'
    opts_dict['casejson'] = ''
    opts_dict['npick'] = 10
    opts_dict['pepsi_gm'] = False
    opts_dict['test_failure'] = True
    opts_dict['pop_tol'] = 3.0
    opts_dict['pop_threshold'] = 0.90
    opts_dict['prn_std_mean'] = False
    opts_dict['lev'] = 0
    opts_dict['eet'] = 0
    opts_dict['json_case'] = ''
    opts_dict['mach'] = 'gaea'
    opts_dict['histfolder'] = 1
    opts_dict['usetiles'] = 1
    opts_dict['outputfreq'] = 1

    # Call utility library getopt_parseconfig to parse the option keys
    # and save to the dictionary
    caller = 'GFDL'
    gmonly = False
    opts_dict = pyEnsLib.getopt_parseconfig(opts,optkeys,caller,opts_dict)
    popens = opts_dict['popens']

    # Create a mpi simplecomm object
    if opts_dict['mpi_enable']:
        me=simplecomm.create_comm()
    else:
        me=simplecomm.create_comm(not opts_dict['mpi_enable'])

    # Print out timestamp, input ensemble file and new run directory
    dt=datetime.now()
    verbose = opts_dict['verbose']

    if me.get_rank()==0:
        print '--------pyCECT--------'
        print ' '
        print dt.strftime("%A, %d. %B %Y %I:%M%p")
        print ' '
        print 'Ensemble summary file root = '+opts_dict['sumfile']
        print ' '
        print 'Testcase file directory = '+opts_dict['indir']    
        print ' '
        print ' '

    # Ensure sensible EET value
    if opts_dict['eet'] and opts_dict['numRunFile'] > opts_dict['eet']:
        pyEnsLib.CECT_usage()
        sys.exit(2)

    # only require tag, compset, and res if cesm
    print opts_dict['mach']
    if opts_dict['mach'] == 'gaea':
       # Put output history folder options (histfolder=0,1,2) in a dictionary
       output_dict=dict()
       output_dict[0] = ['1x12m0d']
       output_dict[1] = ['1x1m0d']   
       output_dict[2] = ['1x0m2d','1x0m8d','2x0m1d']
       output_allowed_opts = [0,1,4,8]
       if opts_dict['outputfreq'] not in output_allowed_opts:
          print 'Error: output frequency must be 0 (monthly)\n 1 (annual)\n 4 (4xdaily)\n 8 (8xdaily)'
          sys.exit()
       histfolder = output_dict[opts_dict['histfolder']]

       outputfreq_dict=dict()
       outputfreq_dict[0] = ['month']
       outputfreq_dict[1] = ['year']   
       outputfreq_dict[4] = ['4xdaily']
       outputfreq_dict[8] = ['8xdaily']
                    
       outputfreq = outputfreq_dict[opts_dict['outputfreq']][0]
       use_tiles = opts_dict['usetiles']
       tiles = [int(i) for i in use_tiles.split(',')]
       print tiles   
    #output files
    # output directory tree is root/platform.compiler-debug/prod/repro(-openmp)/duration_pelayout/history/
    #input_dir='/net2/Jessica.Liptak/'
    #
    # Form ensembles, each missing one member; compute RMSZs and global means
    #    #for each variable, we also do max norm also (currently done in pyStats)
    tslice = opts_dict['tslice']
    input_dir = opts_dict['indir']      

    if (os.path.exists(input_dir)):
       # Get the list of files
       in_dirs_temp = os.listdir(input_dir)
       in_files_temp =[]
       # get directories with output for desired time average
       for d in in_dirs_temp:
           for root,subdirs,files in os.walk(os.path.join(input_dir,d)):
               if len(files)>0:
                  # remove input_dir from root, since input_dir is joined
                  # to file name path during read/write
                  new_root = root.replace(input_dir,'')
                  for f in files:
                      for tx in use_tiles:
                          if '.tile' + str(tx) in f: 
                             file_path=os.path.join(new_root,f)
                             #print file_path
                             for h in histfolder:
                                 if h in file_path and outputfreq in file_path:
                                    in_files_temp.append(file_path)
      
    # loop through individual tiles
    for tx in tiles:
        print 'Tile is ' + str(tx)
        in_files=[]
        tstring = 'tile' + str(tx)
        in_files = [f for f in in_files_temp if tstring in f]
        
        num_file = len(in_files)
        print num_file
        if opts_dict['numRunFile'] > num_file:
           print "You requested more numRunFile than it is available at the indir, please change"
           sys.exit()

        in_files.sort()

        if popens:
        #Partition the input file list 
           in_files_list=me.partition(in_files,func=EqualStride(),involved=True)

        else:
        # Random pick non pop files
           in_files_list=pyEnsLib.Random_pickup(in_files,opts_dict)

        ifiles=[]
        for frun_file in in_files_list:
            if frun_file.find(opts_dict['indir']) != -1:
               frun_temp=frun_file
            else:
               frun_temp=opts_dict['indir']+'/'+frun_file
            if os.path.isfile(frun_temp):
               ifiles.append(Nio.open_file(frun_temp,"r"))
            else:
               print "COULD NOT LOCATE FILE " +frun_temp+" EXISTING"
               sys.exit()


if __name__ == "__main__":
    main(sys.argv[1:])
