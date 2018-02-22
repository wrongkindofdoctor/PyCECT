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
def main(argv):


    # Get command line stuff and store in a dictionary
    s = 'tag= compset= esize= tslice= res= sumfile= indir= sumfiledir= mach= verbose jsonfile= mpi_enable maxnorm gmonly popens cumul regx= startMon= endMon= fIndex='
    optkeys = s.split()
    try: 
        opts, args = getopt.getopt(argv, "h", optkeys)
    except getopt.GetoptError:
        pyEnsLib.EnsSum_usage()
        sys.exit(2)

    # Put command line options in a dictionary - also set defaults
    opts_dict={}
    
    # Defaults
    opts_dict['tag'] = ''
    opts_dict['compset'] = ''
    opts_dict['mach'] = 'gaea'
    opts_dict['esize'] = 151
    opts_dict['tslice'] = 0
    opts_dict['res'] = ''
    opts_dict['sumfile'] = 'ens.summary.nc'
    opts_dict['indir'] = './'
    opts_dict['sumfiledir'] = './'
    opts_dict['jsonfile'] = 'exclude_empty.json'
    opts_dict['verbose'] = False
    opts_dict['mpi_enable'] = False
    opts_dict['maxnorm'] = False
    opts_dict['gmonly'] = False
    opts_dict['popens'] = False
    opts_dict['cumul'] = False
    opts_dict['regx'] = 'test'
    opts_dict['startMon'] = 1
    opts_dict['endMon'] = 1
    opts_dict['fIndex'] = 151

    # This creates the dictionary of input arguments 
    opts_dict = pyEnsLib.getopt_parseconfig(opts,optkeys,'ES',opts_dict)

    verbose = opts_dict['verbose']

    st = opts_dict['esize']
    esize = int(st)
    # only require tag, compset, and res if cesm
    if opts_dict['mach'] == 'gaea':
       # Put output time frequency options (tslice=0,1,2) in a dictionary
       output_dict=dict()
       output_dict[0] = ['1x12m0d']
       output_dict[1] = ['1x1m0d','1x1m0d']   
       output_dict[2] = ['1x0m2d','1x0m8d','2x0m1d']

    else:
       if not (opts_dict['tag'] and opts_dict['compset'] and opts_dict['mach'] or opts_dict['res']):
          print 'Please specify --tag, --compset, --mach and --res options'
          sys.exit()
    #output files
    # output directory tree is root/platform.compiler-debug/prod/repro(-openmp)/duration_pelayout/history/
    #input_dir='/net2/Jessica.Liptak/'
    in_files=[]
    #
    # Form ensembles, each missing one member; compute RMSZs and global means
    #    #for each variable, we also do max norm also (currently done in pyStats)
    tslice = opts_dict['tslice']
    time_out=output_dict[tslice]
    input_dir = opts_dict['indir']

    # The var list that will be excluded
    ex_varlist=[]
    inc_varlist=[]

    # Create a mpi simplecomm object
    if opts_dict['mpi_enable']:
        me=simplecomm.create_comm()
    else:
        me=simplecomm.create_comm(not opts_dict['mpi_enable'])
    
    if me.get_rank() == 0:
       print 'Running pyEnsSum!'

    if me.get_rank() ==0 and (verbose == True):
        print opts_dict
        print 'Ensemble size for summary = ', esize

    exclude=False
    if me.get_rank() == 0:
        if opts_dict['jsonfile']:
            inc_varlist=[]
            # Read in the excluded or included var list
            ex_varlist,exclude=pyEnsLib.read_jsonlist(opts_dict['jsonfile'],'ES')
            if exclude == False:
               inc_varlist=ex_varlist
               ex_varlist=[]
            # Read in the included var list
            #inc_varlist=pyEnsLib.read_jsonlist(opts_dict['jsonfile'],'ES')

    # Broadcast the excluded var list to each processor
    # if opts_dict['mpi_enable']:
    #   ex_varlist=me.partition(ex_varlist,func=Duplicate(),involved=True)
    # Broadcast the excluded var list to each processor
    if opts_dict['mpi_enable']:
        exclude=me.partition(exclude,func=Duplicate(),involved=True)
        if exclude:
           ex_varlist=me.partition(ex_varlist,func=Duplicate(),involved=True)
        else:
           inc_varlist=me.partition(inc_varlist,func=Duplicate(),involved=True)

    if (os.path.exists(input_dir)):
       # Get the list of files
       in_dirs_temp = os.listdir(input_dir)
       in_files_temp =[]
       # get directories with output for desired time average
       for d in in_dirs_temp:
           for root,subdirs,files in os.walk(os.path.join(input_dir,d)):
               if len(files)>0:
                  for f in files:
                      file_path=os.path.join(root,f)
                      for t in time_out:
                          if t in file_path and '19790101' in file_path:
                             in_files_temp.append(file_path)

       in_files=sorted(in_files_temp)
       #print 'in_files=',in_files
       # Make sure we have enough
       num_files = len(in_files)
       if me.get_rank()==0 and (verbose == True):
         print 'Number of files in input directory = ', num_files
       if (num_files < esize):
          if me.get_rank()==0 and (verbose == True):
             print 'Number of files in input directory (',num_files,\
             ') is less than specified ensemble size of ', esize
             sys.exit(2)



if __name__ == "__main__":
    main(sys.argv[1:])
