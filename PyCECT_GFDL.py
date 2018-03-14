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
    s='indir= input_globs= sumfile= tslice= nPC= sigMul= verbose minPCFail= minRunFail= numRunFile= printVarTest popens jsonfile= mpi_enable nbin= minrange= maxrange= outfile= casejson= npick= pepsi_gm test_failure pop_tol= pop_threshold= prn_std_mean lev= eet= json_case= mach= histfolder= usetiles= outputfreq= model= landmask'
    optkeys = s.split()
    try: 
        opts, args = getopt.getopt(argv, "h", optkeys)
    except getopt.GetoptError:
        pyEnsLib.EnsSum_usage()
        sys.exit(2)

    print opts
    print args
    # Put command line options in a dictionary - also set defaults
    opts_dict={}
    
    # Defaults
    opts_dict['input_globs'] = ''
    opts_dict['indir'] = ''
    opts_dict['sumfile'] = 'ens_summary.nc'
    opts_dict['tslice'] = 0
    opts_dict['nPC'] = 5
    opts_dict['sigMul'] = 2
    opts_dict['verbose'] = False
    opts_dict['minPCFail'] = 3
    opts_dict['minRunFail'] = 2
    opts_dict['numRunFile'] = 10
    opts_dict['printVarTest'] = False
    opts_dict['popens'] = False
    opts_dict['jsonfile'] = ''
    opts_dict['mpi_enable'] = False
    opts_dict['nbin'] = 40
    opts_dict['minrange'] = 0.0
    opts_dict['maxrange'] = 4.0
    opts_dict['outfile'] = 'testcase.result'
    opts_dict['casejson'] = ''
    opts_dict['npick'] = 5
    opts_dict['pepsi_gm'] = False
    opts_dict['test_failure'] = True
    opts_dict['pop_tol'] = 3.0
    opts_dict['pop_threshold'] = 0.90
    opts_dict['prn_std_mean'] = False
    opts_dict['lev'] = 0
    opts_dict['eet'] = 20
    opts_dict['json_case'] = ''
    opts_dict['mach'] = 'gaea'
    opts_dict['histfolder'] = 1
    opts_dict['usetiles'] = 1
    opts_dict['outputfreq'] = 0
    opts_dict['model'] = 'gfdl'
    opts_dict['landmask'] = False

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

    # Ensure sensible EET value
    if opts_dict['eet'] and opts_dict['numRunFile'] > opts_dict['eet']:
       'print error: EET must be greater than numRunFile (default = 10)'
       sys.exit(2)

    # only require tag, compset, and res if cesm
    this_sumfile = []
   
    if opts_dict['model'] == 'gfdl':
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
       print opts_dict['outputfreq']
       use_tiles = [opts_dict['usetiles']]
      
       if len(use_tiles) > 1:
          tiles = [int(i) for i in use_tiles.split(',')]
       else:
          tiles = [int(use_tiles[0])]

    # summary ensemble file
       
       for tx in tiles: 
           tstring = 'tile' + str(tx)
           this_sumfile.append(opts_dict['sumfile'].strip('.nc') + '_' + outputfreq + '_' + tstring + '.nc')
       print 'Ensemble summary file root = ', this_sumfile
    else:
       this_sumfile.append(opts_dict['sumfile'])

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
    #output files
    # output directory tree is root/platform.compiler-debug/prod/repro(-openmp)/duration_pelayout/history/
    #input_dir='/net2/Jessica.Liptak/'
    #
    # Form ensembles, each missing one member; compute RMSZs and global means
    #    #for each variable, we also do max norm also (currently done in pyStats)
   
    input_dir = opts_dict['indir']      
    if opts_dict['model'] == 'gfdl':
       if os.path.exists(input_dir):
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
                         for tx in tiles:
                             if '.tile' + str(tx) in f and f[len(f)-3:len(f)] == '.nc': 
                                file_path=os.path.join(new_root,f)
                                #print file_path
                                for h in histfolder:
                                    if h in file_path and outputfreq in file_path:
                                       in_files_temp.append(file_path)
    else:
       tslice = opts_dict['tslice']
      
    # loop through individual tiles
    for tx in tiles:
        ifiles=[]
        in_files=[]
        # Random pick pop files from not_pick_files list
        if opts_dict['casejson']:
           sumfile = opts_dict['sumfile']
           with open(opts_dict['casejson']) as fin:
                result=json.load(fin)
                in_files_first=result['not_pick_files']
                in_files=random.sample(in_files_first,opts_dict['npick'])
                print 'Testcase files:'
                print '\n'.join(in_files)
           
        elif opts_dict['json_case']: 
           json_file=opts_dict['json_case']
           sumfile = opts_dict['sumfile']
           if os.path.exists(json_file):
              fd=open(json_file)
              metainfo=json.load(fd)
              if 'CaseName' in metainfo:
                 casename=metainfo['CaseName']
                 if os.path.exists(opts_dict['indir']):
                    for name in casename: 
                        wildname='*.'+name+'.*'
                        full_glob_str=os.path.join(opts_dict['indir'],wildname)
                        glob_file=glob.glob(full_glob_str)
                        in_files.extend(glob_file)
           else:
              print "Error: "+opts_dict['json_case']+" does not exist"
              sys.exit()

        elif len(opts_dict['input_globs']) > 0: 
           wildname='*'+opts_dict['input_globs']+'*'
           # Open all input files
           if os.path.exists(opts_dict['indir']):
              full_glob_str=os.path.join(opts_dict['indir'],wildname)
              glob_files=glob.glob(full_glob_str)
              in_files.extend(glob_files)
              sumfile = opts_dict['sumfile']
        elif opts_dict['model'] == 'gfdl':
           print 'Tile is ' + str(tx)
           tstring = 'tile' + str(tx)
           in_files = [f for f in in_files_temp if tstring in f]
           sumfile = [f for f in this_sumfile if tstring in f][0]

           tslice = pyEnsLib.get_timestep_number(input_dir,in_files[0])
       
           print 'tslice is ', tslice
        
        num_file = len(in_files)
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
        if popens:    
           # Read in the included var list
           Var2d,Var3d=pyEnsLib.read_jsonlist(opts_dict['jsonfile'],'ESP')
           print ' '
           print 'Z-score tolerance = '+'{:3.2f}'.format(opts_dict['pop_tol'])
           print 'ZPR = '+'{:.2%}'.format(opts_dict['pop_threshold'])
           zmall,n_timeslice=pyEnsLib.compare_raw_score(opts_dict,ifiles,me.get_rank(),Var3d,Var2d)  
           #zmall = np.concatenate((Zscore3d,Zscore2d),axis=0)
           np.set_printoptions(threshold=np.nan)

           if opts_dict['mpi_enable']:
              zmall = pyEnsLib.gather_npArray_pop(zmall,me,(me.get_size(),len(Var3d)+len(Var2d),len(ifiles),opts_dict['nbin'])) 
              if me.get_rank()==0:
                 fout = open(opts_dict['outfile'],"w")
                 for i in range(me.get_size()):
                     for j in zmall[i]:
                         np.savetxt(fout,j,fmt='%-7.2e')
        else:
           # Read all variables from the ensemble summary file
           if opts_dict['landmask']:
              ens_var_name,ens_avg,ens_stddev,ens_rmsz,ens_gm,num_var3d,mu_gm,sigma_gm,loadings_gm,sigma_scores_gm,is_SE_sum,std_gm, tile_mean, land_mean,std_land_mean, mu_land_mean,sigma_land_mean,loadings_land_mean,sigma_scores_land_mean = pyEnsLib.read_ensemble_summary(sumfile, opts_dict) 
           else:    
              ens_var_name,ens_avg,ens_stddev,ens_rmsz,ens_gm,num_var3d,mu_gm,sigma_gm,loadings_gm,sigma_scores_gm,is_SE_sum,std_gm, tile_mean = pyEnsLib.read_ensemble_summary(sumfile, opts_dict)
        if len(ens_rmsz) == 0:
           gmonly = True
        # Add ensemble rmsz and global mean to the dictionary "variables"
        variables={}
        if not gmonly:
           for k,v in ens_rmsz.iteritems():
               pyEnsLib.addvariables(variables,k,'zscoreRange',v)

        for k,v in ens_gm.iteritems():
            pyEnsLib.addvariables(variables,k,'gmRange',v)

        if opts_dict['landmask']:
           for k,v in ens_gm.iteritems():
               pyEnsLib.addvariables(variables,k,'landRange',v)
        # Get 3d variable name list and 2d variable name list seperately
        var_name3d=[]
        var_name2d=[]
        for vcount,v in enumerate(ens_var_name):
            if vcount < num_var3d:
               var_name3d.append(v)
            else:
               var_name2d.append(v)

        n3d = len(var_name3d)
        n2d = len(var_name2d)
       
        # Get ncol and nlev value
        npts3d,npts2d,is_SE=pyEnsLib.get_ncol_nlev(ifiles[0])
        #print 'n3d: ', n3d
        #print 'n2d: ', n2d
        if is_SE ^ is_SE_sum:
           print 'Warning: please note the ensemble summary file is different from the testing files, they use different grids'

        # Compare the new run and the ensemble summary file to get rmsz score
        results={}
        countzscore=np.zeros(len(ifiles),dtype=np.int32)
        countgm=np.zeros(len(ifiles),dtype=np.int32)
        if opts_dict['landmask']:
           countzscore_land=np.zeros(len(ifiles),dtype=np.int32)
           count_landmean=np.zeros(len(ifiles),dtype=np.int32)
        if not gmonly:
           for fcount,fid in enumerate(ifiles): 
               otimeSeries = fid.variables 
               for var_name in ens_var_name: 
                   orig=otimeSeries[var_name]
                   Zscore,has_zscore=pyEnsLib.calculate_raw_score(var_name,orig[opts_dict['tslice']],npts3d,npts2d,ens_avg,ens_stddev,is_SE,opts_dict,0,0,0) 
                   if opts_dict['landmask']:
                      Zscore_land,has_zscore_land=pyEnsLib.calculate_raw_score(var_name,orig[opts_dict['tslice']],npts3d,npts2d,land_mean,ens_stddev,is_SE,opts_dict,0,0,0) 
                   if has_zscore:
                      # Add the new run rmsz zscore to the dictionary "results"
                      pyEnsLib.addresults(results,'zscore',Zscore,var_name,'f'+str(fcount))
           # Evaluate the new run rmsz score if is in the range of the ensemble summary rmsz zscore range
           for fcount,fid in enumerate(ifiles):
               countzscore[fcount]=pyEnsLib.evaluatestatus('zscore','zscoreRange',variables,'ens',results,'f'+str(fcount))

        # Calculate the new run global mean
        #mean3d,mean2d,varlist=pyEnsLib.generate_global_mean_for_summary(ifiles,var_name3d,var_name2d,is_SE,opts_dict['pepsi_gm'],opts_dict)
        
        if not opts_dict['landmask']:
           mean3d,mean2d,var_list = pyEnsLib.generate_global_mean_for_summary(ifiles,var_name3d,var_name2d, is_SE, False,opts_dict,tx)
        elif opts_dict['landmask'] and n3d > 0 and n2d > 0:
           mean3d,mean2d,mean3d_land,var_list = pyEnsLib.generate_global_mean_for_summary(ifiles,var_name3d,var_name2d, is_SE, False,opts_dict,tx)  
           land_means = np.concatenate((mean3d_land,mean2d_land),axis=0)
        elif opts_dict['landmask'] and n3d == 0 and n2d > 0:
           mean3d,mean2d,mean2d_land,var_list = pyEnsLib.generate_global_mean_for_summary(ifiles,var_name3d,var_name2d, is_SE, False,opts_dict,tx)  
           land_means = mean2d_land
        elif opts_dict['landmask'] and n3d > 0 and n2d == 0:
           mean3d,mean2d,mean3d_land,mean2d_land,var_list = pyEnsLib.generate_global_mean_for_summary(ifiles,var_name3d,var_name2d, is_SE, False,opts_dict,tx)  
           land_means = mean3d_land
 
        if n3d == 0 and n2d > 0:
           means = mean2d
        elif n3d > 0 and n2d == 0:
           means = mean3d
        elif n3d and n2d > 0:
           means=np.concatenate((mean3d,mean2d),axis=0)
   
        # Add the new run global mean to the dictionary "results"
        for i in range(means.shape[1]): # n files
            for j in range(means.shape[0]): # n variables
                pyEnsLib.addresults(results,'means',means[j][i],ens_var_name[j],'f'+str(i))
        if opts_dict['landmask']:
           for i in range(land_means.shape[1]): # n files
               for j in range(land_means.shape[0]): # n variables
                   pyEnsLib.addresults(results,'land_means',land_means[j][i],ens_var_name[j],'f'+str(i))
                
        # Evaluate the new run global mean (CESM) or tile mean (GFDL) if it is in the range of the ensemble summary global/tile mean range
        for fcount,fid in enumerate(ifiles):
            countgm[fcount]=pyEnsLib.evaluatestatus('means','gmRange',variables,'gm',results,'f'+str(fcount))
        print countgm[2]

        if opts_dict['landmask']:
           for fcount,fid in enumerate(ifiles):
               count_land[fcount]=pyEnsLib.evaluatestatus('land_means','landRange',variables,'land_mean',results,'f'+str(fcount))

        # Calculate the PCA scores of the new run
        new_scores,var_list,comp_std_gm=pyEnsLib.standardized(means,mu_gm,sigma_gm,loadings_gm,ens_var_name,opts_dict,ens_avg,me)
        print new_scores.shape
        print sigma_scores_gm.shape
        run_index,decision=pyEnsLib.comparePCAscores(ifiles,new_scores,sigma_scores_gm,opts_dict,me)

        if opts_dict['landmask']:
          # Calculate the PCA scores of the new run
          new_scores_land,var_list_land,comp_std_land=pyEnsLib.standardized(land_means,mu_land_mean,sigma_land_mean,loadings_land_mean,ens_var_name,opts_dict,ens_avg,me)
         
          run_index,decision=pyEnsLib.comparePCAscores(ifiles,new_scores_land,sigma_scores_land_mean,opts_dict,me)
      
        # If there is failure, plot out standardized mean and compared standardized mean in box plots
        if opts_dict['prn_std_mean'] and decision == 'FAILED':
            import seaborn as sns
            category={"all_outside99":[],"two_outside99":[],"one_outside99":[],"all_oneside_outside1QR":[]}
            b=list(pyEnsLib.chunk(ens_var_name,10))
            for f,alist in enumerate(b):
                for fc,avar in enumerate(alist):
                    dist_995=np.percentile(std_gm[avar],99.5)
                    dist_75=np.percentile(std_gm[avar],75)
                    dist_25=np.percentile(std_gm[avar],25)
                    dist_05=np.percentile(std_gm[avar],0.5)
                    c=0
                    d=0
                    p=0
                    q=0
                    for i in range(comp_std_gm[f+fc].size):
                        if comp_std_gm[f+fc][i]>dist_995:
                           c=c+1
                        elif comp_std_gm[f+fc][i]<dist_05:
                           d=d+1
                        elif (comp_std_gm[f+fc][i]<dist_995 and comp_std_gm[f+fc][i]>dist_75):
                           p=p+1
                        elif (comp_std_gm[f+fc][i]>dist_05 and comp_std_gm[f+fc][i]<dist_25):
                           q=q+1
                    if c == 3 or d == 3:
                       category["all_outside99"].append((avar,f+fc))
                    elif c == 2 or d == 2:    
                       category["two_outside99"].append((avar,f+fc))
                    elif c == 1 or d == 1:
                       category["one_outside99"].append((avar,f+fc))
                    if p == 3 or q == 3:
                       category["all_oneside_outside1QR"].append((avar,f+fc))
            part_name=opts_dict['indir'].split('/')[-1]
            if not part_name:
                part_name=opts_dict['indir'].split('/')[-2]
            for key in sorted(category):
                list_array=[]
                list_array2=[]
                list_var=[]
                value=category[key]
                print "value len=",key,len(value)
                for each_var in value:
                    list_array.append(std_gm[each_var[0]])
                    list_array2.append(comp_std_gm[each_var[1]])
                    list_var.append(each_var[0])
                if len(value) !=0 :
                    ax=sns.boxplot(data=list_array,whis=[0.5,99.5],fliersize=0.0)
                    sns.stripplot(data=list_array2,jitter=True,color="r")
                    sns.plt.xticks(range(len(list_array)),list_var,fontsize=8,rotation=-45)
                    if decision == 'FAILED':
                       sns.plt.savefig(part_name+"_"+key+"_fail.png")
                    else:
                       sns.plt.savefig(part_name+"_"+key+"_pass.png")
                    sns.plt.clf()
                
            '''
            if len(run_index)>0:
               json_file=opts_dict['json_case']
               if (os.path.exists(json_file)):
                  fd=open(json_file)
                  metainfo=json.load(fd)
                  caseindex=metainfo['CaseIndex']
                  enspath=str(metainfo['EnsPath'][0])
                  #print caseindex
                  if (os.path.exists(enspath)):
                     i=0
                     comp_file=[]
                     search = '\.[0-9]{3}\.'
                     for name in in_files_list:
                        s=re.search(search,name)
                        in_files_index=s.group(0)
                        if in_files_index[1:4] in caseindex:
                           ens_index=str(caseindex[in_files_index[1:4]])
                           wildname='*.'+ens_index+'.*'
                           full_glob_str=os.path.join(enspath,wildname)
                           glob_file=glob.glob(full_glob_str)
                           comp_file.extend(glob_file)
                     print "comp_file=",comp_file                
                     pyEnsLib.plot_variable(in_files_list,comp_file,opts_dict,var_list,run_index,me)
            '''
        # Print out 
        if opts_dict['printVarTest']:
            print '*********************************************** '
            print 'Variable-based testing (for reference only - not used to determine pass/fail)'
            print '*********************************************** '
            for fcount,fid in enumerate(ifiles):
                print ' '
                print 'Run '+str(fcount+1)+":"
                print ' '
                if not gmonly:
                    print '***'+str(countzscore[fcount])," of "+str(len(ens_var_name))+' variables are outside of ensemble RMSZ distribution***'
                    pyEnsLib.printsummary(results,'ens','zscore','zscoreRange',(fcount),variables,'RMSZ')
                    print ' '
                print '***'+str(countgm[fcount])," of "+str(len(ens_var_name))+' variables are outside of ensemble global mean distribution***'
                pyEnsLib.printsummary(results,'gm','means','gmRange',fcount,variables,'global mean')
                print ' '
                print '----------------------------------------------------------------------------'
    if me.get_rank() == 0:
        print ' '
        print "Testing complete."
print ' '
        
if __name__ == "__main__":
    main(sys.argv[1:])
