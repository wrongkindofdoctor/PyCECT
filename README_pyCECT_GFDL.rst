README.pyCECT_GFDL
=========================
CESM Ensemble Consistency Test modfied for compatibility with NOAA
GFDL model: pyGECT
------------------------------
The modified CESM-ECT package is used to compare the results of a set of new
GFDL simulations against the accepted ensemble.  
For details on CAM-ECT, see: 

A.H. Baker, D.M. Hammerling, M.N. Levy, H. Xu, J.M. Dennis, B.E. Eaton, J. Edwards, 
C. Hannay, S.A. Mickelson, R.B. Neale, D. Nychka, J. Shollenberger, J. Tribbia, 
M. Vertenstein, and D. Williamson,"A new ensemble-based consistency test for the 
community earth system model."Geoscientific Model Development, 8, pp. 2829-2840, 
doi:10.5194/gmd-8-2829-2015, 2015.

For details on POP-ECT, see:  

A.H.Baker, Y.Hu, D.M. Hammerling, H. Xu, Y.Tseng, 
X. Huang, and F.O. Byan "Evaluating Consistency in the Ocean Model Component of 
the Community Earth System Model", submitted to Geoscientific Model Development 
Discussions, Dec. 2015

:AUTHORS: Haiying Xu, Allison Baker, J. L.
:VERSION: 3.0.0
:COPYRIGHT: See the document entitled LICENSE.txt

pyCECT_GFDL, pyEnsSumGFDL:
:AUTHOR: Jessica Liptak
:CONTACT email: jessica.liptak@noaa.gov

This package includes:  
----------------------
     	pyEnsSumGFDL.py             
                            A script that generates an ensemble summary file 
     		            from a collection of files. This script has been modified
                            from the pyEnsSum.py script to read GFDL history files generated
                            by simulations run with fre. The history files contain output 
                            from individual tiles on the cubed sphere grid. Output variables are specified
                            using namelists, and may vary depending on the frequency of output 
                            (e.g., monthly mean files may contain different variables than 6-hourly or 
                             3-hourly files).

        pyEnsLib.py     
                            Library python script used by pyEnsSum.py. This script has been modfied 
                            to process either 2D or 3D variables (i.e., it doesn't assume that history files 
                            contain both 2D and 3D variables) on cubed sphere tiles, and has additional routines
                            for computing weighted averages of GFDL output. 


        ens_excluded_varlist.json
                            The variable list that will excluded from
                            reading and processing

        included_varlist.json
                            The variable list that will be included for
                            reading and processing

	exclude_empty.json
	                   An empty exclude variable list, useful for 
			   determining from scratch which variables to exclude
 

        pyCECT_GFDL.py
                            A script which compares the new AM results to the 
                            accepted ensemble and issues an overall "pass" or "fail".
                            This script is a modfied version of pyCECT.py that reads in 
                            GFDL history files.

Software requirements:
--------------------------------------------
To use the package on GFDL workstations, I recommend creating an anaconda
environment and installing PyNIO, PyNGL, SciPy, mpi4py, and ASAPTools packages,
as there are several required libraries that may not be compatible with 
your system defaults. Anaconda should be available on gaea and 
workstations running RHEL7. Workstations running RHEL6 or earlier 
will require a local install of the 'miniconda' package, which is 
a stripped-down version of anaconda that has a greater likelihood of 
fitting in the user space constraints. 

MINICONDA installation instructions:
--------------------------------------------
To install miniconda, go to https://conda.io/miniconda.html
and download the Python2.7 installation shell script for 
the Linux 64-bit architecture.

cd to the directory with the downloaded script and
give yourself permission to run the executable: 
> chmod u+x Miniconda2-latest-Linux-x86_64.sh

Run the installation script, accept the terms and conditions,
and follow instructions to select a build directory
> sh ./Miniconda2-latest-Linux-x86_64.sh

when the installation is complete, either alias the 'conda' 
command in your .cshrc file: 
alias conda  '/<miniconda_install_directory>/bin/conda'
or create a soft link:
> ln -s /<miniconda_install_directory>/bin/conda conda

Installing required Python packages
--------------------------------------------
After verifying that anaconda or miniconda is installed on your machine,
switch to a bash shell:
: exec sh

then, create a new environment: 
sh-4.2$ conda create -n <my_environment>

activate the environment:
sh-4.2$: source activate <my_environment>

install the following packages from conda-forge:
sh-4.2$: conda install --channel conda-forge pynio
sh-4.2$: conda install --channel conda-forge pyngl
sh-4.2$: conda install --channel conda-forge scipy
sh-4.2$: conda install --channel conda-forge mpi4py
sh-4.2$: pip install --user ASAPTools

OR: recreate the environment with the file pycect_package-list.txt
sh-4.2$: conda create -n <my_environment> --file  pycect_package-list.txt

then install ASAPTools:
sh-4.2$: pip install --user ASAPTools

Notes and examples:
--------------------------------------------
(1) Activate the anaconda environment with the 
required packages:
sh-4.2$: source activate <my_environment>

(2) Create a summary file with pyEnsSumGFDL.py (example doesn't use mpi)

export indir='/net2/$USER/'
export tag=warsaw_201710
export sumfile='ens.summary.nc'

python pyEnsSumGFDL.py --indir $indir --sumfile $sumfile --esize 50 --outputfreq 4 --tag $tag  --mach gaea --jsonfile included_varlist_AM4_4xdaily.json --histfolder 1 --usetiles 1,2,3,4,5,6

Options used in this example are:
--indir: name of root directory with history files. The script has been modified to
loop through the default subdirectory structure created by FRE.

--sumfile: root name of summary file(s) that will be generated (the tiles specified by --usetiles will be appended)

--esize: Number of ensemble members to use in analysis, must be <= the number of history files analyzed; default is 50

--outputfreq: option for history file output frequency may be monthly (0, default), annual (1), 4xdaily (4), 8xdaily (8) 

--tag: optional name of model code tag

--mach: machine used to run simulations that generated the history files; default is 'gfdl'

--jsonfile: file of variables to include or exclude

--histfolder: index of history folder(s) to process [N segments x number of months x number of days]
         
        0 = '1x12m0d' [1 segmentt of 12-month 0-day output aka 1 year of output] 
        1 = '1x1m0d'  [1 segment of 1-month 0-day output aka 1 month of output]
        2 = '1x0m2d','1x0m8d','2x0m1d' [1 or 2 segments of 0-month, multi-day output] 

--tiles: tiles to generate summary stats for; default analyzes only tile 1 history files

Other options and their descriptions are listed in PyEnsLib.py

(3) run the pyCECT_GFDL script

   (A) Options for all GFDL-ECT approaches:

     Required:

         To specify the summary file generated by pyEnsSum.py:
	    --sumfile  ens.summary.nc

     	 To specify the directory path that contains the run(s) to be evaluated:
	    --indir /home/$USER/output_directory
         By default, the script assumes a FRE-generated output structure of 
         <platform>-<target>/<Nxnmxnd>_<NprocessorsxNa>/history/: 
         e.g., the paths to the files in 'indir' are:

         ncrc4.pgi-prod/1x1m0d_432x1a/history/
         ncrc4.intel17-repro-openmp-avx/1x1m0d_432x1a/history/
         ncrc3.cce-debug/1x1m0d_432x1a/history/ ...
         

    Optional:
	 Verbose information:
	     --verbose

    (B) AM-ECT and UF-AM-ECT specific options (and summary file generated by pyEnsSumGFDL.py)

        Note that AM-ECT is the default test.

    Note that the parameters setting the pass/fail criteria are all set by 
    default (ie. sigMul, minPCFail, minRunFail, numRunFile, and nPC).  But 
    if the specified indir contains more files than the number (num) specified by 
    "--numRunFile <num>"  (default= 3), then <num> files will be chosen at random 
    from that directory. Ensemble Exhaustive Test (EET) is specified by --eet <num>. 
    This tool computes the failure rate of <num> tests taken <numRunFile> at a time.
    Therefore, when specifying --eet <num>, <num> must be greater than or equal to
    <numRunFile>. 

    To enable printing of extra variable information:
       --printVarTest

    By default, AM-ECT looks at monthly averages

       --outputfreq 0
    For other intervals (e.g, 4xdaily), AM-ECT analyzes the last time step in each history file
    unless a different value is specified with --tslice. 
 
Examples:
--------------------------------------
    export nprocs=2
    export indir='/net2/$USER/'
    export tag=warsaw_201710

    Example using the default settings:
    
    python pyGECT_GFDL.py --sumfile  AM4p0.ens.summary.nc --indir  $indir

    Example using EET:

    python pyCECT_GFDL.py --sumfile  AM4p0.ens.summary.nc --indir $indir --eet 10
         

    Example run in parallel:
         
    mpirun -n $nprocs python pyCECT_GFDL.py --indir $indir --esize 41--outputfreq 4 --tag $tag --sumfile am4p0.summary.nc --mpi_enable --mach gaea --jsonfile included_varlist_AM4_monthly.json --histfolder 1 --verbose


