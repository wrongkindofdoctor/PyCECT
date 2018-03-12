#! /bin/sh -f
# load platform environment module
source activate pycect_env

export nprocs=2
export indir='/net2/$USER/'
export tag=warsaw_201710
export sumfile='ens.summary.nc'

# 4x daily files 
python pyEnsSumGFDL.py --indir $indir --sumfile $sumfile --esize 50 --outputfreq 4 --tag $tag --mach gaea --jsonfile included_varlist_AM4_4xdaily.json --histfolder 1 --usetiles 1,2,3,4,5,6

python PyCECT_GFDL.py --indir $indir --sumfile $sumfile -outputfreq 4 --mach gaea --jsonfile included_varlist_AM4_4xdaily.json --histfolder 1 --usetiles 1,2,3,4,5,6

# monthly files
python pyEnsSumGFDL.py --indir $indir --sumfile $sumfile --esize 60 --outputfreq 0 --tag $tag  --mach gaea --jsonfile included_varlist_AM4_monthly.json --histfolder 1 --usetiles 1,2,3,4,5,6

python PyCECT_GFDL.py --indir $indir --sumfile $sumfile --eet 20 --outputfreq 0 --mach gaea --jsonfile included_varlist_AM4_monthly.json --histfolder 1 --usetiles 1,2,3,4,5,6
