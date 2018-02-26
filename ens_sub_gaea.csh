#! /bin/tcsh -f

#FRE scheduler-options
#PBS -d /lustre/f1/Jessica.Liptak/awg/pycect/ens.%J.stdout
#PBS -o /lustre/f1/Jessica.Liptak/awg/pycect/ens.%J.stdout
#PBS -N pycect_c4
#PBS -l walltime=1:00:00
#PBS -W umask=026
#PBS -S /bin/tcsh
#PBS -r y
#PBS -q batch
#PBS -l partition=c4
#PBS -m a
#PBS -j oe
#PBS -E
#PBS -l nodes=6
#PBS -l qos=windfall
#PBS -A gfdl_f

# Platform environment defaults from /ncrc/home2/fms/local/opt/fre-commands/bronx-12/site/ncrc3/env.defaults.intel
module unload PrgEnv-pgi PrgEnv-intel PrgEnv-gnu PrgEnv-cray
module unload cray-netcdf cray-hdf5 fre
module load PrgEnv-intel/6.0.3
module swap intel intel/17.0.4.196
module load fre/bronx-12
module load cray-hdf5/1.8.16
module load python/2.7.12

set -r htopt = -j1
set -r executable = /lustre/f1/unswept/Jessica.Liptak/awg/warsaw_201710/pycect
alias runCommand /usr/bin/time -p `which aprun` $htopt -n 1 -d 1

$runCommand python pyEnsSum.py --indir /glade/p/tdd/asap/verification/cesm1_3_beta11/sz151-yellowstone-intel --esize 151   --tslice 1  --tag cesm1_3_beta11 --sumfile /glade/scratch/haiyingx/cesm1.3.b11.ne30_ne30.FC5_V6.nc --mpi_enable --mach yellowstone --compset FC5 --res ne30_ne30 --jsonfile included_varlist.json --gmonly
