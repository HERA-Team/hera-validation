#!/bin/bash
#PBS -N validation_abscal_model_prep
#PBS -q hera
#PBS -l nodes=1:ppn=1
#PBS -l walltime=48:00:00
#PBS -l vmem=256g
#PBS -j oe
#PBS -o /lustre/aoc/projects/hera/Validation/test-4.0.0/logs/abscal_file_prep.out
#PBS -m be
#PBS -M r.pascua+nrao@berkeley.edu

# This is intended to be run on the heramgr account.
date
source ~/.bashrc
conda activate validation

# Choose one of the days to use as a reference for antenna adjustment.
# This is a little roundabout, but it works.
declare obsfiles=()
for f in $(ls /lustre/aoc/projects/hera/H1C_IDR2/IDR2_2/2458098/*.HH.uvh5); do
    obsfiles+=($f)
done
declare obsfile=${obsfiles[0]}

# Finish setup.
declare testdir=/lustre/aoc/projects/hera/Validation/test-4.0.0
declare simfile=$testdir/data/visibilities/sum.uvh5
declare savedir=$testdir/data/visibilities/abscal_model
declare Nint_per_file=60
declare lst_min=0
declare lst_max=12
declare Trx=100 # Make sure this is the same as used in the file prep!
declare Nchunks=4 # Process in batches of 3 hours each

# Make sure the tools are up to date and switch to the right directory.
cd ~/hera_software/hera-validation
declare branches=$(git branch)
declare remotes=$(git branch -r)
if [[ " ${branches[@]} " =~ "test-4.0.0" ]]; then
    git checkout test-4.0.0
    git pull origin test-4.0.0
elif [[ " ${remotes[@]} " =~ "test-4.0.0" ]]; then
    git checkout --track origin/test-4.0.0
fi
cd test-series/4

echo "python -m tools.prepare_abscal_model_files ${simfile} ${obsfile} ${savedir} ${Nint_per_file} ${lst_min} ${lst_max} --Trx ${Trx} --Nchunks ${Nchunks} --add_noise --clobber --verbose"
python -m tools.prepare_abscal_model_files ${simfile} ${obsfile} ${savedir} ${Nint_per_file} ${lst_min} ${lst_max} --Trx ${Trx} --Nchunks ${Nchunks} --add_noise --clobber --verbose

# Update permissions.
chmod ug+r $savedir/*
date
