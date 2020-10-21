#!/bin/bash -l

#### select resources

#PBS -l place=scatter:excl

#PBS -l select=1:ncpus=36

#PBS -l walltime=48:00:00

#### Replace [budget code] below with your project code (e.g. t01)      
#PBS -A sc004

# Change to the directory that the job was submitted from
#cd $PBS_O_WORKDIR
cd $SLURM_SUBMIT_DIR

#### load matlab module (setup environment) 
module load matlab/R2018a

cp -R ../nufft ../nufft_${SLURM_JOB_ID}

rm -rf /tmp/dot_matlab_${SLURM_JOB_ID}
mkdir /tmp/dot_matlab_${SLURM_JOB_ID}
cp -pr ~/.matlab /tmp/dot_matlab_${SLURM_JOB_ID}
export MATLAB_PREFDIR=/tmp/dot_matlab_${SLURM_JOB_ID}/.matlab/R2018a
matlab_wrapper -nodisplay -nodisplay -nodesktop -r "main_simulated_data_FINAL_clean("$alg",'"$ver"',"$gamma",'"$PBS_JOBID"')" > ${alg}_${ver}_${gamma}.txt

rm -rf /tmp/dot_matlab_${SLURM_JOB_ID}
rm -r ../nufft_${SLURM_JOB_ID}
