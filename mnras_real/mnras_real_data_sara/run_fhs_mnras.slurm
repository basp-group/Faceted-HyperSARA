#!/bin/bash -l
#### select resources
#SBATCH --time=96:00:00
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --account=ec110-guest
#SBATCH --partition=standard
#SBATCH --qos=standard

# Change to the directory that the job was submitted from
cd $SLURM_SUBMIT_DIR
# mkdir pref_$SLURM_JOB_ID
mkdir results/

# module purge             # remove all previously loaded modules
module load matlab/R2019a  # version to be checked on cirrus

echo  "main_real_data("$Qx", "$Qy",'"$algoversion"', "$ind", "$nreweights", "$extractdata", "$computenorm", "$solve", "$ch","$gam","$rw","$alpha")"  > temporary_file.m
matlab -nodisplay < temporary_file.m > ind=${ind}_ch=${ch}_Qx=${Qx}_Qy=${Qy}_gam=${gam}_alpha=${alpha}_rw=${rw}.txt
