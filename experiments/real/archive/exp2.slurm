#!/bin/bash -l

#### select resources
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=18
#SBATCH --cpus-per-task=1
####SBATCH --exclusive
#SBATCH --account=sc004
#SBATCH --partition=standard
#SBATCH --qos=standard


#### redirect error & output files
#### change to current working directory


#### create tmp dir

#### load matlab module (setup environment)

module load matlab/R2020b

#### create a  copy of file 

matlab -nodisplay < fixed_main_input_exp2.m>> exp2_$SLURM_JOB_ID.log

