#!/bin/bash -l

#### select resources

#SBATCH --time=2:00:00
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=36


# Replace [budget code] below with your budget code (e.g. t01)
#SBATCH --account=ec110-guest
# We use the "standard" partition as we are running on CPU nodes
#SBATCH --partition=standard
# We use the "standard" QoS as our runtime is less than 4 days
#SBATCH --qos=standard



#### load matlab module (setup environment)
module load matlab/R2020b # version to be checked on cirrus

#### dir. 
echo "we are here: "$SLURM_SUBMIT_DIR

# dir. of MATLAB temp files
mkdir $SLURM_SUBMIT_DIR/tmp_jobs/
JOB_TMP=$SLURM_SUBMIT_DIR/tmp_jobs/pref_$SLURM_JOB_ID
mkdir $JOB_TMP
echo " tmp files will be saved at: "$JOB_TMP


# dir. of MATLAB output
mkdir $SLURM_SUBMIT_DIR/logs/
echo "logs are in: "$SLURM_SUBMIT_DIR/logs/
MLOG=$SLURM_SUBMIT_DIR/logs/run_$SLURM_JOB_ID.log

# set path of temp files in matlab
export MATLAB_PREFDIR=$JOB_TMP

# run imaging
matlab -nosplash -nodesktop -nodisplay < main_input_imaging.m > $MLOG

# delete tmp dir. 
rm -rf $JOB_TMP