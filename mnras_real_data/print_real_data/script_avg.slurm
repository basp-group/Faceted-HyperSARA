#!/bin/bash
#SBATCH --job-name=avg
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
# #SBATCH --cpus-per-task=36

#SBATCH --mail-type=ALL # email notifications about the job: start, end, ...
#SBATCH --mail-user=pierre-antoine.thouvenin@centralelille.fr
#SBATCH --account=ec110-guest
# We use the "standard" partition as we are running on CPU nodes
#SBATCH --partition=standard
# We use the "standard" QoS as our runtime is less than 4 days
#SBATCH --qos=standard

# Change to the directory that the job was submitted from
cd $SLURM_SUBMIT_DIR
# mkdir pref_$SLURM_JOB_ID

# module purge             # remove all previously loaded modules
module load matlab/R2019a  # version to be checked on cirrus

# Matlab option to deactivate threading
# -singleCompThread

# matlab_wrapper_2019 -nosplash -nodisplay -nodesktop -batch \
# "main_synth_data("$numworkers", "$Nx", "$L", "$p", "$snr", \
# "$gencube", "$gencov", "$genvis")" \
# > logs/hypersara_N=${Nx}_L=${L}_p=${p}_snr=${snr}.log

$HOME/matlab_wrapper_2019_fix \
"generate_avg_cube \
> log_avg.txt
