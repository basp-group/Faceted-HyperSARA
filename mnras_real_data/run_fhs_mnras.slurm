#!/bin/bash
# #SBATCH --job-name=facet-hypersara
#SBATCH --time=96:00:00
# #SBATCH --exclusive
#SBATCH --nodes=1
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
mkdir results/

# module purge             # remove all previously loaded modules
module load matlab/R2019a  # version to be checked on cirrus

mkdir $SLURM_SUBMIT_DIR/pref_$SLURM_JOB_ID
cp -R ../lib/measurement-operator/nufft ../lib/measurement-operator/nufft_$SLURM_JOB_ID
export MATLAB_PREFDIR=$SLURM_SUBMIT_DIR/pref_$SLURM_JOB_ID

# Matlab option to deactivate threading
# -singleCompThread

# matlab_wrapper_2019 -nosplash -nodisplay -nodesktop -batch \
# "main_synth_data("$numworkers", "$Nx", "$L", "$p", "$snr", \
# "$gencube", "$gencov", "$genvis")" \
# > logs/hypersara_N=${Nx}_L=${L}_p=${p}_snr=${snr}.log

$HOME/matlab_wrapper_2019_fix \
"main_real_data("$Qx", "$Qy", "$Qc", \
'"$algoversion"', '"$wintype"', "$ncdata", "$ind", "$overlapsize", "$nreweights", \
"$extractdata", "$computenorm", "$solve", \
"$gam0","$gam","$rw","$alpha")" \
> ${ind}_Qx=${Qx}_Qy=${Qy}_Qc=${Qc}_overlap=${overlapsize}_gam0=${gam0}_gam=${gam}_alpha=${alpha}_rw=${rw}.txt

rm -r $SLURM_SUBMIT_DIR/pref_$SLURM_JOB_ID
rm -r ../lib/measurement-operator/nufft_$SLURM_JOB_ID
