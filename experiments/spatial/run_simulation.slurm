#!/bin/bash -l
#### select resources
#SBATCH --time=96:00:00
# #SBATCH --exclusive
#SBATCH --nodes=1
# #SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --account=ec110-guest
#SBATCH --partition=standard
#SBATCH --qos=standard

# # option indicated in the cirrus documentation
# #SBATCH --tasks-per-node=36
# #SBATCH --cpus-per-task=1

# Change to the directory that the job was submitted from
cd $SLURM_SUBMIT_DIR
mkdir results/
mkdir $SLURM_SUBMIT_DIR/pref_$SLURM_JOB_ID
export MATLAB_PREFDIR=$SLURM_SUBMIT_DIR/pref_$SLURM_JOB_ID
mkdir -p /lustre/home/sc004/$USER/$SLURM_JOB_ID

# log_path = results/${imagename}_${exptype}/${algoversion}/logs
# mkdir results/${imagename}_${exptype}/${algoversion}/logs
# $SLURM_SUBMIT_DIR/$log_path

# module purge             # remove all previously loaded modules
module load matlab/R2020b  # version to be checked on cirrus

echo "fprintf('Job id: %s \n', num2str("$SLURM_JOB_ID")); main_simulation('"$imagename"', "$nchannels", "$Qx", "$Qy", "$Qc", '"$algoversion"',  '"$wintype"', "$ncdata", "$ind", ["$overlapy","$overlapx"], "$nreweights", '"$covpath"', "$gam", "$gambar", "$rw", '"$exptype"', "$superresolution", "$isnr", "$genvis", "$computenorm", "$solve", 0, 1, "$flaghomotopy")" > temporary_file_${SLURM_JOB_ID}.m

matlab -nosplash -nodesktop -nodisplay < temporary_file_${SLURM_JOB_ID}.m > ${logpath}/spatial_${imagename}_${algoversion}_${wintype}_ind=${ind}_Qx=${Qx}_Qy=${Qy}_Qc=${Qc}_overlap=${overlapx}_gam=${gam}_gambar=${gambar}_rw=${rw}_homotopy=${flaghomotopy}_exptype=${exptype}_srf=${superresolution}_snr=${isnr}.txt

# Cleanup local work directory
rm temporary_file_${SLURM_JOB_ID}.m
rm -rf $SLURM_SUBMIT_DIR/pref_$SLURM_JOB_ID
rm -rf /lustre/home/sc004/$USER/$SLURM_JOB_ID
