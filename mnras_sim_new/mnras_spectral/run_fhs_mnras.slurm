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
# mkdir pref_$SLURM_JOB_ID
mkdir results/
mkdir $SLURM_SUBMIT_DIR/pref_$SLURM_JOB_ID
export MATLAB_PREFDIR=$SLURM_SUBMIT_DIR/pref_$SLURM_JOB_ID
mkdir -p /lustre/home/sc004/$USER/$SLURM_JOB_ID

# module purge             # remove all previously loaded modules
module load matlab/R2020b  # version to be checked on cirrus

echo  "main_simulated_data_mnras('"$imagename"',"$nchannels","$Qx","$Qy","$Qc",'"$algoversion"','"$wintype"',"$ncdata","$ind",["$overlapy","$overlapx"],"$nreweights",'"$gencube"',"$genvis","$computenorm","$solve",strcat('"$imagename"','_',num2str("$nchannels")),'"$covpath"',"$gam","$rw","$flaghomotopy", "$lowerbounds", '"$rwtype"', "$gambar", '"$exptype"', "$superresolution","$isnr")" > temporary_file_${SLURM_JOB_ID}.m
matlab -nodisplay < temporary_file_${SLURM_JOB_ID}.m > spectral__${imagename}_${algoversion}_${wintype}_id=${ind}_Qx=${Qx}_Qy=${Qy}_Qc=${Qc}_overlap=${overlapx}_gam=${gam}_gambar=${gambar}_rw=${rw}_homotopy=${flaghomotopy}_rwt=${rwtype}_exptype=${exptype}_srf=${superresolution}_snr=${isnr}.txt

# Cleanup local work directory
rm temporary_file_${SLURM_JOB_ID}.m
rm -r $SLURM_SUBMIT_DIR/pref_$SLURM_JOB_ID
rm -rf /lustre/home/sc004/$USER/$SLURM_JOB_ID
