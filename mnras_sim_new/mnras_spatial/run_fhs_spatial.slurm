#!/bin/bash
# #SBATCH --job-name=facet-hypersara
#SBATCH --time=96:00:00
# #SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --account=ec110-guest
#SBATCH --partition=standard
# #SBATCH --cpus-per-task=36

# Change to the directory that the job was submitted from
cd $SLURM_SUBMIT_DIR
# mkdir pref_$SLURM_JOB_ID
mkdir results/

# module purge             # remove all previously loaded modules
module load matlab/R2020b  # version to be checked on cirrus

# mkdir $SLURM_SUBMIT_DIR/pref_$SLURM_JOB_ID
# cp -R ../lib/measurement-operator/nufft ../lib/measurement-operator/nufft_$SLURM_JOB_ID
# export MATLAB_PREFDIR=$SLURM_SUBMIT_DIR/pref_$SLURM_JOB_ID

# Matlab option to deactivate threading
# -singleCompThread

# matlab_wrapper_2019 -nosplash -nodisplay -nodesktop -batch \
# "main_synth_data("$numworkers", "$Nx", "$L", "$p", "$snr", \
# "$gencube", "$gencov", "$genvis")" \
# > logs/hypersara_N=${Nx}_L=${L}_p=${p}_snr=${snr}.log

echo  "main_simulated_data_mnras('"$imagename"',"$nchannels","$Qx","$Qy","$Qc",'"$algoversion"','"$wintype"',"$ncdata","$ind",["$overlapy","$overlapx"],"$nreweights",'"$gencube"',"$genvis","$computenorm","$solve",strcat('"$imagename"',num2str("$nchannels")),'"$covpath"',"$gam","$rw","$flag_homotopy", "$lowerbounds")" > temporary_file_${SLURM_JOB_ID}.m

matlab -nodisplay < temporary_file_${SLURM_JOB_ID}.m > spatial_${wintype}_id=${ind}_Qx=${Qx}_Qy=${Qy}_Qc=${Qc}_overlap=${overlapx}_gam=${gam}_rw=${rw}_homotopy=${homotopy}.txt

# rm -r $SLURM_SUBMIT_DIR/pref_$SLURM_JOB_ID
# rm -r ../lib/measurement-operator/nufft_$SLURM_JOB_ID
