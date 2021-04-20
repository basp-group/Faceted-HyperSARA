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

# module purge             # remove all previously loaded modules
module load matlab/R2020b  # version to be checked on cirrus

echo  "main_simulated_data_mnras('"$imagename"',"$nchannels","$Qx","$Qy","$Qc",'"$algoversion"','"$wintype"',"$ncdata","$ind",["$overlapy","$overlapx"],"$nreweights",'"$gencube"',"$genvis","$computenorm","$solve",strcat('"$imagename"','_',num2str("$nchannels")),'"$covpath"',"$gam","$rw","$flaghomotopy", "$lowerbounds", '"$rwtype"')" > temporary_file_${SLURM_JOB_ID}.m
matlab -nodisplay < temporary_file_${SLURM_JOB_ID}.m > spectral_${wintype}_id=${ind}_Qx=${Qx}_Qy=${Qy}_Qc=${Qc}_overlap=${overlapx}_gam=${gam}_rw=${rw}_homotopy=${flaghomotopy}_rwt=${rwtype}.txt

rm temporary_file_${SLURM_JOB_ID}.m
