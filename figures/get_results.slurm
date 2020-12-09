#!/bin/bash
#SBATCH --job-name=results
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=4
# #SBATCH --export=ALL
#SBATCH --mail-type=ALL # email notifications about the job: start, end, ...
#SBATCH --mail-user=pierre-antoine.thouvenin@centralelille.fr
#SBATCH --account=ec110-guest
# We use the "standard" partition as we are running on CPU nodes
#SBATCH --partition=standard
# We use the "standard" QoS as our runtime is less than 4 days
#SBATCH --qos=standard
# #SBATCH --cpus-per-task=36
# #SBATCH --exclusive

# Change to the directory that the job was submitted from
cd $SLURM_SUBMIT_DIR
# mkdir pref_$SLURM_JOB_ID

# module purge             # remove all previously loaded modules
module load matlab/R2019a  # version to be checked on cirrus

mkdir $SLURM_SUBMIT_DIR/pref_$SLURM_JOB_ID
# cp -R $SLURM_SUBMIT_DIR/lib/nufft $SLURM_SUBMIT_DIR/lib/nufft_$SLURM_JOB_ID
# export MATLAB_PREFDIR=$SLURM_SUBMIT_DIR/pref_$SLURM_JOB_ID

# Matlab option to deactivate threading
# -singleCompThread

# matlab_wrapper_2019 -nosplash -nodisplay -nodesktop -batch \
# "main_synth_data("$numworkers", "$Nx", "$L", "$p", "$snr", \
# "$gencube", "$gencov", "$genvis")" \
# > logs/hypersara_N=${Nx}_L=${L}_p=${p}_snr=${snr}.log

$HOME/matlab_wrapper_2019_fix \
"get_ground_truth('/lustre/home/shared/sc004/mnras_faceted_corrected'); \
print_images_spectral('/lustre/home/shared/sc004/'); \
print_images_facets('/lustre/home/shared/sc004/'); \
print_images_overlap('/lustre/home/shared/sc004/')" > results.log

rm -r $SLURM_SUBMIT_DIR/pref_$SLURM_JOB_ID
# rm -r $SLURM_SUBMIT_DIR/lib/nufft_$SLURM_JOB_ID
