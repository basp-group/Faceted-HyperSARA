#!/bin/bash
#SBATCH --job-name=results
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
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

# Matlab option to deactivate threading
# -singleCompThread

# matlab_wrapper_2019 -nosplash -nodisplay -nodesktop -batch \
# "main_synth_data("$numworkers", "$Nx", "$L", "$p", "$snr", \
# "$gencube", "$gencov", "$genvis")" \
# > logs/hypersara_N=${Nx}_L=${L}_p=${p}_snr=${snr}.log

module load matlab/R2020b  # version to be checked on cirrus

Ny = '512'
Nx = '1024'
nChannels = '20'
gam = '1'
flaghomotopy = '1'
rwtype = 'dirty'
exptype = 'test'
srf = '2'
isnr = '40'
pathToGrountruth = '../../data/cygASband_Cube_512_1024_20.fits'

echo  "addpath('../../figures'); get_sara_cube('./results/cygASband_Cube_512_1024_20_test/',"$nChannels","$Ny","$Nx","$nChannels","$gam","$flaghomotopy", '"$rwtype"', '"$pathToGrountruth"')" > temporary_file_${SLURM_JOB_ID}.m
matlab -nosplash -nodisplay -nodesktop < temporary_file_${SLURM_JOB_ID}.m > cygASband_Cube_512_1024_20_sara_gam=${gam}_homotopy=${flaghomotopy}_rwt=${rwtype}_exptype=${exptype}_srf=${superresolution}_snr=${isnr}.txt

rm temporary_file_${SLURM_JOB_ID}.m
