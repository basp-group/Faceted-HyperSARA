#!/bin/bash -l

#### select resources

#SBATCH --time=1:00:00
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



#### redirect error & output files
#### change to current working directory

#### load matlab module (setup environment)
module unload anaconda
module use /lustre/home/shared/sc004/software/modules
module load fftw/3.3.8-intel19

module load owlcat
module load pyxis/latest
module load py-casacore
module load tigger
module load meqtrees-cattery
module load meqtrees-timba
export MEQTREES_CATTERY_PATH=/lustre/home/shared/sc004/software/meqtrees-cattery/1.7.0
module load pyfits
module load mpt
module load wsclean/2.10.1


## DATASET/CONFIG 1
pyxis MSBW1=CYG-A-5-8M10S.MS MSBW2=CYG-A-7-8M10S.MS CHF=0 CHL=15 FIELDID=2 OUTPUT=CYGA-ConfigA getdata_ms_concat_bandwidth

## DATASET/CONFIG  2
pyxis MSBW1=CYG-C-5-8M10S.MS MSBW2=CYG-C-7-8M10S.MS CHF=0 CHL=15 FIELDID=2 OUTPUT=CYGA-ConfigC getdata_ms_concat_bandwidth

