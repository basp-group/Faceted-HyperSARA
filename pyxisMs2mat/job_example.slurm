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
module use /work/sc004/shared/software/modules
module load fftw/3.3.8-intel19

module load owlcat
module load kittens
module load pyxis/latest
module load py-casacore
module load tigger
module load meqtrees-cattery
module load meqtrees-timba
export MEQTREES_CATTERY_PATH=/work/sc004/shared/software/meqtrees-cattery/1.7.0
module load astropy
module load mpt
module load wsclean/2.10.1


#### Example 1: 4 MSs spanning 2 frequency bands and associated with 2 configurations of the VLA (data imaged in Thouvenin et al. 2022)
## DATASETS: VLA configuration A at two frequency bands
pyxis MSLOW=CYG-A-5-8M10S.MS MSHIGH=CYG-A-7-8M10S.MS FRISTCH=0 LASTCH=15 FIELDID=2 SRCNAME=CYGA MSTAG=CYGA-ConfigA getdata_ms_concat_bandwidth
## DATASETS: VLA configuration C at two frequency bands
pyxis MSLOW=CYG-C-5-8M10S.MS MSHIGH=CYG-C-7-8M10S.MS FRISTCH=0 LASTCH=15 FIELDID=2 SRCNAME=CYGA MSTAG=CYGA-ConfigC getdata_ms_concat_bandwidth

#### Example 2: single MS 
#pyxis MS=CYG-C-5-8M10S.MS FRISTCH=0 LASTCH=15 FIELDID=2 SRCNAME=CYGA  getdata_ms

#### Example 3: 2 MSs spanning the same frequency band and associated with 2 configurations of the VLA
#pyxis MS=CYG-C-5-8M10S.MS FRISTCH=0 LASTCH=15 FIELDID=2 SRCNAME=CYGA MSTAG=CYGA-ConfigC  getdata_ms
#pyxis MS=CYG-A-5-8M10S.MS FRISTCH=0 LASTCH=15 FIELDID=2 SRCNAME=CYGA MSTAG=CYGA-ConfigA  getdata_ms
