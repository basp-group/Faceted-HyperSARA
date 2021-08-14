#!/bin/bash -l
  
#### select resources

#SBATCH --time=3:00:00
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

module load matlab


matlab -nodisplay < 'script_extract_subcubes.m' > subcubes.log 
