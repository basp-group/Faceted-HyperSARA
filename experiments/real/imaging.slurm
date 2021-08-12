#!/bin/bash -l

#### select resources
#SBATCH --time=24:00:00
#SBATCH --nodes=2
#SBATCH --tasks-per-node=36
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --account=ec110-guest
#SBATCH --partition=standard
#SBATCH --qos=standard


#### redirect error & output files
#### change to current working directory

cd $PDIR

#### load matlab module (setup environment)

module load matlab/R2020b

#### create a  copy of file 
RUNFILENAME='./runs/cube-'$SUBCUBE'_job-'$RUNID'.m'
cp fixed_main_input.m    $RUNFILENAME 

echo "subcube_ind="$SUBCUBE  >>$RUNFILENAME
echo "param_reg.gam="$GAM  >>$RUNFILENAME
echo "param_reg.gam_bar="$GAMBAR  >>$RUNFILENAME
echo "param_reg.rw="$RW  >>$RUNFILENAME
echo "param_reg.nReweights="$nRW  >>$RUNFILENAME

echo "%% main script"  >>$RUNFILENAME
echo "main_real_data_dev_ad(imagename,cube_filename,subcube_ind,channels2image,algoversion,ncoredata,param_global,param_reg,input_flags,cirrus) ;" >>$RUNFILENAME

cat $RUNFILENAME
### run Matlab job
##matlab -nodisplay < $RUNFILENAME  >> $LOG

