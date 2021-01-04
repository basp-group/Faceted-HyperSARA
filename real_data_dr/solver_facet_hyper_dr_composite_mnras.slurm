#!/bin/bash 

#### select resources
#SBATCH --time=4:00:00
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36

#SBATCH --account=ec110
#SBATCH --partition=standard
#SBATCH --qos=standard

#### change to current working directory
#SBATCH -D /lustre/home/sc004/mjiang/workspace/matlab/mnras/Dimensionality-reduced-hyper-SARA-sdwt/real_data_dr

#### Redirect error & output files
#SBATCH -e /lustre/home/sc004/mjiang/workspace/matlab/mnras/Dimensionality-reduced-hyper-SARA-sdwt/real_data_dr/solver_facethypersara_dr_err.log
#SBATCH -o /lustre/home/sc004/mjiang/workspace/matlab/mnras/Dimensionality-reduced-hyper-SARA-sdwt/real_data_dr/solver_facethypersara_dr_out.log

#### load matlab module (setup environment)
module load matlab/R2019b

### execute program
mkdir ${SLURM_SUBMIT_DIR}/pref_${SLURM_JOB_ID}
cp -R ../lib/measurement-operator/nufft ../lib/measurement-operator/nufft_${SLURM_JOB_ID}
export MATLAB_PREFDIR=${SLURM_SUBMIT_DIR}/pref_${SLURM_JOB_ID}
# $HOME/matlab_wrapper_2019_fix "func_solver_fouRed_real_data_composite('${datadir}', '${name}', ${Qx}, ${Qy}, ${Qc}, ${gamma0}, ${gamma}, [1:17,19,21:32], ${subInd}, 2, 1, 2, ${fouRed_gamma}, 1, ${adpteps})">./logs/solver_facet_hyper_dr_mnras_perc${fouRed_gamma}_gamma${gamma}_gamma0_${gamma0}_adapteps${adpteps}_Qx${Qx}_Qy${Qy}_Qc${Qc}_run.txt
echo "func_solver_fouRed_real_data_composite('${datadir}', '${name}', ${Qx}, ${Qy}, ${Qc}, ${gamma0}, ${gamma}, [1:17,19,21:32], ${subInd}, 2, 1, 2, ${fouRed_gamma}, 1, ${adpteps})" > dummy_test.m
matlab -nodisplay < dummy_test.m > ./logs/solver_facet_hyper_dr_mnras_perc${fouRed_gamma}_gamma${gamma}_gamma0_${gamma0}_adapteps${adpteps}_Qx${Qx}_Qy${Qy}_Qc${Qc}_run.txt

rm dummy_test.m
rm -rf ${SLURM_SUBMIT_DIR}/pref_${SLURM_JOB_ID}
rm -r ../lib/measurement-operator/nufft_${SLURM_JOB_ID}

module unload matlab/R2019a
