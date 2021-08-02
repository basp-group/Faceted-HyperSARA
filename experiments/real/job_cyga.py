#!/usr/bin/env python3

import os
import pathlib
import subprocess

import numpy as np

Qx = 2
Qy = 2
overlapx = 0.5
overlapy = 0.5
computenorm = 0
lowerbounds = 1
solve = 1

imagename = "cygASband_Cube_1024_2048_20"
algoversion = 'cw'  # solver to be run: 'sara', 'cw' or 'hypersara'
nchannels = 20  # total number of channels
ind = 1  # spectral subcube index
Qc = 1  # number of spectral subcubes
rw = -1  # reweighting step to restart from
gam = ["1"]  # multiplicative factor affecting the l21-norm
gambar = ["1"]  # multiplicative factor affecting the nuclear-norm
nreweights = 5  # total number of reweighting steps
wintype = "triangular"  # type of the apodization window
covpath = "../../data/msSpecs.mat"  # path to coverage file
ncdata = 20  # number of cores for the data fidelity terms
flaghomotopy = 0  # flag to activate homotopy (or not)
superresolution_factor = 2  # superresolution factor


if Qx * Qy <= 1:
    algoversion = "hypersara"
else:
    algoversion = "cw"
ncores = np.minimum(
    Qx * Qy + int(ncdata) + 1, 36
)  # max number of cpus = 36
print("Total number of cpus: {0}".format(ncores))

# path to log and slurm log files
slurmlogpath = os.path.join(
    os.getcwd(),
    "results",
    imagename,
    algoversion,
    "slurm_logs",
)
logpath = os.path.join(
    os.getcwd(), "results", imagename, algoversion, "logs"
)
pathlib.Path(slurmlogpath).mkdir(parents=True, exist_ok=True)
pathlib.Path(logpath).mkdir(parents=True, exist_ok=True)

for g in gam:

    for g_bar in gambar:

        slurm_command = r"""sbatch --job-name=cyga-id={7}- --ntasks-per-node={21} \
        -e {22}/cyga_{16}_L={1}_Qx={2}_Qy={3}_Qc={6}_id={7}_overlapx={4}_overlapy={5}_gamma={13}_gammabar={14}_rw={12}.err \
        -o {22}/cyga_{16}_L={1}_Qx={2}_Qy={3}_Qc={6}_id={7}_overlapx={4}_overlapy={5}_gamma={13}_gammabar={14}_rw={12}.out \
        -v --export=ALL,imagename={0},nchannels={1},Qx={2},Qy={3},overlapx={4},overlapy={5},Qc={6},ind={7},computenorm={8},lowerbounds={9},solve={10},algoversion,={11},rw={12},gam={13},gambar={14},nreweights={15},wintype={16},covpath={17},ncdata={18},flaghomotopy={19},superresolution_factor={20},ncores={21},slurmlogpath={22},logpath={23} \
        run_fhs_mnras.slurm""".format(
            imagename,
            nchannels,
            Qx,
            Qy,
            overlapx,
            overlapy,
            Qc,
            ind,
            computenorm,
            lowerbounds,
            solve,
            algoversion,
            rw,
            gam,
            gambar,
            nreweights,
            wintype,
            covpath,
            ncdata,
            flaghomotopy,
            superresolution_factor,
            ncores,
            slurmlogpath,
            logpath,
        )

        # print(slurm_command) # Uncomment this line when testing to view the
        # sbatch command

        # Comment the following 3 lines when testing to prevent jobs from 
        # being submitted
        exit_status = subprocess.call(slurm_command, shell=True)
        if exit_status == 1:  # Check to make sure the job submitted
            print("Job {0} failed to submit".format(slurm_command))

print("Submission complete.")
