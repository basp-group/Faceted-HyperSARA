#!/usr/bin/env python3

import csv
import os
import pathlib
import subprocess
import time

import numpy as np

parameter_file_full_path = "job_spatial.csv"

Qx = 1
Qy = 1
overlapx = 0
overlapy = 0

# to be activated only for the first run (generating the data), and can systematically deactivated afertwards)
genvis = 0
computenorm = 1
solve = 1

# kept fixed throughout all the simulations from this folder
imagename = "cygASband_Cube_1024_2048_20"  # cygASband_Cube_512_1024_20
algoversion = "sara"
nchannels = 20
Qc = nchannels
rw = -1
gam = ["1"]
nreweights = 5
wintype = "none"
covpath = "../../data/msSpecs.mat"  # '../../data/vla_7.95h_dt10s.uvw.mat'
ncdata = 9  # number of workers in this case (one per dictionary)
flaghomotopy = 0
exp_type = "spatial"  # 'test'
superresolution_factor = 2
isnr = 40
updatereg = 0

ncores = ncdata + 3

slurm_log_path = os.path.join(
    os.getcwd(), "results", imagename + "_" + exp_type, algoversion, "slurm_logs"
)
log_path = os.path.join(
    os.getcwd(), "results", imagename + "_" + exp_type, algoversion, "logs"
)
pathlib.Path(slurm_log_path).mkdir(parents=True, exist_ok=True)
pathlib.Path(log_path).mkdir(parents=True, exist_ok=True)

for g in gam:

    params = [
        imagename,
        algoversion,
        nchannels,
        Qc,
        rw,
        g,
        nreweights,
        wintype,
        covpath,
        ncdata,
        flaghomotopy,
        genvis,
        computenorm,
        solve,
        Qx,
        Qy,
        overlapx,
        overlapy,
    ]

    for cubeid in range(1, nchannels + 1):

        print("Total number of cpus: {0}".format(ncores))

        slurm_command = r"""sbatch --job-name=spa_{1}_h{10}_a{5} --ntasks-per-node={19} \
        -e {23}/{0}_{1}_L={2}_Qx={14}_Qy={15}_Qc={3}_ind={18}_overlapx={16}_overlapy={17}_gamma={5}_rw={4}_exptype={20}_srf={21}_isnr={22}_homotopy={10}.err \
        -o {23}/{0}_{1}_L={2}_Qx={14}_Qy={15}_Qc={3}_ind={18}_overlapx={16}_overlapy={17}_gamma={5}_rw={4}_exptype={20}_srf={21}_issnr={22}_homotopy={10}.out \
        -v --export=ALL,imagename={0},algoversion={1},nchannels={2},ind={18},Qx={14},Qy={15},Qc={3},wintype={7},overlapx={16},overlapy={17},gam={5},nreweights={6},gencube={11},genvis={11},computenorm={12},solve={13},covpath={8},ncdata={9},rw={4},flaghomotopy={10},gambar=1,exptype={20},superresolution={21},isnr={22},logpath={24}\
        run_simulation.slurm""".format(
            *params,  # +18
            cubeid,  # 18
            ncores,  # 19
            exp_type,  # 20
            superresolution_factor,  # 21
            isnr,  # 22
            slurm_log_path,  # 23
            log_path,  # 24
        )

        # Uncomment this line when testing to view the sbatch command
        # print(slurm_command)

        # Comment the following 3 lines when testing to prevent jobs from being
        # submitted
        exit_status = subprocess.call(slurm_command, shell=True)
        if exit_status == 1:
            print("Job {0} failed to submit".format(slurm_command))

print("Submission complete.")
