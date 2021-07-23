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
gencube = 0
genvis = 0
computenorm = 1
lowerbounds = 1
solve = 1

# kept fixed throughout all the simulations from this folder
imagename = "cygASband_Cube_1024_2048_20"  # cygASband_Cube_512_1024_20
algoversion = "sara"
nchannels = 20
Qc = nchannels
rw = -1
gam = ["3"]
nreweights = 5
wintype = "none"
covpath = "../../data/msSpecs.mat"  # '../../data/vla_7.95h_dt10s.uvw.mat'
ncdata = 9  # number of workers in this case (one per dictionary)
flaghomotopy = 0
rw_type = "heuristic"  # 'ground_truth' 'dirty' 'heuristic' 'heuristic2'
exp_type = "spatial"  # 'test'
superresolution_factor = 2
isnr = 40
updatereg = 0
regtype = "heuristic"  # 'inv' 'log' 'heuristic' 'heuristic2'
xapprox = "none"  # 'none' 'precond'
noise_transfer = "none"  # 'none' 'precond'
reg_option = "none"  # 'none' 'dirty'

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
        gencube,
        genvis,
        computenorm,
        lowerbounds,
        solve,
        Qx,
        Qy,
        overlapx,
        overlapy,
    ]

    for cubeid in range(1, nchannels + 1):

        print("Total number of cpus: {0}".format(ncores))

        slurm_command = r"""sbatch --job-name=spa_{1}_h{10}_reg{26}_rt={28}_xapprox={29}_nt={30}_a{5}_ropt{32} --ntasks-per-node={21} \
        -e {27}/{0}_{1}_L={2}_Qx={16}_Qy={17}_Qc={3}_id={20}_overlapx={18}_overlapy={19}_gamma={5}_rw={4}_rwt={23}_exptype={22}_srf={24}_snr={25}_homotopy={10}_updatereg={26}_regtype={28}_xapprox={29}_nt={30}_ropt={32}.err \
        -o {27}/{0}_{1}_L={2}_Qx={16}_Qy={17}_Qc={3}_id={20}_overlapx={18}_overlapy={19}_gamma={5}_rw={4}_rwt={23}_exptype={22}_srf={24}_snr={25}_homotopy={10}_updatereg={26}_regtype={28}_xapprox={29}_nt={30}_ropt={32}.out \
        -v --export=ALL,imagename={0},algoversion={1},nchannels={2},ind={20},Qx={16},Qy={17},Qc={3},wintype={7},overlapx={18},overlapy={19},gam={5},nreweights={6},gencube={11},genvis={12},computenorm={13},solve={15},covpath={8},ncdata={9},rw={4},flaghomotopy={10},lowerbounds={14},gambar=1,exptype={22},rwtype={23},superresolution={24},isnr={25},updatereg={26},regtype={28},xapprox={29},noisetransfer={30},logpath={31},regoption={32} \
        run_fhs_mnras.slurm""".format(
            *params,
            cubeid,
            ncores,
            exp_type,
            rw_type,
            superresolution_factor,
            isnr,
            updatereg,
            slurm_log_path,
            regtype,
            xapprox,
            noise_transfer,
            log_path,
            reg_option
        )

        # print(slurm_command) # Uncomment this line when testing to view the sbatch command

        # # Comment the following 3 lines when testing to prevent jobs from being submitted
        exit_status = subprocess.call(slurm_command, shell=True)
        if exit_status is 1:  # Check to make sure the job submitted
            print("Job {0} failed to submit".format(slurm_command))

        # time.sleep(0.5)

print("Submission complete.")
