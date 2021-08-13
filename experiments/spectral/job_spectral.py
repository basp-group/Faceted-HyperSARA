#!/usr/bin/env python3

import csv
import os
import pathlib
import subprocess

import numpy as np

parameter_file_full_path = "job_spectral.csv"

# # to be set in the csv file
# Qc = 1, 2, 3, 7, 10, 16 # check in the paper to be sure
# algoversion = 'hypersara' or 'sara'
# sara, 100

# to be activated only for the first run (generating the data), and can 
# systematically deactivated afertwards)
genvis = 1
computenorm = 1
solve = 1

# kept fixed throughout all the simulations from this folder
imagename = "cygASband_Cube_256_512_100"
nchannels = 100
Qx = 1
Qy = 1
overlapx = 0
overlapy = 0
rw = -1
gam = "1"  # of the order of '1e-3'
gam_bar = "1"
nreweights = 5
wintype = "none"
covpath = "../../data/msSpecs.mat"  # '../../data/vla_7.95h_dt10s.uvw256.mat'
ncdata = 10
flaghomotopy = 0
exp_type = "spectral"
superresolution_factor = 2
isnr = 40

params = [
    imagename,
    nchannels,
    Qx,
    Qy,
    overlapx,
    overlapy,
    rw,
    gam,
    nreweights,
    wintype,
    covpath,
    ncdata,
    flaghomotopy,
    genvis,
    computenorm,
    solve,
]  # 16 params

with open(parameter_file_full_path, "r") as csvfile:

    reader = csv.reader(csvfile)

    for job in reader:

        ncores = np.minimum(int(ncdata) + 5, 36)  # max number of cpus = 36
        print("Total number of cpus: {0}".format(ncores))

        slurm_log_path = os.path.join(
            os.getcwd(), "results", imagename + "_" + exp_type, job[0], "slurm_logs"
        )
        pathlib.Path(slurm_log_path).mkdir(parents=True, exist_ok=True)
        log_path = os.path.join(
            os.getcwd(), "results", imagename + "_" + exp_type, job[0], "logs"
        )
        pathlib.Path(log_path).mkdir(parents=True, exist_ok=True)

        for cubeid in range(1, int(job[1]) + 1):

            slurm_command = r"""sbatch --job-name=spectral_{1}_{17}_{18} --ntasks-per-node={19} \
            -e {24}/{0}_{16}_L={1}_Qx={2}_Qy={3}_Qc={17}_ind={18}_overlapx={4}_overlapy={5}_gamma={7}_gammabar={20}_rw={6}_exptype={21}_srf={22}_isnr={23}_homotopy={12}.err \
            -o {24}/{0}_{16}_L={1}_Qx={2}_Qy={3}_Qc={17}_ind={18}_overlapx={4}_overlapy={5}_gamma={7}_gammabar={20}_rw={6}_exptype={21}_srf={22}_isnr={23}_homotopy={12}.out \
            -v --export=ALL,imagename={0},algoversion={16},nchannels={1},ind={18},Qx={2},Qy={3},Qc={17},wintype={9},overlapx={4},overlapy={5},gam={7},nreweights={8},genvis={13},computenorm={14},solve={15},covpath={10},ncdata={11},rw={6},flaghomotopy={12},gambar={20},exptype={21},superresolution={22},isnr={23},logpath={25}\
            run_fhs_mnras.slurm""".format(
                *params,  # +16
                *job,  # + 2
                cubeid,  # 18
                ncores,  # 19
                gam_bar,  # 20
                exp_type,  # 21
                superresolution_factor,  # 22
                isnr,  # 23
                slurm_log_path,  # 24
                log_path,  # 25
            )

            # Uncomment this line when testing to view the sbatch command
            # print(slurm_command)

            # Comment the following 3 lines when testing to prevent jobs from
            # being submitted
            exit_status = subprocess.call(slurm_command, shell=True)
            if exit_status == 1:
                print("Job {0} failed to submit".format(slurm_command))

print("Submission complete.")
