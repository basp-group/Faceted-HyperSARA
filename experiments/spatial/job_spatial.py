#!/usr/bin/env python3

import csv
import os
import pathlib
import subprocess
import time

import numpy as np

parameter_file_full_path = "job_spatial.csv"

# # to be set in the csv file
# Qx = 2
# Qy = 2
# # compute overlap from Qx, Qy right from the Matlab script?
# overlapx = 0.5
# overlapy = 0.5

# to be activated only for the first run (generating the data), and can systematically deactivated afertwards)
genvis = 0
computenorm = 0
solve = 1

# kept fixed throughout all the simulations from this folder
imagename = "cygASband_Cube_1024_2048_20"  # cygASband_Cube_512_1024_20
nchannels = 20
ind = 1
Qc = 1
rw = -1
gam = ["1"]  # multiplicative factor affecting the ratio -> '1e-5' order of magnitude
gam_bar = ["1"]
nreweights = 5
wintype = "triangular"
covpath = "../../data/msSpecs.mat"  #'../../data/vla_7.95h_dt10s.uvw.mat'
ncdata = 20
flaghomotopy = 0
exp_type = "spatial"  # 'test'
superresolution_factor = 2
isnr = 40

# os.path.join(dir_name, base_filename + suffix)

for g in gam:

    params = [
        imagename,
        nchannels,
        ind,
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
    ]

    for g_bar in gam_bar:
        with open(parameter_file_full_path, "r") as csvfile:

            reader = csv.reader(csvfile)

            for job in reader:

                if int(job[1]) * int(job[2]) <= 1:
                    job[0] = "hs"
                else:
                    job[0] = "fhs"
                ncores = np.minimum(
                    int(job[1]) * int(job[2]) + int(ncdata) + 1, 36
                )  # max number of cpus = 36
                print("Total number of cpus: {0}".format(ncores))

                slurm_log_path = os.path.join(
                    os.getcwd(),
                    "results",
                    imagename + "_" + exp_type,
                    job[0],
                    "slurm_logs",
                )

                log_path = os.path.join(
                    os.getcwd(), "results", imagename + "_" + exp_type, job[0], "logs"
                )

                pathlib.Path(slurm_log_path).mkdir(parents=True, exist_ok=True)
                pathlib.Path(log_path).mkdir(parents=True, exist_ok=True)

                slurm_command = r"""sbatch --job-name=spa_{14}_h{10}_a{5}_ab{20} --ntasks-per-node={19} \
                -e {24}/{0}_{14}_L={1}_Qx={15}_Qy={16}_Qc={3}_ind={2}_overlapx={17}_overlapy={18}_gamma={5}_gammabar={20}_rw={4}_exptype={21}_srf={22}_isnr={23}_homotopy={10}.err \
                -o {24}/{0}_{14}_L={1}_Qx={15}_Qy={16}_Qc={3}_ind={2}_overlapx={17}_overlapy={18}_gamma={5}_gammabar={20}_rw={4}_exptype={21}_srf={22}_isnr={23}_homotopy={10}.out \
                -v --export=ALL,imagename={0},algoversion={14},nchannels={1},ind={2},Qx={15},Qy={16},Qc={3},wintype={7},overlapx={17},overlapy={18},gam={5},nreweights={6},gencube={11},genvis={11},computenorm={12},solve={13},covpath={8},ncdata={9},rw={4},flaghomotopy={10},gambar={20},exptype={21},superresolution={22},isnr={23},logpath={25}\
                run_simulation.slurm""".format(
                    *params,  # +14
                    *job,  # +5 hs,1,1,0,0
                    ncores,  # 19
                    g_bar,  # 20
                    exp_type,  # 21
                    superresolution_factor,  # 22
                    isnr,  # 23
                    slurm_log_path,  # 24
                    log_path,  # 25
                )

                # Uncomment this line when testing to view the sbatch command
                # print(slurm_command)

                # Comment the following 3 lines when testing to prevent jobs 
                # from being submitted
                exit_status = subprocess.call(slurm_command, shell=True)
                if exit_status == 1:
                    print("Job {0} failed to submit".format(slurm_command))

print("Submission complete.")
