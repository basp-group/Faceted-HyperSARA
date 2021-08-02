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
gencube = 0
genvis = 0
computenorm = 0
lowerbounds = 1
solve = 1

# kept fixed throughout all the simulations from this folder
imagename = "cygASband_Cube_1024_2048_20"  # cygASband_Cube_512_1024_20
# algoversion = 'cw'
nchannels = 20
ind = 1
Qc = 1
rw = -1
gam = ["1"]  # multiplicative factor affecting the ratio -> '1e-5' order of magnitude
gam_bar = ["3"]
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
        gencube,
        genvis,
        computenorm,
        lowerbounds,
        solve,
    ]

    for g_bar in gam_bar:
        with open(parameter_file_full_path, "r") as csvfile:

            reader = csv.reader(csvfile)

            for job in reader:

                if int(job[1]) * int(job[2]) <= 1:
                    job[0] = "hypersara"
                else:
                    job[0] = "cw"
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

                slurm_command = r"""sbatch --job-name=spa_{16}_h{10}_a{5}_ab{21} --ntasks-per-node={20} \
                -e {28}/{0}_{16}_L={1}_Qx={17}_Qy={18}_Qc={3}_id={2}_overlapx={19}_overlapy={20}_gamma={5}_gammabar={21}_rw={4}_exptype={22}_srf={23}_snr={24}_homotopy={10}.err \
                -o {28}/{0}_{16}_L={1}_Qx={17}_Qy={18}_Qc={3}_id={2}_overlapx={19}_overlapy={20}_gamma={5}_gammabar={21}_rw={4}_exptype={22}_srf={23}_snr={24}_homotopy={10}.out \
                -v --export=ALL,imagename={0},algoversion={16},nchannels={1},ind={2},Qx={17},Qy={18},Qc={3},wintype={7},overlapx={19},overlapy={20},gam={5},nreweights={6},gencube={11},genvis={12},computenorm={13},solve={15},covpath={8},ncdata={9},rw={4},flaghomotopy={10},lowerbounds={14},gambar={23},exptype={22},superresolution={23},isnr={24},logpath={26}\
                run_simulation.slurm""".format(
                    *params,  # +16
                    *job,  # +5
                    ncores,  # 20
                    g_bar,  # 21
                    exp_type,  # 22
                    superresolution_factor,  # 23
                    isnr,  # 24
                    slurm_log_path,  # 25
                    log_path,  # 26
                )

                # print(slurm_command) # Uncomment this line when testing to view the sbatch command

                # Comment the following 3 lines when testing to prevent jobs from being submitted
                exit_status = subprocess.call(slurm_command, shell=True)
                if exit_status == 1:  # Check to make sure the job submitted
                    print("Job {0} failed to submit".format(slurm_command))

                # time.sleep(0.5)

print("Submission complete.")
