#!/usr/bin/env python3

import csv, subprocess, os
import numpy as np

parameter_file_full_path = "job_params_mnras.csv"

with open(parameter_file_full_path, "r") as csvfile:

    reader = csv.reader(csvfile)

    for job in reader:

        os.makedirs("results/{0}/logs".format(job[0]), exist_ok=True) # if the directory already exists, do nothing
        ncores = np.minimum(int(job[2])*int(job[3]) + int(job[10]) + 1, 36) # max number of cpus = 36
        # ncores = max(Qx*Qy + n_data_cores + 1, 36) if solver limited to a single node

        slurm_command = """sbatch --job-name={0} --ntasks-per-node={23} \
            -e /lustre/home/sc004/pthouven/Faceted-Hyper-SARA/results/{0}/logs/{0}_L={1}_Qx={2}_Qy={3}_Qc={4}_id={11}_overlap={12}_p={6}_snr={7}.err \
            -o /lustre/home/sc004/pthouven/Faceted-Hyper-SARA/results/{0}/logs/{0}_L={1}_Qx={2}_Qy={3}_Qc={4}_id={11}_overlap={12}_p={6}_snr={7}.out \
            -v --export=ALL,imgname={0},Qx={2},Qy={3},Qc={4},nchannels={1},p={6},snr={7},algoversion={8},wintype={9},ncdata={10},ind={11},overlapsize={12},gencube={13},gencov={14},genvis={15},genundersampledcube={16},computenorm={17},solve={18},cubepath={19},covpath={20},nreweights={21},gamma={22} \
            run_fhs_mnras.slurm""".format(*job,ncores)

        # print(slurm_command) # Uncomment this line when testing to view the sbatch command

        # Comment the following 3 lines when testing to prevent jobs from being submitted
        exit_status = subprocess.call(slurm_command, shell=True)
        if exit_status is 1:  # Check to make sure the job submitted
            print("Job {0} failed to submit".format(slurm_command))

print("Submission complete.")

# name, nchannels, Qx, Qy, Qc, ndatacores, p, isnr, algo_version, window_type, ndataCores, ind, overlapsize, gencube, gencov, genvis, undersampledCube, computeNorm, solve, cubepath, covpath, reweights, gamma
