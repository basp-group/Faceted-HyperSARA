#!/usr/bin/env python3

import csv, subprocess, os
import time

parameter_file_full_path = "job_test.csv" #"job_params_mnras.csv"

with open(parameter_file_full_path, "r") as csvfile:

    reader = csv.reader(csvfile)

    for job in reader:

        slurm_command = """sbatch --job-name={1} \
            -e {0}_{1}_L={2}_Qx={4}_Qy={5}_Qc={6}_id={3}_overlap={8}_gamma={10}_rw={11}.err \
            -o {0}_{1}_L={2}_Qx={4}_Qy={5}_Qc={6}_id={3}_overlap={8}_gamma={10}_rw={11}.out \
            -v --export=ALL,imgname={0},algoversion={1},nchannels={2},ind={3},Qx={4},Qy={5},Qc={6},wintype={7},overlapx={8},overlapy={9},gam={10},
            nreweights={11},gencube={12},gencov={13},genvis={14},computenorm={15},solve={16},covpath={17},ncdata={18},rw={19},flaghomotopy={20},lowerbounds={21},ncpus={22} \
            run_fhs_mnras.slurm""".format(*job)

        # print(slurm_command) # Uncomment this line when testing to view the sbatch command

        # Comment the following 3 lines when testing to prevent jobs from being submitted
        exit_status = subprocess.call(slurm_command, shell=True)
        if exit_status is 1:  # Check to make sure the job submitted
            print("Job {0} failed to submit".format(slurm_command))

        time.sleep(10)

print("Submission complete.")
