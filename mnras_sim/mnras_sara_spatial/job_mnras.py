#!/usr/bin/env python3

import csv, subprocess, os

parameter_file_full_path = "job_params_mnras.csv"

with open(parameter_file_full_path, "r") as csvfile:

    reader = csv.reader(csvfile)

    for job in reader:

        slurm_command = """sbatch --job-name={5}_{1} --ntasks-per-node={19} \
            -e {5}_{0}_L={1}_ch={21}_gamma={18}_rw={20}.err \
            -o {5}_{0}_L={1}_ch={21}_gamma={18}_rw={20}.out \
            -v --export=ALL,imgname={0},nchannels={1},Qx={2},Qy={3},Qc={4},algoversion={5},wintype={6},ncdata={7},ind={8},overlapsize={9},nreweights={10},gencube={11},gencov={12},genvis={13},genundersampledcube={14},computenorm={15},solve={16},covpath={17},gam={18},ncpus={19},rw={20},ch={21} \
            run_fhs_mnras.slurm""".format(*job)

        # print(slurm_command) # Uncomment this line when testing to view the sbatch command

        # Comment the following 3 lines when testing to prevent jobs from being submitted
        exit_status = subprocess.call(slurm_command, shell=True)
        if exit_status is 1:  # Check to make sure the job submitted
            print("Job {0} failed to submit".format(slurm_command))

print("Submission complete.")
