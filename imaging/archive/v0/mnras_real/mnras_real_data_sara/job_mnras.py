#!/usr/bin/env python3

import csv
import os
import subprocess
import time

parameter_file_full_path = "job_params_mnras.csv"

with open(parameter_file_full_path, "r") as csvfile:

    reader = csv.reader(csvfile)

    for job in reader:

        slurm_command = """sbatch --job-name={3}_{8} \
            -e ind={3}_ch={8}_Qx={0}_Qy={1}_gam={9}_alpha={11}_rw={10}.err \
            -o ind={3}_ch={8}_Qx={0}_Qy={1}_gam={9}_alpha={11}_rw={10}.out \
            -v --export=ALL,Qx={0},Qy={1},algoversion={2},ind={3},nreweights={4},extractdata={5},computenorm={6},solve={7},ch={8},gam={9},rw={10},alpha={11} \
            run_fhs_mnras.slurm""".format(
            *job
        )

        # print(slurm_command) # Uncomment this line when testing to view the sbatch command

        # Comment the following 3 lines when testing to prevent jobs from being submitted
        exit_status = subprocess.call(slurm_command, shell=True)
        if exit_status is 1:  # Check to make sure the job submitted
            print("Job {0} failed to submit".format(slurm_command))

        # time.sleep(60)

print("Submission complete.")
