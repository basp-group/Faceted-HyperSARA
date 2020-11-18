#!/usr/bin/env python3

import csv, subprocess, os

parameter_file_full_path = "job_params_mnras.csv"

with open(parameter_file_full_path, "r") as csvfile:

    reader = csv.reader(csvfile)

    for job in reader:

        os.makedirs("results/{0}/logs".format(job[0]), exist_ok=True) # if the directory already exists, do nothing

        slurm_command = """sbatch --job-name={0} --ntasks-per-node={1} \
            -e mnras_{0}_L={1}_Qx={2}_Qy={3}_Qc={4}_id={10}_overlap={11}.err \
            -o mnras_{0}_L={1}_Qx={2}_Qy={3}_Qc={4}_id={10}_overlap={11}.out \
            -v --export=ALL,imgname={0},nchannels={1},Qx={2},Qy={3},Qc={4},p={5},snr={6},algoversion={7},wintype={8},ncdata={9},ind={10},overlapsize={11},nreweights={12},gencube={13},gencov={14},genvis={15},genundersampledcube={16},computenorm={17},solve={18},cubepath={19},covpath={20},gam={21} \
            run_fhs_mnras.slurm""".format(*job)

        # print(slurm_command) # Uncomment this line when testing to view the sbatch command

        # Comment the following 3 lines when testing to prevent jobs from being submitted
        exit_status = subprocess.call(slurm_command, shell=True)
        if exit_status is 1:  # Check to make sure the job submitted
            print("Job {0} failed to submit".format(slurm_command))

print("Submission complete.")
