#!/usr/bin/env python3

import csv, subprocess, os

parameter_file_full_path = "job_params.csv"

with open(parameter_file_full_path, "r") as csvfile:

    reader = csv.reader(csvfile)

    for job in reader:

        os.makedirs("results/{0}/logs".format(job[0]), exist_ok=True) # if the directory already exists, do nothing

        slurm_command = """sbatch --job-name={3} --ntasks-per-node={4} \
        -e /lustre/home/sc004/aa61/HyperSARA_dist_Full/test7_20_triangular_d256/{0}_{1}_{2}.err \
        -o /lustre/home/sc004/aa61/HyperSARA_dist_Full/test7_20_triangular_d256/{0}_{1}_{2}.out \
        -v --export=ALL,alg={0},ver={1},gamma={2} run_batch.slurm""".format(*job)

        # print(slurm_command) # Uncomment this line when testing to view the sbatch command

        # # Comment the following 3 lines when testing to prevent jobs from being submitted
        exit_status = subprocess.call(slurm_command, shell=True)
        if exit_status is 1:  # Check to make sure the job submitted
            print("Job {0} failed to submit".format(slurm_command))

print("Submission complete.")
