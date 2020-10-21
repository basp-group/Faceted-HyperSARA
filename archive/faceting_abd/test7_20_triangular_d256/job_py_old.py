#!/usr/bin/env python

import csv, subprocess, os

parameter_file_full_path = "/lustre/home/sc004/aa61/HyperSARA_dist_Full/test7_20_triangular_d256/job_params.csv"

with open(parameter_file_full_path, "r") as csvfile:

    reader = csv.reader(csvfile)

    for job in reader:

        slurm_command = """sbatch --job-name={0} --ntasks-per-node=36 \
        -e /lustre/home/sc004/aa61/HyperSARA_dist_Full/test7_20_triangular_d256/{0}_{1}_{2}.err \
        -o /lustre/home/sc004/aa61/HyperSARA_dist_Full/test7_20_triangular_d256/{0}_{1}_{2}.out \
        -v --export=ALL,alg={0},ver={1},gamma={2} run_batch.slurm""".format(*job)

        # print qsub_command # Uncomment this line when testing to view the qsub command

        # Comment the following 3 lines when testing to prevent jobs from being submitted
        exit_status = subprocess.call(qsub_command, shell=True)
        if exit_status is 1:  # Check to make sure the job submitted
            print "Job {0} failed to submit".format(slurm_command)

print "Submission complete."
