#!/usr/bin/env python

import csv, subprocess

parameter_file_full_path = "/lustre/home/sc004/aa61/HyperSARA_dist_Full/test88_20_zero_padded/job_params.csv"

with open(parameter_file_full_path, "rb") as csvfile:

    reader = csv.reader(csvfile)

    for job in reader:

        qsub_command = """qsub -q R1005972 -N {0}_{1}_{2} \
        -e /lustre/home/sc004/aa61/HyperSARA_dist_Full/test88_20_zero_padded/{0}_{1}_{2}.err \
        -o /lustre/home/sc004/aa61/HyperSARA_dist_Full/test88_20_zero_padded/{0}_{1}_{2}.out \
        -v alg={0},ver={1},gamma={2},w={3} run_batch.pbs""".format(*job)

        # print qsub_command # Uncomment this line when testing to view the qsub command

        # Comment the following 3 lines when testing to prevent jobs from being submitted
        exit_status = subprocess.call(qsub_command, shell=True)
        if exit_status is 1:  # Check to make sure the job submitted
            print "Job {0} failed to submit".format(qsub_command)

print "Submission complete."
