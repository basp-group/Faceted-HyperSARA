#!/usr/bin/env python

import csv, subprocess

parameter_file_full_path = "/lustre/home/sc004/aa61/HyperSARA_dist_Full/real_data33/job_params.csv"

with open(parameter_file_full_path, "rb") as csvfile:

    reader = csv.reader(csvfile)

    for job in reader:

        qsub_command = """qsub -N {0}_{1} \
        -e /lustre/home/sc004/aa61/HyperSARA_dist_Full/real_data33/{0}_{1}.err \
        -o /lustre/home/sc004/aa61/HyperSARA_dist_Full/real_data33/{0}_{1}.out \
        -v ind={0},gamma={1},rw={2} run_batch.pbs""".format(*job)

        # print qsub_command # Uncomment this line when testing to view the qsub command

        # Comment the following 3 lines when testing to prevent jobs from being submitted
        exit_status = subprocess.call(qsub_command, shell=True)
        if exit_status is 1:  # Check to make sure the job submitted
            print "Job {0} failed to submit".format(qsub_command)

print "Submission complete."
