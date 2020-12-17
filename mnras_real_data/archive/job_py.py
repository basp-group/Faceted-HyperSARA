#!/usr/bin/env python

import csv, subprocess

parameter_file_full_path = "/lustre/home/shared/sc004/FacetedHyperSARA/real_data_full_30_all/job_params.csv"

with open(parameter_file_full_path, "rb") as csvfile:

    reader = csv.reader(csvfile)

    for job in reader:

        qsub_command = """qsub -N {0}_{1}_{2}_{4} \
        -e /lustre/home/shared/sc004/FacetedHyperSARA/real_data_full_30_all/{0}_{1}_{2}_{4}.err \
        -o /lustre/home/shared/sc004/FacetedHyperSARA/real_data_full_30_all/{0}_{1}_{2}_{4}.out \
        -v ind={0},gamma0={1},gamma={2},rw={3},alpha={4} run_batch.pbs""".format(*job)

        # print qsub_command # Uncomment this line when testing to view the qsub command

        # Comment the following 3 lines when testing to prevent jobs from being submitted
        exit_status = subprocess.call(qsub_command, shell=True)
        if exit_status is 1:  # Check to make sure the job submitted
            print "Job {0} failed to submit".format(qsub_command)

print "Submission complete."
