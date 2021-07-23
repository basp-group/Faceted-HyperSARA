#!/usr/bin/env python

import csv
import subprocess

parameter_file_full_path = "job_params.csv"

with open(parameter_file_full_path, "rb") as csvfile:

    reader = csv.reader(csvfile)

    for job in reader:

        qsub_command = """qsub -N {1}_xy{2}_c{4} -l select=1:ncpus={7} \
        -e /lustre/home/sc004/pthouven/Dimensionality-reduced-hyper-SARA-sdwt/logs/hyperSARA_{1}_Qx={2}_Qy={3}_Qc={4}.err \
        -o /lustre/home/sc004/pthouven/Dimensionality-reduced-hyper-SARA-sdwt/logs/hyperSARA_{1}_Qx={2}_Qy={3}_Qc={4}.out \
        -v flag_algo={0},parallel_version={1},Qx={2},Qy={3},Qc={4},tot={5},num_chunk={6} run_batch.pbs""".format(
            *job
        )

        # print qsub_command # Uncomment this line when testing to view the qsub command

        # Comment the following 3 lines when testing to prevent jobs from being submitted
        exit_status = subprocess.call(qsub_command, shell=True)
        if exit_status is 1:  # Check to make sure the job submitted
            print "Job {0} failed to submit".format(qsub_command)

print "Submission complete."
