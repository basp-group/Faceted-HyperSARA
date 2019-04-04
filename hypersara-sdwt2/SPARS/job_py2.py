#!/usr/bin/env python
import csv, subprocess

parameter_file_full_path = "/lustre/home/sc004/aa61/HyperSARA_dist/hypersara-sdwt2/SPARS/job_params2.csv"

with open(parameter_file_full_path, "rb") as csvfile:
    reader = csv.reader(csvfile)
    for job in reader:
	qsub_command = """qsub -N {0}_{1}_{2} -e /home/sc004/aa61/HyperSARA_dist/hypersara-sdwt2/SPARS/{0}_{1}_{2}.err -o /home/sc004/aa61/HyperSARA_dist/hypersara-sdwt2/SPARS/{0}_{1}_{2}.out -v flag_algo={0},tot={1},num_chunk={2} test2.pbs""".format(*job)

        print qsub_command # Uncomment this line when testing to view the qsub command

        # Comment the following 3 lines when testing to prevent jobs from being submitted
        exit_status = subprocess.call(qsub_command, shell=True)
        if exit_status is 1:  # Check to make sure the job submitted
            print "Job {0} failed to submit".format(qsub_command)
	print "Done submitting jobs!"
