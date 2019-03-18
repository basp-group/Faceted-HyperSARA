#!/usr/bin/env python
import csv, subprocess

parameter_file_full_path = "/lustre/home/sc004/aa61/HyperSARA_dist/hypersara-sdwt2/full_test/job_params.csv"

with open(parameter_file_full_path, "rb") as csvfile:
    reader = csv.reader(csvfile)
    for job in reader:
	qsub_command = """qsub -N {1}_{2}_{3}_{6} -e /home/sc004/aa61/HyperSARA_dist/hypersara-sdwt2/full_test/{1}_{2}_{3}_{6}.err -o /home/sc004/aa61/HyperSARA_dist/hypersara-sdwt2/full_test/{1}_{2}_{3}_{6}.out -v num_workers={0},tot={1},chunk_width={2},step={3},Qx={4},Qy={5},flag_algo={6} test.pbs""".format(*job)

        print qsub_command # Uncomment this line when testing to view the qsub command

        # Comment the following 3 lines when testing to prevent jobs from being submitted
        exit_status = subprocess.call(qsub_command, shell=True)
        if exit_status is 1:  # Check to make sure the job submitted
            print "Job {0} failed to submit".format(qsub_command)
	print "Done submitting jobs!"
