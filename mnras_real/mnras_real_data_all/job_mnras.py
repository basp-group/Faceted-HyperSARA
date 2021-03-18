#!/usr/bin/env python3

import csv, subprocess, os
import time

parameter_file_full_path = "job_test.csv" #"job_params_mnras.csv"

with open(parameter_file_full_path, "r") as csvfile:

    reader = csv.reader(csvfile)

    for job in reader:

        # previous version
        # slurm_command = """sbatch --job-name={6} \
        #     -e {6}_Qx={0}_Qy={1}_Qc={2}_overlap={7}_gam0={12}_gam={13}_alpha={15}_rw={14}.err \
        #     -o {6}_Qx={0}_Qy={1}_Qc={2}_overlap={7}_gam0={12}_gam={13}_alpha={15}_rw={14}.out \
        #     -v --export=ALL,Qx={0},Qy={1},Qc={2},algoversion={3},wintype={4},ncdata={5},ind={6},overlapsize={7},nreweights={8},extractdata={9},computenorm={10},solve={11},gam0={12},gam={13},rw={14},alpha={15},ncpus={16} \
        #     run_fhs_mnras.slurm""".format(*job)

        slurm_command = """sbatch --job-name={6} \
            -e {6}_Qx={0}_Qy={1}_Qc={2}_overlap={7}_gam0={12}_gam={13}_alpha={15}_rw={14}.err \
            -o {6}_Qx={0}_Qy={1}_Qc={2}_overlap={7}_gam0={12}_gam={13}_alpha={15}_rw={14}.out \
            -v --export=ALL,Qx={0},Qy={1},Qc={2},algoversion={3},wintype={4},ncdata={5},ind={6},overlapsize={7},nreweights={8},extractdata={9},computenorm={10},solve={11},gam0={12},gam={13},rw={14},alpha={15},ncpus={16},primal={17},homotopy={18} \
            run_fhs_mnras.slurm""".format(*job)

        # print(slurm_command) # Uncomment this line when testing to view the sbatch command

        # Comment the following 3 lines when testing to prevent jobs from being submitted
        exit_status = subprocess.call(slurm_command, shell=True)
        if exit_status is 1:  # Check to make sure the job submitted
            print("Job {0} failed to submit".format(slurm_command))

        time.sleep(10)

print("Submission complete.")
