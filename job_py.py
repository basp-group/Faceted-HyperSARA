#!/usr/bin/env python3

import csv, subprocess

parameter_file_full_path = "job_params.csv"

with open(parameter_file_full_path, "r") as csvfile:

    reader = csv.reader(csvfile)

    for job in reader:

        slurm_command = """sbatch --job-name={0} --ntasks-per-node={1} --time=10:00:00 \
            --mail-type=ALL --mail-user=pierre-antoine.thouvenin@centralelille.fr \
            -e /lustre/home/sc004/pthouven/Faceted-Hyper-SARA/logs/{0}_Qx={2}_Qy={3}_Qc={4}_L={5}_p={6}_snr={7}.err \
            -o /lustre/home/sc004/pthouven/Faceted-Hyper-SARA/logs/{0}_Qx={2}_Qy={3}_Qc={4}_L={5}_p={6}_snr={7}.out \
            -v --export=ALL,imgname={0},Qx={2},Qy={3},Qc={4},nchannels={5},p={6},snr={7},\
            algoversion={8},wintype={9},ncdata={10},ind={11},overlapsize={12},\
            gencube={13},gencov={14},genvis={15},genundersampledcube={16},\
            computenorm={16},solve={17},\
            cubepath={18},covpath={19} \
            run_faceted_hypersara.slurm""".format(*job)

        # print(slurm_command) # Uncomment this line when testing to view the sbatch command

        # # Comment the following 3 lines when testing to prevent jobs from being submitted
        exit_status = subprocess.call(slurm_command, shell=True)
        if exit_status is 1:  # Check to make sure the job submitted
            print("Job {0} failed to submit".format(slurm_command))

print("Submission complete.")
