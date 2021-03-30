#!/usr/bin/env python3

import csv, subprocess, os
import numpy as np
import time

parameter_file_full_path = "job_spatial.csv"

# # to be set in the csv file
# Qx = 2
# Qy = 2
# # compute overlap from Qx, Qy right from the Matlab script?
# overlapx = 512
# overlapy = 512

# to be activated only for the first run (generating the data), and can systematically deactivated afertwards)
gencube = 1
genvis = 1
computenorm = 1
lowerbounds = 1
solve = 1

# kept fixed throughout all the simulations from this folder
imagename = 'W28_1024'
algoversion = 'cw'
nchannels = 20
ind = 1
Qc = 1
rw = 1
gam = 1e-5
nreweights = 30
wintype = 'triangular'
covpath = '../../data/vla_7.95h_dt10s.uvw.mat'
ncdata = 20
flaghomotopy = 1

params = [imagename,algoversion,nchannels,ind,Qc,rw,gam,nreweights,wintype,covpath,ncdata,flaghomotopy,gencube,genvis,computenorm,lowerbounds,solve]

with open(parameter_file_full_path, "r") as csvfile:

    reader = csv.reader(csvfile)

    for job in reader:

        ncores = np.minimum(int(job[0])*int(job[1]) + int(ncdata) + 1, 36) # max number of cpus = 36
        print("Total number of cpus: {0}".format(ncores))

        slurm_command = r"""sbatch --job-name=spatial_{1} --cpus-per-task={21} \
        -e {0}_{1}_L={2}_Qx={17}_Qy={18}_Qc={4}_id={3}_overlapx={19}_overlapy={20}_gamma={6}_rw={5}.err \
        -o {0}_{1}_L={2}_Qx={17}_Qy={18}_Qc={4}_id={3}_overlapx={19}_overlapy={20}_gamma={6}_rw={5}.out \
        -v --export=ALL,imagename={0},algoversion={1},nchannels={2},ind={3},Qx={17},Qy={18},Qc={4},wintype={8},overlapx={19},overlapy={20},gam={10},nreweights={7},gencube={12},genvis={13},computenorm={14},solve={16},covpath={9},ncdata={10},rw={5},flaghomotopy={11},lowerbounds={15} \
        run_fhs_mnras.slurm""".format(*params,*job,ncores)

        print(slurm_command) # Uncomment this line when testing to view the sbatch command

        # Comment the following 3 lines when testing to prevent jobs from being submitted
        exit_status = subprocess.call(slurm_command, shell=True)
        if exit_status is 1:  # Check to make sure the job submitted
            print("Job {0} failed to submit".format(slurm_command))

        time.sleep(5)

print("Submission complete.")
