#!/usr/bin/env python3

import csv, subprocess, os
import numpy as np
import time

parameter_file_full_path = "job_spatial.csv"

Qx = 1
Qy = 1
overlapx = 0
overlapy = 0 

# to be activated only for the first run (generating the data), and can systematically deactivated afertwards)
gencube = 0
genvis = 0
computenorm = 1
lowerbounds = 1
solve = 1

# kept fixed throughout all the simulations from this folder
imagename = 'W28_1024'
algoversion = 'sara'
nchannels = 20
Qc = 1
rw = -1
gam = '1e-4'
nreweights = 30
wintype = 'none'
covpath = '../../data/vla_7.95h_dt10s.uvw256.mat'
ncdata = 9 # number of workers in this case (one per dictionary)
flaghomotopy = 1

params = [imagename,algoversion,nchannels,Qc,rw,gam,nreweights,wintype,covpath,ncdata,flaghomotopy,gencube,genvis,computenorm,lowerbounds,solve,Qx,Qy,overlapx,overlapy]

ncores = ncdata + 3

for cubeid in range(1,nchannels+1):

    print("Total number of cpus: {0}".format(ncores))

    slurm_command = r"""sbatch --job-name=spatial_{1} --ntasks-per-node={21} \
    -e {0}_{1}_L={2}_Qx={16}_Qy={17}_Qc={3}_id={20}_overlapx={18}_overlapy={19}_gamma={5}_rw={4}.err \
    -o {0}_{1}_L={2}_Qx={16}_Qy={17}_Qc={4}_id={3}_overlapx={18}_overlapy={19}_gamma={5}_rw={4}.out \
    -v --export=ALL,imagename={0},algoversion={1},nchannels={2},ind={20},Qx={16},Qy={17},Qc={3},wintype={7},overlapx={18},overlapy={19},gam={5},nreweights={6},gencube={11},genvis={12},computenorm={13},solve={15},covpath={8},ncdata={9},rw={4},flaghomotopy={10},lowerbounds={14} \
    run_fhs_mnras.slurm""".format(*params,cubeid,ncores)

    # print(slurm_command) # Uncomment this line when testing to view the sbatch command

    # Comment the following 3 lines when testing to prevent jobs from being submitted
    exit_status = subprocess.call(slurm_command, shell=True)
    if exit_status is 1:  # Check to make sure the job submitted
        print("Job {0} failed to submit".format(slurm_command))

    time.sleep(5)

print("Submission complete.")
