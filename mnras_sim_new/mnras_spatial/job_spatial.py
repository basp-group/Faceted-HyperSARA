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
gencube = 0
genvis = 0
computenorm = 0
lowerbounds = 1
solve = 1

# kept fixed throughout all the simulations from this folder
imagename = 'W28_1024'
algoversion = 'cw'
nchannels = 20
ind = 1
Qc = 1
rw = -1
gam = '1'  # multiplicative factor affecting the ratio -> '1e-5' order of magnitude
gam_bar = '1'
nreweights = 30
wintype = 'triangular'
covpath = '../../data/vla_7.95h_dt10s.uvw.mat'
ncdata = 20
flaghomotopy = 1
rw_type = 'dirty' # heuristic, dirty

params = [imagename,algoversion,nchannels,ind,Qc,rw,gam,nreweights,wintype,covpath,ncdata,flaghomotopy,gencube,genvis,computenorm,lowerbounds,solve]

with open(parameter_file_full_path, "r") as csvfile:

    reader = csv.reader(csvfile)

    for job in reader:

        if int(job[0])*int(job[1]) <= 1:
            params[1] = 'hypersara'
        else:
            params[1] = 'cw'
        ncores = np.minimum(int(job[0])*int(job[1]) + int(ncdata) + 1, 36) # max number of cpus = 36
        print("Total number of cpus: {0}".format(ncores))

        slurm_command = r"""sbatch --job-name=spatial_{1} --ntasks-per-node={21} \
        -e {0}_{1}_L={2}_Qx={17}_Qy={18}_Qc={4}_id={3}_overlapx={19}_overlapy={20}_gamma={6}_gammabar={23}_rw={5}_rwt={22}.err \
        -o {0}_{1}_L={2}_Qx={17}_Qy={18}_Qc={4}_id={3}_overlapx={19}_overlapy={20}_gamma={6}_gammabar={23}_rw={5}_rwt={22}.out \
        -v --export=ALL,imagename={0},algoversion={1},nchannels={2},ind={3},Qx={17},Qy={18},Qc={4},wintype={8},overlapx={19},overlapy={20},gam={6},nreweights={7},gencube={12},genvis={13},computenorm={14},solve={16},covpath={9},ncdata={10},rw={5},flaghomotopy={11},lowerbounds={15},rwtype={22},gambar={23} \
        run_fhs_mnras.slurm""".format(*params,*job,ncores,rw_type,gam_bar)

        print(slurm_command) # Uncomment this line when testing to view the sbatch command

        # Comment the following 3 lines when testing to prevent jobs from being submitted
        exit_status = subprocess.call(slurm_command, shell=True)
        if exit_status is 1:  # Check to make sure the job submitted
            print("Job {0} failed to submit".format(slurm_command))

        time.sleep(5)

print("Submission complete.")
