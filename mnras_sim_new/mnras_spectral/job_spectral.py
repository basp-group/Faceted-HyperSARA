#!/usr/bin/env python3

import csv, subprocess, os
import numpy as np
import time
import pathlib
import os

parameter_file_full_path = "job_spectral.csv"

# # to be set in the csv file
# Qc = 1, 2, 3, 7, 10, 16 # check in the paper to be sure
# algoversion = 'hypersara' or 'sara'
# sara, 100

# to be activated only for the first run (generating the data), and can systematically deactivated afertwards)
gencube = 0
genvis = 1
computenorm = 1
lowerbounds = 1
solve = 1

# kept fixed throughout all the simulations from this folder
imagename = 'cygASband_Cube_256_512_100'
nchannels = 100
Qx = 1
Qy = 1
overlapx = 0
overlapy = 0
rw = -1
gam = 1 # '1e-3'
gam_bar = '1'
nreweights = 1
wintype = 'none'
covpath = '../../data/vla_7.95h_dt10s.uvw256.mat'
ncdata = 5
flaghomotopy = 0 # 13 elements
rw_type = 'dirty' # heuristic, dirty
exp_type = 'spectral'
superresolution_factor = 2
isnr = 40
updatereg = 0
regtype = 'log'
xapprox = 'precond'
noise_transfer = 'precond'

params = [imagename,nchannels,Qx,Qy,overlapx,overlapy,rw,gam,nreweights,wintype,covpath,ncdata,flaghomotopy,gencube,genvis,computenorm,lowerbounds,solve] # 18 params

with open(parameter_file_full_path, "r") as csvfile:

    reader = csv.reader(csvfile)

    for job in reader:

        ncores = np.minimum(int(ncdata) + 2, 36) # max number of cpus = 36
        print("Total number of cpus: {0}".format(ncores))

        slurm_log_path = os.path.join(os.getcwd(), 'results', imagename + '_' + exp_type, job[0], 'slurm_logs')
        log_path = os.path.join(os.getcwd(), 'results', imagename + '_' + exp_type, job[0], 'logs')
        pathlib.Path(log_path).mkdir(parents=True, exist_ok=True)

        for cubeid in range(1,int(job[1])+1):

            slurm_command = r"""sbatch --job-name=spectral_{18}_{19}_{20} --ntasks-per-node={21} \
            -e {28}/{0}_{18}_L={1}_Qx={2}_Qy={3}_Qc={19}_id={20}_overlapx={4}_overlapy={5}_gamma={7}_gammabar={23}_rw={6}_rwt={22}_exptype={24}_srf={25}_snr={26}_homotopy={12}_updatereg={27}.err \
            -o {28}/{0}_{18}_L={1}_Qx={2}_Qy={3}_Qc={19}_id={20}_overlapx={4}_overlapy={5}_gamma={7}_gammabar={23}_rw={6}_rwt={22}_exptype={24}_srf={25}_snr={26}_homotopy={12}_updatereg={27}.out \
            -v --export=ALL,imagename={0},algoversion={18},nchannels={1},ind={20},Qx={2},Qy={3},Qc={19},wintype={9},overlapx={4},overlapy={5},gam={7},nreweights={8},gencube={13},genvis={14},computenorm={15},solve={17},covpath={10},ncdata={11},rw={6},flaghomotopy={12},lowerbounds={16},rwtype={22},gambar={23},exptype={24},superresolution={25},isnr={26},updatereg={27},regtype={29},xapprox={30},noisetransfer={31},logpath={32} \
            run_fhs_mnras.slurm""".format(*params,*job,cubeid,ncores,rw_type,gam_bar,exp_type,superresolution_factor,isnr,updatereg,slurm_log_path,regtype,xapprox,noise_transfer,log_path)

            # print(slurm_command) # Uncomment this line when testing to view the sbatch command

            # Comment the following 3 lines when testing to prevent jobs from being submitted
            exit_status = subprocess.call(slurm_command, shell=True)
            if exit_status is 1:  # Check to make sure the job submitted
                print("Job {0} failed to submit".format(slurm_command))

            # time.sleep(0.5)

print("Submission complete.")
