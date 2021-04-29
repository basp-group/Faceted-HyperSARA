#!/usr/bin/env python3

import csv, subprocess, os
import numpy as np
import time

parameter_file_full_path = "job_spatial.csv"

# # to be set in the csv file
# Qx = 2
# Qy = 2
# # compute overlap from Qx, Qy right from the Matlab script?
# overlapx = 0.5
# overlapy = 0.5

# to be activated only for the first run (generating the data), and can systematically deactivated afertwards)
gencube = 0
genvis = 1
computenorm = 1
lowerbounds = 1
solve = 1

# kept fixed throughout all the simulations from this folder
imagename = 'cygASband_Cube_512_1024_20' # cygASband_Cube_1024_2048_20
# algoversion = 'cw'
nchannels = 20
ind = 1
Qc = 1
rw = -1
gam = ['1']  # multiplicative factor affecting the ratio -> '1e-5' order of magnitude
gam_bar = ['1e-1','1','10']
nreweights = 30
wintype = 'triangular'
covpath = '../../data/msSpecs.mat' #'../../data/vla_7.95h_dt10s.uvw.mat'
ncdata = 20
flaghomotopy = 1
rw_type = 'dirty' # 'ground_truth' 'dirty'
exp_type = 'test'
superresolution_factor = 2


for g in gam:

    params = [imagename,nchannels,ind,Qc,rw,g,nreweights,wintype,covpath,ncdata,flaghomotopy,gencube,genvis,computenorm,lowerbounds,solve]

    for g_bar in gam_bar:
        with open(parameter_file_full_path, "r") as csvfile:

            reader = csv.reader(csvfile)

            for job in reader:

                if int(job[1])*int(job[2]) <= 1:
                    job[0] = 'hypersara'
                else:
                    job[0] = 'cw'
                ncores = np.minimum(int(job[1])*int(job[2]) + int(ncdata) + 1, 36) # max number of cpus = 36
                print("Total number of cpus: {0}".format(ncores))

                slurm_command = r"""sbatch --job-name=spatial_{16} --ntasks-per-node={21} \
                -e {0}_{16}_L={1}_Qx={17}_Qy={18}_Qc={3}_id={2}_overlapx={19}_overlapy={20}_gamma={5}_gammabar={23}_rw={4}_rwt={22}_exptype={24}_srf={25}.err \
                -o {0}_{16}_L={1}_Qx={17}_Qy={18}_Qc={3}_id={2}_overlapx={19}_overlapy={20}_gamma={5}_gammabar={23}_rw={4}_rwt={22}_exptype={24}_srf={25}.out \
                -v --export=ALL,imagename={0},algoversion={16},nchannels={1},ind={2},Qx={17},Qy={18},Qc={3},wintype={7},overlapx={19},overlapy={20},gam={5},nreweights={6},gencube={11},genvis={12},computenorm={13},solve={15},covpath={8},ncdata={9},rw={4},flaghomotopy={10},lowerbounds={14},rwtype={22},gambar={23},exptype={24},superresolution={25} \
                run_fhs_mnras.slurm""".format(*params,*job,ncores,rw_type,g_bar,exp_type,superresolution_factor)

                print(slurm_command) # Uncomment this line when testing to view the sbatch command

                # Comment the following 3 lines when testing to prevent jobs from being submitted
                # exit_status = subprocess.call(slurm_command, shell=True)
                # if exit_status is 1:  # Check to make sure the job submitted
                #     print("Job {0} failed to submit".format(slurm_command))

                # time.sleep(1)

print("Submission complete.")
