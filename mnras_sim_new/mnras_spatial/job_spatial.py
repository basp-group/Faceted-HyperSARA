#!/usr/bin/env python3

import csv, subprocess, os
import numpy as np
import time
import pathlib
import os

parameter_file_full_path = "job_spatial.csv"

# # to be set in the csv file
# Qx = 2
# Qy = 2
# # compute overlap from Qx, Qy right from the Matlab script?
# overlapx = 0.5
# overlapy = 0.5

# to be activated only for the first run (generating the data), and can systematically deactivated afertwards)
gencube = 0
genvis = 0
computenorm = 0
lowerbounds = 1
solve = 1

# kept fixed throughout all the simulations from this folder
imagename = 'cygASband_Cube_1024_2048_20' # cygASband_Cube_512_1024_20
# algoversion = 'cw'
nchannels = 20
ind = 1
Qc = 1
rw = -1
gam = ['1']  # multiplicative factor affecting the ratio -> '1e-5' order of magnitude
gam_bar = ['3']
nreweights = 5
wintype = 'triangular'
covpath = '../../data/msSpecs.mat' #'../../data/vla_7.95h_dt10s.uvw.mat'
ncdata = 20
flaghomotopy = 0
rw_type = 'heuristic2' # 'ground_truth' 'dirty' 'heuristic' 'heuristic2'
exp_type = 'spatial' # 'test'
superresolution_factor = 2
isnr = 40
updatereg = 0
regtype = 'heuristic2' # 'inv' 'log' 'heuristic' 'heuristic2'
xapprox = 'none' # 'none' 'precond'
noise_transfer = 'none' # 'none' 'precond'
reg_option = 'none' # 'none' 'dirty'

# os.path.join(dir_name, base_filename + suffix)

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
            
                slurm_log_path = os.path.join(os.getcwd(), 'results', imagename + '_' + exp_type, job[0], 'slurm_logs') 

                log_path = os.path.join(os.getcwd(), 'results', imagename + '_' + exp_type, job[0], 'logs') 

                pathlib.Path(slurm_log_path).mkdir(parents=True, exist_ok=True)
                pathlib.Path(log_path).mkdir(parents=True, exist_ok=True)

                slurm_command = r"""sbatch --job-name=spa_{16}_h{10}_reg{27}_rt={29}_xa={30}_nt={31}_a{5}_ab{23}_ropt{33} --ntasks-per-node={21} \
                -e {28}/{0}_{16}_L={1}_Qx={17}_Qy={18}_Qc={3}_id={2}_overlapx={19}_overlapy={20}_gamma={5}_gammabar={23}_rw={4}_rwt={22}_exptype={24}_srf={25}_snr={26}_homotopy={10}_updatereg={27}_regtype={29}_xapprox={30}_nt={31}_ropt={33}.err \
                -o {28}/{0}_{16}_L={1}_Qx={17}_Qy={18}_Qc={3}_id={2}_overlapx={19}_overlapy={20}_gamma={5}_gammabar={23}_rw={4}_rwt={22}_exptype={24}_srf={25}_snr={26}_homotopy={10}_updatereg={27}_regtype={29}_xapprox={30}_nt={31}_ropt={33}.out \
                -v --export=ALL,imagename={0},algoversion={16},nchannels={1},ind={2},Qx={17},Qy={18},Qc={3},wintype={7},overlapx={19},overlapy={20},gam={5},nreweights={6},gencube={11},genvis={12},computenorm={13},solve={15},covpath={8},ncdata={9},rw={4},flaghomotopy={10},lowerbounds={14},rwtype={22},gambar={23},exptype={24},superresolution={25},isnr={26},updatereg={27},regtype={29},xapprox={30},noisetransfer={31},logpath={32},regoption={33} \
                run_fhs_mnras.slurm""".format(*params,*job,ncores,rw_type,g_bar,exp_type,superresolution_factor,isnr,updatereg,slurm_log_path,regtype,xapprox,noise_transfer,log_path,reg_option)

                # print(slurm_command) # Uncomment this line when testing to view the sbatch command

                # Comment the following 3 lines when testing to prevent jobs from being submitted
                exit_status = subprocess.call(slurm_command, shell=True)
                if exit_status is 1:  # Check to make sure the job submitted
                    print("Job {0} failed to submit".format(slurm_command))

                # time.sleep(0.5)

print("Submission complete.")
