#!/usr/bin/env python

import csv, subprocess
import numpy as np

datadir = '/lustre/home/shared/sc004/dr_2b_result_real_data/'
name = 'facethyperSARA'
Qx = 5
Qy = 3
Qc = 15
gamma0 = 0.01
gamma = 5e-6
#chInd = "[1:18,19,21:33]"
subInd = "[1:16]"
fouRed_gamma = 15
adpteps = 0
computelowerbounds = 1 # only needed the first time (saved once it has been computed, loaded whenever the flag is 0)

primal = 0
homotopy = 0

# experiments to be carried out:
# primal = 0, homotopy = 0
# primal = 1, homotopy = 0
# primal = 1, homotopy = 1

# qsub_command = "sbatch -J solver_facet_hyper_dr_mnras_gamma0{5}_gamma{6} --export datadir='{0}',name='{1}',Qx={2},Qy={3},Qc={4},gamma0={5},gamma={6},subInd={7},fouRed_gamma={8},adpteps={9} solver_facet_hyper_dr_composite_mnras.slurm".format(datadir, name, Qx, Qy, Qc, gamma0, gamma, subInd, fouRed_gamma, adpteps)

# mnras current option
# qsub_command = "sbatch -J solver_facet_hyper_dr_mnras_gamma0{5}_gamma{6} --export datadir='{0}',name='{1}',Qx={2},Qy={3},Qc={4},gamma0={5},gamma={6},subInd={7},fouRed_gamma={8},adpteps={9},initoption={10} solver_facet_hyper_dr_composite_mnras.slurm".format(datadir, name, Qx, Qy, Qc, gamma0, gamma, subInd, fouRed_gamma, adpteps,initoption)

# new option
qsub_command = "sbatch -J solver_facet_hyper_dr_mnras_gamma0{5}_gamma{6} --export datadir='{0}',name='{1}',Qx={2},Qy={3},Qc={4},gamma0={5},gamma={6},subInd={7},fouRed_gamma={8},adpteps={9},primal={10},homotopy={11},computelowerbounds={12} solver_facet_hyper_dr_composite_mnras.slurm".format(datadir, name, Qx, Qy, Qc, gamma0, gamma, subInd, fouRed_gamma, adpteps, primal, homotopy, computelowerbounds)

print(qsub_command) # Uncomment this line when testing to view the qsub command

# Comment the following 3 lines when testing to prevent jobs from being submitted
exit_status = subprocess.call(qsub_command, shell=True)
if exit_status is 1:  # Check to make sure the job submitted
    print("Job {0} failed to submit".format(qsub_command))

print("Submission complete.")
