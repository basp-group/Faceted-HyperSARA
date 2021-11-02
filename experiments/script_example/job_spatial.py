#!/usr/bin/env python3
import json
import os
import pathlib
import subprocess

import numpy as np

if __name__ == "__main__":

    jsonfilename = "setup_matlab.json"

    # Reading elements from json configuration file
    with open(jsonfilename) as f:
        json_config = json.load(f)

    algo_version = json_config[0]["mandatory"]["algo_version"]
    Qx = int(json_config[0]["mandatory"]["Qx"])
    Qy = int(json_config[0]["mandatory"]["Qy"])
    Qc = int(json_config[0]["mandatory"]["Qc"])

    # Specify options, override when necessary (print warning message whenever
    # this happens)
    log_path = os.path.join(
        os.getcwd(),
        "results",
        json_config[0]["general"]["image_name"]
        + "_"
        + json_config[0]["synth"]["exp_type"],
        algo_version,
        "logs",
    )

    slurm_log_path = os.path.join(
        os.getcwd(),
        "results",
        json_config[0]["general"]["image_name"]
        + "_"
        + json_config[0]["synth"]["exp_type"],
        algo_version,
        "slurm_logs",
    )

    pathlib.Path(slurm_log_path).mkdir(parents=True, exist_ok=True)
    pathlib.Path(log_path).mkdir(parents=True, exist_ok=True)

    ncdata = int(json_config[0]["mandatory"]["ncores_data"])
    if int(Qx) * int(Qy) <= 1:
        algo_version = "hs"
    else:
        algo_version = "fhs"
    ncores = np.minimum(
        int(Qx) * int(Qy) + int(ncdata) + 1, 36
    )  # max number of cpus = 36
    print("Total number of cpus: {0}".format(ncores))

    # TODO: add options to override default specifications
    # sara should always run with Qc = nchannels
    # Qx, Qy, overlap not taken into account for hs or sara
    for ind in range(1, Qc + 1):

        logfile = "{0}/spatial_{1}_{2}_{3}_ind={4}_Qx={5}_Qy={6}_Qc={7}_overlap={8}_gam={9}_gambar={10}_rw={11}_exptype={12}_srf={13}_snr={14}.txt".format(
            log_path,
            json_config[0]["general"]["image_name"],
            algo_version,
            json_config[0]["mandatory"]["window_type"],
            ind,
            Qx,
            Qy,
            Qc,
            json_config[0]["mandatory"]["overlap_fraction"],
            json_config[0]["mandatory"]["gamma"],
            json_config[0]["mandatory"]["gamma_bar"],
            json_config[0]["mandatory"]["warmstart_iteration"],
            json_config[0]["synth"]["exp_type"],
            json_config[0]["synth"]["superresolution_factor"],
            json_config[0]["synth"]["isnr"],
        )

        slurm_logname = "{0}/spatial_{1}_{2}_{3}_ind={4}_Qx={5}_Qy={6}_Qc={7}_overlap={8}_gam={9}_gambar={10}_rw={11}_exptype={12}_srf={13}_snr={14}".format(
            slurm_log_path,
            json_config[0]["general"]["image_name"],
            algo_version,
            json_config[0]["mandatory"]["window_type"],
            ind,
            Qx,
            Qy,
            Qc,
            json_config[0]["mandatory"]["overlap_fraction"],
            json_config[0]["mandatory"]["gamma"],
            json_config[0]["mandatory"]["gamma_bar"],
            json_config[0]["mandatory"]["warmstart_iteration"],
            json_config[0]["synth"]["exp_type"],
            json_config[0]["synth"]["superresolution_factor"],
            json_config[0]["synth"]["isnr"],
        )

        slurm_command = r"""sbatch --job-name=spa_{0}_g{1}_gb{2} --ntasks-per-node={3} \
        -e {4}.err \
        -o {4}.out \
        -v --export=ALL,jsonfilename={5},ind={6} \
        run_simulation.slurm""".format(
            algo_version,
            json_config[0]["mandatory"]["gamma"],
            json_config[0]["mandatory"]["gamma_bar"],
            ncores,
            slurm_logname,
            jsonfilename,
            ind,
        )

        # Uncomment this line when testing to view the sbatch command
        print(slurm_command)

        # Comment the following 3 lines when testing to prevent jobs
        # from being submitted
        exit_status = subprocess.call(slurm_command, shell=True)
        if exit_status == 1:
            print("Job {0} failed to submit".format(slurm_command))

    print("Submission complete.")
