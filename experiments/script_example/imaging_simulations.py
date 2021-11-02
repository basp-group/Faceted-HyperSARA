#!/usr/bin/env python3
import json
import os
import pathlib
import subprocess

import numpy as np

if __name__ == "__main__":

    # spatial faceting experiment
    jsonfilename = "setup_spatial.json"
    Qc_ = [1]
    Qx_ = [1, 2, 3, 4]
    overlap_ = [0, 0.1, 0.25, 0.4, 0.5]

    # spectral faceting experiment
    # jsonfilename = "setup_spectral.json"
    # Qc_ = [1, 2, 5, 10, 20]
    # Qx_ = [1]
    # overlap_ = [0.]

    slurm_filename = "run_simulation.slurm"

    # Reading elements from json configuration file
    with open(jsonfilename) as f:
        json_config = json.load(f)

    algo_version = json_config[0]["mandatory"]["algo_version"]
    # Qx = int(json_config[0]["mandatory"]["Qx"])
    # Qy = int(json_config[0]["mandatory"]["Qy"])
    # Qc = int(json_config[0]["mandatory"]["Qc"])
    # overlap = json_config[0]["mandatory"]["overlap_fraction"]

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

    if algo_version == "sara":
        ncdata = 9
        Qc_ = [int(json_config[0]["mandatory"]["nchannels"])]
    else:
        ncdata = int(json_config[0]["mandatory"]["ncores_data"])

    for qc in Qc_:

        Qc = qc

        for q in Qx_:

            Qx = q
            Qy = q

            if int(Qx) * int(Qy) <= 1:
                algo_version = "hs"
            else:
                algo_version = "fhs"

            # max number of cpus = 36
            if json_config[0]["synth"]["exp_type"] == "spatial":
                ncores = np.minimum(int(Qx) * int(Qy) + int(ncdata) + 1, 36)
            elif json_config[0]["synth"]["exp_type"] == "spectral":
                ncores = np.minimum(int(ncdata) + 5, 36)
            else:
                raise ValueError(r"Unknown experiment type.")
            print("Total number of cpus: {0}".format(ncores))

            for overlap in overlap_:

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
                        overlap,
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
                        overlap,
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
                    -v --export=ALL,jsonfilename={5},Qx={6},Qy={7},Qc={8},\
                    overlap={9},ind={10} \
                    {11}""".format(
                        algo_version,
                        json_config[0]["mandatory"]["gamma"],
                        json_config[0]["mandatory"]["gamma_bar"],
                        ncores,
                        slurm_logname,
                        jsonfilename,
                        Qx,
                        Qy,
                        Qc,
                        overlap,
                        ind,
                        slurm_filename,
                    )

                    # Uncomment this line when testing to view the sbatch
                    # command
                    # print(slurm_command)

                    # Comment the following 3 lines when testing to prevent
                    # jobs from being submitted
                    exit_status = subprocess.call(slurm_command, shell=True)
                    if exit_status == 1:
                        print("Job {0} failed to submit".format(slurm_command))

    print("Submission complete.")
