#!/usr/bin/env python
import argparse
import os
import subprocess
from pathlib import Path


def main():

    # Creating auxiliary folders and setting paths
    project_path_run = os.path.dirname(os.path.realpath(__file__))
    log_path = "%s/logs" % project_path_run
    run_path = "%s/runs" % project_path_run

    Path(log_path).mkdir(parents=True, exist_ok=True)
    Path(run_path).mkdir(parents=True, exist_ok=True)

    slurm_file = "%s/imaging.slurm" % project_path_run
    print(slurm_file)

    parser = argparse.ArgumentParser()
    print("Current experiment Dir %s" % project_path_run)

    parser.add_argument(
        "-gam",
        "--gamma",
        help="Multiplicative constant to adjust the value of the heuristic low-rankness regularization parameter.",
        default=1.0,
        type=float,
    )

    parser.add_argument(
        "-gambar",
        "--gammabar",
        help="Multiplicative constant to adjust the value of the heuristic average joint-sparsity regularization parameter.",
        default=1.0,
        type=float,
    )

    parser.add_argument(
        "-img",
        "--image_name",
        help="Name of the image considered (only used to set the name of auxiliary files).",
        type=str,
    )

    parser.add_argument(
        "--solver", help="Sovler used (sara, hs or fhs).", default="fhs", type=str,
    )

    parser.add_argument(
        "--job_id", help="Integer identifier of the current run.", default=1, type=int,
    )

    parser.add_argument(
        "--ncubes", help="Number of spectral facets considered.", default=1, type=int,
    )

    parser.add_argument(
        "-rw",
        "--warmstart",
        help="Index of the last reweighting iteration to restart from.",
        default=0,
        type=int,
    )

    parser.add_argument(
        "-nRW",
        "--nreweight",
        help="Maximum number of reweighting iterations.",
        default=15,
        type=int,
    )

    args = parser.parse_args()

    # Sanity checks
    # assert condition, error_message"

    # Job submission
    for cube_id in range(1, args.ncubes + 1):

        job_name = """{0}_{1}_ID{2}_cube-{3}_gam-{4}_gambar-{5}""".format(
            args.solver,
            args.image_name,
            args.job_id,
            cube_id,
            args.gamma,
            args.gammabar,
        )

        # log file
        log_file = """{logpath}/{jobname}""".format(jobname=job_name, logpath=log_path)

        # command to launch
        sbatch_command = """sbatch --job-name {jobname} --error {log}.err --output {log}.out  --export=LOG={log}.log,PDIR={projectpath},IMAGENAME={0},ALGO={1},GAM={2},GAMBAR={3},SUBCUBE={4},RUNID={5},RW={6},nRW={7},LOG={log}.log  {slurmFile}""".format(
            args.image_name,
            args.solver,
            args.gamma,
            args.gammabar,
            cube_id,
            args.job_id,
            args.warmstart,
            args.nreweight,
            jobname=job_name,
            slurmFile=slurm_file,
            log=log_file,
            projectpath=project_path_run,
        )

        print(" ")
        print(
            sbatch_command
        )  # Uncomment this line when testing to view the sbatch command

        # launch job
        exit_status = subprocess.call(sbatch_command, shell=True)

        if exit_status == 1:  # Check to make sure the job submitted
            print("Job failed to submit: {}".format(sbatch_command))


if __name__ == "__main__":
    main()
