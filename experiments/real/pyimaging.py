#!/usr/bin/env python
import csv, subprocess
import sys
import os
# global vars
project_path_run=os.path.dirname(os.path.realpath(__file__))
#project_path_run="%s/experiments/real/"%project_path
#setting paths 
log_path  = "%s/logs"%project_path_run
try:
    os.mkdir(log_path)
except OSError as error:
    print(error) 

try:
   run_path= "%s/runs"%project_path_run
   os.mkdir(run_path)
except OSError as error:
    print(error)

slurm_file = "%s/imaging.slurm"%project_path_run

def main():

   parameter_file ="%s/%s"%(project_path_run, sys.argv[1]);

   print("Current experiment Dir %s"%project_path_run)
   print("Param file: %s"%parameter_file)

   with open(parameter_file, "rt") as csvfile:
        reader = csv.reader(csvfile)
        for job in reader:
            job_name="""{0}_ID{1}_Subcube{2}Algo{3}""".format(*job)
            #log file
            log_file = """{logpath}/{jobname}""".format(jobname=job_name,logpath=log_path)
            
            # command to launch
            sbatch_command = """sbatch --job-name {jobname} --error {log}.err --output {log}.out \
                    --export=LOG={log}.log,PDIR={projectpath},GAM={4},GAMBAR={5},  \
                    SUBCUBE={2},ALGO={3},RUNID={1},IMAGENAME={0},RW={6},nRW={7},LOG={log}.log \ 
                    {slurmFile}""".format(*job,jobname=job_name,slurmFile=slurm_file,log=log_file,projectpath=project_path_run) 
            
            print(sbatch_command) # Uncomment this line when testing to view the sbatch command
            
            #launch job
            exit_status = subprocess.call(sbatch_command, shell=True)
           
            if exit_status is 1:  # Check to make sure the job submitted
              print("Job {11} failed to submit".format(*job))



if __name__ == "__main__":
   main()
