#!/bin/bash

algoversion=cw
rwtype=dirty
windowtype=triangular
overlap=0.5
Qx=2
Qy=2
Qc=1
rw=-1
# updatereg=0

logpath=$(pwd)/results/cygASband_Cube_512_1024_20_test/${algoversion}/slurm_logs
outputfile=$(pwd)/results/cygASband_Cube_512_1024_20_test/${algoversion}/slurm_errors.txt

for updatereg in 0 1; do
    for homotopy in 0 1; do # {0..1}
        for rwtype in dirty ground_truth; do
            for alph in 1e-1 1 10; do
                for alphbar in 1e-1 1 10; do

                    echo h=${homotopy}, updatereg=${updatereg}, alpha=${alph}, alpha_bar=${alphbar} >> ${outputfile}
                    echo ----- >> ${outputfile}

                    # filename=${logpath}/spatial_cygASband_Cube_512_1024_20_${algoversion}_L=${nchannels}_Qx=${Qx}_Qy=${Qy}_Qc=${Qc}_id=1_overlapx=${overlap}_overlapy=${overlap}_gamma=${alph}_gammabar=${alphbar}_rw=${rw}_rwt=${rwtype}_exptype=test_srf=2_snr=40_homotopy=${homotopy}.err

                    filename=${logpath}/spatial_cygASband_Cube_512_1024_20_${algoversion}_L=${nchannels}_Qx=${Qx}_Qy=${Qy}_Qc=${Qc}_id=1_overlapx=${overlap}_overlapy=${overlap}_gamma=${alph}_gammabar=${alphbar}_rw=${rw}_rwt=${rwtype}_exptype=test_srf=2_snr=40_homotopy=${homotopy}_updatereg=${updatereg}.err

                    grep -A 7 Error ${filename} >> ${outputfile}
                    echo ===== >> ${outputfile}
                done
            done
        done
    done
done