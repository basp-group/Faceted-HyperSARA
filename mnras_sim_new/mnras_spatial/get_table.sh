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

logpath=$(pwd)/results/cygASband_Cube_512_1024_20_test/${algoversion}/logs
outputfile=$(pwd)/results/cygASband_Cube_512_1024_20_test/${algoversion}/output_parameters.txt

# grep pattern + sourrounding lines
# https://stackoverflow.com/questions/9081/grep-a-file-but-show-several-surrounding-lines

for updatereg in 0 1; do
    for homotopy in 0 1; do # {0..1}
        for rwtype in dirty ground_truth; do
            for alph in 1e-1 1 10; do
                for alphbar in 1e-1 1 10; do

                    echo h=${homotopy}, updatereg=${updatereg}, alpha=${alph}, alpha_bar=${alphbar} >> ${outputfile}
                    echo ----- >> ${outputfile}

                    # filename=${logpath}/spatial_cygASband_Cube_512_1024_20_${algoversion}_${windowtype}_id=1_Qx=${Qx}_Qy=${Qy}_Qc=${Qc}_overlap=${overlap}_gam=${alph}_gambar=${alphbar}_rw=-1_homotopy=${homotopy}_rwt=${rwtype}_exptype=test_srf=2_snr=40.txt 
                    # grep -A 7 'Reweighting type' ${filename} >> ${outputfile} # for older log files

                    filename=${logpath}/spatial_cygASband_Cube_512_1024_20_${algoversion}_${windowtype}_id=1_Qx=${Qx}_Qy=${Qy}_Qc=${Qc}_overlap=${overlap}_gam=${alph}_gambar=${alphbar}_rw=${rw}_homotopy=${homotopy}_rwt=${rwtype}_exptype=test_srf=2_snr=40_updatereg=${updatereg}.txt 
                    
                    grep 'Job id:' ${filename} >> ${outputfile}
                    grep 'Rwt:' ${filename} >> ${outputfile}
                    grep 'Reweighting:' ${filename} | tail -1 >> ${outputfile}
                    grep 'epsilon ='  ${filename} | tail -1 >> ${outputfile}
                    grep 'SNR ='  ${filename} | tail -1 >> ${outputfile}
                    echo ===== >> ${outputfile}
                done
            done
        done
    done
done