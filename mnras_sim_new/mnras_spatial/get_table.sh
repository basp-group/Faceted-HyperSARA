#!/bin/bash

wintype=triangular
overlap=0.5
rw=-1
alph=1
alphbar=1
updatereg=0
homotopy=0
isnr=40
superresolution=3

# grep pattern + sourrounding lines
# https://stackoverflow.com/questions/9081/grep-a-file-but-show-several-surrounding-lines

imagename=cygASband_Cube_512_1024_20
exptype=spatial
algoversion=cw
Qx=2
Qy=2
Qc=1
nchannels=20
rwtype=heuristic # dirty heuristic
regtype=heuristic # inv log heuristic
xapprox=none # dirty precond
noisetransfer=none # precond none
regoption=none # dirty none

logpath=$(pwd)/results/${imagename}_${exptype}/${algoversion}/logs
outputfile=$(pwd)/${algoversion}_rw=${rwtype}_reg=${regtype}.txt

for l in {1..${Qc}}; do
    # for updatereg in 0 1; do
    #     for homotopy in 0; do # {0..1} 0 1
    #         for rwtype in dirty; do # dirty ground_truth
    #             for alph in 1e-1 1 10; do
    #                 for alphbar in 1e-1 1 10; do

    echo h=${homotopy}, updatereg=${updatereg}, regtype=${regtype}, alpha=${alph}, alpha_bar=${alphbar} >> ${outputfile}
    echo ----- >> ${outputfile}

    filename=${logpath}/${exptype}_${imagename}_${algoversion}_${wintype}_id=${ind}_Qx=${Qx}_Qy=${Qy}_Qc=${Qc}_overlap=${overlapx}_gam=${alph}_gambar=${alphbar}_rw=${rw}_homotopy=${homotopy}_rwt=${rwtype}_exptype=${exptype}_srf=${superresolution}_snr=${isnr}_updatereg=${updatereg}_regtype=${regtype}_xapprox=${xapprox}_nt=${noisetransfer}_ropt=${regoption}.txt
    
    grep 'Job id:' ${filename} >> ${outputfile}
    grep 'Rwt:' ${filename} >> ${outputfile}
    grep -A 7 'Updated reg' ${filename} >> ${outputfile} 
    grep 'Reweighting:' ${filename} | tail -1 >> ${outputfile}
    grep 'epsilon ='  ${filename} | tail -1 >> ${outputfile}
    grep 'SNR ='  ${filename} | tail -1 >> ${outputfile}
    echo ===== >> ${outputfile}
    #                 done
    #             done
    #         done
    #     done
    # done
done