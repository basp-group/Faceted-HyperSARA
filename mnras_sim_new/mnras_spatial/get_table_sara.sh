#!/bin/bash

algoversion=sara
wintype=none
overlap=0
Qx=1
Qy=1
Qc=1
rw=-1
superresolution=2
isnr=40
updatereg=0
homotopy=0
alph=1
alphbar=1

# grep pattern + sourrounding lines
# https://stackoverflow.com/questions/9081/grep-a-file-but-show-several-surrounding-lines

imagename=cygASband_Cube_512_1024_20
exptype=spatial
nchannels=20
rwtype=heuristic # dirty heuristic
regtype=heuristic # inv log heuristic
xapprox=none # dirty precond
noisetransfer=none # precond none
regoption=none # dirty none

logpath=$(pwd)/results/cygASband_Cube_512_1024_20_test/${algoversion}/logs
outputfile=$(pwd)/${algoversion}_rw=${rwtype}_reg=${regtype}.txt

for l in {1..${nchannels}}; do
    echo h=${homotopy}, updatereg=${updatereg}, regtype=${regtype}, alpha=${alph}, alpha_bar=${alphbar} >> ${outputfile}
    echo ----- >> ${outputfile}

    filename=${logpath}/${exptype}_${imagename}_${algoversion}_${wintype}_id=${l}_Qx=${Qx}_Qy=${Qy}_Qc=${Qc}_overlap=${overlap}_gam=${alph}_gambar=${alphbar}_rw=${rw}_homotopy=${homotopy}_rwt=${rwtype}_exptype=${exptype}_srf=${superresolution}_snr=${isnr}_updatereg=${updatereg}_regtype=${regtype}_xapprox=${xapprox}_nt=${noisetransfer}_ropt=${regoption}.txt
    
    grep 'Job id:' ${filename} >> ${outputfile}
    grep 'Rwt:' ${filename} >> ${outputfile}
    grep -A 7 'Updated reg' ${filename} >> ${outputfile} 
    grep 'Reweighting:' ${filename} | tail -1 >> ${outputfile}
    grep 'epsilon ='  ${filename} | tail -1 >> ${outputfile}
    grep 'SNR ='  ${filename} | tail -1 >> ${outputfile}
    echo ===== >> ${outputfile}
done
