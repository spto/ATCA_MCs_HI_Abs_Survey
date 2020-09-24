#!/bin/bash
# to activate the script, excute under script directory: chmod +x scriptname.sh
# to call the script: calib180126.sh

cd /mnt/nau/bliu/ATCA/$1

band=$2.1419.5.2
phs=$3.1419.5.2
targ=$4.1419.5.2

mfcal vis=$band options=interpolate
gpcal vis=$band
gpcopy vis=$band out=$phs
gpcal vis=$phs options=nopass
gpboot vis=$phs cal=$band
gpcopy vis=$phs out=$targ

echo "finished"
