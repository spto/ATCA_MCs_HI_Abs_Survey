#!/bin/bash
# to activate the script, excute under script directory: chmod +x scriptname.sh
# to call the script: loadnsplit.sh 170210

cd /mnt/nau/bliu/ATCA/$1
atlod in=* out=$1.uv restfreq=1.4204058 options=birdie,rfiflag,noauto,xycorr,bary
uvsplit vis=$1.uv
