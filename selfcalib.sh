#!/bin/bash
# to activate the script, excute under script directory: chmod +x scriptname.sh
# to call the script: selfcalib.sh source date start end imsize
# eg. selfcalib.sh 0527-6549 3700 ${datearr[@]}

rm -rf /mnt/science1/bliu/ATCA/$1_selfcal_0601/
mkdir /mnt/science1/bliu/ATCA/$1_selfcal_0601/
cd /mnt/science1/bliu/ATCA/$1_selfcal_0601

comma=,

if [ ! -n "$3" ]; then
  echo "Date1 IS NULL"
else
  visstring=../$3/$1.1419.5.2
fi

if [ ! -n "$4" ]; then
  echo "Date2 IS NULL"
else
  visstringadd=../$4/$1.1419.5.2
  visstring="$visstring$comma$visstringadd"
fi

if [ ! -n "$5" ]; then
  echo "Date3 IS NULL"
else
  visstringadd=../$5/$1.1419.5.2
  visstring="$visstring$comma$visstringadd"
fi

if [ ! -n "$6" ]; then
  echo "Date4 IS NULL"
else
  visstringadd=../$6/$1.1419.5.2
  visstring="$visstring$comma$visstringadd"
fi

if [ ! -n "$7" ]; then
  echo "Date5 IS NULL"
else
  visstringadd=../$7/$1.1419.5.2
  visstring="$visstring$comma$visstringadd"
fi

echo $visstring
uvcat vis=$visstring stokes=i,v out=$1.1419.5.2.cat     #generate the visibility with the calibration table

invert vis=$1.1419.5.2.cat map=$1.1419.5.2.mfs.map.i,$1.1419.5.2.mfs.map.v beam=$1.1419.5.2.mfs.beam imsize=2048,2048 cell=2 stokes=i,v robust=-2 options=double,mfs
imhist in=$1.1419.5.2.mfs.map.v > imhist_$1.1419.5.2.mfs.map.v.txt

text=$(sed -n 11p imhist_$1.1419.5.2.mfs.map.v.txt)
arr=($text)
rms=${arr[2]}
rmsdec=$(echo $rms | awk '{printf("%f\n",$0)}')
factor1=3
factor2=7
threshold1=$(echo "$factor1 * $rmsdec" | bc)
threshold2=$(echo "$factor2 * $rmsdec" | bc)
echo $threshold1
echo $threshold2

clean map=$1.1419.5.2.mfs.map.i beam=$1.1419.5.2.mfs.beam out=$1.1419.5.2.mfs.clean cutoff=$threshold1
restor model=$1.1419.5.2.mfs.clean beam=$1.1419.5.2.mfs.beam map=$1.1419.5.2.mfs.map.i out=$1.1419.5.2.mfs.restor
fits in=$1.1419.5.2.mfs.clean op=xyout out=$1.1419.5.2.mfs.clean_b4slfcl.fits
fits in=$1.1419.5.2.mfs.restor op=xyout out=$1.1419.5.2.mfs.restor_b4slfcl.fits

imhist in=$1.1419.5.2.mfs.restor > $1.1419.5.2.mfs.restor_b4slfcl.txt
text=$(sed -n 11p $1.1419.5.2.mfs.restor_b4slfcl.txt)
arr=($text)
rms=${arr[2]}
rms_b4slfcl=$(echo $rms | awk '{printf("%f\n",$0)}')
echo $rms_b4slfcl

uvspec vis=$1.1419.5.2.cat stokes=i interval=1000 hann=9 options=avall,nobase nxy=1,1 log=$1.1419.5.2-hann9_b4slfcl.log
selfcal vis=$1.1419.5.2.cat model=$1.1419.5.2.mfs.clean clip=$threshold2 options=mfs
uvspec vis=$1.1419.5.2.cat stokes=i interval=1000 hann=9 options=avall,nobase nxy=1,1 log=$1.1419.5.2-hann9_aftslfcl.log

rm -rf $1.1419.5.2.mfs.beam
rm -rf $1.1419.5.2.mfs.clean
rm -rf $1.1419.5.2.mfs.restor
invert vis=$1.1419.5.2.cat map=$1.1419.5.2.mfs.map beam=$1.1419.5.2.mfs.beam imsize=2048,2048 cell=2 stokes=i robust=-2 line=chan,1000,$2 options=double,mfs
clean map=$1.1419.5.2.mfs.map beam=$1.1419.5.2.mfs.beam out=$1.1419.5.2.mfs.clean cutoff=$threshold1
restor model=$1.1419.5.2.mfs.clean beam=$1.1419.5.2.mfs.beam map=$1.1419.5.2.mfs.map out=$1.1419.5.2.mfs.restor
fits in=$1.1419.5.2.mfs.restor op=xyout out=$1.1419.5.2.mfs.restor_aftslfcl.fits

imhist in=$1.1419.5.2.mfs.restor > $1.1419.5.2.mfs.restor_aftslfcl.txt
text=$(sed -n 11p $1.1419.5.2.mfs.restor_aftslfcl.txt)
arr=($text)
rms=${arr[2]}
rms_b4slfcl=$(echo $rms | awk '{printf("%f\n",$0)}')
echo $rms_aftslfcl

cp ../$1_selfcal_1/cgcurs.region cgcurs.region
if [ ! -f "cgcurs.region" ]; then
  cgcurs in=$1.1419.5.2.mfs.restor range=0,0,log device=/xs nxy=1,1 options=cursor,mark,region,logfile
fi

invert vis=$1.1419.5.2.cat map=$1.1419.5.2.map beam=$1.1419.5.2.beam imsize=2048,2048 cell=2 stokes=i robust=-2 slop=0.5 line=chan,1000,$2 options=double
clean map=$1.1419.5.2.map beam=$1.1419.5.2.beam out=$1.1419.5.2.clean cutoff=$threshold1 region=@cgcurs.region
restor model=$1.1419.5.2.clean beam=$1.1419.5.2.beam map=$1.1419.5.2.map out=$1.1419.5.2.restor
fits in=$1.1419.5.2.restor op=xyout out=$1.1419.5.2.restor.fits

mkdir /mnt/science1/bliu/ATCA/results20190601/
mkdir /mnt/science1/bliu/ATCA/results20190601/$1
mv -b -f $1.1419.5.2.mfs.restor_b4slfcl.fits /mnt/science1/bliu/ATCA/results20190601/$1/
mv -b -f $1.1419.5.2.mfs.restor_aftslfcl.fits /mnt/science1/bliu/ATCA/results20190601/$1/
mv -b -f $1.1419.5.2-hann9_b4slfcl.log /mnt/science1/bliu/ATCA/results20190601/$1/
mv -b -f $1.1419.5.2-hann9_aftslfcl.log /mnt/science1/bliu/ATCA/results20190601/$1/
mv -b -f $1.1419.5.2.restor.fits /mnt/science1/bliu/ATCA/results20190601/$1/
