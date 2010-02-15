#!/bin/ksh
#  ice.ksh : build a line of matlab file for ice output
#set -x
dir=$( basename `pwd` )
CONFIG=${dir%-DIAGS}

if [  ! -d ../${CONFIG}-MONITOR ] ; then mkdir ../${CONFIG}-MONITOR ; fi

\rm -f  ${CONFIG}_ice.mtl

# if no ice.txt files in the dir, skip
 ls *ice.txt 1> /dev/null 2>&1
 if [ $? != 0 ] ; then echo no ice to deal with ... ; exit ; fi
n=0

for file in *ice.txt
do
n=$(( $n +1 ))
year=$( head -1 $file | awk '{ print $2}' ) 
nvol=$(  cat $file | grep -e 'NVolume' | grep -v NVolumet | awk '{ printf "%.0f  ", $4}' ) 
svol=$(  cat $file | grep -e 'SVolume' | grep -v SVolumet | awk '{ printf "%.0f  ", $4}' )
narea=$(  cat $file | grep -e 'NArea' | awk '{ printf "%.0f  ", $4}' ) 
sarea=$(  cat $file | grep -e 'SArea' | awk '{ printf "%.0f  ", $4}' )
nextent=$(  cat $file | grep -e 'NExtend' | awk '{ printf "%.0f  ", $4}' )
sextent=$(  cat $file | grep -e 'SExtend' | awk '{ printf "%.0f  ", $4}' )

if [ $n == 1 ] ; then
   echo 0000 02 03 08 09  02 03 08 09 02 03 08 09  02 03 08 09 02 03 08 09  02 03 08 09 > ${CONFIG}_ice.mtl
fi

echo $year $nvol $svol $narea $sarea $nextent $sextent >> ${CONFIG}_ice.mtl

done

mv ${CONFIG}_ice.mtl ../${CONFIG}-MONITOR/


