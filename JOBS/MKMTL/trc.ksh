#!/bin/ksh
# This script scan the years of TRC*.dat diags
#set -x
dir=$( basename `pwd` )
CONFIG=${dir%-DIAGS}

firstline() {  cat $1 | grep -n J | awk -F: '{print $1+1}' ; }
# if no TRCmean.dat files in the dir, skip
 ls *TRCmean.dat 1> /dev/null 2>&1
 if [ $? != 0 ] ; then echo no TRCmean.dat with ... ; exit ; fi

 \rm -f trc_cfc.mtl
 \rm -f trc_c14.mtl

n=0
for f in ${CONFIG}_y????_TRCmean.dat ; do
 n=$(( n + 1 ))
 zonalmean=$( echo $f | sed -e 's/TRCmean/TRCzonalmean/' )
 zonalsum=$( echo $f | sed -e 's/TRCmean/TRCzonalsum/' )
 zonalsurf=$( echo $f | sed -e 's/TRCmean/TRCzonalsurf/' )

 year=$( cat $f | awk '{ print  $1}' )
 cfcmean=$( cat $f | awk '{ print  $2}' ) 
 c14mean=$( cat $f | awk '{ print  $3}' )

if [ $n == 1 ] ; then
 echo "% CFC monitoring file for " $CONFIG >> trc_cfc.mtl
 echo "% Each year is represented by  3 raws of data (latitude) zonal mean/int : CFC inv , CFC inv (integral), CFC conc (surf)" >> trc_cfc.mtl
 echo "% yr   Total  <-----------  latitude   ----------------------------------------------- .... ----------> " >> trc_cfc.mtl

 echo "% B-C14 monitoring file for " $CONFIG >> trc_c14.mtl
 echo "% Each year is represented by  3 raws of data (latitude) zonal mean/int : C14 inv , C14 inv (integral), C14 conc (surf)" >> trc_c14.mtl
 echo "% yr   Total  <-----------  latitude   ----------------------------------------------- .... ----------> " >> trc_c14.mtl

 deb=$(firstline $zonalmean) 
 lat=$( cat $zonalmean | awk '{ if ( NR >= deb ) {print $2 }}  END {printf "\n" }' deb=$deb )

 echo   0 0 | awk '{ printf "%04d  % 13.6e " , $1 ,$2 }'>>  trc_cfc.mtl
 echo   0 0 | awk '{ printf "%04d  % 13.6e " , $1 ,$2 }'>>  trc_c14.mtl
 echo $lat  >> trc_cfc.mtl
 echo $lat  >> trc_c14.mtl
fi

 deb=$(firstline $zonalmean)
 cfcinv=$( cat $zonalmean | awk '{ if ( NR >= deb ) {print $3 }}  END {printf "\n" }' deb=$deb )
 c14inv=$( cat $zonalmean | awk '{ if ( NR >= deb ) {print $4 }}  END {printf "\n" }' deb=$deb )

 deb=$(firstline $zonalsum)
 cfcint=$( cat $zonalsum | awk '{ if ( NR >= deb ) {print $3 }}  END {printf "\n" }' deb=$deb )
 c14int=$( cat $zonalsum | awk '{ if ( NR >= deb ) {print $4 }}  END {printf "\n" }' deb=$deb )

 deb=$(firstline $zonalsurf)
 cfcsurf=$( cat $zonalsurf | awk '{ if ( NR >= deb ) {print $3 }}  END {printf "\n" }' deb=$deb )
 c14surf=$( cat $zonalsurf | awk '{ if ( NR >= deb ) {print $4 }}  END {printf "\n" }' deb=$deb )

 echo $year $cfcmean | awk '{ printf "%04d  % 13.6e " , $1 ,$2 }'>>  trc_cfc.mtl
 echo $cfcinv  >> trc_cfc.mtl
 echo $year $cfcmean | awk '{ printf "%04d  % 13.6e " , $1 ,$2 }'>>  trc_cfc.mtl
 echo $cfcint  >> trc_cfc.mtl
 echo $year $cfcmean | awk '{ printf "%04d  % 13.6e " , $1 ,$2 }'>>  trc_cfc.mtl
 echo $cfcsurf  >> trc_cfc.mtl

 echo $year $c14mean | awk '{ printf "%04d  % 13.6e " , $1 ,$2 }'>>  trc_c14.mtl
 echo $c14inv >>  trc_c14.mtl
 echo $year $c14mean | awk '{ printf "%04d  % 13.6e " , $1 ,$2 }'>>  trc_c14.mtl
 echo $c14int  >> trc_c14.mtl
 echo $year $c14mean | awk '{ printf "%04d  % 13.6e " , $1 ,$2 }'>>  trc_c14.mtl
 echo $c14surf  >> trc_c14.mtl


done

mv trc_cfc.mtl ../${CONFIG}-MONITOR/${CONFIG}_trc_cfc.mtl
mv trc_c14.mtl ../${CONFIG}-MONITOR/${CONFIG}_trc_c14.mtl

