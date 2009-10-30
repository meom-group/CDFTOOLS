#!/bin/ksh
#set -x

dir=$( basename `pwd` )
CONFCASE=${dir%-DIAGS}

CONFIG=${CONFCASE%-*}
CASE=${CONFCASE#*-}


if [  ! -d ../${CONFCASE}-MONITOR ] ; then mkdir ../${CONFCASE}-MONITOR ; fi

n=0
# mini and maxi of sigma0
  mini=25.2 ; maxi=28.5

cd TRPSIG

for file in ${CONFCASE}*_01_Denmark_strait_trpsig.txt ; do
   tmp=${file#*_y}
   year=${tmp%%_*.txt}
   fil2=$( echo $file | sed -e 's/01_Denmark_strait/02_Faoes_Bank_Channel/' )
   n=$(( $n +1 ))
   #awk '{ if ( FNR > 10  && FNR < 20 ) { print } }'
   sig=$( cat $file | grep -v -e '^#' | awk '{ if ( $1 >   25.2  && $1 < 28.5 ) {printf "%8.3f" ,  $1 }}' mini=$mini  maxi=$maxi )
   trp01=$( cat $file | grep -v -e '^#' | awk '{ if ( $1 > 25.2  && $1 < 28.5 ) {printf "%13.4e" ,  $2 } }' mini=$mini  maxi=$maxi)
   trp02=$( cat $fil2 | grep -v -e '^#' | awk '{ if ( $1 > 25.2  && $1 < 28.5 ) {printf "%13.4e" ,  $2 } }' mini=$mini  maxi=$maxi)

   if [ $n = 1 ] ; then
     echo 000000  $sig > ../${CONFCASE}_TRPSIG.mtl
   fi
   echo $year $trp01 >>  ../${CONFCASE}_TRPSIG.mtl
   echo $year $trp02 >>  ../${CONFCASE}_TRPSIG.mtl

done

cd ../
mv ${CONFCASE}_TRPSIG.mtl ../${CONFCASE}-MONITOR


