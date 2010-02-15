#!/bin/ksh

dir=$( basename `pwd` )
CONFIGCASE=${dir%-DIAGS}

if [  ! -d ../${CONFIGCASE}-MONITOR ] ; then mkdir ../${CONFIGCASE}-MONITOR ; fi

CONFIG=${CONFIGCASE%-*}
CASE=${CONFIGCASE#*-}

# if no maxmoc40.txt files in the dir, skip
 ls *maxmoc40.txt 1> /dev/null 2>&1
 if [ $? != 0 ] ; then echo no maxmoc40 to deal with ... ; exit ; fi
 ls *heattrp.dat 1> /dev/null 2>&1
 if [ $? != 0 ] ; then echo no heattrp to deal with ... ; exit ; fi

\rm -f ${CONFIGCASE}_maxmoc40.mtl
touch  ${CONFIGCASE}_maxmoc40.mtl

for file in *_maxmoc40.txt 
do
  year=$( head -1 $file )
  mht=${CONFIGCASE}_y${year}_heattrp.dat

  # GLO  max a 40 et -30
  tmp=$( cat $file    | grep -e '^Glo' | grep Max | awk '{ printf "%8.3f" ,  $3 }' )
  maxglo40n=$( echo $tmp | awk '{print $1}' )
  maxglo30s=$( echo $tmp | awk '{print $2}' )

  # ATL  max a 40 et -30
  tmp=$( cat $file    | grep -e '^Atl' | grep Max | awk '{ printf "%8.3f" ,  $3 }' )
  maxatl40n=$( echo $tmp | awk '{print $1}' )
  maxatl30s=$( echo $tmp | awk '{print $2}' )

  # INP : 1 Min at -30 S
  mininp30s=$( cat $file    | grep -e '^Inp' | grep Min | awk '{ printf "%8.3f" ,  $3 }' )

  # AUS 1 max at -50 S
  maxaus50s=$( cat $file    | grep -e '^Aus' | grep Max | awk '{ printf "%8.3f" ,  $3 }' )

  # heattrp at 20 N
#  heattrp=$(  cat $mht | awk '{  if ( $2 >= 20 ) { atlmht=$4 ; glomht=$3 } } END { printf "%6.3f %6.3f ", glomht, atlmht }' )
#  mhtglo=$( echo $heattrp | awk '{print $1}' )
#  mhtatl=$( echo $heattrp | awk '{print $2}' )
 
 echo $year $maxglo40n $maxglo30s  $maxatl40n $maxatl30s  $mininp30s $maxaus50s >>  ${CONFIGCASE}_maxmoc40.mtl

done
mv ${CONFIGCASE}_maxmoc40.mtl   ../${CONFIGCASE}-MONITOR
