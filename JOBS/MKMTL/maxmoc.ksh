#!/bin/ksh

dir=$( basename `pwd` )
CONFIGCASE=${dir%-DIAGS}

if [  ! -d ../${CONFIGCASE}-MONITOR ] ; then mkdir ../${CONFIGCASE}-MONITOR ; fi

CONFIG=${CONFIGCASE%-*}
CASE=${CONFIGCASE#*-}

\rm -f -r ${CONFIGCASE}_maxmoc.mtl
# if no minmaxmoc.txt files in the dir, skip
 ls *minmaxmoc.txt 1> /dev/null 2>&1
 if [ $? != 0 ] ; then echo no minmaxmoc to deal with ... ; exit ; fi
 ls *heattrp.dat 1> /dev/null 2>&1
 if [ $? != 0 ] ; then echo no heattrp to deal with ... ; exit ; fi
touch  ${CONFIGCASE}_maxmoc.mtl

for file in *_minmaxmoc.txt 
do
  year=$( head -1 $file )
  mht=${CONFIGCASE}_y${year}_heattrp.dat

  # GLO
  maxglo=$( cat $file    | grep -e '^Glo' | grep Max | awk '{ printf "%8.3f" ,  $3 }' )
# maxglolat=$( cat $file | grep -e '^Glo' | grep Max | awk '{ printf "%8.2f" ,  $6 }' )
# maxglodep=$( cat $file | grep -e '^Glo' | grep Max | awk '{ printf "%8.2f" ,  $9 }' )

  minglo=$( cat $file    | grep -e '^Glo' | grep Min | awk '{ printf "%8.3f" ,  $3 }' )
# minglolat=$( cat $file | grep -e '^Glo' | grep Min | awk '{ printf "%8.2f" ,  $6 }' )
# minglodep=$( cat $file | grep -e '^Glo' | grep Min | awk '{ printf "%8.2f" ,  $9 }' )

  # ATL
  maxatl=$( cat $file    | grep -e '^Atl' | grep Max | awk '{ printf "%8.3f" ,  $3 }' )
# maxatllat=$( cat $file | grep -e '^Atl' | grep Max | awk '{ printf "%8.2f" ,  $6 }' )
# maxatldep=$( cat $file | grep -e '^Atl' | grep Max | awk '{ printf "%8.2f" ,  $9 }' )

  minatl=$( cat $file    | grep -e '^Atl' | grep Min | awk '{ printf "%8.3f" ,  $3 }' )
# minatllat=$( cat $file | grep -e '^Atl' | grep Min | awk '{ printf "%8.2f" ,  $6 }' )
# minatldep=$( cat $file | grep -e '^Atl' | grep Min | awk '{ printf "%8.2f" ,  $9 }' )

  # INP : attention we have 2 Minimum for INP mininp1 and mininp2
  tmp=$( cat $file    | grep -e '^Inp' | grep Min | awk '{ printf "%8.3f" ,  $3 }' )
# mininplat=$( cat $file | grep -e '^Inp' | grep Min | awk '{ printf "%8.2f" ,  $6 }' )
# mininpdep=$( cat $file | grep -e '^Inp' | grep Min | awk '{ printf "%8.2f" ,  $9 }' )
  mininp1=$( echo $tmp | awk '{print $1}' )
  mininp2=$( echo $tmp | awk '{print $2}' )

  # AUS
  maxaus=$( cat $file    | grep -e '^Aus' | grep Max | awk '{ printf "%8.3f" ,  $3 }' )
# maxauslat=$( cat $file | grep -e '^Aus' | grep Max | awk '{ printf "%8.2f" ,  $6 }' )
# maxausdep=$( cat $file | grep -e '^Aus' | grep Max | awk '{ printf "%8.2f" ,  $9 }' )

  minaus=$( cat $file    | grep -e '^Aus' | grep Min | awk '{ printf "%8.3f" ,  $3 }' )
# minauslat=$( cat $file | grep -e '^Aus' | grep Min | awk '{ printf "%8.2f" ,  $6 }' )
# minausdep=$( cat $file | grep -e '^Aus' | grep Min | awk '{ printf "%8.2f" ,  $9 }' )

  # heattrp at 20 N
  heattrp=$(  cat $mht | awk '{  if ( $2 >= 20 ) { atlmht=$4 ; glomht=$3 } } END { printf "%6.3f %6.3f ", glomht, atlmht }' )
  mhtglo=$( echo $heattrp | awk '{print $1}' )
  mhtatl=$( echo $heattrp | awk '{print $2}' )
 
# echo $year $maxglo $maxglolat $maxglodep $minglo $minglolat $minglodep \
#            $maxatl $maxatllat $maxatldep $minatl $minatllat $minatldep \
#            $maxinp $maxinplat $maxinpdep $mininp $mininplat $mininpdep \
#            $maxaus $maxauslat $maxausdep $minaus $minauslat $minausdep  >>  ${CONFIGCASE}_maxmoc.mtl

 echo $year $maxglo $minglo $mhtglo $maxatl $minatl $mhtatl  $mininp1 $mininp2 0000  \
            $maxaus $minaus 0000 >>  ${CONFIGCASE}_maxmoc.mtl

done
mv ${CONFIGCASE}_maxmoc.mtl   ../${CONFIGCASE}-MONITOR
