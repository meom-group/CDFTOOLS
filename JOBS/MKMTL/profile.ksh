#!/bin/ksh
#  section.ksh : build a line of matlab file for section output
#set -x
dir=$( basename `pwd` )
CONFIG=${dir%-DIAGS}

if [  ! -d ../${CONFIG}-MONITOR ] ; then mkdir ../${CONFIG}-MONITOR ; fi

# if no TMEAN.txt files in the dir, skip
 ls *TMEAN.txt 1> /dev/null 2>&1
 if [ $? != 0 ] ; then echo no TMEAN to deal with ... ; exit ; fi
 ls *SMEAN.txt 1> /dev/null 2>&1
 if [ $? != 0 ] ; then echo no SMEAN to deal with ... ; exit ; fi
 ls *SSHMEAN.txt 1> /dev/null 2>&1
 if [ $? != 0 ] ; then echo no SSHMEAN to deal with ... ; exit ; fi
n=0

for file in ${CONFIG}*TMEAN.txt ; do
   n=$(( $n +1 ))
   year=$( head -1 $file )
   dep=$( cat $file | grep -e 'Mean value at level' | awk '{ printf "%10.1f" ,  $7 }' )
   mean=$( cat $file | grep -e 'over' | awk '{ printf "%10.4f" ,  $6 }' )
   tem=$( cat $file | grep  -e 'Mean value at level' | awk '{ printf "%10.4f" ,  $9 }' )

   f=${CONFIG}_y${year}_SSHMEAN.txt
   sshmean=$( cat $f | grep ocean | awk '{ printf " %8.4f " , $6 }' )

   if [ $n = 1 ] ; then
     echo "% T  SSH mean diags for " $CONFIG >  ${CONFIG}_TMEAN.mtl
     echo "% yr     ssh     Tmean     <----------------------- depth ------ ...... "  >>  ${CONFIG}_TMEAN.mtl
     echo 0 0 0 | awk '{ printf "%04d % 8.4f % 8.4f ",$1,$2,$3 }'  >> ${CONFIG}_TMEAN.mtl
     echo $dep | awk '{  for ( i=1 ; i <= NF ; i++ ) printf "% 8.1f ", $i }'  >>  ${CONFIG}_TMEAN.mtl
     printf "\n" >> ${CONFIG}_TMEAN.mtl
   fi
   printf "%04d " $year  >>  ${CONFIG}_TMEAN.mtl
   echo $sshmean $mean | awk '{ printf "% 8.4f % 8.4f ", $1, $2 }' >>  ${CONFIG}_TMEAN.mtl
   echo $tem |  awk '{  for ( i=1 ; i <= NF ; i++ )  printf "% 8.4f " ,$i }' >> ${CONFIG}_TMEAN.mtl
   printf "\n" >> ${CONFIG}_TMEAN.mtl
done


mv ${CONFIG}_TMEAN.mtl ../${CONFIG}-MONITOR/

n=0
for file in ${CONFIG}*SMEAN.txt ; do 
   n=$(( $n +1 ))
   year=$( head -1 $file )
   dep=$( cat $file | grep -e 'Mean value at level' | awk '{ printf "%10.1f" ,  $7 }' )
   mean=$( cat $file | grep -e 'over' | awk '{ printf "%10.4f" ,  $6 }' )
   tem=$( cat $file | grep  -e 'Mean value at level' | awk '{ printf "%10.4f" ,  $9 }' )

   f=${CONFIG}_y${year}_SSHMEAN.txt
   sshmean=$( cat $f | grep ocean | awk '{ printf " %8.4f " , $6 }' )

   if [ $n = 1 ] ; then
     echo "% S  SSH mean diags for " $CONFIG >  ${CONFIG}_SMEAN.mtl
     echo "% yr     ssh     Smean     <----------------------- depth ------ ...... "  >>  ${CONFIG}_SMEAN.mtl
     echo 0 0 0 | awk '{ printf "%04d % 8.4f % 8.4f ",$1,$2,$3 }'  >> ${CONFIG}_SMEAN.mtl
     echo $dep | awk '{  for ( i=1 ; i <= NF ; i++ ) printf "% 8.1f ", $i }'  >>  ${CONFIG}_SMEAN.mtl
     printf "\n" >> ${CONFIG}_SMEAN.mtl
   fi
   printf "%04d " $year  >>  ${CONFIG}_SMEAN.mtl
   echo $sshmean $mean | awk '{ printf "% 8.4f % 8.4f ", $1, $2 }' >>  ${CONFIG}_SMEAN.mtl
   echo $tem |  awk '{  for ( i=1 ; i <= NF ; i++ )  printf "% 8.4f " ,$i }' >> ${CONFIG}_SMEAN.mtl
   printf "\n" >> ${CONFIG}_SMEAN.mtl
done

mv ${CONFIG}_SMEAN.mtl ../${CONFIG}-MONITOR/

