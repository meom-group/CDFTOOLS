#!/bin/ksh

# gib.ksh : create a matlab file with
#   1rst line = depth
#   2nd line  = Tlevitus
#   3rd line  = Slevitus
#   then a pair of lines (T S ) foreach year. First column gives the year
dir=$( basename `pwd` )
CONFIG=${dir%-DIAGS}
if [  ! -d ../${CONFIG}-MONITOR ] ; then mkdir ../${CONFIG}-MONITOR ; fi

\rm -f ${CONFIG}_gib.mtl
# if no *TGIB.txt files in the dir, skip
 ls *TGIB.txt 1> /dev/null 2>&1
 if [ $? != 0 ] ; then echo no gib files to deal with ... ; exit ; fi

dep=$( cat  LEVITUS_y0000_TGIB.txt | grep 'Mean value at level' | awk '{ printf "%8.1f",$7 }' )
Tlev=$( cat LEVITUS_y0000_TGIB.txt | grep 'Mean value at level' | awk '{ printf "%8.4f", $9 }' )
Slev=$( cat LEVITUS_y0000_SGIB.txt | grep 'Mean value at level' | awk '{ printf "%8.4f", $9 }' )

echo 0000 $dep > ${CONFIG}_gib.mtl
echo 0000 $Tlev >> ${CONFIG}_gib.mtl
echo 0000 $Slev >> ${CONFIG}_gib.mtl

for tfile  in ${CONFIG}_y????_TGIB.txt ; do
   sfile=$( echo $tfile | sed -e 's/TGIB/SGIB/' )
   year=$( head -1 $tfile )
   Tcur=$( cat $tfile | grep 'Mean value at level' | awk '{ printf "%8.4f", $9 }' )
   Scur=$( cat $sfile | grep 'Mean value at level' | awk '{ printf "%8.4f", $9 }' )
  echo $year $Tcur >> ${CONFIG}_gib.mtl
  echo $year $Scur >> ${CONFIG}_gib.mtl
done

mv ${CONFIG}_gib.mtl ../${CONFIG}-MONITOR/
