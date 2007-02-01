#!/bin/ksh

cd $tmpdir
for f in *MOC* ; do
  CONFIG=${f%-*}
  tmp=${f#*-} ; CASE=${tmp%%_*}
  tmp=${f%_*} ; year=${tmp#*_y}
  echo $CONFIG $CASE $year
  outfile=${CONFIG}-${CASE}_y${year}_maxmoc40.txt
  \rm -f $outfile

    echo $year > $outfile
# GLO
printf "%s" 'Glo ' >>  $outfile ; cdfmaxmoc $f glo 40 40 500 2000 | grep Maximum >> $outfile
printf "%s" 'Glo ' >>  $outfile ; cdfmaxmoc $f glo -30 -30 2000 5500 | grep Maximum >> $outfile
# ATL
printf "%s" 'Atl ' >>  $outfile ; cdfmaxmoc $f atl 40 40 500 2000 | grep Maximum >> $outfile
printf "%s" 'Atl ' >>  $outfile ; cdfmaxmoc $f atl -30 -30  500 2000 | grep Maximum >> $outfile
#INP
printf "%s" 'Inp ' >>  $outfile ; cdfmaxmoc $f inp -30 -30 1000 5500  | grep Minimum >> $outfile
#AUS
printf "%s" 'Aus ' >>  $outfile ; cdfmaxmoc $f glo -50 -50 0 2000   | grep Maximum >> $outfile

mfput $outfile ${CONFIG}/${CONFIG}-${CASE}-DIAGS/
done


