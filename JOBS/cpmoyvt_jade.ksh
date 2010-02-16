#!/bin/ksh
. config_def.ksh

run=1
while (( $run==1 )); do
sleep 3
qstat -u $USER | grep -q metamoyvt || run=0
echo "run = $run"
done

for year in `ls -d ???? `
do
echo $year
nfilenc=$( ls -l $year/*.nc | wc -l )
echo "nfilenc = $nfilenc"
nfileok=$( ls -l $year/OK? | wc -l )
echo "nfileok = $nfileok"
if [ $nfilenc -eq 130 -a $nfileok -eq 6 ]; then
echo "cp file for year  "
if [ ! -d /scratch/$USER/${CONFIG}-${CASE}-MEAN/ ] ; then mkdir /scratch/$USER/${CONFIG}-${CASE}-MEAN/ ; fi;
if [ ! -d /data/$USER/$CONFIG/${CONFIG}-${CASE}-MEAN/$year ] ; then mkdir /data/$USER/$CONFIG/${CONFIG}-${CASE}-MEAN/$year ; fi;
if [ ! -d /scratch/$USER/${CONFIG}-${CASE}-MEAN/$year ] ; then mkdir /scratch/$USER/${CONFIG}-${CASE}-MEAN/$year ; fi;
cp -f $year/*.nc /data/${USER}/$CONFIG/${CONFIG}-${CASE}-MEAN/$year/.
mv -f $year/*.nc /scratch/${USER}/${CONFIG}-${CASE}-MEAN/$year/. 
fi
done
echo "end"
