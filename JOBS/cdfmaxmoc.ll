#!/bin/csh
# @ cpu_limit  = 3600
# @ data_limit = 1gb
# Nom du travail LoadLeveler
# @ job_name   = maxmoc
# Fichier de sortie standard du travail
# @ output     = $(job_name).$(jobid)
# Fichier de sortie d'erreur du travail
# @ error      =  $(job_name).$(jobid)
# @ queue


set CONFIG=ORCA025
set CASE=G32

set year=0001
set yrfin=0010

#-----------------------------------------------------------------------------
set CONFCASE=${CONFIG}-${CASE}
set CDFTOOLS=$HOME/CDFTOOLS-2.0/

cd $TMPDIR

cp $CDFTOOLS/cdfmaxmoc .
while ( $year < $yrfin ) 
set year=`printf "%04d" $year


# MAX and MIN of MOC
#-------------------
set f=${CONFCASE}_y${year}_MOC.nc
mfget ${CONFIG}/${CONFIG}-${CASE}-DIAGS/$f ./
set outfile=${CONFCASE}_y${year}_minmaxmoc.txt
    echo $year > $outfile
# GLO
printf "%s" 'Glo ' >>  $outfile ; ./cdfmaxmoc $f glo 20 60 500 2000 | grep Maximum >> $outfile
printf "%s" 'Glo ' >>  $outfile ; ./cdfmaxmoc $f glo -40 30 2000 5500 | grep Minimum >> $outfile
# ATL
printf "%s" 'Atl ' >>  $outfile ; ./cdfmaxmoc $f atl 0 60 500 2000 | grep Maximum >> $outfile
printf "%s" 'Atl ' >>  $outfile ; ./cdfmaxmoc $f atl -20 40 2000 5500 | grep Minimum  >> $outfile
#INP
printf "%s" 'Inp ' >>  $outfile ; ./cdfmaxmoc $f inp 15 50 100 1000 | grep Minimum >> $outfile
printf "%s" 'Inp ' >>  $outfile ; ./cdfmaxmoc $f inp -30 20 1000 5500  | grep Minimum >> $outfile
#AUS
printf "%s" 'Aus ' >>  $outfile ; ./cdfmaxmoc $f glo -70 0 0 2000   | grep Maximum >> $outfile
printf "%s" 'Aus ' >>  $outfile ; ./cdfmaxmoc $f glo -70 0 2000 5500  | grep Minimum >> $outfile

mfput $outfile ${CONFIG}/${CONFCASE}-DIAGS/

# Max and Min of MOC at some specific latitudes
set f=${CONFCASE}_y${year}_MOC.nc
set outfile=${CONFIG}-${CASE}_y${year}_maxmoc40.txt
  \rm -f $outfile

    echo $year > $outfile
# GLO  MAX at 40 N and 30S
printf "%s" 'Glo ' >>  $outfile ; cdfmaxmoc $f glo 40 40 500 2000 | grep Maximum >> $outfile
printf "%s" 'Glo ' >>  $outfile ; cdfmaxmoc $f glo -30 -30 500  5500 | grep Maximum >> $outfile
# ATL  MAX at 40N and 30S
printf "%s" 'Atl ' >>  $outfile ; cdfmaxmoc $f atl 40 40 500 2000 | grep Maximum >> $outfile
printf "%s" 'Atl ' >>  $outfile ; cdfmaxmoc $f atl -30 -30  500 5000 | grep Maximum >> $outfile
#INP  Min at 30 S
printf "%s" 'Inp ' >>  $outfile ; cdfmaxmoc $f inp -30 -30 1000 5500  | grep Minimum >> $outfile
#AUS  MAX at 50 S
printf "%s" 'Aus ' >>  $outfile ; cdfmaxmoc $f glo -50 -50 0 2000   | grep Maximum >> $outfile

mfput $outfile ${CONFIG}/${CONFIG}-${CASE}-DIAGS/

@ year ++
end

