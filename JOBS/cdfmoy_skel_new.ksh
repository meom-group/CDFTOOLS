#!/bin/ksh
# @ cpu_limit  = 7200
# @ data_limit = 1gb
# Nom du travail LoadLeveler
# @ job_name   = moy-YYYY
# Fichier de sortie standard du travail       
# @ output     = $(job_name).$(jobid)
# Fichier de sortie d'erreur du travail
# @ error      =  $(job_name).$(jobid)
# @ queue                   

# Example of script for time-average computation
# This script produce monthly mean and annual mean. The annual mean is now a weighted mean of the monthly means.
# This works for DRAKKAR output ( 5day average, no leap year).
# $Rev$
# $Date$
# $Id$

set -x

CONFIG=NATL025
CASE=GH01

YEARS=YYYY
YEARE=YYYE
#
YEARLST=""
y=$YEARS
while (( $y <= $YEARE )) ; do 
  YEARLST="$YEARLST $y "
  y=$(( y + 1 ))
done
# define some usefull functions
chkdirg() { rsh gaya " if [ ! -d $1 ] ; then mkdir $1 ; fi " ; }
getmonth() { for f  in $( rsh gaya ls $SDIR/${CONFCASE}_y${YEAR}m${1}\*_$2.nc )  ; do
              mfget $f ./
             done  ; }
putmonth() { mfput cdfmoy.nc $MDIR/${CONFCASE}_y${YEAR}m${1}_$2.nc ; \rm ${CONFCASE}_y${YEAR}m${1}d??_$2.nc ; }
putmonth2() { mfput cdfmoy2.nc $MDIR/${CONFCASE}_y${YEAR}m${1}_$2.nc ; }

putannual() { mfput cdfmoy_annual.nc $MDIR/${CONFCASE}_y${YEAR}_ANNUAL_$1.nc ; \rm cdfmoy_annual.nc ;}


CONFCASE=${CONFIG}-${CASE}
CDFTOOLS=~rcli002/CDFTOOLS-2.1/
cd $TMPDIR
mkdir MONTHLY

for YEAR in $YEARLST ; do 
   SDIR=${CONFIG}/${CONFCASE}-S/$YEAR
   MDIR=${CONFIG}/${CONFCASE}-MEAN/$YEAR
   chkdirg $MDIR

 # Monthly mean
 #
 for month in 01 02 03 04 05 06 07 08 09 10 11 12  ; do 
   for grid in gridT gridU gridV gridW icemod  ; do 
    getmonth $month $grid
    $CDFTOOLS/cdfmoy ${CONFCASE}_y${YEAR}m${month}d??_$grid.nc
    case $grid in 
      icemod) putmonth $month $grid ;;
           *) putmonth $month $grid ;
              putmonth2 $month ${grid}2 ;;
    esac
   done
 done
 \rm -f cdfmoy*nc
  # ANNUAL
  # at this point all monthly averages are in MONTHLY
  cd MONTHLY
  for grid in gridT gridU gridV gridW icemod gridT2 gridU2 gridV2 gridW2 ; do
    $CDFTOOLS/cdfmoy_annual ${CONFCASE}_y${YEAR}m??_$grid.nc
    putannual  $grid ;
    \rm ${CONFCASE}_y${YEAR}m??_$grid.nc
  done
  cd $TMPDIR
done
