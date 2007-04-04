#!/bin/ksh
# @ cpu_limit  = 3600
# @ data_limit = 1gb
# Nom du travail LoadLeveler
# @ job_name   = vt-YYYY
# Fichier de sortie standard du travail       
# @ output     = $(job_name).$(jobid)
# Fichier de sortie d'erreur du travail
# @ error      =  $(job_name).$(jobid)
# @ queue                   

# Example of script for time-average computation for VT files
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
putmonth() { mfput vt.nc $MDIR/${CONFCASE}_y${YEAR}m${1}_VT.nc ; \rm ${CONFCASE}_y${YEAR}m${month}d??_grid[UVT].nc ; }

putannual() { mfput cdfmoy_annual.nc $MDIR/${CONFCASE}_y${YEAR}_ANNUAL_VT.nc ; \rm cdfmoy_annual.nc ; }
#
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
   getmonth $month gridT
   getmonth $month gridU
   getmonth $month gridV

   list=''
   for f in ${CONFCASE}_y${YEAR}m${month}d??_gridT.nc ; do
     tag=$( echo $f | awk -F_ '{print $2}' )
     list="$list $tag"
   done

   $CDFTOOLS/cdfvT $CONFCASE $list
   putmonth $month
   mv vt.nc MONTHLY/${CONFCASE}_y${YEAR}m${$month}_VT.nc
 done

 # annual mean  (uses a ponderation to compute the exact annual mean ).
 cd $TMPDIR/MONTHLY
 $CDFTOOLS/cdfmoy_annual ${CONFCASE}_y${YEAR}m??_VT.nc
 putannual
 cd $TMPDIR

done
