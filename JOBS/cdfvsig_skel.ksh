#!/bin/ksh
# @ cpu_limit  = 3600
# @ data_limit = 1gb
# Nom du travail LoadLeveler
# @ job_name   = vsig-YYYY
# Fichier de sortie standard du travail       
# @ output     = $(job_name).$(jobid)
# Fichier de sortie d'erreur du travail
# @ error      =  $(job_name).$(jobid)
# @ queue                   

# Example of script for time-average computation for U x sigma0 terms
# This script produce  annual mean.
# This works for DRAKKAR output ( 5day average, no leap year).

# $Rev$
# $Date$
# $Id$


set -x

CONFIG=ORCA025
CASE=G70
MESH_MASK_ID=ORCA025-G70

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
putannual() { mfput $1 $MDIR/$1 ; }

getmask () { mfget ${CONFIG}/${CONFIG}-I/${MESH_MASK_ID}_byte_mask.nc mask.nc ; }
#
CONFCASE=${CONFIG}-${CASE}
CDFTOOLS=~rcli002/CDFTOOLS-2.1/
cd $TMPDIR
mkdir MONTHLY
getmask

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
   getmonth $month gridW
 done
 
 # perform cdfvsig on the whole year

   list=''
   for f in ${CONFCASE}_y${YEAR}m??d??_gridT.nc ; do
     tag=$( echo $f | awk -F_ '{print $2}' )
     list="$list $tag"
   done

   $CDFTOOLS/cdfvsig $CONFCASE $list
   mv usig.nc ${CONFCASE}_y${YEAR}_USIG.nc
   mv vsig.nc ${CONFCASE}_y${YEAR}_VSIG.nc
   mv wsig.nc ${CONFCASE}_y${YEAR}_WSIG.nc

   putannual ${CONFCASE}_y${YEAR}_USIG.nc
   putannual ${CONFCASE}_y${YEAR}_VSIG.nc
   putannual ${CONFCASE}_y${YEAR}_WSIG.nc

   \rm *gridT*
   \rm *gridU*
   \rm *gridV*
   \rm *gridW*
   \rm *SIG.nc

 cd $TMPDIR

done
