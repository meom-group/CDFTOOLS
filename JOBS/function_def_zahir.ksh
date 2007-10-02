#!/bin/ksh
# function_def.ksh file for zahir and gaya. To be used with cdfmoy and cdfvT jobs


#  $Rev$
#  $Date$
#  $Id$


# chkdirg  path : check existence of directory path  on (remote) archiving machine. If it does not exist, create it.
chkdirg() { rsh gaya -l $REMOTE_USER " if [ ! -d $1 ] ; then mkdir $1 ; fi " ; }

# getmonth  mm type : retrieve all 5 days average files for month mm and grid type 'type', corresponding to current year, current confcase.A
#                     ex: getmonth 04 gridU  : retrieve all april files for gridU
getmonth() { for f  in $( rsh gaya -l $REMOTE_USER ls $SDIR/${CONFCASE}_y${YEAR}m${1}\*_$2.nc )  ; do
              mfget -u $REMOTE_USER $f ./
             done  ; }
# putmonth mm type : write back monthly mean for month mm type 'type' on remote machine in -MEAN/YEAR/ directory.
#                    also move the localfile to local MONTHLY dir for further annual mean computing

putmonth() { mfput -u $REMOTE_USER cdfmoy.nc $MDIR/${CONFCASE}_y${YEAR}m${1}_$2.nc ;\
             mv cdfmoy.nc MONTHLY/${CONFCASE}_y${YEAR}m${1}_$2.nc ; \rm ${CONFCASE}_y${YEAR}m${1}d??_$2.nc ; }

# putmonth2 mm type : write back monthly quadratic  mean for month mm type 'type' on remote machine in -MEAN/YEAR/ directory.
#                    also move the localfile to local MONTHLY dir for further annual mean computing
putmonth2() { mfput -u $REMOTE_USER cdfmoy2.nc $MDIR/${CONFCASE}_y${YEAR}m${1}_$2.nc ; \
             mv cdfmoy2.nc MONTHLY/${CONFCASE}_y${YEAR}m${1}_$2.nc ; }

# putannual type : write annual MEAN to remote -MEAN dir, in the corresponding year. Clean local files
putannual() { mfput -u $REMOTE_USER cdfmoy_annual.nc $MDIR/${CONFCASE}_y${YEAR}_ANNUAL_$1.nc ; \rm cdfmoy_annual.nc ;}

# putvtmonth mm : write back monthly mean for month mm type 'VT' on remote machine in -MEAN/YEAR/ directory.
#                    also move the localfile to local MONTHLY dir for further annual mean computing
putvtmonth() { mfput -u $REMOTE_USER vt.nc $MDIR/${CONFCASE}_y${YEAR}m${1}_VT.nc ; \
             mv vt.nc MONTHLY/${CONFCASE}_y${YEAR}m${1}_VT.nc ; \rm ${CONFCASE}_y${YEAR}m${month}d??_grid[UVT].nc ; }

# putvtannual type : write annual MEAN to remote -MEAN dir, in the corresponding year. Clean local files
putvtannual() { mfput -u $REMOTE_USER cdfmoy_annual.nc $MDIR/${CONFCASE}_y${YEAR}_ANNUAL_VT.nc ; \rm cdfmoy_annual.nc ; }
#

