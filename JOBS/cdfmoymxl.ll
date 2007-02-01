#!/bin/csh
# @ cpu_limit  = 4200
# @ data_limit = 1gb
# Nom du travail LoadLeveler
# @ job_name   = cdfmoymxl
# Fichier de sortie standard du travail       
# @ output     = $(job_name).$(jobid)
# Fichier de sortie d'erreur du travail
# @ error      =  $(job_name).$(jobid)
# @ queue                   


set echo


set CONFIG=ORCA025
set CASE=G43b

set YEAR=0004
#

set CONFCASE=${CONFIG}-${CASE}
set CDFTOOLS=~rcli002/CDFTOOLS-2.0/


cd $TMPDIR
mkdir MONTHLY

cp $CDFTOOLS/att.txt .
rsh gaya mkdir ${CONFIG}/${CONFCASE}-DIAGS/$YEAR/

#goto annual
#goto quarterly
# Monthly mean
#
foreach month (01 02 03 04 05 06 07 08 09 10 11 12 )
  foreach f ( `rsh gaya ls ${CONFIG}/${CONFCASE}-DIAGS/$YEAR/${CONFCASE}_y${YEAR}m${month}\*_MXL.nc `)
   mfget $f ./
  end 

  set list=''
  foreach f ( ${CONFCASE}_y${YEAR}m${month}d??_MXL.nc )
    set list=($list $f )
  end
   $CDFTOOLS/cdfmoy $list
   mfput cdfmoy.nc ${CONFIG}/${CONFCASE}-DIAGS/$YEAR/${CONFCASE}_y${YEAR}m${month}_MXL.nc
   \rm $list cdfmoy.nc cdfmoy2.nc
end
