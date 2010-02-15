#!/bin/csh
# @ cpu_limit  = 3500
# @ data_limit = 1gb
# Nom du travail LoadLeveler
# @ job_name   = cdfmoy-ets
# Fichier de sortie standard du travail       
# @ output     = $(job_name).$(jobid)
# Fichier de sortie d'erreur du travail
# @ error      =  $(job_name).$(jobid)
# @ queue                   


set echo


set CONFIG=ORCA025
set CASE=G22

set YEAR=0010
#

set CONFCASE=${CONFIG}-${CASE}

set CDFTOOLS=~rcli002/CDFTOOLS-2.0


cd $TMPDIR
mkdir MONTHLY

cp $CDFTOOLS/att.txt .
rsh gaya mkdir ${CONFIG}/${CONFCASE}-MEAN/$YEAR/

foreach month (01 02 03 04 05 06 07 08 09 10 11 12 )
  foreach f ( `rsh gaya ls ${CONFIG}/${CONFCASE}-DIAGS/$YEAR/${CONFCASE}_y${YEAR}m${month}\*_ETS.nc `)
   mfget $f ./
  end 

  set list=''
  foreach f ( ${CONFCASE}_y${YEAR}m${month}d??_ETS.nc )
    set list=($list $f )
  end
   $CDFTOOLS/cdfmoy_sp $list
   mfput cdfmoy.nc ${CONFIG}/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}m${month}_ETS.nc
   \rm $list cdfmoy.nc cdfmoy2.nc

end
