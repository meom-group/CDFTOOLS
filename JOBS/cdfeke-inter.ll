#!/bin/csh
# @ cpu_limit  = 3500
# @ data_limit = 1gb
# Nom du travail LoadLeveler
# @ job_name   = EKE_int
# Fichier de sortie standard du travail       
# @ output     = $(job_name).$(jobid)
# Fichier de sortie d'erreur du travail
# @ error      =  $(job_name).$(jobid)
# @ queue                   


set echo


set CONFIG=ORCA025

set YEAR=0008-0010
#

set CDFTOOLS=~rcli002/CDFTOOLS-2.0
foreach CASE ( G42 )
  set CONFCASE=${CONFIG}-${CASE}

  cd $TMPDIR
  cp $CDFTOOLS/att.txt .
  cp $CDFTOOLS/cdfeke .
  chmod 755 cdfeke

  foreach f ( `rsh gaya ls $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_grid\[UV\]\*nc ` )
   mfget $f ./
  end
   mfget $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_gridT2.nc  ./

  ./cdfeke ${CONFCASE}_y${YEAR}_gridU.nc ${CONFCASE}_y${YEAR}_gridU2.nc ${CONFCASE}_y${YEAR}_gridV.nc ${CONFCASE}_y${YEAR}_gridV2.nc ${CONFCASE}_y${YEAR}_gridT2.nc

  mv eke.nc ${CONFCASE}_y${YEAR}_EKE.nc
  mfput ${CONFCASE}_y${YEAR}_EKE.nc $CONFIG/${CONFCASE}-MEAN/$YEAR/

  \rm *${YEAR}* *EKE*
end

