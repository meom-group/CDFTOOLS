#!/bin/csh
# @ cpu_limit  = 3500
# @ data_limit = 1gb
# Nom du travail LoadLeveler
# @ job_name   = cdfeke
# Fichier de sortie standard du travail       
# @ output     = $(job_name).$(jobid)
# Fichier de sortie d'erreur du travail
# @ error      =  $(job_name).$(jobid)
# @ queue                   


set echo


set CONFIG=ORCA025

set YEAR=0005
#
set CASE=G32
set CDFTOOLS=~rcli002/CDFTOOLS-2.0
#foreach CASE ( G22 G23 G03 )
foreach YEAR (   0010 )
  set CONFCASE=${CONFIG}-${CASE}

  cd $TMPDIR
  cp $CDFTOOLS/att.txt .
  cp $CDFTOOLS/cdfeke .
  chmod 755 cdfeke

  foreach f ( `rsh gaya ls $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_ANNUAL_grid\[UV\]\*nc ` )
   mfget $f ./
  end
   mfget $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_ANNUAL_gridT2.nc  ./

  ./cdfeke ${CONFCASE}_y${YEAR}_ANNUAL_gridU.nc ${CONFCASE}_y${YEAR}_ANNUAL_gridU2.nc ${CONFCASE}_y${YEAR}_ANNUAL_gridV.nc ${CONFCASE}_y${YEAR}_ANNUAL_gridV2.nc ${CONFCASE}_y${YEAR}_ANNUAL_gridT2.nc

  mv eke.nc ${CONFCASE}_y${YEAR}_EKE.nc
  mfput ${CONFCASE}_y${YEAR}_EKE.nc $CONFIG/${CONFCASE}-MEAN/$YEAR/

  \rm *ANNUAL* *EKE*
end

