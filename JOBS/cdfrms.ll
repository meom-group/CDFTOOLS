#!/bin/csh
# @ cpu_limit  = 3500
# @ data_limit = 1gb
# Nom du travail LoadLeveler
# @ job_name   = cdfrms
# Fichier de sortie standard du travail       
# @ output     = $(job_name).$(jobid)
# Fichier de sortie d'erreur du travail
# @ error      =  $(job_name).$(jobid)
# @ queue                   


set echo


set CONFIG=ORCA05

set YEAR=0008-0010
#

set CDFTOOLS=~rcli002/CDFTOOLS-2.0
#foreach CASE ( G22 G23 G03 )
foreach CASE ( G32  )
  set CONFCASE=${CONFIG}-${CASE}

  cd $TMPDIR
  cp $CDFTOOLS/att.txt .
  cp $CDFTOOLS/cdfrmsssh ./
  chmod 755 cdfrmsssh
  cp $CDFTOOLS/cdfstdevw ./
  chmod 755 cdfstdevw

   mfget $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_gridT.nc ./
   mfget $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_gridT2.nc ./

  ./cdfrmsssh  ${CONFCASE}_y${YEAR}_gridT.nc ${CONFCASE}_y${YEAR}_gridT2.nc

  mfput rms.nc $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_RMSSSH.nc
  \rm rms.nc

   mfget $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_gridW.nc ./
   mfget $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_gridW2.nc ./

  ./cdfstdevw  ${CONFCASE}_y${YEAR}_gridW.nc ${CONFCASE}_y${YEAR}_gridW2.nc

  mfput rmsw.nc $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_STDEVW.nc
  \rm rmsw.nc


end

