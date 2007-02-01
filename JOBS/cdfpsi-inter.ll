#!/bin/csh
# @ cpu_limit  = 3500
# @ data_limit = 1gb
# Nom du travail LoadLeveler
# @ job_name   = cdfpsi-inter
# Fichier de sortie standard du travail       
# @ output     = $(job_name).$(jobid)
# Fichier de sortie d'erreur du travail
# @ error      =  $(job_name).$(jobid)
# @ queue                   


set echo

set CONFIG=ORCA025
set CASE=G42
set YEAR=0008-0010

#
  set CDFTOOLS=~rcli002/CDFTOOLS-2.0

cd $TMPDIR

mfget $CONFIG/${CONFIG}-I/${CONFIG}-${CASE}_mesh_hgr.nc mesh_hgr.nc
mfget $CONFIG/${CONFIG}-I/${CONFIG}-${CASE}_mesh_zgr.nc mesh_zgr.nc
mfget $CONFIG/${CONFIG}-I/${CONFIG}-${CASE}_byte_mask.nc mask.nc

  set CONFCASE=${CONFIG}-${CASE}

  cp $CDFTOOLS/att.txt .
  cp $CDFTOOLS/cdfpsi .
  chmod 755 cdfpsi

  foreach f ( `rsh gaya ls $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_grid\[UV\].nc ` )
   mfget $f ./
  end

  ./cdfpsi ${CONFCASE}_y${YEAR}_gridU.nc ${CONFCASE}_y${YEAR}_gridV.nc 

  mv psi.nc ${CONFCASE}_y${YEAR}_PSI.nc
  mfput ${CONFCASE}_y${YEAR}_PSI.nc $CONFIG/${CONFCASE}-MEAN/$YEAR/

  \rm  *PSI*

