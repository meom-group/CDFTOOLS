#!/bin/csh
# @ cpu_limit  = 3500
# @ data_limit = 1gb
# Nom du travail LoadLeveler
# @ job_name   = cdfpsi
# Fichier de sortie standard du travail       
# @ output     = $(job_name).$(jobid)
# Fichier de sortie d'erreur du travail
# @ error      =  $(job_name).$(jobid)
# @ queue                   


set echo

set CONFIG=ORCA025
set CASE=G50

set YEARS=(1948 1949 1950 1951 )
set MESH_MASK_ID='ORCA025-G45b'  


#
set CONFCASE=${CONFIG}-${CASE}
set CDFTOOLS=~rcli002/CDFTOOLS-2.0

cd $TMPDIR
cp $CDFTOOLS/att.txt .
cp $CDFTOOLS/cdfpsi .
chmod 755 cdfpsi

mfget ${CONFIG}/${CONFIG}-I/${MESH_MASK_ID}_byte_mask.nc mask.nc
mfget ${CONFIG}/${CONFIG}-I/${MESH_MASK_ID}_mesh_hgr.nc mesh_hgr.nc
mfget ${CONFIG}/${CONFIG}-I/${MESH_MASK_ID}_mesh_zgr.nc mesh_zgr.nc

foreach YEAR ( $YEARS )


  foreach f ( `rsh gaya ls $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_ANNUAL_grid\[UV\].nc ` )
   mfget $f ./
  end

  ./cdfpsi ${CONFCASE}_y${YEAR}_ANNUAL_gridU.nc ${CONFCASE}_y${YEAR}_ANNUAL_gridV.nc 

  mv psi.nc ${CONFCASE}_y${YEAR}_PSI.nc
  mfput ${CONFCASE}_y${YEAR}_PSI.nc $CONFIG/${CONFCASE}-MEAN/$YEAR/

  \rm *ANNUAL* *PSI*
end

