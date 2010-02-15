#!/bin/csh
# @ cpu_limit  = 3500
# @ data_limit = 1gb
# Nom du travail LoadLeveler
# @ job_name   = cdfmoc
# Fichier de sortie standard du travail       
# @ output     = $(job_name).$(jobid)
# Fichier de sortie d'erreur du travail
# @ error      =  $(job_name).$(jobid)
# @ queue                   


set echo


set CONFIG=ORCA025
set CASE=G44

set YEARS=(0002 0003 0004 0005)
set MESH_MASK_ID='ORCA025-G44'

#

set CONFCASE=${CONFIG}-${CASE}

set CDFTOOLS=~rcli002/CDFTOOLS-2.0


cd $TMPDIR
cp $CDFTOOLS/cdfmoc .
cp $CDFTOOLS/att.txt .
mfget ${CONFIG}/${CONFIG}-I/new_maskglo.nc 
mfget ${CONFIG}/${CONFIG}-I/${MESH_MASK_ID}_byte_mask.nc mask.nc
mfget ${CONFIG}/${CONFIG}-I/${MESH_MASK_ID}_mesh_hgr.nc mesh_hgr.nc
mfget ${CONFIG}/${CONFIG}-I/${MESH_MASK_ID}_mesh_zgr.nc mesh_zgr.nc


set CONFCASE=${CONFIG}-${CASE}
rsh gaya mkdir ${CONFIG}/${CONFCASE}-DIAGS/

foreach year ( $YEARS  )

mfget ${CONFIG}/${CONFCASE}-MEAN/$year/${CONFCASE}_y${year}_ANNUAL_gridV.nc .

./cdfmoc ${CONFCASE}_y${year}_ANNUAL_gridV.nc
mv moc.nc  ${CONFCASE}_y${year}_MOC.nc
mfput  ${CONFCASE}_y${year}_MOC.nc  ${CONFIG}/${CONFCASE}-DIAGS/

end
