#!/bin/csh
# @ cpu_limit  = 3500
# @ data_limit = 1gb
# Nom du travail LoadLeveler
# @ job_name   = cdfmhst
# Fichier de sortie standard du travail       
# @ output     = $(job_name).$(jobid)
# Fichier de sortie d'erreur du travail
# @ error      =  $(job_name).$(jobid)
# @ queue                   


set echo

set CONFIG=ORCA025
set CASE=G70

set YEARS=(1958)
set MESH_MASK_ID='ORCA025-G70'


set CONFCASE=${CONFIG}-${CASE}

set CDFTOOLS=~rcli858/DEV/CDFTOOLS-2.1-DEV
cd $TMPDIR
cp $CDFTOOLS/cdfmhst .
mfget ../../cache1/rcli002/${CONFIG}/${CONFIG}-I/new_maskglo.nc new_maskglo.nc
mfget ../../cache1/rcli002/${CONFIG}/${CONFIG}-I/${MESH_MASK_ID}_byte_mask.nc mask.nc
mfget ../../cache1/rcli002/${CONFIG}/${CONFIG}-I/${MESH_MASK_ID}_mesh_hgr.nc mesh_hgr.nc
mfget ../../cache1/rcli002/${CONFIG}/${CONFIG}-I/${MESH_MASK_ID}_mesh_zgr.nc mesh_zgr.nc

foreach year ( $YEARS )
set CONFCASE=${CONFIG}-${CASE}
rsh gaya mkdir ${CONFIG}/${CONFCASE}-DIAGS/

mfget ../../cache1/rcli002/${CONFIG}/${CONFCASE}-MEAN/$year/${CONFCASE}_y${year}_ANNUAL_VT.nc .

./cdfmhst  ${CONFCASE}_y${year}_ANNUAL_VT.nc
mv zonal_heat_trp.dat ${CONFCASE}_y${year}_heattrp.dat
mv zonal_salt_trp.dat ${CONFCASE}_y${year}_salttrp.dat

mfput ${CONFCASE}_y${year}_heattrp.dat ${CONFIG}/${CONFCASE}-DIAGS/
cp ${CONFCASE}_y${year}_heattrp.dat $WORKDIR/${CONFCASE}-S/

mfput ${CONFCASE}_y${year}_salttrp.dat ${CONFIG}/${CONFCASE}-DIAGS/
cp ${CONFCASE}_y${year}_salttrp.dat $WORKDIR/${CONFCASE}-S/

\rm -f *.dat ${CONFCASE}_y${year}_ANNUAL_VT.nc

end
