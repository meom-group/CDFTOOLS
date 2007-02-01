#!/bin/csh
# @ cpu_limit  = 3500
# @ data_limit = 1gb
# Nom du travail LoadLeveler
# @ job_name   = cdfmhst-full
# Fichier de sortie standard du travail       
# @ output     = $(job_name).$(jobid)
# Fichier de sortie d'erreur du travail
# @ error      =  $(job_name).$(jobid)
# @ queue                   


set echo


set CONFIG=ORCA025
set CASE=G03

#

set CONFCASE=${CONFIG}-${CASE}

set CDFTOOLS=~rcli002/CDFTOOLS-2.0
cd $TMPDIR
cp $CDFTOOLS/cdfmhst-full .
mfget ${CONFIG}/${CONFIG}-I/ORCA025-G30_mesh_hgr.nc mesh_hgr.nc
mfget ${CONFIG}/${CONFIG}-I/ORCA025-G30_mesh_zgr.nc mesh_zgr.nc
mfget ${CONFIG}/${CONFIG}-I/ORCA025-G30_byte_mask.nc mask.nc
mfget ${CONFIG}/${CONFIG}-I/new_maskglo.nc .

foreach year ( 0001 0002 0003 0004 )
set CONFCASE=${CONFIG}-${CASE}
rsh gaya mkdir ${CONFIG}/${CONFCASE}-DIAGS/

mfget ${CONFIG}/${CONFCASE}-MEAN/$year/${CONFCASE}_y${year}_ANNUAL_VT.nc .

./cdfmhst-full  ${CONFCASE}_y${year}_ANNUAL_VT.nc
mv zonal_heat_trp.dat ${CONFCASE}_y${year}_heattrp.dat
mv zonal_salt_trp.dat ${CONFCASE}_y${year}_salttrp.dat

mfput ${CONFCASE}_y${year}_heattrp.dat ${CONFIG}/${CONFCASE}-DIAGS/
cp ${CONFCASE}_y${year}_heattrp.dat $WORKDIR/${CONFCASE}-S/

mfput ${CONFCASE}_y${year}_salttrp.dat ${CONFIG}/${CONFCASE}-DIAGS/
cp ${CONFCASE}_y${year}_salttrp.dat $WORKDIR/${CONFCASE}-S/

\rm -f *.dat ${CONFCASE}_y${year}_ANNUAL_VT.nc

end

