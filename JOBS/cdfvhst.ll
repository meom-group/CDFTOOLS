#!/bin/csh
# @ cpu_limit  = 3500
# @ data_limit = 1gb
# Nom du travail LoadLeveler
# @ job_name   = cdfvhst
# Fichier de sortie standard du travail       
# @ output     = $(job_name).$(jobid)
# Fichier de sortie d'erreur du travail
# @ error      =  $(job_name).$(jobid)
# @ queue                   


set echo


set CONFIG=ORCA025
set CASE=G30

set YEAR=0008-0010
#

set CONFCASE=${CONFIG}-${CASE}
set CDFTOOLS=~rcli002/CDFTOOLS-2.0


cd $TMPDIR
cp $CDFTOOLS/cdfvhst .
cp $CDFTOOLS/att.txt .
mfget ${CONFIG}/${CONFIG}-I/ORCA025-G30_mesh_hgr.nc mesh_hgr.nc
mfget ${CONFIG}/${CONFIG}-I/ORCA025-G30_mesh_zgr.nc mesh_zgr.nc
mfget ${CONFIG}/${CONFIG}-I/ORCA025-G30_byte_mask.nc mask.nc
mfget ${CONFIG}/${CONFIG}-I/new_maskglo.nc .

foreach year ( 0008-0010 )
set CONFCASE=${CONFIG}-${CASE}
rsh gaya mkdir ${CONFIG}/${CONFCASE}-DIAGS/

mfget ${CONFIG}/${CONFCASE}-MEAN/$year/${CONFCASE}_y${year}_VT.nc .

./cdfvhst ${CONFCASE}_y${year}_VT.nc
mv zonal_heat_trp.dat ${CONFCASE}_y${year}_heattrp.dat
mv zonal_salt_trp.dat ${CONFCASE}_y${year}_salttrp.dat
mv trp.nc  ${CONFCASE}_y${year}_trp.nc

mfput  ${CONFCASE}_y${year}_trp.nc ${CONFIG}/${CONFCASE}-DIAGS/

mfput ${CONFCASE}_y${year}_heattrp.dat ${CONFIG}/${CONFCASE}-DIAGS/
cp ${CONFCASE}_y${year}_heattrp.dat $WORKDIR/${CONFCASE}-S/

mfput ${CONFCASE}_y${year}_salttrp.dat ${CONFIG}/${CONFCASE}-DIAGS/
cp ${CONFCASE}_y${year}_salttrp.dat $WORKDIR/${CONFCASE}-S/

\rm -f *.dat ${CONFCASE}_y${year}_VT.nc  ${CONFCASE}_y${year}_trp.nc


end

