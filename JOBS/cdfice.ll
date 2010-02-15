#!/bin/csh
# @ cpu_limit  = 3500
# @ data_limit = 1gb
# Nom du travail LoadLeveler
# @ job_name   = cdfice
# Fichier de sortie standard du travail       
# @ output     = $(job_name).$(jobid)
# Fichier de sortie d'erreur du travail
# @ error      =  $(job_name).$(jobid)
# @ queue                   


set echo

set CONFIG=ORCA025
set CASE=G50

set YEARS=(1949 1950 1951 )
set MESH_MASK_ID='ORCA025-G45b'  


#
set CONFCASE=${CONFIG}-${CASE}
set CDFTOOLS=~rcli002/CDFTOOLS-2.0

cd $TMPDIR

mfget ${CONFIG}/${CONFIG}-I/${MESH_MASK_ID}_byte_mask.nc mask.nc
mfget ${CONFIG}/${CONFIG}-I/${MESH_MASK_ID}_mesh_hgr.nc mesh_hgr.nc

cp $CDFTOOLS/cdficediags .
chmod 755 cdficediags

foreach YEAR ( $YEARS )
# ice control for m02 m03   m08 m09

mfget $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}m02_icemod.nc ./
mfget $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}m03_icemod.nc ./
mfget $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}m08_icemod.nc ./
mfget $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}m09_icemod.nc ./

echo '###' $YEAR 02 > ${CONFCASE}_y${YEAR}_ice.txt
./cdficediags ${CONFCASE}_y${YEAR}m02_icemod.nc  >> ${CONFCASE}_y${YEAR}_ice.txt
echo '###' $YEAR 03 >> ${CONFCASE}_y${YEAR}_ice.txt
./cdficediags ${CONFCASE}_y${YEAR}m03_icemod.nc  >> ${CONFCASE}_y${YEAR}_ice.txt
echo '###' $YEAR 08 >> ${CONFCASE}_y${YEAR}_ice.txt
./cdficediags ${CONFCASE}_y${YEAR}m08_icemod.nc  >> ${CONFCASE}_y${YEAR}_ice.txt
echo '###' $YEAR 09 >> ${CONFCASE}_y${YEAR}_ice.txt
./cdficediags ${CONFCASE}_y${YEAR}m09_icemod.nc  >> ${CONFCASE}_y${YEAR}_ice.txt

mfput ${CONFCASE}_y${YEAR}_ice.txt ${CONFIG}/${CONFCASE}-DIAGS


end
