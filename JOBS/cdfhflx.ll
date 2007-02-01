#!/bin/csh
# @ cpu_limit  = 3500
# @ data_limit = 1gb
# Nom du travail LoadLeveler
# @ job_name   = cdfhflx
# Fichier de sortie standard du travail       
# @ output     = $(job_name).$(jobid)
# Fichier de sortie d'erreur du travail
# @ error      =  $(job_name).$(jobid)
# @ queue                   


set echo

set CONFIG=ORCA025
set CASE=G50

set YEARS=(1949 1950 1951 1952 1953 1954 1955 1956 1957 1958 1959  )
set MESH_MASK_ID='ORCA025-G50'


set CONFCASE=${CONFIG}-${CASE}

set CDFTOOLS=~rcli002/CDFTOOLS-2.0
cd $TMPDIR
cp $CDFTOOLS/cdfhflx .
mfget ${CONFIG}/${CONFIG}-I/new_maskglo.nc new_maskglo.nc
mfget ${CONFIG}/${CONFIG}-I/${MESH_MASK_ID}_byte_mask.nc mask.nc
mfget ${CONFIG}/${CONFIG}-I/${MESH_MASK_ID}_mesh_hgr.nc mesh_hgr.nc
mfget ${CONFIG}/${CONFIG}-I/${MESH_MASK_ID}_mesh_zgr.nc mesh_zgr.nc

foreach year ( $YEARS )
set CONFCASE=${CONFIG}-${CASE}
rsh gaya mkdir ${CONFIG}/${CONFCASE}-DIAGS/

mfget ${CONFIG}/${CONFCASE}-MEAN/$year/${CONFCASE}_y${year}_ANNUAL_T.nc .

./cdfhflx  ${CONFCASE}_y${year}_ANNUAL_T.nc
mv hflx.out ${CONFCASE}_y${year}_hflx.dat

mfput ${CONFCASE}_y${year}_hflx.dat ${CONFIG}/${CONFCASE}-DIAGS/

\rm -f *.dat ${CONFCASE}_y${year}_ANNUAL_T.nc

end
