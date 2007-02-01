#!/bin/csh
# @ cpu_limit  = 3600
# @ data_limit = 1gb
# Nom du travail LoadLeveler
# @ job_name   = gibraltar
# Fichier de sortie standard du travail       
# @ output     = $(job_name).$(jobid)
# Fichier de sortie d'erreur du travail
# @ error      =  $(job_name).$(jobid)
# @ queue                   


set echo


set CONFIG=ORCA025

set CASE=G70
set MESH_MASK_ID='ORCA025-G70'  
#

set CDFTOOLS=~rcli002/CDFTOOLS-2.0
set CONFCASE=${CONFIG}-${CASE}

# limit for ORCA025
set GIB=(1094 1109 653 674 )


cd $TMPDIR

  cp $CDFTOOLS/att.txt .
  cp $CDFTOOLS/cdfmean .

mfget ${CONFIG}/${CONFIG}-I/${MESH_MASK_ID}_byte_mask.nc mask.nc
mfget ${CONFIG}/${CONFIG}-I/${MESH_MASK_ID}_mesh_hgr.nc mesh_hgr.nc
mfget ${CONFIG}/${CONFIG}-I/${MESH_MASK_ID}_mesh_zgr.nc mesh_zgr.nc

foreach YEAR ( 1958 1959 1960 1961 1962 1963 1964 1965 )
 mfget ${CONFIG}/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_ANNUAL_gridT.nc ./

 \rm -f ${CONFCASE}_y${YEAR}_TGIB.txt
     echo $YEAR > ${CONFCASE}_y${YEAR}_TGIB.txt
    ./cdfmean ${CONFCASE}_y${YEAR}_ANNUAL_gridT.nc  votemper T $GIB  0 0  >>  ${CONFCASE}_y${YEAR}_TGIB.txt
 \rm -f ${CONFCASE}_y${YEAR}_SGIB.txt
     echo $YEAR > ${CONFCASE}_y${YEAR}_SGIB.txt
    ./cdfmean ${CONFCASE}_y${YEAR}_ANNUAL_gridT.nc  vosaline T $GIB  0 0 >>  ${CONFCASE}_y${YEAR}_SGIB.txt

mfput ${CONFCASE}_y${YEAR}_TGIB.txt ${CONFIG}/${CONFCASE}-DIAGS
mfput ${CONFCASE}_y${YEAR}_SGIB.txt ${CONFIG}/${CONFCASE}-DIAGS

# clean the space
 \rm ${CONFCASE}_y${YEAR}_ANNUAL_gridT.nc

end

