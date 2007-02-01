#!/bin/csh
# @ cpu_limit  = 3500
# @ data_limit = 1gb
# Nom du travail LoadLeveler
# @ job_name   = cdftransportiz-full
# Fichier de sortie standard du travail       
# @ output     = $(job_name).$(jobid)
# Fichier de sortie d'erreur du travail
# @ error      =  $(job_name).$(jobid)
# @ queue                   


set echo


set CONFIG=ORCA025

set year=0008-0010
#


set CDFTOOLS=~rcli002/CDFTOOLS-2.0


cd $TMPDIR
cp $CDFTOOLS/cdftransportiz-full  .
cp $CDFTOOLS/att.txt .
cp $CDFTOOLS/section.dat .
mfget ${CONFIG}/${CONFIG}-I/ORCA025_PS_mesh_hgr.nc mesh_hgr.nc
mfget ${CONFIG}/${CONFIG}-I/ORCA025_PS_mesh_zgr.nc mesh_zgr.nc

foreach CASE (G03 G04 )

set CONFCASE=${CONFIG}-${CASE}
rsh gaya mkdir ${CONFIG}/${CONFCASE}-DIAGS/

mfget ${CONFIG}/${CONFCASE}-MEAN/$year/${CONFCASE}_y${year}_VT.nc .
mfget ${CONFIG}/${CONFCASE}-MEAN/$year/${CONFCASE}_y${year}_gridU.nc .
mfget ${CONFIG}/${CONFCASE}-MEAN/$year/${CONFCASE}_y${year}_gridV.nc .

./cdftransportiz-full  ${CONFCASE}_y${year}_VT.nc \
   ${CONFCASE}_y${year}_gridU.nc \
   ${CONFCASE}_y${year}_gridV.nc  1250 3500 < section.dat > ${CONFCASE}_sections.txt

 grep -v Give ${CONFCASE}_sections.txt > tmp
mv -f tmp ${CONFCASE}_sections.txt
mfput  ${CONFCASE}_sections.txt  ${CONFIG}/${CONFCASE}-DIAGS/


end


