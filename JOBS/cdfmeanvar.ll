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

#set echo


set CONFIG=ORCA2
set CASE=G70

set YEARS=(1990 1991 1992 1993 1994 1995 1996 1997 1998 )
#set YEARS=(1990 )
set MESH_MASK_ID='ORCA2-G70'

set usergaya=/u/rech/cli/rcli544

#

set CONFCASE=${CONFIG}-${CASE}

set CDFTOOLS=~rcli002/CDFTOOLS-2.1


cd $TMPDIR
mfget $usergaya/${CONFIG}/${CONFIG}-I/new_maskglo.nc 
mfget $usergaya/${CONFIG}/${CONFIG}-I/${MESH_MASK_ID}_byte_mask.nc mask.nc
mfget $usergaya/${CONFIG}/${CONFIG}-I/${MESH_MASK_ID}_mesh_hgr.nc mesh_hgr.nc
mfget $usergaya/${CONFIG}/${CONFIG}-I/${MESH_MASK_ID}_mesh_zgr.nc mesh_zgr.nc


set CONFCASE=${CONFIG}-${CASE}
rsh gaya mkdir ${CONFIG}
rsh gaya mkdir ${CONFIG}/${CONFCASE}-DIAGS/

foreach year ( $YEARS  )
     echo $year > ${CONFCASE}_y${year}_meanvar.txt
    foreach vfile ( ` rsh  gaya " ls $usergaya/${CONFIG}/${CONFCASE}-S/$year/${CONFCASE}_y${year}m??d??_gridV.nc"  ` )
       echo $vfile
       mfget $vfile ./
       set vf=`basename $vfile `
       ncdump -v  time_counter $vf | grep time | tail -1 | awk '{ printf "%f ", $3/86400. }' >>  ${CONFCASE}_y${year}_meanvar.txt

       cdfmeanvar  $vf vomecrty V | grep over | awk '{ printf "%s ", $NF}'  >> ${CONFCASE}_y${year}_meanvar.txt
       printf "\n" >> ${CONFCASE}_y${year}_meanvar.txt

       mfput  ${CONFCASE}_y${year}_meanvar.txt  ${CONFIG}/${CONFCASE}-DIAGS/
    \rm $vf
    end

end
