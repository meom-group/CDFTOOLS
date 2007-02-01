#!/bin/csh
# @ cpu_limit  = 3500
# @ data_limit = 1gb
# Nom du travail LoadLeveler
# @ job_name   = cdfmxl
# Fichier de sortie standard du travail       
# @ output     = $(job_name).$(jobid)
# Fichier de sortie d'erreur du travail
# @ error      =  $(job_name).$(jobid)
# @ queue                   


set echo


set CONFIG=ORCA025
set CASE=G42
set YEAR=0005
set MESH_ZGR=ORCA025-G42_mesh_zgr.nc
#
set CDFTOOLS=~rcli002/CDFTOOLS-2.0

  set CONFCASE=${CONFIG}-${CASE}

  cd $TMPDIR
  cp $CDFTOOLS/att.txt .
  cp $CDFTOOLS/cdfmxl .
  chmod 755 cdfsig0
  mfget $CONFIG/${CONFIG}-I/$MESH_ZGR mesh_zgr.nc

  rsh gaya mkdir $CONFIG/${CONFCASE}-DIAGS/
  rsh gaya mkdir $CONFIG/${CONFCASE}-DIAGS/$YEAR

foreach f ( `rsh gaya ls $CONFIG/${CONFCASE}-S/$YEAR/${CONFCASE}_y${YEAR}\*gridT.nc ` )
   mfget $f ./
   set g=`basename $f | sed -e 's/gridT/MXL/' `

  ./cdfmxl `basename $f`

  mfput mxl.nc $CONFIG/${CONFCASE}-DIAGS/$YEAR/$g
  \rm `basename $f` mxl.nc

end

