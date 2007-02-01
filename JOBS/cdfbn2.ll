#!/bin/csh
# @ cpu_limit  = 3500
# @ data_limit = 1gb
# Nom du travail LoadLeveler
# @ job_name   = cdfbn2
# Fichier de sortie standard du travail       
# @ output     = $(job_name).$(jobid)
# Fichier de sortie d'erreur du travail
# @ error      =  $(job_name).$(jobid)
# @ queue                   


set echo


set CONFIG=ORCA025

set YEAR=0010
#
set CDFTOOLS=~rcli002/CDFTOOLS-2.0

set CASE=G22
  set CONFCASE=${CONFIG}-${CASE}

  cd $TMPDIR
  cp $CDFTOOLS/att.txt .
  cp $CDFTOOLS/cdfbn2 .
  chmod 755 cdfbn2
  mfget $CONFIG/${CONFIG}-I/${CONFIG}_PS_mesh_hgr.nc mesh_hgr.nc
  mfget $CONFIG/${CONFIG}-I/${CONFIG}_PS_mesh_zgr.nc mesh_zgr.nc

  rsh gaya mkdir $CONFIG/${CONFCASE}-DIAGS/
  rsh gaya mkdir $CONFIG/${CONFCASE}-DIAGS/$YEAR

foreach f ( `rsh gaya ls $CONFIG/${CONFCASE}-S/$YEAR/${CONFCASE}_y${YEAR}\*gridT.nc ` )
   mfget $f ./
   set g=`basename $f | sed -e 's/gridT/N2/' `

  ./cdfbn2 `basename $f`

  mfput bn2.nc $CONFIG/${CONFCASE}-DIAGS/$YEAR/$g
  \rm `basename $f` bn2.nc

end

