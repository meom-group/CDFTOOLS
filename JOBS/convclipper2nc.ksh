#!/bin/ksh
# @ cpu_limit  = 36000
# @ data_limit = 1gb
# Nom du travail LoadLeveler
# @ job_name   = convclipp
# Fichier de sortie standard du travail
# @ output     = $(job_name).$(jobid)
# Fichier de sortie d'erreur du travail
# @ error      =  $(job_name).$(jobid)
# @ queue


GAYA=/cache2/rost011
CONFIG=ATL6
CASE=V6

CONFCASE=${CONFIG}-${CASE}

MEANDIR=${GAYA}/${CONFIG}/${CONFCASE}-MEAN
MEANNCDIR=${GAYA}/${CONFIG}/${CONFCASE}NC-MEAN
IDIR=${GAYA}/${CONFIG}/${CONFIG}-V0-I/

cd $TMPDIR

y1=1980
y2=2000

y=$y1

mfget $IDIR/${CONFCASE}_mesh_hgr.nc  mesh_hgr.nc
mfget $IDIR/${CONFCASE}_mesh_zgr.nc  mesh_zgr.nc

while (( $y <= $y2 )) ; do
  # list 2D files in MEANDIR/y
   lst=$( rsh gaya ls $MEANDIR/$y/\*2D\*dimg )
    for f in $lst ; do
     g=$( basename $f )
     tag=$( echo $g | awk -F_ '{ print $3}' | sed -e 's/.dimg//' )
     u=$( echo $f | sed -e 's/_2D_/_U_/' )
     v=$( echo $f | sed -e 's/_2D_/_V_/' )
     t=$( echo $f | sed -e 's/_2D_/_T_/' )
     s=$( echo $f | sed -e 's/_2D_/_S_/' )
     ssh=$( echo $f | sed -e 's/_2D_/_SSH_/' )
     uu=$( echo $f | sed -e 's/_2D_/_UU_/' )
     vv=$( echo $f | sed -e 's/_2D_/_VV_/' )
     mfget $f ; mfget $u ; mfget $v ; mfget $t ; mfget $s ; mfget $ssh  ; mfget $uu ; mfget $vv
     cdfconvert $tag $CONFCASE
     rcp ${CONFCASE}_${tag}_*.nc rost011@gaya:$MEANNCDIR/
     \rm  ${CONFCASE}_${tag}_*.nc *.dimg
    done
  y=$(( y + 1 ))
done

