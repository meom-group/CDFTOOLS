# Ice Volume area and extent for all months: input file : icemod, and mesh_mask
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  if [ $ICEMONTH == 1 ] ; then
   # get icemod files
   m=1
   while (( $m <= 12 )) ; do
    mm=$( printf "%02d" $m )
    rapatrie  ${CONFCASE}_y${YEAR}m${mm}_icemod.nc $MEANY ${CONFCASE}_y${YEAR}m${mm}_icemod.nc
    m=$(( m + 1 ))
   done
 
   # get mesh mask files
   rapatrie  ${MESH_MASK_ID}_byte_mask.nc $IDIR mask.nc
   rapatrie  ${MESH_MASK_ID}_mesh_hgr.nc $IDIR mesh_hgr.nc
   rapatrie  ${MESH_MASK_ID}_mesh_zgr.nc $IDIR mesh_zgr.nc
 
   # Ascii output file:
   fice=${CONFCASE}_y${YEAR}_icemonth.txt
 
   m=1
   while (( $m <= 12 )) ; do
    mm=$( printf "%02d" $m )
    case $mm in 
    01) echo '###' $YEAR $mm > $fice ;;
    *)  echo '###' $YEAR $mm >> $fice ;;
    esac
    cdficediags ${CONFCASE}_y${YEAR}m${mm}_icemod.nc  >> $fice
    m=$(( m + 1 ))
   done
 
   expatrie $fice $DIAGS $fice 
 
#### Append corresponding lines to matlab file for time series
   #ice
   month='01 02 03 04 05 06 07 08 09 10 11'
   if [ $(chkfile $MONITOR/${CONFCASE}_icemonth.mtl ) == present ] ; then
     rapatrie ${CONFCASE}_icemonth.mtl $MONITOR ${CONFCASE}_icemonth.mtl
   else
    # first time: create file and add header
    echo 0000 $month  $month $month $month $month $month > ${CONFCASE}_icemonth.mtl
   fi
 
   year=$( head -1 $fice | awk '{ print $2}' )
   nvol=$(  cat $fice | grep -e 'NVolume' | grep -v NVolumet | awk '{ printf "%.0f  ", $4}' )
   svol=$(  cat $fice | grep -e 'SVolume' | grep -v SVolumet | awk '{ printf "%.0f  ", $4}' )
   narea=$(  cat $fice | grep -e 'NArea' | awk '{ printf "%.0f  ", $4}' )
   sarea=$(  cat $fice | grep -e 'SArea' | awk '{ printf "%.0f  ", $4}' )
   nextent=$(  cat $fice | grep -e 'NExtend' | awk '{ printf "%.0f  ", $4}' )
   sextent=$(  cat $fice | grep -e 'SExtend' | awk '{ printf "%.0f  ", $4}' )
 
   echo $year $nvol $svol $narea $sarea $nextent $sextent >> ${CONFCASE}_icemonth.mtl
 
   expatrie ${CONFCASE}_icemonth.mtl  $MONITOR ${CONFCASE}_icemonth.mtl
 
#### cp to web site
   cptoweb ${CONFCASE}_icemonth.mtl
 
   # clean up a little bit
   \rm ${CONFCASE}_icemonth.mtl
  fi

