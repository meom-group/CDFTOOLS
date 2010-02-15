#!/bin/ksh
  ##  $Rev$
  ##  $Date$
  ##  $Id$

CONFIG=ORCA025
CASE=G70
TAG=y1980-2004
CDFTOOLS=~molines/CDFTOOLS-2.1

CONFCASE=${CONFIG}-${CASE}
m=1
while (( m <= 12 )) ; do
  mm=$( printf "%02d" $m )
  f=${CONFCASE}_${TAG}m${mm}_gridT.nc
  r=runoff_m${mm}.nc
  wflx=$(echo $f | sed -e "s/gridT/wflx/" )
  $CDFTOOLS/cdfwflx $f $r
  mv wflx.nc $wflx

  m=$(( m + 1 ))
done

# ANNUAL
  f=${CONFCASE}_${TAG}_gridT.nc
  r=runoff_ANNUAL.nc
  wflx=$(echo $f | sed -e "s/gridT/wflx/" )

  $CDFTOOLS/cdfwflx $f $r
  mv wflx.nc $wflx
