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
  buoyflx=$(echo $f | sed -e "s/gridT/buoyflx/" )
  $CDFTOOLS/cdfbuoyflx $f $r
  mv buoyflx.nc $buoyflx

  m=$(( m + 1 ))
done
# concatenations of monthly files

ncrcat -h -a ${CONFCASE}_${TAG}m??_buoyflx.nc ${CONFCASE}_${TAG}_1m_buoyflx.nc

# ANNUAL
  f=${CONFCASE}_${TAG}_gridT.nc
  r=runoff_ANNUAL.nc
  buoyflx=$(echo $f | sed -e "s/gridT/buoyflx/" )

  $CDFTOOLS/cdfbuoyflx $f $r
  mv buoyflx.nc $buoyflx
