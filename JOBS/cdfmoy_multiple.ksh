#!/bin/ksh

CONFCASE=$1
YEAR=$2

CONFIG=${CONFCASE%-*}
CASE=${CONFCASE#*-}

# main script cp cdfmoy to the tmpdir. So, CDFTOOLS path  is  now  set  to ../
CDFTOOLS=../

 rsh gaya mkdir ${CONFIG}/${CONFCASE}-MEAN0/$YEAR/
# gridT files
  for f in $( rsh gaya ls ${CONFIG}/${CONFCASE}-S/$YEAR/${CONFCASE}_y${YEAR}m??d??\*_gridT.nc  ) ; do
   mfget $f ./
  done

   $CDFTOOLS/cdfmoy *gridT.nc
   mfput cdfmoy.nc ${CONFIG}/${CONFCASE}-MEAN0/$YEAR/${CONFCASE}_y${YEAR}_gridT.nc
   mfput cdfmoy2.nc ${CONFIG}/${CONFCASE}-MEAN0/$YEAR/${CONFCASE}_y${YEAR}_gridT2.nc
   \rm $list cdfmoy.nc cdfmoy2.nc *gridT.nc

# gridU files
  for f in $( rsh gaya ls ${CONFIG}/${CONFCASE}-S/$YEAR/${CONFCASE}_y${YEAR}m??d??\*_gridU.nc  ) ; do
   mfget $f ./
  done

   $CDFTOOLS/cdfmoy *gridU.nc
   mfput cdfmoy.nc ${CONFIG}/${CONFCASE}-MEAN0/$YEAR/${CONFCASE}_y${YEAR}_gridU.nc
   mfput cdfmoy2.nc ${CONFIG}/${CONFCASE}-MEAN0/$YEAR/${CONFCASE}_y${YEAR}_gridU2.nc
   \rm $list cdfmoy.nc cdfmoy2.nc *gridU.nc

# gridV files
  for f in $( rsh gaya ls ${CONFIG}/${CONFCASE}-S/$YEAR/${CONFCASE}_y${YEAR}m??d??\*_gridV.nc  ) ; do
   mfget $f ./
  done

   $CDFTOOLS/cdfmoy *gridV.nc
   mfput cdfmoy.nc ${CONFIG}/${CONFCASE}-MEAN0/$YEAR/${CONFCASE}_y${YEAR}_gridV.nc
   mfput cdfmoy2.nc ${CONFIG}/${CONFCASE}-MEAN0/$YEAR/${CONFCASE}_y${YEAR}_gridV2.nc
   \rm $list cdfmoy.nc cdfmoy2.nc *gridV.nc



# gridW files
  for f in $( rsh gaya ls ${CONFIG}/${CONFCASE}-S/$YEAR/${CONFCASE}_y${YEAR}m??d??\*_gridW.nc  ) ; do
   mfget $f ./
  done

   $CDFTOOLS/cdfmoy *gridW.nc
   mfput cdfmoy.nc ${CONFIG}/${CONFCASE}-MEAN0/$YEAR/${CONFCASE}_y${YEAR}_gridW.nc
   mfput cdfmoy2.nc ${CONFIG}/${CONFCASE}-MEAN0/$YEAR/${CONFCASE}_y${YEAR}_gridW2.nc
   \rm $list cdfmoy.nc cdfmoy2.nc *gridW.nc


# icemod files
  for f in $( rsh gaya ls ${CONFIG}/${CONFCASE}-S/$YEAR/${CONFCASE}_y${YEAR}m??d??\*_icemod.nc  ) ; do
   mfget $f ./
  done

   $CDFTOOLS/cdfmoy *icemod.nc
   mfput cdfmoy.nc ${CONFIG}/${CONFCASE}-MEAN0/$YEAR/${CONFCASE}_y${YEAR}_icemod.nc
   \rm $list cdfmoy.nc cdfmoy2.nc *icemod.nc

# TRC files
  for f in $( rsh gaya ls ${CONFIG}/${CONFCASE}-S/$YEAR/${CONFCASE}_y${YEAR}m??d??\*_ptrcT.nc  ) ; do
   mfget $f ./
  done

   $CDFTOOLS/cdfmoy *ptrcT.nc
   mfput cdfmoy.nc ${CONFIG}/${CONFCASE}-MEAN0/$YEAR/${CONFCASE}_y${YEAR}_ptrcT.nc
   \rm $list cdfmoy.nc cdfmoy2.nc *ptrcT.nc

# All is finished : touch a done file in ../

touch ../$YEAR.done


