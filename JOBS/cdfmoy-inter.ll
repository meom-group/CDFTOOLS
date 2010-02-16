#!/bin/csh
# @ wall_clock_limit = 20:00:00
# @ as_limit = 1gb
# @ job_name   = cdfmoy-inter
# @ output     = $(job_name).$(jobid)
# @ error      =  $(job_name).$(jobid)
# @ notify_user = molines@hmg.inpg.fr
# @ notification = error
# @ queue



set echo

set INTER=1980-2004
set CONFIG=ORCA025.L75
set CASELIST=( G85 )

set CDFTOOLS=~/CDFTOOLS-2.1

#######################
set tmp=`echo $INTER | sed -e 's/-/ /'`
set year1=`echo $tmp[1] | awk '{printf "%04d", $1 }'`
set year2=`echo $tmp[2] | awk '{printf "%04d", $1 }'`

set nyear=`expr  $year2 - $year1 + 1 `
set year=$year1
set n=1
set YEARLIST=''
while ( $n   <= $nyear )
  set YEARLIST=($YEARLIST $year)
  set year=`expr $year + 1 `
  set year=`echo $year  | awk '{printf "%04d", $1 }'`
  @ n ++
end

   cd $TMPDIR

foreach CASE ( $CASELIST )
   set CONFCASE=${CONFIG}-${CASE}
   rsh gaya mkdir $CONFIG/${CONFCASE}-MEAN/$INTER


###  GRID T ###
###############
  set list=''
  if ( ! -f $HOMEGAYA/${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}_gridT.nc )  then
  foreach YEAR ( $YEARLIST )
      mfget  $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_ANNUAL_gridT.nc
      set list=($list ${CONFCASE}_y${YEAR}_ANNUAL_gridT.nc )
  end 

  $CDFTOOLS/cdfmoy $list
   mfput cdfmoy.nc ${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}_gridT.nc
   \rm -f cdfmoy.nc cdfmoy2.nc $list
  endif

### GRID T2 ###
###############
  set list=''
  if ( ! -f $HOMEGAYA/${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}_gridT2.nc )  then
  foreach YEAR ( $YEARLIST )
      mfget  $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_ANNUAL_gridT2.nc
     set list=($list ${CONFCASE}_y${YEAR}_ANNUAL_gridT2.nc )
  end

  $CDFTOOLS/cdfmoy $list
   mfput cdfmoy.nc ${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}_gridT2.nc
   \rm -f cdfmoy.nc cdfmoy2.nc $list
  endif

###  GRID U ###
###############
  set list=''
  if ( ! -f $HOMEGAYA/${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}_gridU.nc )  then
  foreach YEAR ( $YEARLIST )
      mfget  $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_ANNUAL_gridU.nc
     set list=($list ${CONFCASE}_y${YEAR}_ANNUAL_gridU.nc )
  end

  $CDFTOOLS/cdfmoy $list
   mfput cdfmoy.nc ${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}_gridU.nc
   \rm -f cdfmoy.nc cdfmoy2.nc $list
  endif
###  GRID U2 ###
################
  set list=''
  if ( ! -f $HOMEGAYA/${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}_gridU2.nc )  then
  foreach YEAR ( $YEARLIST )
      mfget  $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_ANNUAL_gridU2.nc
     set list=($list ${CONFCASE}_y${YEAR}_ANNUAL_gridU2.nc )
  end

  $CDFTOOLS/cdfmoy $list
   mfput cdfmoy.nc  ${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}_gridU2.nc
   \rm -f cdfmoy.nc cdfmoy2.nc $list
  endif

###  GRID V ###
###############
  set list=''
  if ( ! -f $HOMEGAYA/${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}_gridV.nc )  then
  foreach YEAR ( $YEARLIST )
      mfget  $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_ANNUAL_gridV.nc
     set list=($list ${CONFCASE}_y${YEAR}_ANNUAL_gridV.nc )
  end

  $CDFTOOLS/cdfmoy $list
   mfput cdfmoy.nc ${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}_gridV.nc
   \rm -f cdfmoy.nc cdfmoy2.nc $list
  endif
###  GRID V2 ###
################
  set list=''
  if ( ! -f $HOMEGAYA/${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}_gridV2.nc )  then
  foreach YEAR ( $YEARLIST )
      mfget  $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_ANNUAL_gridV2.nc
     set list=($list ${CONFCASE}_y${YEAR}_ANNUAL_gridV2.nc )
  end

  $CDFTOOLS/cdfmoy $list
   mfput cdfmoy.nc ${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}_gridV2.nc
   \rm -f cdfmoy.nc cdfmoy2.nc $list
   endif

###  GRID W ###
###############
  set list=''
  if ( ! -f $HOMEGAYA/${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}_gridW.nc )  then
  foreach YEAR ( $YEARLIST )
      mfget  $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_ANNUAL_gridW.nc
     set list=($list ${CONFCASE}_y${YEAR}_ANNUAL_gridW.nc )
  end

  $CDFTOOLS/cdfmoy $list
   mfput cdfmoy.nc ${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}_gridW.nc
   \rm -f cdfmoy.nc cdfmoy2.nc $list
  endif
###  GRID W2 ###
################
  set list=''
  if ( ! -f $HOMEGAYA/${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}_gridW2.nc )  then
  foreach YEAR ( $YEARLIST )
      mfget  $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_ANNUAL_gridW2.nc
     set list=($list ${CONFCASE}_y${YEAR}_ANNUAL_gridW2.nc )
  end

  $CDFTOOLS/cdfmoy $list
   mfput cdfmoy.nc ${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}_gridW2.nc
   \rm -f cdfmoy.nc cdfmoy2.nc $list
   endif
###  ICEMOD ###
###############

  set list=''
  if ( ! -f $HOMEGAYA/${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}_icemod.nc )  then
  foreach YEAR ( $YEARLIST )
      mfget  $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_ANNUAL_icemod.nc
     set list=($list ${CONFCASE}_y${YEAR}_ANNUAL_icemod.nc )
  end

  $CDFTOOLS/cdfmoy $list
   mfput cdfmoy.nc ${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}_icemod.nc
   \rm -f cdfmoy.nc cdfmoy2.nc $list
  endif
  ## m03 ##
  set list=''
  if ( ! -f $HOMEGAYA/${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}m03_icemod.nc )  then
  foreach YEAR ( $YEARLIST )
      mfget  $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}m03_icemod.nc
     set list=($list ${CONFCASE}_y${YEAR}m03_icemod.nc )
  end

  $CDFTOOLS/cdfmoy $list
   mfput cdfmoy.nc ${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}m03_icemod.nc
   \rm -f cdfmoy.nc cdfmoy2.nc $list
  endif

  ## m09 ##
  set list=''
  if ( ! -f $HOMEGAYA/${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}m09_icemod.nc )  then
  foreach YEAR ( $YEARLIST )
      mfget  $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}m09_icemod.nc
     set list=($list ${CONFCASE}_y${YEAR}m09_icemod.nc )
  end

  $CDFTOOLS/cdfmoy $list
   mfput cdfmoy.nc ${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}m09_icemod.nc
   \rm -f cdfmoy.nc cdfmoy2.nc $list
  endif

### EKE ###
  set list=''
  if ( ! -f $HOMEGAYA/${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}_EKE.nc )  then
  foreach YEAR ( $YEARLIST )
      mfget  $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_EKE.nc
     set list=($list ${CONFCASE}_y${YEAR}_EKE.nc )
  end

  $CDFTOOLS/cdfmoy $list
   mfput cdfmoy.nc ${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}_EKE.nc
   \rm -f cdfmoy.nc cdfmoy2.nc $list
  endif
### MOC ### 
  set list=''
  if ( ! -f $HOMEGAYA/${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}_MOC.nc )  then
  foreach YEAR ( $YEARLIST )
      mfget  $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_MOC.nc
     set list=($list ${CONFCASE}_y${YEAR}_MOC.nc )
  end

  $CDFTOOLS/cdfmoy $list
   mfput cdfmoy.nc ${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}_MOC.nc
   \rm -f cdfmoy.nc cdfmoy2.nc $list
  endif
### PSI ###
  set list=''
  if ( ! -f $HOMEGAYA/${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}_PSI.nc )  then
  foreach YEAR ( $YEARLIST )
      mfget  $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_PSI.nc
     set list=($list ${CONFCASE}_y${YEAR}_PSI.nc )
  end

  $CDFTOOLS/cdfmoy $list
   mfput cdfmoy.nc ${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}_PSI.nc
   \rm -f cdfmoy.nc cdfmoy2.nc $list
  endif

### MXL ###

  ## m03 ##
  set list=''
  if ( ! -f $HOMEGAYA/${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}m03_MXL.nc )  then
  foreach YEAR ( $YEARLIST )
      mfget  $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}m03_MXL.nc
     set list=($list ${CONFCASE}_y${YEAR}m03_MXL.nc )
  end

  $CDFTOOLS/cdfmoy $list
   mfput cdfmoy.nc ${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}m03_MXL.nc
   \rm -f cdfmoy.nc cdfmoy2.nc $list
  endif

  ## m09 ##
  set list=''
  if ( ! -f $HOMEGAYA/${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}m09_MXL.nc )  then
  foreach YEAR ( $YEARLIST )
      mfget  $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}m09_MXL.nc
     set list=($list ${CONFCASE}_y${YEAR}m09_MXL.nc )
  end

  $CDFTOOLS/cdfmoy $list
   mfput cdfmoy.nc ${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}m09_MXL.nc
   \rm -f cdfmoy.nc cdfmoy2.nc $list
  endif

end
