#!/bin/csh
# @ cpu_limit  = 3500
# @ data_limit = 1gb
# Nom du travail LoadLeveler
# @ job_name   = cdfmoy-inter
# Fichier de sortie standard du travail       
# @ output     = $(job_name).$(jobid)
# Fichier de sortie d'erreur du travail
# @ error      =  $(job_name).$(jobid)
# @ queue                   


set echo

set INTER=0008-0010
set CONFIG=ORCA05
set CASELIST=( GERA GMAR)

set CDFTOOLS=~rcli002/CDFTOOLS-2.0

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
   cp $CDFTOOLS/att.txt .

foreach CASE ( $CASELIST )
   set CONFCASE=${CONFIG}-${CASE}
   rsh gaya mkdir $CONFIG/${CONFCASE}-MEAN/$INTER


###  GRID T ###
###############
  set list=''
  foreach YEAR ( $YEARLIST )
      mfget  $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_ANNUAL_gridT.nc
      set list=($list ${CONFCASE}_y${YEAR}_ANNUAL_gridT.nc )
  end 

  $CDFTOOLS/cdfmoy $list
   mfput cdfmoy.nc ${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}_gridT.nc
   \rm -f cdfmoy.nc cdfmoy2.nc $list

### GRID T2 ###
###############
  set list=''
  foreach YEAR ( $YEARLIST )
      mfget  $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_ANNUAL_gridT2.nc
     set list=($list ${CONFCASE}_y${YEAR}_ANNUAL_gridT2.nc )
  end

  $CDFTOOLS/cdfmoy $list
   mfput cdfmoy.nc ${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}_gridT2.nc
   \rm -f cdfmoy.nc cdfmoy2.nc $list

###  GRID U ###
###############
  set list=''
  foreach YEAR ( $YEARLIST )
      mfget  $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_ANNUAL_gridU.nc
     set list=($list ${CONFCASE}_y${YEAR}_ANNUAL_gridU.nc )
  end

  $CDFTOOLS/cdfmoy $list
   mfput cdfmoy.nc ${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}_gridU.nc
   \rm -f cdfmoy.nc cdfmoy2.nc $list

###  GRID U2 ###
################
  set list=''
  foreach YEAR ( $YEARLIST )
      mfget  $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_ANNUAL_gridU2.nc
     set list=($list ${CONFCASE}_y${YEAR}_ANNUAL_gridU2.nc )
  end

  $CDFTOOLS/cdfmoy $list
   mfput cdfmoy.nc  ${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}_gridU2.nc
   \rm -f cdfmoy.nc cdfmoy2.nc $list


###  GRID V ###
###############
  set list=''
  foreach YEAR ( $YEARLIST )
      mfget  $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_ANNUAL_gridV.nc
     set list=($list ${CONFCASE}_y${YEAR}_ANNUAL_gridV.nc )
  end

  $CDFTOOLS/cdfmoy $list
   mfput cdfmoy.nc ${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}_gridV.nc
   \rm -f cdfmoy.nc cdfmoy2.nc $list

###  GRID V2 ###
################
  set list=''
  foreach YEAR ( $YEARLIST )
      mfget  $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_ANNUAL_gridV2.nc
     set list=($list ${CONFCASE}_y${YEAR}_ANNUAL_gridV2.nc )
  end

  $CDFTOOLS/cdfmoy $list
   mfput cdfmoy.nc ${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}_gridV2.nc
   \rm -f cdfmoy.nc cdfmoy2.nc $list


###  GRID W ###
###############
  set list=''
  foreach YEAR ( $YEARLIST )
      mfget  $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_ANNUAL_gridW.nc
     set list=($list ${CONFCASE}_y${YEAR}_ANNUAL_gridW.nc )
  end

  $CDFTOOLS/cdfmoy $list
   mfput cdfmoy.nc ${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}_gridW.nc
   \rm -f cdfmoy.nc cdfmoy2.nc $list

###  GRID W2 ###
################
  set list=''
  foreach YEAR ( $YEARLIST )
      mfget  $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_ANNUAL_gridW2.nc
     set list=($list ${CONFCASE}_y${YEAR}_ANNUAL_gridW2.nc )
  end

  $CDFTOOLS/cdfmoy $list
   mfput cdfmoy.nc ${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}_gridW2.nc
   \rm -f cdfmoy.nc cdfmoy2.nc $list

###  ICEMOD ###
###############

  set list=''
  foreach YEAR ( $YEARLIST )
      mfget  $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_ANNUAL_icemod.nc
     set list=($list ${CONFCASE}_y${YEAR}_ANNUAL_icemod.nc )
  end

  $CDFTOOLS/cdfmoy $list
   mfput cdfmoy.nc ${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}_icemod.nc
   \rm -f cdfmoy.nc cdfmoy2.nc $list

end
