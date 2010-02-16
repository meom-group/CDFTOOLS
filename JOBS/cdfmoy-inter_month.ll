#!/bin/csh
# @ cpu_limit  = 7200
# @ data_limit = 1gb
# Nom du travail LoadLeveler
# @ job_name   = cdfmoy-inter
# Fichier de sortie standard du travail       
# @ output     = $(job_name).$(jobid)
# Fichier de sortie d'erreur du travail
# @ error      =  $(job_name).$(jobid)
# @ queue                   


set echo

set INTER=1990-2000
set CONFIG=ORCA025
set CASELIST=( G70 )

set CDFTOOLS=~rcli002/CDFTOOLS-2.1

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
  foreach month ( 01 02 03 04 05 06 07 08 09 10 11 12 )
  set list=''
  foreach YEAR ( $YEARLIST )
      mfget  $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}m${month}_gridT.nc
      set list=($list ${CONFCASE}_y${YEAR}m${month}_gridT.nc  )
  end 

  $CDFTOOLS/cdfmoy $list
   mfput cdfmoy.nc ${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}m${month}_gridT.nc
   \rm -f cdfmoy.nc cdfmoy2.nc $list

### GRID T2 ###
###############
  set list=''
  foreach YEAR ( $YEARLIST )
      mfget  $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}m${month}_gridT2.nc
     set list=($list ${CONFCASE}_y${YEAR}m${month}_gridT2.nc )
  end

  $CDFTOOLS/cdfmoy $list
   mfput cdfmoy.nc ${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}m${month}_gridT2.nc
   \rm -f cdfmoy.nc cdfmoy2.nc $list

###  GRID U ###
###############
  set list=''
  foreach YEAR ( $YEARLIST )
      mfget  $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}m${month}_gridU.nc
     set list=($list ${CONFCASE}_y${YEAR}m${month}_gridU.nc )
  end

  $CDFTOOLS/cdfmoy $list
   mfput cdfmoy.nc ${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}m${month}_gridU.nc
   \rm -f cdfmoy.nc cdfmoy2.nc $list

###  GRID U2 ###
################
  set list=''
  foreach YEAR ( $YEARLIST )
      mfget  $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}m${month}_gridU2.nc
     set list=($list ${CONFCASE}_y${YEAR}m${month}_gridU2.nc )
  end

  $CDFTOOLS/cdfmoy $list
   mfput cdfmoy.nc  ${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}m${month}_gridU2.nc
   \rm -f cdfmoy.nc cdfmoy2.nc $list


###  GRID V ###
###############
  set list=''
  foreach YEAR ( $YEARLIST )
      mfget  $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}m${month}_gridV.nc
     set list=($list ${CONFCASE}_y${YEAR}m${month}_gridV.nc )
  end

  $CDFTOOLS/cdfmoy $list
   mfput cdfmoy.nc ${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}m${month}_gridV.nc
   \rm -f cdfmoy.nc cdfmoy2.nc $list

###  GRID V2 ###
################
  set list=''
  foreach YEAR ( $YEARLIST )
      mfget  $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}m${month}_gridV2.nc
     set list=($list ${CONFCASE}_y${YEAR}m${month}_gridV2.nc )
  end

  $CDFTOOLS/cdfmoy $list
   mfput cdfmoy.nc ${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}m${month}_gridV2.nc
   \rm -f cdfmoy.nc cdfmoy2.nc $list


###  GRID W ###
###############
  set list=''
  foreach YEAR ( $YEARLIST )
      mfget  $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}m${month}_gridW.nc
     set list=($list ${CONFCASE}_y${YEAR}m${month}_gridW.nc )
  end

  $CDFTOOLS/cdfmoy $list
   mfput cdfmoy.nc ${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}m${month}_gridW.nc
   \rm -f cdfmoy.nc cdfmoy2.nc $list

###  GRID W2 ###
################
  set list=''
  foreach YEAR ( $YEARLIST )
      mfget  $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}m${month}_gridW2.nc
     set list=($list ${CONFCASE}_y${YEAR}m${month}_gridW2.nc )
  end

  $CDFTOOLS/cdfmoy $list
   mfput cdfmoy.nc ${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}m${month}_gridW2.nc
   \rm -f cdfmoy.nc cdfmoy2.nc $list

###  ICEMOD ###
###############

  set list=''
  foreach YEAR ( $YEARLIST )
      mfget  $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}m${month}_icemod.nc
     set list=($list ${CONFCASE}_y${YEAR}m${month}_icemod.nc )
  end

  $CDFTOOLS/cdfmoy $list
   mfput cdfmoy.nc ${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}m${month}_icemod.nc
   \rm -f cdfmoy.nc cdfmoy2.nc $list

end

end
