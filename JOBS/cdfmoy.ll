#!/bin/csh
# @ cpu_limit  = 4200
# @ data_limit = 1gb
# Nom du travail LoadLeveler
# @ job_name   = cdfmoy
# Fichier de sortie standard du travail       
# @ output     = $(job_name).$(jobid)
# Fichier de sortie d'erreur du travail
# @ error      =  $(job_name).$(jobid)
# @ queue                   


set echo


set CONFIG=ORCA025
set CASE=G33

set YEAR=0001
#

set CONFCASE=${CONFIG}-${CASE}
set CDFTOOLS=~rcli002/CDFTOOLS-2.0/


cd $TMPDIR
mkdir MONTHLY

cp $CDFTOOLS/att.txt .
rsh gaya mkdir ${CONFIG}/${CONFCASE}-MEAN/$YEAR/

#goto annual
#goto quarterly
# Monthly mean
#
foreach month (01 02 03 04 05 06 07 08 09 10 11 12 )
  foreach f ( `rsh gaya ls ${CONFIG}/${CONFCASE}-S/$YEAR/${CONFCASE}_y${YEAR}m${month}\*_gridT.nc `)
   mfget $f ./
  end 

  set list=''
  foreach f ( ${CONFCASE}_y${YEAR}m${month}d??_gridT.nc )
    set list=($list $f )
  end
   $CDFTOOLS/cdfmoy $list
   mfput cdfmoy.nc ${CONFIG}/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}m${month}_gridT.nc
   mfput cdfmoy2.nc ${CONFIG}/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}m${month}_gridT2.nc
   \rm $list cdfmoy.nc cdfmoy2.nc

  foreach f ( `rsh gaya ls ${CONFIG}/${CONFCASE}-S/$YEAR/${CONFCASE}_y${YEAR}m${month}\*_gridU.nc `)
   mfget $f ./
  end 

  set list=''
  foreach f ( ${CONFCASE}_y${YEAR}m${month}d??_gridU.nc )
    set list=($list $f )
  end
   $CDFTOOLS/cdfmoy $list
   mfput cdfmoy.nc ${CONFIG}/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}m${month}_gridU.nc
   mfput cdfmoy2.nc ${CONFIG}/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}m${month}_gridU2.nc
   \rm $list cdfmoy.nc cdfmoy2.nc

  foreach f ( `rsh gaya ls ${CONFIG}/${CONFCASE}-S/$YEAR/${CONFCASE}_y${YEAR}m${month}\*_gridV.nc `)
   mfget $f ./
  end 

  set list=''
  foreach f ( ${CONFCASE}_y${YEAR}m${month}d??_gridV.nc )
    set list=($list $f )
  end
   $CDFTOOLS/cdfmoy $list
   mfput cdfmoy.nc ${CONFIG}/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}m${month}_gridV.nc
   mfput cdfmoy2.nc ${CONFIG}/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}m${month}_gridV2.nc
   \rm $list cdfmoy.nc cdfmoy2.nc

  foreach f ( `rsh gaya ls ${CONFIG}/${CONFCASE}-S/$YEAR/${CONFCASE}_y${YEAR}m${month}\*_gridW.nc `)
   mfget $f ./
  end 

  set list=''
  foreach f ( ${CONFCASE}_y${YEAR}m${month}d??_gridW.nc )
    set list=($list $f )
  end
   $CDFTOOLS/cdfmoy $list
   mfput cdfmoy.nc ${CONFIG}/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}m${month}_gridW.nc
   mfput cdfmoy2.nc ${CONFIG}/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}m${month}_gridW2.nc
   \rm $list cdfmoy.nc cdfmoy2.nc

  foreach f ( `rsh gaya ls ${CONFIG}/${CONFCASE}-S/$YEAR/${CONFCASE}_y${YEAR}m${month}\*_icemod.nc `)
   mfget $f ./
  end 

  set list=''
  foreach f ( ${CONFCASE}_y${YEAR}m${month}d??_icemod.nc )
    set list=($list $f )
  end
   $CDFTOOLS/cdfmoy $list
   mfput cdfmoy.nc ${CONFIG}/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}m${month}_icemod.nc
   \rm $list cdfmoy.nc cdfmoy2.nc

 end

quarterly:
# Quarterly mean 
  cd $TMPDIR
  mkdir QUARTERLY
  cd MONTHLY

 foreach season ( WINT SPRI SUMM FALL )
   switch ($season)
    case WINT: 
     set smon=(01 02 03 )
     breaksw
    case SPRI: 
     set smon=(04 05 06 )
     breaksw
    case SUMM: 
     set smon=(07 08 09 )
     breaksw
    case FALL: 
     set smon=(10 11 12 )
     breaksw
    default: 
     echo error ; exit 1
     breaksw
   endsw

  ln -s ../att.txt .

  set list=''
  foreach month ( $smon )
    mfget $CONFIG/${CONFCASE}-MEAN/${YEAR}/${CONFCASE}_y${YEAR}m${month}_gridT.nc ./
    set list=($list ${CONFCASE}_y${YEAR}m${month}_gridT.nc )
  end

   $CDFTOOLS/cdfmoy $list
   mv  cdfmoy.nc ${CONFCASE}_y${YEAR}_${season}_gridT.nc
   mfput ${CONFCASE}_y${YEAR}_${season}_gridT.nc ${CONFIG}/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_${season}_gridT.nc
   mv ${CONFCASE}_y${YEAR}_${season}_gridT.nc ../QUARTERLY
   \rm -f cdfmoy.nc cdfmoy2.nc $list

  set list=''
  foreach month ( $smon )
    mfget $CONFIG/${CONFCASE}-MEAN/${YEAR}/${CONFCASE}_y${YEAR}m${month}_gridT2.nc ./
    set list=($list ${CONFCASE}_y${YEAR}m${month}_gridT2.nc )
  end
   $CDFTOOLS/cdfmoy $list
   mv  cdfmoy.nc ${CONFCASE}_y${YEAR}_${season}_gridT2.nc
   mfput ${CONFCASE}_y${YEAR}_${season}_gridT2.nc ${CONFIG}/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_${season}_gridT2.nc
   mv ${CONFCASE}_y${YEAR}_${season}_gridT2.nc ../QUARTERLY
   \rm -f cdfmoy.nc cdfmoy2.nc $list

  set list=''
  foreach month ( $smon )
    mfget $CONFIG/${CONFCASE}-MEAN/${YEAR}/${CONFCASE}_y${YEAR}m${month}_gridU.nc ./
    set list=($list ${CONFCASE}_y${YEAR}m${month}_gridU.nc )
  end
   $CDFTOOLS/cdfmoy $list
   mv  cdfmoy.nc ${CONFCASE}_y${YEAR}_${season}_gridU.nc
   mfput ${CONFCASE}_y${YEAR}_${season}_gridU.nc ${CONFIG}/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_${season}_gridU.nc
   mv ${CONFCASE}_y${YEAR}_${season}_gridU.nc ../QUARTERLY
   \rm -f cdfmoy.nc cdfmoy2.nc $list

  set list=''
  foreach month ( $smon )
    mfget $CONFIG/${CONFCASE}-MEAN/${YEAR}/${CONFCASE}_y${YEAR}m${month}_gridU2.nc ./
    set list=($list ${CONFCASE}_y${YEAR}m${month}_gridU2.nc )
  end
   $CDFTOOLS/cdfmoy $list
   mv  cdfmoy.nc ${CONFCASE}_y${YEAR}_${season}_gridU2.nc
   mfput ${CONFCASE}_y${YEAR}_${season}_gridU2.nc ${CONFIG}/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_${season}_gridU2.nc
   mv ${CONFCASE}_y${YEAR}_${season}_gridU2.nc ../QUARTERLY
   \rm -f cdfmoy.nc cdfmoy2.nc $list

  set list=''
  foreach month ( $smon )
    mfget $CONFIG/${CONFCASE}-MEAN/${YEAR}/${CONFCASE}_y${YEAR}m${month}_gridV.nc ./
    set list=($list ${CONFCASE}_y${YEAR}m${month}_gridV.nc )
  end
   $CDFTOOLS/cdfmoy $list
   mv  cdfmoy.nc ${CONFCASE}_y${YEAR}_${season}_gridV.nc
   mfput ${CONFCASE}_y${YEAR}_${season}_gridV.nc ${CONFIG}/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_${season}_gridV.nc
   mv ${CONFCASE}_y${YEAR}_${season}_gridV.nc ../QUARTERLY
   \rm -f cdfmoy.nc cdfmoy2.nc $list

  set list=''
  foreach month ( $smon )
    mfget $CONFIG/${CONFCASE}-MEAN/${YEAR}/${CONFCASE}_y${YEAR}m${month}_gridV2.nc ./
    set list=($list ${CONFCASE}_y${YEAR}m${month}_gridV2.nc )
  end
   $CDFTOOLS/cdfmoy $list
   mv  cdfmoy.nc ${CONFCASE}_y${YEAR}_${season}_gridV2.nc
   mfput ${CONFCASE}_y${YEAR}_${season}_gridV2.nc ${CONFIG}/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_${season}_gridV2.nc
   mv ${CONFCASE}_y${YEAR}_${season}_gridV2.nc ../QUARTERLY
   \rm -f cdfmoy.nc cdfmoy2.nc $list

  set list=''
  foreach month ( $smon )
    mfget $CONFIG/${CONFCASE}-MEAN/${YEAR}/${CONFCASE}_y${YEAR}m${month}_gridW.nc ./
    set list=($list ${CONFCASE}_y${YEAR}m${month}_gridW.nc )
  end
   $CDFTOOLS/cdfmoy $list
   mv  cdfmoy.nc ${CONFCASE}_y${YEAR}_${season}_gridW.nc
   mfput ${CONFCASE}_y${YEAR}_${season}_gridW.nc ${CONFIG}/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_${season}_gridW.nc
   mv ${CONFCASE}_y${YEAR}_${season}_gridW.nc ../QUARTERLY
   \rm -f cdfmoy.nc cdfmoy2.nc $list

  set list=''
  foreach month ( $smon )
    mfget $CONFIG/${CONFCASE}-MEAN/${YEAR}/${CONFCASE}_y${YEAR}m${month}_gridW2.nc ./
    set list=($list ${CONFCASE}_y${YEAR}m${month}_gridW2.nc )
  end
   $CDFTOOLS/cdfmoy $list
   mv  cdfmoy.nc ${CONFCASE}_y${YEAR}_${season}_gridW2.nc
   mfput ${CONFCASE}_y${YEAR}_${season}_gridW2.nc ${CONFIG}/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_${season}_gridW2.nc
   mv ${CONFCASE}_y${YEAR}_${season}_gridW2.nc ../QUARTERLY
   \rm -f cdfmoy.nc cdfmoy2.nc $list

  set list=''
  foreach month ( $smon )
    mfget $CONFIG/${CONFCASE}-MEAN/${YEAR}/${CONFCASE}_y${YEAR}m${month}_icemod.nc ./
    set list=($list ${CONFCASE}_y${YEAR}m${month}_icemod.nc )
  end
   $CDFTOOLS/cdfmoy $list
   mv  cdfmoy.nc ${CONFCASE}_y${YEAR}_${season}_icemod.nc
   mfput ${CONFCASE}_y${YEAR}_${season}_icemod.nc ${CONFIG}/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_${season}_icemod.nc
   mv ${CONFCASE}_y${YEAR}_${season}_icemod.nc ../QUARTERLY
   \rm -f cdfmoy.nc cdfmoy2.nc $list

 end


annual:
# ANNUAL
  cd $TMPDIR
  cd QUARTERLY
  ln -s ../att.txt .
  
   set list=''
  foreach season ( WINT SPRI SUMM FALL )
   set list=($list ${CONFCASE}_y${YEAR}_${season}_gridT.nc)
  end
   $CDFTOOLS/cdfmoy $list
   mv cdfmoy.nc ${CONFCASE}_y${YEAR}_ANNUAL_gridT.nc
   mfput ${CONFCASE}_y${YEAR}_ANNUAL_gridT.nc ${CONFIG}/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_ANNUAL_gridT.nc
   \rm -f cdfmoy.nc cdfmoy2.nc $list

   set list=''
  foreach season ( WINT SPRI SUMM FALL )
   set list=($list ${CONFCASE}_y${YEAR}_${season}_gridT2.nc)
  end
   $CDFTOOLS/cdfmoy $list
   mv cdfmoy.nc ${CONFCASE}_y${YEAR}_ANNUAL_gridT2.nc
   mfput ${CONFCASE}_y${YEAR}_ANNUAL_gridT2.nc ${CONFIG}/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_ANNUAL_gridT2.nc
   \rm -f cdfmoy.nc cdfmoy2.nc $list


   set list=''
  foreach season ( WINT SPRI SUMM FALL )
   set list=($list ${CONFCASE}_y${YEAR}_${season}_gridU.nc)
  end
   $CDFTOOLS/cdfmoy $list
   mv cdfmoy.nc ${CONFCASE}_y${YEAR}_ANNUAL_gridU.nc
   mfput ${CONFCASE}_y${YEAR}_ANNUAL_gridU.nc ${CONFIG}/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_ANNUAL_gridU.nc
   \rm -f cdfmoy.nc cdfmoy2.nc $list

   set list=''
  foreach season ( WINT SPRI SUMM FALL )
   set list=($list ${CONFCASE}_y${YEAR}_${season}_gridU2.nc)
  end
   $CDFTOOLS/cdfmoy $list
   mv cdfmoy.nc ${CONFCASE}_y${YEAR}_ANNUAL_gridU2.nc
   mfput ${CONFCASE}_y${YEAR}_ANNUAL_gridU2.nc ${CONFIG}/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_ANNUAL_gridU2.nc
   \rm -f cdfmoy.nc cdfmoy2.nc $list


   set list=''
  foreach season ( WINT SPRI SUMM FALL )
   set list=($list ${CONFCASE}_y${YEAR}_${season}_gridV.nc)
  end
   $CDFTOOLS/cdfmoy $list
   mv cdfmoy.nc ${CONFCASE}_y${YEAR}_ANNUAL_gridV.nc
   mfput ${CONFCASE}_y${YEAR}_ANNUAL_gridV.nc ${CONFIG}/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_ANNUAL_gridV.nc
   \rm -f cdfmoy.nc cdfmoy2.nc $list

   set list=''
  foreach season ( WINT SPRI SUMM FALL )
   set list=($list ${CONFCASE}_y${YEAR}_${season}_gridV2.nc)
  end
   $CDFTOOLS/cdfmoy $list
   mv cdfmoy.nc ${CONFCASE}_y${YEAR}_ANNUAL_gridV2.nc
   mfput ${CONFCASE}_y${YEAR}_ANNUAL_gridV2.nc ${CONFIG}/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_ANNUAL_gridV2.nc
   \rm -f cdfmoy.nc cdfmoy2.nc $list


   set list=''
  foreach season ( WINT SPRI SUMM FALL )
   set list=($list ${CONFCASE}_y${YEAR}_${season}_gridW.nc)
  end
   $CDFTOOLS/cdfmoy $list
   mv cdfmoy.nc ${CONFCASE}_y${YEAR}_ANNUAL_gridW.nc
   mfput ${CONFCASE}_y${YEAR}_ANNUAL_gridW.nc ${CONFIG}/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_ANNUAL_gridW.nc
   \rm -f cdfmoy.nc cdfmoy2.nc $list

   set list=''
  foreach season ( WINT SPRI SUMM FALL )
   set list=($list ${CONFCASE}_y${YEAR}_${season}_gridW2.nc)
  end
   $CDFTOOLS/cdfmoy $list
   mv cdfmoy.nc ${CONFCASE}_y${YEAR}_ANNUAL_gridW2.nc
   mfput ${CONFCASE}_y${YEAR}_ANNUAL_gridW2.nc ${CONFIG}/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_ANNUAL_gridW2.nc
   \rm -f cdfmoy.nc cdfmoy2.nc $list


   set list=''
  foreach season ( WINT SPRI SUMM FALL )
   set list=($list ${CONFCASE}_y${YEAR}_${season}_icemod.nc)
  end
   $CDFTOOLS/cdfmoy $list
   mv cdfmoy.nc ${CONFCASE}_y${YEAR}_ANNUAL_icemod.nc
   mfput ${CONFCASE}_y${YEAR}_ANNUAL_icemod.nc ${CONFIG}/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_ANNUAL_icemod.nc
   \rm -f cdfmoy.nc cdfmoy2.nc $list

