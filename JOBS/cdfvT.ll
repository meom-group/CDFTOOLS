#!/bin/csh
# @ cpu_limit  = 3600
# @ data_limit = 1gb
# Nom du travail LoadLeveler
# @ job_name   = cdfvt
# Fichier de sortie standard du travail       
# @ output     = $(job_name).$(jobid)
# Fichier de sortie d'erreur du travail
# @ error      =  $(job_name).$(jobid)
# @ queue                   


set echo


set CONFIG=ORCA025
set CASE=G32

set YEAR=0010
#

set CONFCASE=${CONFIG}-${CASE}

set CDFTOOLS=~rcli002/CDFTOOLS-2.0


cd $TMPDIR
mkdir MONTHLY

cp $CDFTOOLS/att.txt .
rsh gaya mkdir ${CONFIG}/${CONFCASE}-MEAN/$YEAR/

#goto annual
#goto quarterly
# Monthly mean

foreach month (01 02 03 04 05 06 07 08 09 10 11 12 )
  foreach f ( `rsh gaya ls    ${CONFIG}/${CONFCASE}-S/$YEAR/${CONFCASE}_y${YEAR}m${month}\*_gridT.nc ${CONFIG}/${CONFCASE}-S/$YEAR/${CONFCASE}_y${YEAR}m${month}\*_gridU.nc ${CONFIG}/${CONFCASE}-S/$YEAR/${CONFCASE}_y${YEAR}m${month}\*_gridV.nc `)
   mfget $f ./
  end 

  set list=''
  foreach f ( ${CONFCASE}_y${YEAR}m${month}d??_gridT.nc )
    set tag=`echo $f | sed -e "s/${CONFCASE}_//" -e 's/_gridT.nc//'`
    set list=($list $tag )
  end
   $CDFTOOLS/cdfvT $CONFCASE $list
   mfput vt.nc ${CONFIG}/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}m${month}_VT.nc
   \rm  vt.nc ${CONFCASE}_y${YEAR}m${month}*_grid[TUV].nc 
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
    mfget $CONFIG/${CONFCASE}-MEAN/${YEAR}/${CONFCASE}_y${YEAR}m${month}_VT.nc ./
    set list=($list ${CONFCASE}_y${YEAR}m${month}_VT.nc )
  end

   $CDFTOOLS/cdfmoy  $list
   mv  cdfmoy.nc ${CONFCASE}_y${YEAR}_${season}_VT.nc
   mfput ${CONFCASE}_y${YEAR}_${season}_VT.nc ${CONFIG}/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_${season}_VT.nc
   mv ${CONFCASE}_y${YEAR}_${season}_VT.nc ../QUARTERLY
   \rm -f cdfmoy.nc cdfmoy2.nc $list

 end


annual:
# ANNUAL
  cd $TMPDIR
  cd QUARTERLY
  ln -s ../att.txt .
  
   set list=''
  foreach season ( WINT SPRI SUMM FALL )
   set list=($list ${CONFCASE}_y${YEAR}_${season}_VT.nc)
  end
  $CDFTOOLS/cdfmoy $list
   mv cdfmoy.nc ${CONFCASE}_y${YEAR}_ANNUAL_VT.nc
   mfput ${CONFCASE}_y${YEAR}_ANNUAL_VT.nc ${CONFIG}/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_ANNUAL_VT.nc
   \rm -f cdfmoy.nc cdfmoy2.nc $list

