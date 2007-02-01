#!/bin/csh
# @ cpu_limit  = 3500
# @ data_limit = 1gb
# Nom du travail LoadLeveler
# @ job_name   = cdfvt-inter
# Fichier de sortie standard du travail       
# @ output     = $(job_name).$(jobid)
# Fichier de sortie d'erreur du travail
# @ error      =  $(job_name).$(jobid)
# @ queue                   


set echo

set INTER=0008-0010
set CONFIG=ORCA025
set CASELIST=( G32 )

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


###  VT ###
###############
  set list=''
  foreach YEAR ( $YEARLIST )
      mfget  $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_ANNUAL_VT.nc
      set list=($list ${CONFCASE}_y${YEAR}_ANNUAL_VT.nc )
  end 

  $CDFTOOLS/cdfmoy $list
   mfput cdfmoy.nc ${CONFIG}/${CONFCASE}-MEAN/$INTER/${CONFCASE}_y${INTER}_VT.nc
   \rm -f cdfmoy.nc cdfmoy2.nc $list

end
