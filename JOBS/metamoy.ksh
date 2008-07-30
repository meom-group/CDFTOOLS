#!/bin/ksh

# metamoy.ksh script : to launch both cdfmoy and cdfvT on 2 separated jobs for a given year (argument)

#  $Rev$
#  $Date$
#  $Id$

set -x

if [ $# == 0 ] ; then
    echo USAGE: metamoy.ksh  year
    exit 0
fi

year=$1

. ./config_def.ksh   # CDFTOOLS is set in this script (which is sourced now)

cat $CDFTOOLS/JOBS/cdfmoy_skel_new.ksh | sed -e "s/YYYY/$year/g" -e "s/YYYE/$year/g" -e "s/CCOONNFF/$CONFIG/g" -e "s/CCAASSEE/$CASE/g" \
   > cdfmoytmp.$$.ll

chmod u+x cdfmoytmp.$$.ll
$SUB  cdfmoytmp.$$.ll

cat $CDFTOOLS/JOBS/cdfvT_skel_new.ksh | sed -e "s/YYYY/$year/g" -e "s/YYYE/$year/g" -e "s/CCOONNFF/$CONFIG/g" -e "s/CCAASSEE/$CASE/g" \
  > cdfvTtmp.$$.ll
chmod u+x cdfvTtmp.$$.ll
$SUB cdfvTtmp.$$.ll

#if (( $year > 1958 )) ; then
# cat $CDFTOOLS/JOBS/cdfmoy_trc_skel_new.ksh | sed -e "s/YYYY/$year/g" -e "s/YYYE/$year/g" -e "s/CCOONNFF/$CONFIG/g" -e "s/CCAASSEE/$CASE/g" \
#  > cdfTRCtmp.$$.ll
#chmod u+x cdfTRCtmp.$$.ll
#  $SUB cdfTRCtmp.$$.ll
#fi
