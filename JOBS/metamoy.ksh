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

cat $CDFTOOLS/JOBS/cdfmoy_skel_new.ksh | sed -e "s/YYYY/$year/g" -e "s/YYYE/$year/g" > cdfmoytmp.$$.ll
chmod u+x cdfmoytmp.$$.ll
 oarsub ./cdfmoytmp.$$.ll

cat $CDFTOOLS/JOBS/cdfvT_skel_new.ksh | sed -e "s/YYYY/$year/g" -e "s/YYYE/$year/g" > cdfvTtmp.$$.ll
chmod u+x cdfvTtmp.$$.ll
 oarsub ./cdfvTtmp.$$.ll

#if (( $year > 1958 )) ; then
# cat $CDFTOOLS/JOBS/cdfmoy_trc_skel.ll | sed -e "s/YYYY/$year/g" -e "s/YYYE/$year/g" > cdfTRCtmp.$$.ll
#  llsubmit cdfTRCtmp.$$.ll
#fi
