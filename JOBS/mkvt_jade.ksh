#!/bin/ksh
set -x
cd ..
. config_def.ksh
. function_def.ksh

CONFCASE=${CONFIG}-${CASE}

login_node=service1

YEAR=YYEEAARR
grid=GGRRIIDD
NP=NNPP

cd $YEAR

MEANY=$CONFIG/${CONFIG}-${CASE}-MEAN/$YEAR
SDIRY=$CONFIG/${CONFIG}-${CASE}-S/$YEAR

mkdir VT; cd VT ;
for month in 01 02 03 04 05 06 07 08 09 10 11 12  ; do
    getmonth $month gridT
    getmonth $month gridU
    getmonth $month gridV
        #####################
    list=''
    for f in ${CONFCASE}_y${YEAR}m${month}d??_gridT.nc ; do
        tag=$( echo $f | awk -F_ '{print $2}' )
        list="$list $tag"
    done
    ../../cdfvT $CONFCASE $list
    mv -f vt.nc ${CONFCASE}_y${YEAR}m${month}_VT.nc
    rm -f ${CONFCASE}_y${YEAR}m${month}d??_grid?.nc
done 
 # annual mean  (uses a ponderation to compute the exact annual mean ). ! suppose 5 day averages when creating monthly mean
../../cdfmoy_annual ${CONFCASE}_y${YEAR}m??_VT.nc
mv -f cdfmoy_annual.nc ${CONFCASE}_y${YEAR}_ANNUAL_VT.nc
mv -f ${CONFCASE}_y${YEAR}*_VT.nc ../. 
cd .. ; rmdir VT ; touch OK$NP 
