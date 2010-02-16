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

mkdir $grid; cd $grid;
for month in 01 02 03 04 05 06 07 08 09 10 11 12  ; do 
    getmonth $month $grid
    
    ../../cdfmoy ${CONFCASE}_y${YEAR}m${month}d??_$grid.nc
    mv -f cdfmoy.nc ${CONFCASE}_y${YEAR}m${month}_$grid.nc
    if [ $grid != icemod ] ; then
        mv -f cdfmoy2.nc ${CONFCASE}_y${YEAR}m${month}_${grid}2.nc
    fi
    rm -f ${CONFCASE}_y${YEAR}m${month}d??_$grid.nc
done
../../cdfmoy_annual ${CONFCASE}_y${YEAR}m??_$grid.nc
mv -f cdfmoy_annual.nc ${CONFCASE}_y${YEAR}_ANNUAL_$grid.nc
if [ $grid != icemod ] ; then
    ../../cdfmoy_annual ${CONFCASE}_y${YEAR}m??_${grid}2.nc 
    mv -f cdfmoy_annual.nc ${CONFCASE}_y${YEAR}_ANNUAL_${grid}2.nc
fi
mv -f ${CONFCASE}_y${YEAR}*_${grid}*.nc ../. 
cd .. ; rm -rf $grid ; touch OK$NP ;
