#!/bin/ksh
# This script (mkmtl.ksh) is supposed to be used in the -DIAGS dir of a config.
# it retrieve the CONFIG CASE name from the directory name
# it build the mtl files from the annual info kept in -DIAGS after the monitoring
# it put the mtl files in the -MONITOR directory for the given CONFIG-CASE
#----------------------------------------------------------------------------------
# in order to be launched from the production machine, the absolute path of the 
# formatting scripts is used.
#----------------------------------------------------------------------------------
dir=$( basename `pwd` )
CONFCASE=${dir%-DIAGS}
CONFIG=${CONFCASE%-*}
CASE=${CONFCASE#*-}

bindir=~/MKMTL

echo Transport across section  \.\.\.
$bindir/section.ksh

echo T S Profiles  \.\.\.
$bindir/profile.ksh

echo T S Profiles + Levitus  \.\.\.
$bindir/profile_lev.ksh

echo Ice Diags  \.\.\.
$bindir/ice.ksh

echo Ice Month Diags  \.\.\.
$bindir/ice_month.ksh

echo  Gib Diags  \.\.\.
$bindir/gib.ksh

echo El Nino Diags  \.\.\.
$bindir/nino.ksh

echo Meridional Heat Transport Diags  \.\.\.
$bindir/heat.ksh

echo Max Overturning \.\.\. 
$bindir/maxmoc.ksh

echo Max Overturning at fixed latitude \.\.\. 
$bindir/maxmoc40.ksh

echo DCT monitoring  \.\.\.
$bindir/trpsig.ksh

echo TRACER \.\.\.
$bindir/trc.ksh

cd ../${CONFCASE}-MONITOR/

echo 'Update web site .........'
# check if ad hoc directories exists :
ssh meolipc.hmg.inpg.fr -l drakkar " if [ ! -d DRAKKAR/$CONFIG ] ; then mkdir DRAKKAR/$CONFIG ; fi "
ssh meolipc.hmg.inpg.fr -l drakkar " if [ ! -d DRAKKAR/$CONFIG/$CONFCASE ] ; then mkdir DRAKKAR/$CONFIG/$CONFCASE ; fi "
ssh meolipc.hmg.inpg.fr -l drakkar " if [ ! -d DRAKKAR/$CONFIG/$CONFCASE/DATA ] ; then mkdir DRAKKAR/$CONFIG/$CONFCASE/DATA ; fi "


scp  *mtl  drakkar@meolipc.hmg.inpg.fr:DRAKKAR/$CONFIG/$CONFCASE/DATA/

echo Done.


