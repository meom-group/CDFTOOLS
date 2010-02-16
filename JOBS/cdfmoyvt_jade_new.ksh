#!/bin/ksh
# @ wall_clock_limit = 10:00:00
# @ job_name   = moy-YYYY
# @ as_limit = 1gb
# @ output     = $(job_name).$(jobid)
# @ error      =  $(job_name).$(jobid)
# @ notify_user = molines@hmg.inpg.fr
# @ notification = error
# @ queue                   

### OAR is valid on ZEPHIR
#OAR -n metamoy
#OAR -l /nodes=1/cpu=1,walltime=5:00:00
#OAR -E METAMOY.%jobid%
#OAR -O METAMOY.%jobid%

### PBS is valid on JADE
#PBS -N metamoyvt_jade
#PBS -l select=1:ncpus=8:mpiprocs=8
#PBS -l walltime=4:30:00
#PBS -M molines@hmg.inpg.fr
#PBS -mb -me

#################################################################################
# This script is used to compute time mean averages for DRAKKAR model output.
# It replaces an older script which was also computing quarterly means.
# All customisable variable are set in Part I.
# This script must be launched from metamoy.ksh which edit the years
#
# $Rev: 229 $
# $Date: 2009-03-24 09:34:51 +0100 (mar, 24 mar 2009) $
# $Id: cdfmoy_skel_new.ksh 229 2009-03-24 08:34:51Z rcli002 $
################################################################################

set -x
cd
pwd
. $HOME/.profile
P_CDF_DIR=RUN_CCOONNFF/CCOONNFF-CCAASSEE/CTL/CDF
. $P_CDF_DIR/config_def.ksh
. $P_CDF_DIR/function_def.ksh

chkdir $TMPDIR

#test network
login_node=service1
chknet $login_node
echo "$login_node $NET"
if [ $NET = KO ] ; then
login_node=service2
chknet $login_node
echo "$login_node $NET"
fi;
if [ $NET = KO ] ; then
login_node=service3
chknet $login_node
echo "$login_node $NET"
fi;
if [ $NET = KO ] ; then 
echo network is KO
date
exit
fi

 scp $USER@${login_node}:$P_CDF_DIR/config_def.ksh $TMPDIR/.
 scp $USER@${login_node}:$P_CDF_DIR/function_def.ksh $TMPDIR/.
 scp $USER@${login_node}:$CDFTOOLS/JOBS/mkmoy_jade.ksh $TMPDIR/mkmoy.ksh
 scp $USER@${login_node}:$CDFTOOLS/JOBS/mkvt_jade.ksh $TMPDIR/mkvt.ksh
 scp $USER@${login_node}:$CDFTOOLS/JOBS/testOK.ksh $TMPDIR/testOK.ksh
 chmod 755 $TMPDIR/testOK.ksh

if [ ! -f $TMPDIR/cdfmoy ] ; then scp $USER@${login_node}:$CDFTOOLS/cdfmoy $TMPDIR/. ; fi ;
if [ ! -f $TMPDIR/cdfmoy_annual ] ; then scp $USER@${login_node}:$CDFTOOLS/cdfmoy_annual $TMPDIR/. ; fi ; 
if [ ! -f $TMPDIR/cdfvT ] ; then scp $USER@${login_node}:$CDFTOOLS/cdfvT $TMPDIR/. ; fi ; 

cd $TMPDIR

# Part I : setup config dependent names
#--------------------------------------
. ./config_def.ksh    # this file (or a link) must exist in the current directory
#
# Part II  define some usefull functions
#---------------------------------------
. ./function_def.ksh  # this file (or a link) must exist in the current directory

# Part III : main loops : no more customization below
#-----------------------------------------------------
# set up list of years to process
# Metamoy meta script will subtitute YYYY and YYYE with correct begining and ending years
YEARS=YYYY
YEARE=YYYE

YEARLST=""
y=$YEARS
while (( $y <= $YEARE )) ; do 
  YEARLST="$YEARLST $y "
  y=$(( y + 1 ))
done

# 
CONFCASE=${CONFIG}-${CASE}

# always work in TMPDIR ! not in the data dir as file will be erased at the end of the script !

for YEAR in $YEARLST ; do    
    SDIR=$PREF/${CONFIG}/${CONFCASE}-S/$YEAR
    MDIR=$PREF/${CONFIG}/${CONFCASE}-MEAN/$YEAR
    chkdirg $MDIR
    chkdir $YEAR
    rm -f $YEAR/OK?
 # Monthly mean
 #
    for grid in gridT gridU gridV gridW icemod  ; do 
############### CDFMOY CDFVT #####################################
        if [ $grid = gridT  ]; then NP=0 ; fi ;
        if [ $grid = gridU  ]; then NP=1 ; fi ;
        if [ $grid = gridV  ]; then NP=2 ; fi ;
        if [ $grid = gridW  ]; then NP=3 ; fi ;
        if [ $grid = icemod ]; then NP=4 ; fi ;
##################
        cd $TMPDIR
        cat mkmoy.ksh | sed -e "s/YYEEAARR/$YEAR/g" -e "s/GGRRIIDD/$grid/g" -e "s/NNPP/$NP/g" > $YEAR/tmp_mkmoy_$grid.ksh
        cd $YEAR
        chmod u+x tmp_mkmoy_$grid.ksh
        dplace -c$NP tmp_mkmoy_$grid.ksh > log_$grid 2>&1 &
    done
    NP=5
    cd $TMPDIR
    cat mkvt.ksh | sed -e "s/YYEEAARR/$YEAR/g" -e "s/NNPP/$NP/g" > $YEAR/tmp_mkvt.ksh
    cd $YEAR
    chmod u+x tmp_mkvt.ksh
    dplace -c$NP tmp_mkvt.ksh > log_VT 2>&1 &   #dplace -c$NP 
    dplace -c6 ../testOK.ksh > log_test1 2>&1 
    dplace -c7 ../testOK.ksh > log_test2 2>&1 
cd $TMPDIR
 # clean directory for eventually next year:
################# CDFVT #####################################
done

