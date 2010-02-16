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
#PBS -N metamoy_jade
#PBS -l select=1:ncpus=8:mpiprocs=1
#PBS -l walltime=02:00:00
#PBS -l place=scatter:excl
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
if [ ! -f $TMPDIR/cdfmoy ] ; then scp $USER@${login_node}:$CDFTOOLS/cdfmoy $TMPDIR/. ; fi ;
if [ ! -f $TMPDIR/cdfmoy_annual ] ; then scp $USER@${login_node}:$CDFTOOLS/cdfmoy_annual $TMPDIR/. ; fi ; 

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
cd $TMPDIR   
SDIR=$PREF/${CONFIG}/${CONFCASE}-S/$YEAR
   MDIR=$PREF/${CONFIG}/${CONFCASE}-MEAN/$YEAR
   chkdirg $MDIR
   chkdir $YEAR
   cd $YEAR
 # Monthly mean
 #
 for grid in gridT gridU gridV gridW icemod  ; do 
   for month in 01 02 03 04 05 06 07 08 09 10 11 12  ; do 
    getmonth $month $grid
     
   ../cdfmoy ${CONFCASE}_y${YEAR}m${month}d??_$grid.nc
   mv -f cdfmoy.nc ${CONFCASE}_y${YEAR}m${month}_$grid.nc
   mv -f cdfmoy2.nc ${CONFCASE}_y${YEAR}m${month}_${grid}2.nc
   \rm ${CONFCASE}_y${YEAR}m${month}d??_$grid.nc
   done

   ../cdfmoy_annual ${CONFCASE}_y${YEAR}m??_$grid.nc
   mv -f cdfmoy_annual.nc ${CONFCASE}_y${YEAR}_ANNUAL_$grid.nc
   ../cdfmoy_annual ${CONFCASE}_y${YEAR}m??_${grid}2.nc #; putannual ${grid}2 ;;
   mv -f cdfmoy_annual.nc ${CONFCASE}_y${YEAR}_ANNUAL_${grid}2.nc
 done
done
