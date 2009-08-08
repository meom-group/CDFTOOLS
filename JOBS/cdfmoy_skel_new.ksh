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


#################################################################################
# This script is used to compute time mean averages for DRAKKAR model output.
# It replaces an older script which was also computing quarterly means.
# All customisable variable are set in Part I.
# This script must be launched from metamoy.ksh which edit the years
#
# $Rev$
# $Date$
# $Id$
################################################################################

set -x
. $HOME/.profile
P_CDF_DIR=$PDIR/RUN_CCOONNFF/CCOONNFF-CCAASSEE/CTL/CDF
. $P_CDF_DIR/config_def.ksh

cp $P_CDF_DIR/config_def.ksh $TMPDIR
cp $P_CDF_DIR/function_def.ksh $TMPDIR
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
LOCAL_SAVE=${LOCAL_SAVE:=0}

YEARLST=""
y=$YEARS
while (( $y <= $YEARE )) ; do 
  YEARLST="$YEARLST $y "
  y=$(( y + 1 ))
done

# 
CONFCASE=${CONFIG}-${CASE}

# always work in TMPDIR ! not in the data dir as file will be erased at the end of the script !
cd $TMPDIR
mkdir MONTHLY
   if [ $LOCAL_SAVE = 1 ] ; then
    chkdir $WORKDIR/$CONFIG
    chkdir $WORKDIR/$CONFIG/${CONFCASE}-MEAN
   fi

for YEAR in $YEARLST ; do 
   SDIR=${CONFIG}/${CONFCASE}-S/$YEAR
   MDIR=$PREF/${CONFIG}/${CONFCASE}-MEAN/$YEAR
   chkdirg $MDIR
   if [ $LOCAL_SAVE = 1 ] ; then
    chkdir $WORKDIR/$CONFIG/${CONFCASE}-MEAN/YEAR
   fi

 # Monthly mean
 #
 for grid in gridT gridU gridV gridW icemod  ; do 
   for month in 01 02 03 04 05 06 07 08 09 10 11 12  ; do 
    getmonth $month $grid
    $CDFTOOLS/cdfmoy ${CONFCASE}_y${YEAR}m${month}d??_$grid.nc
    case $grid in 
      icemod) putmonth $month $grid ;;
           *) putmonth $month $grid ;
              putmonth2 $month ${grid}2 ;;
    esac
#  \rm ${CONFCASE}_y${YEAR}m${month}d??_$grid.nc
   done
    

   # all monthes done for given grid, can compute annual mean ...for grid and grid2
   #  suppose 5 day averages when creating monthly mean
   cd MONTHLY
    case $grid in 
      icemod) $CDFTOOLS/cdfmoy_annual ${CONFCASE}_y${YEAR}m??_$grid.nc    ; putannual $grid ;;
           *) $CDFTOOLS/cdfmoy_annual ${CONFCASE}_y${YEAR}m??_$grid.nc    ; putannual $grid ;
              $CDFTOOLS/cdfmoy_annual ${CONFCASE}_y${YEAR}m??_${grid}2.nc ; putannual ${grid}2 ;;
    esac
    # clean MONTHLY from grid and grid2 files
    \rm ${CONFCASE}_y${YEAR}m??_$grid.nc ; \rm ${CONFCASE}_y${YEAR}m??_${grid}2.nc
  cd $TMPDIR
 done
done
