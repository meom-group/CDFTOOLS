#!/bin/ksh

set -x
P_CDF_DIR=$HOME/RUN_CCOONNFF/CCOONNFF-CCAASSEE/CTL/CDF

cd $TMPDIR
if [ ! -d MOY ] ; then mkdir MOY ; fi

MOYDIR=$TMPDIR/MOY
cd $MOYDIR
cp $P_CDF_DIR/config_def.ksh $MOYDIR
cp $P_CDF_DIR/function_def.ksh $MOYDIR

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

# always work in MOYDIR ! not in the data dir as file will be erased at the end of the script !
cd $MOYDIR
mkdir MONTHLY

for YEAR in $YEARLST ; do 
   SDIR=${CONFIG}/${CONFCASE}-S/$YEAR
   MDIR=$PREF/${CONFIG}/${CONFCASE}-MEAN/$YEAR
   chkdirg $MDIR

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
   \rm ${CONFCASE}_y${YEAR}m${month}d??_$grid.nc
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
    mv ${CONFCASE}_y${YEAR}m??_$grid.nc $TMPDIR ; mv ${CONFCASE}_y${YEAR}m??_${grid}2.nc $TMPDIR
  cd $MOYDIR
 done
done

