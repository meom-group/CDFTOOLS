#!/bin/ksh

set -x
P_CDF_DIR=$HOME/RUN_CCOONNFF/CCOONNFF-CCAASSEE/CTL/CDF
cd $P_CDF_DIR

# Part I : setup config dependent names
#--------------------------------------
. ./config_def.ksh   # config_def.ksh may be a link to an existing configuration file

cd $TMPDIR
mkdir VT
VTDIR=$TMPDIR/VT
cd $VTDIR

cp $P_CDF_DIR/config_def.ksh $VTDIR
cp $P_CDF_DIR/function_def.ksh $VTDIR


. ./config_def.ksh   # config_def.ksh may be a link to an existing configuration file

# Part II  define some usefull functions
#---------------------------------------
. ./function_def.ksh # function_def.ksh may be a link to customizable function file

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

# always work in VTDIR ! not in the data dir as file will be erased at the end of the script !
cd $VTDIR
mkdir MONTHLY

for YEAR in $YEARLST ; do
   SDIR=${CONFIG}/${CONFCASE}-S/$YEAR
   MDIR=$PREF/${CONFIG}/${CONFCASE}-MEAN/$YEAR
   chkdirg $MDIR

 # Monthly mean
 #
 for month in 01 02 03 04 05 06 07 08 09 10 11 12  ; do
   getmonth $month gridT
   getmonth $month gridU
   getmonth $month gridV

   list=''
   for f in ${CONFCASE}_y${YEAR}m${month}d??_gridT.nc ; do
     tag=$( echo $f | awk -F_ '{print $2}' )
     list="$list $tag"
   done

   $CDFTOOLS/cdfvT $CONFCASE $list
   putvtmonth $month
   \rm ${CONFCASE}_y${YEAR}m${month}d??_grid[UVT].nc
 done

 # annual mean  (uses a ponderation to compute the exact annual mean ). ! suppose 5 day averages when creating monthly mean
 cd $VTDIR/MONTHLY
 $CDFTOOLS/cdfmoy_annual ${CONFCASE}_y${YEAR}m??_VT.nc
 putvtannual
 
 #  move to TMPDIR for monitoring
 mv ${CONFCASE}_y${YEAR}m??_VT.nc $TMPDIR
 mv ${CONFCASE}_y${YEAR}_ANNUAL_VT.nc $TMPDIR
 cd $VTDIR
done
