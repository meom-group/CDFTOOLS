#!/bin/ksh
# @ wall_clock_limit = 10:00:00
# @ job_name   = moy-YYYY
# @ as_limit = 1gb
# @ output = $(job_name).$(step_name).$(jobid)
# @ error = $(output)
# @ notify_user = molines@hmg.inpg.fr
# @ notification = error

# @ step_name = cdfmoy1
# @ job_type = serial
# @ wall_clock_limit = 7200
# @ data_limit = 0.8Gb
# @ queue

# @ step_name = cdfvt2
# @ job_type = serial
# @ wall_clock_limit = 7200
# @ data_limit = 0.8Gb
# @ queue

# @ step_name = monitor3
# @ dependency = (cdfmoy1 == 0 && cdfvt2 == 0 )
# @ job_type = serial
# @ wall_clock_limit = 7200
# @ data_limit = 0.8Gb
# @ queue

# @ step_name = clean4
# @ dependency = (monitor3 == 0 )
# @ job_type = serial
# @ wall_clock_limit = 7200
# @ data_limit = 0.8Gb
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
# $Rev: 262 $
# $Date: 2009-08-08 11:25:10 +0200 (Sat, 08 Aug 2009) $
# $Id: cdfmoy_skel_new.ksh 262 2009-08-08 09:25:10Z rcli002 $
################################################################################

TMPDIR0=$WORKDIR/METAMOY.$$
if [ ! -d $TMPDIR0 ] ; then mkdir $TMPDIR0 ;fi


set -x

case $LOADL_STEP_NAME in
 cdfmoy1 )
 TMPDIR=$TMPDIR0/CDFMOY
 if [ ! -d $TMPDIR ] ; then mkdir $TMPDIR ;fi
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
    chkdir $WORKDIR/$CONFIG/${CONFCASE}-MEAN/$YEAR
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
 ;;

cdfvt2 ) 

 TMPDIR=$TMPDIR0/CDFVT
 if [ ! -d $TMPDIR ] ; then mkdir $TMPDIR ;fi
set -x
. $HOME/.profile
P_CDF_DIR=$PDIR/RUN_CCOONNFF/CCOONNFF-CCAASSEE/CTL/CDF
. $P_CDF_DIR/config_def.ksh

cp $P_CDF_DIR/config_def.ksh $TMPDIR
cp $P_CDF_DIR/function_def.ksh $TMPDIR
cd $TMPDIR


# Part I : setup config dependent names
#--------------------------------------
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
    chkdir $WORKDIR/$CONFIG/${CONFCASE}-MEAN/$YEAR
   fi

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
 cd $TMPDIR/MONTHLY
 $CDFTOOLS/cdfmoy_annual ${CONFCASE}_y${YEAR}m??_VT.nc
 putvtannual
 
 # clean directory for eventually next year:
 \rm ${CONFCASE}_y${YEAR}m??_VT.nc
 cd $TMPDIR
done
 ;;

monitor3 )
set -x
 TMPDIR=$TMPDIR0/MONITOR
 if [ ! -d $TMPDIR ] ; then mkdir $TMPDIR ;fi

. ./config_def.ksh   # config_def.ksh may be a link to an existing configuration file

# set the list of years you want to monitor 'at once'
yinit=YYYY              # initial year
yend=YYYE               # last year

YEARS=''
while (( $yinit <= $yend )) ; do
 YEARS="$YEARS $yinit "
 yinit=$(( yinit + 1 ))
done

MESH_MASK_ID='ORCA05-G70.112-no-caspian' # root part of the mesh-mask files
                         # (they must be in the -I directory ( $CONFIG/${CONFIG}-I)
                         #  Standard name is thus : ${MESH_MASK_ID}_byte_mask.nc
                         #                          ${MESH_MASK_ID}_mesh_hgr.nc
                         #                          ${MESH_MASK_ID}_mesh_zgr.nc
#
TSCLIM=''    # can be either ''        -> default to Levitus_p2.1
                    #              Gouretski
                    #              Levitus_p2.1

CDFTOOLS=~rcli002/CDFTOOLS-2.1   # PATH for the cdftools executables

# define the I-J window for GIB diags and El NINO DIAG
if [ $CONFIG = 'NATL025' ] ; then
  GIBWIN='338 353 239 260'
  # NOT RELEVANT FOR NATL025. Here for compatibility
  NINO12='790 830 459 499'
  NINO3='550 790 479 519 '
  NINO4='350 550 479 519 '
  NINO34='470 670 479 519 '
elif [ $CONFIG = 'ORCA025.L75' ] ; then
  GIBWIN='1094 1109 653 674 '
  NINO12='790 830 459 499'
  NINO3='550 790 479 519 '
  NINO4='350 550 479 519 '
  NINO34='470 670 479 519 '
elif [ $CONFIG = 'ORCA05' ] ; then
  GIBWIN='547 554 326 337 '
  NINO12='395 415 229 249'
  NINO3='275 395 239 259'
  NINO4='175 275 239 259'
  NINO34='235 335 239 259'
else
  echo GIBWIN and NINO boxes not defined for config $CONFIG
  exit 1
fi


# menu (set to 1 if you want it, to anything else if you do not !)
EKE=1                    # compute EKE
RMSSSH=1                 # compute RMS ssh and w
TSMEAN=1                 # compute TSMEAN and ssh drift
ICE=0                    # compute ice volume, area and extent
ICEMONTH=1               # compute ice volume, area and extent
GIB=1                    # compute Gibraltar diags (restoring zone)
ELNINO=1                 # compute El Nino monitoring SSTs
TRP=1                    # compute barotropic transport accross section as given in section.dat (CTL dir)
MHT=1                    # compute Meridional Heat Transport (advective and from surface fluxes)
MOC=1                    # compute MOC ( need a sub basin mask file called new_maskglo.nc)
MAXMOC=1                 # diagnose the min and max of MOC
BSF=1                    # compute the BSF (psi) from U and V
DCT=1                    # compute density class transports for section given in dens_section.dat (CTL dir)
MXL=1                    # Compute mixed layer depth from 3 criteria for month 03 and 09
TRACER=0                 # Compute passive Tracer statistics

#--------------------- nothing to touch below -----------------------------------------
# copy config and function to the working directory.
cp config_def.ksh $TMPDIR
cp function_def.ksh $TMPDIR

CONFCASE=${CONFIG}-${CASE}
cd $TMPDIR
. ./config_def.ksh

for   YEAR in  $YEARS ; do
   .  $CDFTOOLS/JOBS/monitor_prod.ksh
    cd $WORKDIR/$CONFIG/${CONFCASE}-MEAN/
     \rm -r $YEAR
     cd $TMPDIR
done
 rsh gaya "cd $CONFIG/${CONFCASE}-DIAGS ; ~/bin/mkmtl.ksh "

 ;;
clean4)
  cd $WORKDIR
  \rm -rf $TMPDIR0
 ;;
esac
