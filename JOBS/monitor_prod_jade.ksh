#!/bin/ksh
set -x
# This script is intended to be sourced from a main script. Not Stand Alone
# Basically it runs on the production machine, once the MEAN fields 
# have been computed (monthly, annual) and disposed on the respective 
# CONFIG-CASE-MEAN/YEAR/ directory.

# Each block corresponds to a particular monitoring task. Each block is supposed
# to be independant from the other (in particular, required file are downloaded
# via the rapatrie function, which does the job only if necessary.

# The different tasks are performed with the cdftools programs. BIN is 
# added to the PATH.

#-------------------------------------------------------------------------------
#  $Rev: 231 $
#  $Date: 2009-03-24 11:25:04 +0100 (mar, 24 mar 2009) $
#  $Id: monitor_prod.ksh 231 2009-03-24 10:25:04Z molines $
#-------------------------------------------------------------------------------
YEAR=$1
# define some config dependent variable 
. ./config_def.ksh    # can be a link
#  Define some functions to get/put file from/to gaya (can be easily customized)
. ./function_def.ksh  # can be a link
login_node=service3
CONFCASE=$CONFIG-$CASE
#------------------------------------------------------------------------------
# directory name frequently used:
#------------------------------------------------------------------------------
  # on gaya
  MEANY=$CONFIG/${CONFIG}-${CASE}-MEAN/$YEAR
  SDIRY=$CONFIG/${CONFIG}-${CASE}-S/$YEAR
  DIAGS=${CONFIG}/${CONFIG}-${CASE}-DIAGS
  IDIR=$CONFIG/${CONFIG}-I
  BIN=/scratch/$USER/bin
  # on zahir
  P_CTL=$HOME/RUN_${CONFIG}/${CONFIG}-${CASE}/CTL
  # check existence of some required directories
  # ... on WORKDIR
  
  chkdir ../${CONFIG}-${CASE}-DIAGS
  chkdirg $SDIR/$DIAGS

  R_MONITOR=`pwd`
  cd $YEAR
#-----------------------------------------------------------------------------
# MENU SKEL
#-----------------------------------------------------------------------------
EKE=EEKKEE                    # compute EKE
RMSSSH=RRMMSS                 # compute RMS ssh and w
TSMEAN=TTSSMMEEAANN                 # compute TSMEAN and ssh drift
ICE=IICCEE                    # compute ice volume, area and extent
ICEMONTH=IICCEEMM               # compute ice volume, area and extent
GIB=0                    # compute Gibraltar diags (restoring zone)
ELNINO=0                 # compute El Nino monitoring SSTs
TRP=TTRRPP                    # compute barotropic transport accross section as given in section.dat (CTL dir)
MHT=MMHHTT                    # compute Meridional Heat Transport (advective and from surface fluxes)
MOC=MMOOCC                    # compute MOC ( need a sub basin mask file called new_maskglo.nc)
MAXMOC=MMAAXXMOC                 # diagnose the min and max of MOC
BSF=BBSSFF                    # compute the BSF (psi) from U and V
DCT=DDCCTT                    # compute density class transports for section given in dens_section.dat (CTL dir)
MXL=MMXXLL                    # Compute mixed layer depth from 3 criteria for month 03 and 09
TRACER=TTRRCC                 # Compute passive Tracer statistics
LSPV=LLSSPPVV                   # compute large scale potential vorticity in March and September
#-----------------------------------------------------------------------------
# PATH:
#-----------------------------------------------------------------------------

# check if required cdftools are available, exit if missing
  err=0
  for cdfprog in cdfeke cdfmean cdfrmsssh cdfstdevw cdficediags cdftransportiz\
                  cdfmhst cdfhflx cdfmoc cdfmaxmoc  cdfpsi  cdfsigtrp cdfmxl \
                  cdfzonalmean cdfzonalsum cdfzonalout bimgmoy4 bimgcaltrans cdfmoy ; do
     if [ ! -x /scratch/$USER/bin/$cdfprog ] ; then
       err=$(( err + 1 ))
       echo $cdfprog executable missing. Check your ~/bin for BIN installation
     fi
  done

  if [ $err != 0 ] ; then 
     echo " monitoring cannot be performed, sorry !" ; exit 1 
  fi
#=============================================================================
#  PART I: Derived quantities, to be added to the -MEAN/YEAR directory
#=============================================================================
   # check if we have a NATL config or a ORCA config (to be improved ....)
   atl=$( echo 1 | awk '{ ii=index (config,"NATL") ; print ii  }' config=$CONFIG )
############################################################################################
############################################################################################
# EKE : Eddy Kinetic Energy: Input files gridU, gridV gridU2, gridV2 
#^^^^^^^^^^^^^^^^^^^^^^^^^^
  if [ $EKE == 1 ] ; then
   # retrieve U and V ANNUAL mean files and squared mean
   getannualmean gridU
   getannualmean gridV
   getannualmean gridU2
   getannualmean gridV2
   # retrieve a T file needed for headers only (EKE is computed on the T-point)
   getannualmean  gridT2
   # run cdfeke
   $BIN/cdfeke ${CONFCASE}_y${YEAR}_ANNUAL_gridU.nc \
     ${CONFCASE}_y${YEAR}_ANNUAL_gridU2.nc \
     ${CONFCASE}_y${YEAR}_ANNUAL_gridV.nc \
     ${CONFCASE}_y${YEAR}_ANNUAL_gridV2.nc \
     ${CONFCASE}_y${YEAR}_ANNUAL_gridT2.nc
   mv eke.nc ${CONFCASE}_y${YEAR}_ANNUAL_EKE.nc
   # dispose file on the MEAN directory
   savemeanfile ${CONFCASE}_y${YEAR}_ANNUAL_EKE.nc
  fi
############################################################################################
############################################################################################
# RMS SSH and StdDev W : Input files : gridT, gridT2  gridW, gridW2
#^^^^^^^^^^^^^^^^^^^^^^^
  if [ $RMSSSH == 1 ] ; then 
   # RMSSSH :get gridT gridT2
   getannualmean gridT
   getannualmean gridT2
   # run cdfrmsssh
   $BIN/cdfrmsssh  ${CONFCASE}_y${YEAR}_ANNUAL_gridT.nc ${CONFCASE}_y${YEAR}_ANNUAL_gridT2.nc
   mv rms.nc ${CONFCASE}_y${YEAR}_ANNUAL_RMSSSH.nc
   # dispose file on the MEAN directory
   savemeanfile ${CONFCASE}_y${YEAR}_ANNUAL_RMSSSH.nc
   #####################################################
   # StdDev W :get gridW and gridW2 files
   getannualmean gridW
   getannualmean gridW2   
   # run cdfstdevw
   $BIN/cdfstdevw  ${CONFCASE}_y${YEAR}_ANNUAL_gridW.nc ${CONFCASE}_y${YEAR}_ANNUAL_gridW2.nc
   mv rmsw.nc ${CONFCASE}_y${YEAR}_ANNUAL_STDEVW.nc
   # dispose file on the MEAN directory
   savemeanfile ${CONFCASE}_y${YEAR}_ANNUAL_STDEVW.nc
   fi
############################################################################################
############################################################################################
# Barotropic Transport: Input file: gridU, gridV mesh mask
#^^^^^^^^^^^^^^^^^^^^^
  if [ $BSF == 1 ] ; then
   # get gridU gridV files
   getannualmean gridU
   getannualmean gridV  
   # get mesh mask files 
   getmask      ${MESH_MASK_ID}_byte_mask.nc
   getmesh_hgr  ${MESH_MASK_ID}_mesh_hgr.nc
   getmesh_zgr  ${MESH_MASK_ID}_mesh_zgr.nc
   # run cdfpsi
   $BIN/cdfpsi ${CONFCASE}_y${YEAR}_ANNUAL_gridU.nc ${CONFCASE}_y${YEAR}_ANNUAL_gridV.nc 
   mv psi.nc ${CONFCASE}_y${YEAR}_ANNUAL_PSI.nc
   # dispose and rename on the MEAN directory
   savemeanfile ${CONFCASE}_y${YEAR}_ANNUAL_PSI.nc
  fi
############################################################################################
############################################################################################
# MOC Meridional Overturning Circulation:  Input file: gridV, mesh mask, mask_glo
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  if [ $MOC == 1 ] ; then
   # get gridV  files
   getannualmean gridV 
   # get mesh mask files
   getmask      ${MESH_MASK_ID}_byte_mask.nc
   getmesh_hgr  ${MESH_MASK_ID}_mesh_hgr.nc
   getmesh_zgr  ${MESH_MASK_ID}_mesh_zgr.nc
   if (( $atl == 0 )) ; then rapatrie  new_maskglo.nc $IDIR new_maskglo.nc ; fi
   # run cdfmoc
   $BIN/cdfmoc ${CONFCASE}_y${YEAR}_ANNUAL_gridV.nc
   mv moc.nc ${CONFCASE}_y${YEAR}_ANNUAL_MOC.nc 
   # dispose on gaya MEAN/YEAR directory
   savemeanfile ${CONFCASE}_y${YEAR}_ANNUAL_MOC.nc
  fi
############################################################################################
############################################################################################
# Mixed Layer Diagnostics : Input file : gridT for month 03 and 09 mesh_hgr, mesh_zgr
#^^^^^^^^^^^^^^^^^^^^^^^^
  if [ $MXL == 1 ] ; then
   # get mesh mask files
   getmesh_hgr  ${MESH_MASK_ID}_mesh_hgr.nc
   getmesh_zgr  ${MESH_MASK_ID}_mesh_zgr.nc
   # get mean file and run cdfmxl
   getmonthlymean gridT 03
   $BIN/cdfmxl ${CONFCASE}_y${YEAR}m03_gridT.nc
   mv mxl.nc ${CONFCASE}_y${YEAR}m03_MXL.nc
   getmonthlymean gridT 09
   $BIN/cdfmxl ${CONFCASE}_y${YEAR}m09_gridT.nc
   mv mxl.nc ${CONFCASE}_y${YEAR}m09_MXL.nc
   # dispose on gaya, MEAN/YEAR directory
   savemeanfile ${CONFCASE}_y${YEAR}m03_MXL.nc
   savemeanfile ${CONFCASE}_y${YEAR}m09_MXL.nc
  fi
############################################################################################
############################################################################################
#=============================================================================
#  PART II: Time series: compute some integral quantities relevant for monitor
#           the ocean variability, and the behaviour of the on going run. 
#           Output is basically a small ASCII file, from which a matlab
#           suitable input file  (.mtl) is derived.
#=====================================================================
# Global MEANS: T S SSH Input files: gridT , mesh_hgr, mesh_zgr, mask
#^^^^^^^^^^^^^^
  if [ $TSMEAN == 1 ] ; then
   # get mesh mask files
   getmask      ${MESH_MASK_ID}_byte_mask.nc
   getmesh_hgr  ${MESH_MASK_ID}_mesh_hgr.nc
   getmesh_zgr  ${MESH_MASK_ID}_mesh_zgr.nc
   # get gridT files
   getannualmean gridT 
   # set header on the output file (ASCII)
   fsshmean=${CONFCASE}_y${YEAR}_SSHMEAN.txt
   ftmean=${CONFCASE}_y${YEAR}_TMEAN.txt
   fsmean=${CONFCASE}_y${YEAR}_SMEAN.txt
   echo $YEAR >  $fsshmean ; echo $YEAR >  $ftmean ;  echo $YEAR >  $fsmean
   # 3D means
   $BIN/cdfmean  ${CONFCASE}_y${YEAR}_ANNUAL_gridT.nc sossheig T >> $fsshmean
   $BIN/cdfmean  ${CONFCASE}_y${YEAR}_ANNUAL_gridT.nc votemper T >> $ftmean
   $BIN/cdfmean  ${CONFCASE}_y${YEAR}_ANNUAL_gridT.nc vosaline T >> $fsmean
  
   # dispose ASCII file in the -DIAGS directory 
   savediagfile ${CONFCASE}_y${YEAR}_SSHMEAN.txt
   savediagfile ${CONFCASE}_y${YEAR}_TMEAN.txt
   savediagfile ${CONFCASE}_y${YEAR}_SMEAN.txt
   
   if [ $(chkfile $DIAGS/LEVITUS_y0000_TMEAN.txt ) == absent ] ; then
    # first time : Create header with Levitus equivalent
    # requires  LEVITUS 'same' diags (from the ANNUAL mean )
    T_levitus=Levitus_annual_votemper.nc
    S_levitus=Levitus_annual_vosaline.nc

    getlevitus votemper
    getlevitus vosaline
##    rapatrie $levitus $IDIR $levitus
    $BIN/cdfmean ${T_levitus}  votemper T  >  LEVITUS_y0000_TMEAN.txt
    $BIN/cdfmean ${S_levitus}  vosaline T  >  LEVITUS_y0000_SMEAN.txt
    savediagfile  LEVITUS_y0000_TMEAN.txt
    savediagfile  LEVITUS_y0000_SMEAN.txt
    fi
  fi
############################################################################################
############################################################################################
# Ice Volume area and extent for m02 m03   m08 m09: input file : icemod, and mesh_mask
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  if [ $ICE == 1 ] ; then
   # get icemod file for the month 02 03 and 08  09
   getmonthlymean icemod 02
   getmonthlymean icemod 03
   getmonthlymean icemod 08
   getmonthlymean icemod 09
   # get mesh mask files
   getmask      ${MESH_MASK_ID}_byte_mask.nc
   getmesh_hgr  ${MESH_MASK_ID}_mesh_hgr.nc
   getmesh_zgr  ${MESH_MASK_ID}_mesh_zgr.nc 
   # Ascii output file:
   fice=${CONFCASE}_y${YEAR}_ice.txt
   echo '###' $YEAR 02 > $fice
   $BIN/cdficediags ${CONFCASE}_y${YEAR}m02_icemod.nc  >> $fice
   echo '###' $YEAR 03 >> $fice
   $BIN/cdficediags ${CONFCASE}_y${YEAR}m03_icemod.nc  >> $fice
   echo '###' $YEAR 08 >> $fice
   $BIN/cdficediags ${CONFCASE}_y${YEAR}m08_icemod.nc  >> $fice
   echo '###' $YEAR 09 >> $fice
   $BIN/cdficediags ${CONFCASE}_y${YEAR}m09_icemod.nc  >> $fice
   # save ascii file for plot
   savediagfile ${CONFCASE}_y${YEAR}_ice.txt
  fi
############################################################################################
############################################################################################
# Ice Volume area and extent for all months: input file : icemod, and mesh_mask
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  if [ $ICEMONTH == 1 ] ; then
   # get icemod files
   m=1
   while (( $m <= 12 )) ; do
    mm=$( printf "%02d" $m )
    getmonthlymean icemod $mm
    m=$(( m + 1 ))
   done
   # get mesh mask files
   getmask      ${MESH_MASK_ID}_byte_mask.nc
   getmesh_hgr  ${MESH_MASK_ID}_mesh_hgr.nc
   getmesh_zgr  ${MESH_MASK_ID}_mesh_zgr.nc
   # Ascii output file:
   fice=${CONFCASE}_y${YEAR}_icemonth.txt
   m=1
   while (( $m <= 12 )) ; do
    mm=$( printf "%02d" $m )
    case $mm in 
    01) echo '###' $YEAR $mm > $fice ;;
    *)  echo '###' $YEAR $mm >> $fice ;;
    esac
    $BIN/cdficediags ${CONFCASE}_y${YEAR}m${mm}_icemod.nc  >> $fice
    m=$(( m + 1 ))
   done
   rm -f out.txt 
   # save file txt for plot
   savediagfile ${CONFCASE}_y${YEAR}_icemonth.txt
  fi
############################################################################################
############################################################################################
#
#
#
#
#
#
# ADD GIB PART
# ADD NINO PART
#
#
#
#
#
#
#
############################################################################################
############################################################################################
# Transport: Input files: VT, gridU, gridV, mesh mask, section.dat
#^^^^^^^^^^^
  if [ $TRP == 1 ] ; then
   # section.dat describes the position (I,J) of the sections to monitor
   scp $USER@$login_node:$P_CTL/section.dat . 
   # get VT , gridU, gridV files
   getannualmean VT
   getannualmean gridU
   getannualmean gridV
   # get mesh mask files
   getmask      ${MESH_MASK_ID}_byte_mask.nc
   getmesh_hgr  ${MESH_MASK_ID}_mesh_hgr.nc
   getmesh_zgr  ${MESH_MASK_ID}_mesh_zgr.nc 
   # Ascii output file:
   fsection=${CONFCASE}_y${YEAR}_section_monitor.txt
   echo $YEAR > $fsection
   $BIN/cdftransportiz ${CONFCASE}_y${YEAR}_ANNUAL_VT.nc \
                  ${CONFCASE}_y${YEAR}_ANNUAL_gridU.nc \
                  ${CONFCASE}_y${YEAR}_ANNUAL_gridV.nc  < section.dat >> $fsection
   # eliminate garbage from txt file ...
   grep -v Give $fsection | grep -v level | grep -v IMAX | grep -v FROM > tmp
   mv -f tmp $fsection
   rm -f strp.txt vtrp.txt htrp.txt
   # save txt file for plot
   savediagfile ${CONFCASE}_y${YEAR}_section_monitor.txt
  fi
############################################################################################
############################################################################################ 
# Heat and Salt Meridional Transport : Input files : VT, mesh mask, new_maskglo
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  if [ $MHT == 1 ] ; then
# (a) From advection:
#--------------------
   # get VT  files
   getannualmean VT
   # get mesh mask files
   getmask      ${MESH_MASK_ID}_byte_mask.nc
   getmesh_hgr  ${MESH_MASK_ID}_mesh_hgr.nc
   getmesh_zgr  ${MESH_MASK_ID}_mesh_zgr.nc 
#
#
#
#
#
# need ADD MASKGLO
#
#
#
#
#
#
   # Ascii output file:
   fheat=${CONFCASE}_y${YEAR}_heattrp.dat
   fsalt=${CONFCASE}_y${YEAR}_salttrp.dat
   # run cdfmhst
   $BIN/cdfmhst  ${CONFCASE}_y${YEAR}_ANNUAL_VT.nc
   cp zonal_heat_trp.dat ${CONFCASE}_y${YEAR}_heattrp.dat
   cp zonal_salt_trp.dat ${CONFCASE}_y${YEAR}_salttrp.dat
 
# (b) from Surface Heat fluxes
#-----------------------------
    getannualmean gridT
    # run cdfhflx
    $BIN/cdfhflx  ${CONFCASE}_y${YEAR}_ANNUAL_gridT.nc
    mv hflx.out ${CONFCASE}_y${YEAR}_hflx.dat
    rm -f mhst.nc
    # save dat file for plot
    savediagfile ${CONFCASE}_y${YEAR}_hflx.dat
    savediagfile ${CONFCASE}_y${YEAR}_salttrp.dat
    savediagfile ${CONFCASE}_y${YEAR}_heattrp.dat
  fi
############################################################################################
############################################################################################
# MAX and MIN of MOC: requires that MOC files already exists
#^^^^^^^^^^^^^^^^^^^^
  if [ $MAXMOC == 1  ] ; then
   getannualmean MOC
   f=${CONFCASE}_y${YEAR}_ANNUAL_MOC.nc
   # Ascii output file
   fmaxmoc=${CONFCASE}_y${YEAR}_minmaxmoc.txt
   echo $YEAR > $fmaxmoc
   fmaxmoc40=${CONFIG}-${CASE}_y${YEAR}_maxmoc40.txt
   echo $YEAR > $fmaxmoc40

   if (( atl == 0 )) ; then
   # GLO
   printf "%s" 'Glo ' >>  $fmaxmoc ; $BIN/cdfmaxmoc $f glo 20 60 500 2000 | grep Maximum >> $fmaxmoc
   printf "%s" 'Glo ' >>  $fmaxmoc ; $BIN/cdfmaxmoc $f glo -40 30 2000 5500 | grep Minimum >> $fmaxmoc
   # ATL
   printf "%s" 'Atl ' >>  $fmaxmoc ; $BIN/cdfmaxmoc $f atl 0 60 500 2000 | grep Maximum >> $fmaxmoc
   printf "%s" 'Atl ' >>  $fmaxmoc ; $BIN/cdfmaxmoc $f atl -20 40 2000 5500 | grep Minimum  >> $fmaxmoc
   #INP
   printf "%s" 'Inp ' >>  $fmaxmoc ; $BIN/cdfmaxmoc $f inp 15 50 100 1000 | grep Minimum >> $fmaxmoc
   printf "%s" 'Inp ' >>  $fmaxmoc ; $BIN/cdfmaxmoc $f inp -30 20 1000 5500  | grep Minimum >> $fmaxmoc
   #AUS
   printf "%s" 'Aus ' >>  $fmaxmoc ; $BIN/cdfmaxmoc $f glo -70 0 0 2000   | grep Maximum >> $fmaxmoc
   printf "%s" 'Aus ' >>  $fmaxmoc ; $BIN/cdfmaxmoc $f glo -70 0 2000 5500  | grep Minimum >> $fmaxmoc
   # save file for plot
   savediagfile $fmaxmoc

   # Max and Min of MOC at some specific latitudes
   # GLO  MAX at 40 N and 30S
   printf "%s" 'Glo ' >>  $fmaxmoc40 ; $BIN/cdfmaxmoc $f glo 40 40 500 2000 | grep Maximum >> $fmaxmoc40
   printf "%s" 'Glo ' >>  $fmaxmoc40 ; $BIN/cdfmaxmoc $f glo -30 -30 500  5500 | grep Maximum >> $fmaxmoc40
   # ATL  MAX at 40N and 30S
   printf "%s" 'Atl ' >>  $fmaxmoc40 ; $BIN/cdfmaxmoc $f atl 40 40 500 2000 | grep Maximum >> $fmaxmoc40
   printf "%s" 'Atl ' >>  $fmaxmoc40 ; $BIN/cdfmaxmoc $f atl -30 -30  500 5000 | grep Maximum >> $fmaxmoc40
   #INP  Min at 30 S
   printf "%s" 'Inp ' >>  $fmaxmoc40 ; $BIN/cdfmaxmoc $f inp -30 -30 1000 5500  | grep Minimum >> $fmaxmoc40
   #AUS  MAX at 50 S
   printf "%s" 'Aus ' >>  $fmaxmoc40 ; $BIN/cdfmaxmoc $f glo -50 -50 0 2000   | grep Maximum >> $fmaxmoc40
   savediagfile $fmaxmoc40

   else    # NATL configuration
   # GLO
   printf "%s" 'Glo ' >>  $fmaxmoc ; $BIN/cdfmaxmoc $f glo 20 60 500 2000 | grep Maximum >> $fmaxmoc
   printf "%s" 'Glo ' >>  $fmaxmoc ; $BIN/cdfmaxmoc $f glo -40 30 2000 5500 | grep Minimum >> $fmaxmoc
   savediagfile $fmaxmoc

   # Max and Min of MOC at some specific latitudes
   # GLO  MAX at 40 N and 30S
   printf "%s" 'Glo ' >>  $fmaxmoc40 ; cdfmaxmoc $f glo 40 40 500 2000 | grep Maximum >> $fmaxmoc40
   printf "%s" 'Glo ' >>  $fmaxmoc40 ; cdfmaxmoc $f glo -15 -15 500  5500 | grep Maximum >> $fmaxmoc40
   savediagfile $fmaxmoc40
   fi
  fi

############################################################################################
############################################################################################
# DCT :Density Class transport: Input files : gridT, gridU gridV, mesh mask, dens_section.dat
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  if [ $DCT == 1 ] ; then
  # dens_section.dat describe the sections (either zonal or meridional) where the DCT is computed
  scp $USER@$login_node:$P_CTL/dens_section.dat .
  scp $USER@$login_node:$CDFTOOLS/JOBS/trpsig_postproc_$MACHINE.ksh trpsig_postproc.ksh
   # get mesh mask files
   getmask      ${MESH_MASK_ID}_byte_mask.nc
   getmesh_hgr  ${MESH_MASK_ID}_mesh_hgr.nc
   getmesh_zgr  ${MESH_MASK_ID}_mesh_zgr.nc

  # Required post_processing script: DCT are computed on monthly means, then average is performed
  # for annual values. This process is still done through temporary bimg/dimg files (remnant of the
  # old Clipper times).  By the way, 2 bimgtools are required: bimgmoy4 and bimgcaltrans
  #  In-lining of this script may be confusing. I leave it as an external module.

  # due to the large amount of files that are produced by this diags, we prefer to keep them
  # on a separate directory
  chkdirg $SDIR/$DIAGS/TRPSIG/
  # also need temporary directories in the actual tmpdir:
  for  m in  1 2 3 4 5 6 7 8 9 10 11 12 ; do
    mm=$(printf "%02d" $m)
    tfich=${CONFCASE}_y${YEAR}m${mm}_gridT.nc 
    ufich=$(echo  $tfich | sed -e 's/gridT/gridU/' )
    vfich=$(echo  $tfich | sed -e 's/gridT/gridV/' )
    #get files on gaya
    getmonthlymean gridT $mm
    getmonthlymean gridU $mm
    getmonthlymean gridV $mm
    #retrieve tag time from file name
    tag=$(echo $tfich | sed -e "s/${CONFCASE}_//" -e 's/_gridT.nc//')
    echo $tag > ${CONFCASE}_${tag}_trpsig_monitor.lst
    # run cdfsigtrp
    $BIN/cdfsigtrp $tfich $ufich $vfich 21 30 180 -bimg -print  >>  ${CONFCASE}_${tag}_trpsig_monitor.lst
    # rename file
    for  b in [0-9][0-9]*bimg ; do
        mv  $b ${CONFCASE}_${tag}_$b
    done
    # Idem: for txt files
    mv trpsig.txt ${CONFCASE}_${tag}_trpsig.txt
  done
  # use function def in function_def 
  reorganize
  mean
  # save results on gaya ( as many files as sections in dens_section.dat)
  for f in TRPSIG/${CONFCASE}_y*_trpsig.txt ; do
    savediagfile $f
  done
  fi
# TRACER DIAGS (31/12 of each year) : Input files : ptrcT, mesh mask
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
#
#
# NOT IMPEMENTED FOR JADE
#
#
#
#
#################################################################
# cleaning of directory
rmannualmean gridW
rmannualmean gridW2
rmannualmean gridU
rmannualmean gridU2
rmannualmean gridV
rmannualmean gridV2
rmannualmean gridT
rmannualmean gridT2

rmmonthlymean gridU ??
rmmonthlymean gridU2 ??
rmmonthlymean gridV ??
rmmonthlymean gridV2 ??
rmmonthlymean gridT ??
rmmonthlymean gridT2 ??
rmmonthlymean icemod ??
rmmonthlymean gridW ??

rmmask
rmmesh_hgr
rmmesh_zgr
#######################################################################
#save data on /data and keep it here also for plot
# touch OK file for dplace function
pwd
touch OK_MONITOR

