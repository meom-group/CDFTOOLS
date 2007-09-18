#!/bin/ksh
set -x
# This script is intended to be sourced from a main script. Not Stand Alone
# Basically it runs on the production machine, once the MEAN fields 
# have been computed (monthly, annual) and disposed on the respective 
# CONFIG-CASE-MEAN/YEAR/ directory.

# Each block corresponds to a particular monitoring task. Each block is supposed
# to be independant from the other (in particular, required file are downloaded
# via the rapatrie function, which does the job only if necessary.

# The different tasks are performed with the cdftools programs. CDFTOOLS is 
# added to the PATH.

#-------------------------------------------------------------------------------
#  $Rev$
#  $Date$
#  $Id$
#-------------------------------------------------------------------------------
#  Define some functions to get/put file from/to gaya (can be easily customized)

# rapatrie : Usage: rapatrie  remote_file directory local_file
#   if local_file already here do nothing, else mfget it from gaya,
#   directory/remote_file
rapatrie() { if [ ! -f $3 ] ; then mfget $2/$1 $3 ; else echo $3 is already \
            downloaded ; fi ; }

# expatrie : Usage:  expatrie local_file directory remote_file
#   put local file on gaya in directory/remote_file
#
expatrie() { mfput $1 $2/$3 ; }

# cptoweb : Usage: cptoweb  file.mtl
#    rcp the matlab file to the corresponding DATA dir of the website
cptoweb() { rcp $1 \
       apache@meolipc.hmg.inpg.fr:web/DRAKKAR/$CONFIG/$CONFCASE/DATA/ ; }

# chkfile : Usage: chkfile gaya_file
#    check if a file exists on gaya, return present or absent.
chkfile() { rsh gaya " if [ -f $1 ] ; then echo present ;\
                       else echo absent ; fi " ; }

# chkdirg  : Usage: chkdirg gaya_directory
#    check the existence of a directory on gaya. Create it if not present
chkdirg() { rsh gaya " if [ ! -d $1 ] ; then mkdir $1 ; fi " ; }

# chkdirw  : Usage: chkdirw web_site_directory
#   check the existence of a dir. on the web site. Create it if not present
chkdirw() { rsh meolipc.hmg.inpg.fr -l apache " if [ ! -d web/DRAKKAR/$1 ] ; \
            then mkdir web/DRAKKAR/$1 ; fi " ; }

# chkdir  : Usage: chkdir local_dir 
#   check the existence of a directory. Create it if not present
chkdir() { if [ ! -d $1 ] ; then mkdir $1 ; fi  ; }

#------------------------------------------------------------------------------
# directory name frequently used:
#------------------------------------------------------------------------------
  # on gaya
  MEANY=$CONFIG/${CONFCASE}-MEAN/$YEAR
  SDIRY=$CONFIG/${CONFCASE}-S/$YEAR
  DIAGS=${CONFIG}/${CONFCASE}-DIAGS
  MONITOR=${CONFIG}/${CONFCASE}-MONITORTEST
  IDIR=$CONFIG/${CONFIG}-I

  # on zahir
  P_CTL=$HOME/RUN_${CONFIG}/${CONFCASE}/CTL

  # check existence of some required directories
  # ... on gaya
  chkdirg $MONITOR
  chkdirg $DIAGS

  # ... on the web site
  chkdirw $CONFIG/
  chkdirw $CONFIG/$CONFCASE
  chkdirw $CONFIG/$CONFCASE/DATA/

#------------------------------------------------------------------------------
# PATH:
#-----------------------------------------------------------------------------
  export PATH=$CDFTOOLS/:$PATH

# check if required cdftools are available, exit if missing
  err=0
  for cdfprog in cdfeke cdfmean cdfrmsssh cdfstdevw cdficediags cdftransportiz\
                  cdfmhst cdfhflx cdfmoc cdfmaxmoc  cdfpsi  cdfsigtrp cdfmxl \
                  cdfzonalmean cdfzonalsum cdfzonalout bimgmoy4 bimgcaltrans ; do
     if [ ! -x $CDFTOOLS/$cdfprog ] ; then
       err=$(( err + 1 ))
       echo $cdfprog executable missing. Check your $CDFTOOLS installation
     fi
  done

  if [ $err != 0 ] ; then 
     echo " monitoring cannot be performed, sorry !" ; exit 1 
  fi
#=============================================================================
#  PART I: Derived quantities, to be added to the -MEAN/YEAR directory
#=============================================================================

# EKE : Eddy Kinetic Energy: Input files gridU, gridV gridU2, gridV2 
#^^^^^^^^^^^^^^^^^^^^^^^^^^
  if [ $EKE == 1 ] ; then
   # retrieve U and V ANNUAL mean files and squared mean
     rapatrie  ${CONFCASE}_y${YEAR}_ANNUAL_gridU.nc $MEANY ${CONFCASE}_y${YEAR}_ANNUAL_gridU.nc
     rapatrie  ${CONFCASE}_y${YEAR}_ANNUAL_gridU2.nc $MEANY ${CONFCASE}_y${YEAR}_ANNUAL_gridU2.nc
     rapatrie  ${CONFCASE}_y${YEAR}_ANNUAL_gridV.nc $MEANY ${CONFCASE}_y${YEAR}_ANNUAL_gridV.nc
     rapatrie  ${CONFCASE}_y${YEAR}_ANNUAL_gridV2.nc $MEANY ${CONFCASE}_y${YEAR}_ANNUAL_gridV2.nc
   
   # retrieve a T file needed for headers only (EKE is computed on the T-point)
   rapatrie  ${CONFCASE}_y${YEAR}_ANNUAL_gridT2.nc $MEANY  ${CONFCASE}_y${YEAR}_ANNUAL_gridT2.nc

   cdfeke ${CONFCASE}_y${YEAR}_ANNUAL_gridU.nc \
     ${CONFCASE}_y${YEAR}_ANNUAL_gridU2.nc \
     ${CONFCASE}_y${YEAR}_ANNUAL_gridV.nc \
     ${CONFCASE}_y${YEAR}_ANNUAL_gridV2.nc \
     ${CONFCASE}_y${YEAR}_ANNUAL_gridT2.nc

   # dispose file on the MEAN directory
   expatrie eke.nc $MEANY  ${CONFCASE}_y${YEAR}_EKE.nc
   \rm eke.nc
  fi


# RMS SSH and StdDev W : Input files : gridT, gridT2  gridW, gridW2
#^^^^^^^^^^^^^^^^^^^^^^^
  if [ $RMSSSH == 1 ] ; then 
   # RMSSSH :get gridT gridT2
   rapatrie  ${CONFCASE}_y${YEAR}_ANNUAL_gridT.nc $MEANY  ${CONFCASE}_y${YEAR}_ANNUAL_gridT.nc
   rapatrie  ${CONFCASE}_y${YEAR}_ANNUAL_gridT2.nc $MEANY  ${CONFCASE}_y${YEAR}_ANNUAL_gridT2.nc
   cdfrmsssh  ${CONFCASE}_y${YEAR}_ANNUAL_gridT.nc ${CONFCASE}_y${YEAR}_ANNUAL_gridT2.nc

   # dispose file on the MEAN directory
   expatrie rms.nc $MEANY  ${CONFCASE}_y${YEAR}_ANNUAL_RMSSSH.nc
   \rm rms.nc

   # StdDev W :get gridW and gridW2 files
   rapatrie ${CONFCASE}_y${YEAR}_ANNUAL_gridW.nc $MEANY ${CONFCASE}_y${YEAR}_ANNUAL_gridW.nc
   rapatrie ${CONFCASE}_y${YEAR}_ANNUAL_gridW2.nc $MEANY ${CONFCASE}_y${YEAR}_ANNUAL_gridW2.nc

   cdfstdevw  ${CONFCASE}_y${YEAR}_ANNUAL_gridW.nc ${CONFCASE}_y${YEAR}_ANNUAL_gridW2.nc

   # dispose file on the MEAN directory
   expatrie rmsw.nc $MEANY ${CONFCASE}_y${YEAR}_ANNUAL_STDEVW.nc
  \rm rmsw.nc
  fi

# Barotropic Transport: Input file: gridU, gridV mesh mask
#^^^^^^^^^^^^^^^^^^^^^
  if [ $BSF == 1 ] ; then
   # get gridU gridV files
   rapatrie ${CONFCASE}_y${YEAR}_ANNUAL_gridU.nc $MEANY ${CONFCASE}_y${YEAR}_ANNUAL_gridU.nc
   rapatrie ${CONFCASE}_y${YEAR}_ANNUAL_gridV.nc $MEANY ${CONFCASE}_y${YEAR}_ANNUAL_gridV.nc
 
   # get mesh mask files 
   rapatrie  ${MESH_MASK_ID}_byte_mask.nc $IDIR mask.nc
   rapatrie  ${MESH_MASK_ID}_mesh_hgr.nc $IDIR mesh_hgr.nc
   rapatrie  ${MESH_MASK_ID}_mesh_zgr.nc $IDIR mesh_zgr.nc
 
   cdfpsi ${CONFCASE}_y${YEAR}_ANNUAL_gridU.nc ${CONFCASE}_y${YEAR}_ANNUAL_gridV.nc 
 
   # dispose and rename on the MEAN directory
   expatrie psi.nc  $MEANY ${CONFCASE}_y${YEAR}_PSI.nc
  fi

# MOC Meridional Overturning Circulation:  Input file: gridV, mesh mask, mask_glo
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  if [ $MOC == 1 ] ; then
   # get gridV  files
   rapatrie ${CONFCASE}_y${YEAR}_ANNUAL_gridV.nc $MEANY ${CONFCASE}_y${YEAR}_ANNUAL_gridV.nc
 
   # get mesh mask files + new_maskglo
   rapatrie  ${MESH_MASK_ID}_byte_mask.nc $IDIR mask.nc
   rapatrie  ${MESH_MASK_ID}_mesh_hgr.nc $IDIR mesh_hgr.nc
   rapatrie  ${MESH_MASK_ID}_mesh_zgr.nc $IDIR mesh_zgr.nc
   rapatrie  new_maskglo.nc $IDIR new_maskglo.nc
 
   cdfmoc ${CONFCASE}_y${YEAR}_ANNUAL_gridV.nc
 
   # dispose on gaya MEAN/YEAR directory
   expatrie moc.nc $MEANY ${CONFCASE}_y${YEAR}_MOC.nc
  fi

# Mixed Layer Diagnostics : Input file : gridT for month 03 and 09 mesh_hgr, mesh_zgr
#^^^^^^^^^^^^^^^^^^^^^^^^
  if [ $MXL == 1 ] ; then
   # get mesh mask files
   rapatrie  ${MESH_MASK_ID}_mesh_hgr.nc $IDIR mesh_hgr.nc
   rapatrie  ${MESH_MASK_ID}_mesh_zgr.nc $IDIR mesh_zgr.nc
 
   for m in 3 9  ; do
     f=${CONFCASE}_y${YEAR}m0${m}_gridT.nc 
     g=$(echo $f | sed -e 's/gridT/MXL/')

     rapatrie $f $MEANY $f
 
     cdfmxl  $f
 
     # dispose on gaya, MEAN/YEAR directory
     expatrie mxl.nc $MEANY $g
   done
  fi

#=============================================================================
#  PART II: Time series: compute some integral quantities relevant for monitor
#           the ocean variability, and the behaviour of the on going run. 
#           Output is basically a small ASCII file, from which a matlab
#           suitable input file  (.mtl) is derived.
#=============================================================================
# Global MEANS: T S SSH Input files: gridT , mesh_hgr, mesh_zgr, mask
#^^^^^^^^^^^^^^
  if [ $TSMEAN == 1 ] ; then
   # get mesh mask files
   rapatrie  ${MESH_MASK_ID}_byte_mask.nc $IDIR mask.nc
   rapatrie  ${MESH_MASK_ID}_mesh_hgr.nc $IDIR mesh_hgr.nc
   rapatrie  ${MESH_MASK_ID}_mesh_zgr.nc $IDIR mesh_zgr.nc
  
   # get gridT files
   rapatrie ${CONFCASE}_y${YEAR}_ANNUAL_gridT.nc $MEANY ${CONFCASE}_y${YEAR}_ANNUAL_gridT.nc

   # set header on the output file (ASCII)
   fsshmean=${CONFCASE}_y${YEAR}_SSHMEAN.txt
   ftmean=${CONFCASE}_y${YEAR}_TMEAN.txt
   fsmean=${CONFCASE}_y${YEAR}_SMEAN.txt
   echo $YEAR >  $fsshmean ; echo $YEAR >  $ftmean ;  echo $YEAR >  $fsmean

   # 3D means
   cdfmean  ${CONFCASE}_y${YEAR}_ANNUAL_gridT.nc sossheig T >> $fsshmean
   cdfmean  ${CONFCASE}_y${YEAR}_ANNUAL_gridT.nc votemper T >> $ftmean
   cdfmean  ${CONFCASE}_y${YEAR}_ANNUAL_gridT.nc vosaline T >> $fsmean
 
   # dispose ASCII file in the -DIAGS directory
   expatrie  $fsshmean  $DIAGS $fsshmean
   expatrie  $ftmean  $DIAGS $ftmean
   expatrie  $fsmean  $DIAGS $fsmean
   
####  Append corresponding line to the matlab files. If it does not exists, create it with header
   # TMEAN
   if [ $(chkfile $MONITOR/${CONFCASE}_TMEAN.mtl ) == present ] ; then
     rapatrie ${CONFCASE}_TMEAN.mtl $MONITOR ${CONFCASE}_TMEAN.mtl
   else
    # first time : header
    dep=$( cat $ftmean | grep -e 'Mean value at level' | awk '{ printf "%10.1f" ,  $7 }' )
    echo 0000 0000 0000 $dep > ${CONFCASE}_TMEAN.mtl
   fi
   year=$( head -1 $ftmean )
   mean=$( cat $ftmean | grep -e 'over' | awk '{ printf "%10.4f" ,  $6 }' )
   tem=$( cat $ftmean | grep  -e 'Mean value at level' | awk '{ printf "%10.4f" ,  $9 }' )
   sshmean=$( cat $fsshmean | grep ocean | awk '{ printf " %8.4f " , $6 }' )

   echo $year $sshmean $mean $tem >>  ${CONFCASE}_TMEAN.mtl
   expatrie ${CONFCASE}_TMEAN.mtl  $MONITOR ${CONFCASE}_TMEAN.mtl
 
   # SMEAN
   if [ $(chkfile $MONITOR/${CONFCASE}_SMEAN.mtl ) == present ] ; then
     rapatrie ${CONFCASE}_SMEAN.mtl $MONITOR ${CONFCASE}_SMEAN.mtl
   else
    # first time
    dep=$( cat $fsmean | grep -e 'Mean value at level' | awk '{ printf "%10.1f" ,  $7 }' )
    echo 0000 0000 0000 $dep > ${CONFCASE}_SMEAN.mtl
   fi
   year=$( head -1 $fsmean )
   mean=$( cat $fsmean | grep -e 'over' | awk '{ printf "%10.4f" ,  $6 }' )
   sal=$( cat $fsmean | grep  -e 'Mean value at level' | awk '{ printf "%10.4f" ,  $9 }' )
 
   echo $year $sshmean $mean $sal >>  ${CONFCASE}_SMEAN.mtl
   expatrie ${CONFCASE}_SMEAN.mtl  $MONITOR ${CONFCASE}_SMEAN.mtl
 
#### cp to web site
   cptoweb ${CONFCASE}_TMEAN.mtl
   cptoweb ${CONFCASE}_SMEAN.mtl
 
   # clean up a little bit
   \rm ${CONFCASE}_[TS]MEAN.mtl
  fi


# Ice Volume area and extent for m02 m03   m08 m09: input file : icemod, and mesh_mask
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  if [ $ICE == 1 ] ; then
   # get icemod file for the month 02 03 and 08  09
   rapatrie  ${CONFCASE}_y${YEAR}m02_icemod.nc $MEANY ${CONFCASE}_y${YEAR}m02_icemod.nc
   rapatrie  ${CONFCASE}_y${YEAR}m03_icemod.nc $MEANY ${CONFCASE}_y${YEAR}m03_icemod.nc
   rapatrie  ${CONFCASE}_y${YEAR}m08_icemod.nc $MEANY ${CONFCASE}_y${YEAR}m08_icemod.nc
   rapatrie  ${CONFCASE}_y${YEAR}m09_icemod.nc $MEANY ${CONFCASE}_y${YEAR}m09_icemod.nc
 
   # get mesh mask files
   rapatrie  ${MESH_MASK_ID}_byte_mask.nc $IDIR mask.nc
   rapatrie  ${MESH_MASK_ID}_mesh_hgr.nc $IDIR mesh_hgr.nc
   rapatrie  ${MESH_MASK_ID}_mesh_zgr.nc $IDIR mesh_zgr.nc
 
   # Ascii output file:
   fice=${CONFCASE}_y${YEAR}_ice.txt
 
   echo '###' $YEAR 02 > $fice
   cdficediags ${CONFCASE}_y${YEAR}m02_icemod.nc  >> $fice
   echo '###' $YEAR 03 >> $fice
   cdficediags ${CONFCASE}_y${YEAR}m03_icemod.nc  >> $fice
   echo '###' $YEAR 08 >> $fice
   cdficediags ${CONFCASE}_y${YEAR}m08_icemod.nc  >> $fice
   echo '###' $YEAR 09 >> $fice
   cdficediags ${CONFCASE}_y${YEAR}m09_icemod.nc  >> $fice
 
   expatrie $fice $DIAGS $fice 
 
#### Append corresponding lines to matlab file for time series
   #ice
   month='02 03 08 09 '
   if [ $(chkfile $MONITOR/${CONFCASE}_ice.mtl ) == present ] ; then
     rapatrie ${CONFCASE}_ice.mtl $MONITOR ${CONFCASE}_ice.mtl
   else
    # first time: create file and add header
    echo 0000 $month  $month $month $month $month $month > ${CONFCASE}_ice.mtl
   fi
 
   year=$( head -1 $fice | awk '{ print $2}' )
   nvol=$(  cat $fice | grep -e 'NVolume' | grep -v NVolumet | awk '{ printf "%.0f  ", $4}' )
   svol=$(  cat $fice | grep -e 'SVolume' | grep -v SVolumet | awk '{ printf "%.0f  ", $4}' )
   narea=$(  cat $fice | grep -e 'NArea' | awk '{ printf "%.0f  ", $4}' )
   sarea=$(  cat $fice | grep -e 'SArea' | awk '{ printf "%.0f  ", $4}' )
   nextent=$(  cat $fice | grep -e 'NExtend' | awk '{ printf "%.0f  ", $4}' )
   sextent=$(  cat $fice | grep -e 'SExtend' | awk '{ printf "%.0f  ", $4}' )
 
   echo $year $nvol $svol $narea $sarea $nextent $sextent >> ${CONFCASE}_ice.mtl
 
   expatrie ${CONFCASE}_ice.mtl  $MONITOR ${CONFCASE}_ice.mtl
 
#### cp to web site
   cptoweb ${CONFCASE}_ice.mtl
 
   # clean up a little bit
   \rm ${CONFCASE}_ice.mtl
  fi

# Ice Volume area and extent for all months: input file : icemod, and mesh_mask
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  if [ $ICEMONTH == 1 ] ; then
   # get icemod files
   m=1
   while (( $m <= 12 )) ; do
    mm=$( printf "%02d" $m )
    rapatrie  ${CONFCASE}_y${YEAR}m${mm}_icemod.nc $MEANY ${CONFCASE}_y${YEAR}m${mm}_icemod.nc
    m=$(( m + 1 ))
   done
 
   # get mesh mask files
   rapatrie  ${MESH_MASK_ID}_byte_mask.nc $IDIR mask.nc
   rapatrie  ${MESH_MASK_ID}_mesh_hgr.nc $IDIR mesh_hgr.nc
   rapatrie  ${MESH_MASK_ID}_mesh_zgr.nc $IDIR mesh_zgr.nc
 
   # Ascii output file:
   fice=${CONFCASE}_y${YEAR}_icemonth.txt
 
   m=1
   while (( $m <= 12 )) ; do
    mm=$( printf "%02d" $m )
    case $mm in 
    01) echo '###' $YEAR $mm > $fice ;;
    *)  echo '###' $YEAR $mm >> $fice ;;
    esac
    cdficediags ${CONFCASE}_y${YEAR}m${mm}_icemod.nc  >> $fice
    m=$(( m + 1 ))
   done
 
   expatrie $fice $DIAGS $fice 
 
#### Append corresponding lines to matlab file for time series
   #ice
   month='01 02 03 04 05 06 07 08 09 10 11'
   if [ $(chkfile $MONITOR/${CONFCASE}_icemonth.mtl ) == present ] ; then
     rapatrie ${CONFCASE}_icemonth.mtl $MONITOR ${CONFCASE}_icemonth.mtl
   else
    # first time: create file and add header
    echo 0000 $month  $month $month $month $month $month > ${CONFCASE}_icemonth.mtl
   fi
 
   year=$( head -1 $fice | awk '{ print $2}' )
   nvol=$(  cat $fice | grep -e 'NVolume' | grep -v NVolumet | awk '{ printf "%.0f  ", $4}' )
   svol=$(  cat $fice | grep -e 'SVolume' | grep -v SVolumet | awk '{ printf "%.0f  ", $4}' )
   narea=$(  cat $fice | grep -e 'NArea' | awk '{ printf "%.0f  ", $4}' )
   sarea=$(  cat $fice | grep -e 'SArea' | awk '{ printf "%.0f  ", $4}' )
   nextent=$(  cat $fice | grep -e 'NExtend' | awk '{ printf "%.0f  ", $4}' )
   sextent=$(  cat $fice | grep -e 'SExtend' | awk '{ printf "%.0f  ", $4}' )
 
   echo $year $nvol $svol $narea $sarea $nextent $sextent >> ${CONFCASE}_icemonth.mtl
 
   expatrie ${CONFCASE}_icemonth.mtl  $MONITOR ${CONFCASE}_icemonth.mtl
 
#### cp to web site
   cptoweb ${CONFCASE}_icemonth.mtl
 
   # clean up a little bit
   \rm ${CONFCASE}_icemonth.mtl
  fi

# Vertical T-S profiles off the coast of Portugal for Gib monitoring: input file: gridT, mesh_mask
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  if [ $GIB == 1 ] ; then
   # get gridT file
   rapatrie ${CONFCASE}_y${YEAR}_ANNUAL_gridT.nc $MEANY ${CONFCASE}_y${YEAR}_ANNUAL_gridT.nc
 
   # get mesh mask files
   rapatrie  ${MESH_MASK_ID}_byte_mask.nc $IDIR mask.nc
   rapatrie  ${MESH_MASK_ID}_mesh_hgr.nc $IDIR mesh_hgr.nc
   rapatrie  ${MESH_MASK_ID}_mesh_zgr.nc $IDIR mesh_zgr.nc
 
   # Ascii output files:
   ftgib=${CONFCASE}_y${YEAR}_TGIB.txt
   fsgib=${CONFCASE}_y${YEAR}_SGIB.txt
 
   echo $YEAR > $ftgib
   cdfmean ${CONFCASE}_y${YEAR}_ANNUAL_gridT.nc  votemper T $GIBWIN  0 0  >>  $ftgib
   echo $YEAR > $fsgib
   cdfmean ${CONFCASE}_y${YEAR}_ANNUAL_gridT.nc  vosaline T $GIBWIN  0 0  >>  $fsgib
 
   expatrie $ftgib $DIAGS $ftgib
   expatrie $fsgib $DIAGS $fsgib
 
#### Append corresponding lines to matlab file for time series
   #GIB
   if [ $(chkfile $MONITOR/${CONFCASE}_gib.mtl ) == present ] ; then
     rapatrie ${CONFCASE}_gib.mtl $MONITOR ${CONFCASE}_gib.mtl
   else
    # first time : Create header with Levitus equivalent
    # requires  LEVITUS 'same' diags (from the ANNUAL mean )
    levitus=Levitus_p2.1_ANNUAL_TS_masked_$( echo $CONFIG | tr 'A-Z' 'a-z').nc
    rapatrie $levitus $IDIR $levitus
    cdfmean $levitus  votemper T $GIBWIN  0 0  >  LEVITUS_y0000_TGIB.txt
    cdfmean $levitus  vosaline T $GIBWIN  0 0  >  LEVITUS_y0000_SGIB.txt
    expatrie  LEVITUS_y0000_TGIB.txt $DIAGS  LEVITUS_y0000_TGIB.txt
    expatrie  LEVITUS_y0000_SGIB.txt $DIAGS  LEVITUS_y0000_SGIB.txt

    dep=$( cat  LEVITUS_y0000_TGIB.txt | grep 'Mean value at level' | awk '{ printf "%8.1f",$7 }' )
    Tlev=$( cat LEVITUS_y0000_TGIB.txt | grep 'Mean value at level' | awk '{ printf "%8.4f", $9 }' )
    Slev=$( cat LEVITUS_y0000_SGIB.txt | grep 'Mean value at level' | awk '{ printf "%8.4f", $9 }' )
 
    echo 0000 $dep > ${CONFCASE}_gib.mtl
    echo 0000 $Tlev >> ${CONFCASE}_gib.mtl
    echo 0000 $Slev >> ${CONFCASE}_gib.mtl
   fi
 
    year=$( head -1 $ftgib )
    Tcur=$( cat $ftgib | grep 'Mean value at level' | awk '{ printf "%8.4f", $9 }' )
    Scur=$( cat $fsgib | grep 'Mean value at level' | awk '{ printf "%8.4f", $9 }' )
    echo $year $Tcur >> ${CONFCASE}_gib.mtl
    echo $year $Scur >> ${CONFCASE}_gib.mtl
 
    expatrie  ${CONFCASE}_gib.mtl $MONITOR  ${CONFCASE}_gib.mtl
 
#### cp to web site
   cptoweb ${CONFCASE}_gib.mtl
 
   # clean up a little bit
   \rm ${CONFCASE}_gib.mtl
  fi

# El nino indexes : Input files : monthly gridT,  mesh mask
#^^^^^^^^^^^^^^^^^^
  if [ $ELNINO == 1 ] ; then
   # get mesh mask files
   rapatrie  ${MESH_MASK_ID}_byte_mask.nc $IDIR mask.nc
   rapatrie  ${MESH_MASK_ID}_mesh_hgr.nc $IDIR mesh_hgr.nc
   rapatrie  ${MESH_MASK_ID}_mesh_zgr.nc $IDIR mesh_zgr.nc
 
   # Ascii outputfile
   fnino=${CONFCASE}_y${YEAR}_NINO.txt
 
   # get monthly mean gridT files and compute mean SST on each NINO box
   for  m in  1 2 3 4 5 6 7 8 9 10 11 12 ; do
     mm=$(printf "%02d" $m)
     f=${CONFCASE}_y${YEAR}m${mm}_gridT.nc 

     rapatrie $f $MEANY  $f
 
     #  header
     printf "%04d %02d" $YEAR $m >>   $fnino
 
    # nino 1+2   [ -90 W -- -80 W, -10 S -- 10 N ]
    cdfmean  $f votemper T $NINO12 1 1 | tail -1 | awk '{ printf " %8.5f 0.00", $6 }'  >> $fnino
    # nino 3     [ -150 W -- -90 W, -5 S -- 5 N ]
    cdfmean  $f votemper T $NINO3 1 1  | tail -1 | awk '{ printf " %8.5f 0.00", $6 }'  >> $fnino
    # nino 4     [ -200 W -- -150 W, -5 S -- 5 N ]
    cdfmean  $f votemper T $NINO4 1 1 | tail -1 | awk '{ printf " %8.5f 0.00", $6 }'  >> $fnino
    # nino 3.4   [ -170 W -- -120 W, -% S -- % N ]
    cdfmean  $f votemper T $NINO34 1 1 | tail -1 | awk '{ printf " %8.5f 0.00\n", $6 }'  >> $fnino
 
    \rm $f
   done
 
   expatrie $fnino $DIAGS $fnino
 
#### Append corresponding lines to matlab file for time series
   #NINO
   if [ $(chkfile $MONITOR/${CONFCASE}_nino.mtl ) == present ] ; then
     rapatrie ${CONFCASE}_nino.mtl $MONITOR ${CONFCASE}_nino.mtl
   fi
 
   cat $fnino >> ${CONFCASE}_nino.mtl
   expatrie  ${CONFCASE}_nino.mtl $MONITOR  ${CONFCASE}_nino.mtl
 
#### cp to web site
   cptoweb ${CONFCASE}_nino.mtl
 
   # clean up a little bit
   \rm ${CONFCASE}_nino.mtl
  fi

# Transport: Input files: VT, gridU, gridV, mesh mask, section.dat
#^^^^^^^^^^^
  if [ $TRP == 1 ] ; then
   # section.dat describes the position (I,J) of the sections to monitor
   cp $P_CTL/section.dat .
 
   # get VT , gridU, gridV files
   rapatrie ${CONFCASE}_y${YEAR}_ANNUAL_VT.nc $MEANY ${CONFCASE}_y${YEAR}_ANNUAL_VT.nc
   rapatrie ${CONFCASE}_y${YEAR}_ANNUAL_gridU.nc $MEANY ${CONFCASE}_y${YEAR}_ANNUAL_gridU.nc
   rapatrie ${CONFCASE}_y${YEAR}_ANNUAL_gridV.nc $MEANY ${CONFCASE}_y${YEAR}_ANNUAL_gridV.nc
 
   # get mesh mask files
   rapatrie  ${MESH_MASK_ID}_byte_mask.nc $IDIR mask.nc
   rapatrie  ${MESH_MASK_ID}_mesh_hgr.nc $IDIR mesh_hgr.nc
   rapatrie  ${MESH_MASK_ID}_mesh_zgr.nc $IDIR mesh_zgr.nc
 
   # Ascii output file:
   fsection=${CONFCASE}_y${YEAR}_section_monitor.txt
 
   echo $YEAR > $fsection
 
   cdftransportiz ${CONFCASE}_y${YEAR}_ANNUAL_VT.nc \
                  ${CONFCASE}_y${YEAR}_ANNUAL_gridU.nc \
                  ${CONFCASE}_y${YEAR}_ANNUAL_gridV.nc  < section.dat >> $fsection
 
   # eliminate garbage from txt file ...
   grep -v Give $fsection | grep -v level | grep -v IMAX | grep -v FROM > tmp
   mv -f tmp $fsection
   
   expatrie  $fsection $DIAGS $fsection
 
#### Append corresponding lines to matlab file for time series
   #sections
   if [ $(chkfile $MONITOR/${CONFCASE}_matrix.mtl ) == present ] ; then
     rapatrie ${CONFCASE}_matrix.mtl $MONITOR ${CONFCASE}_matrix.mtl
   fi
 
   year=$( head -1 $fsection )
   mass=$( cat $fsection | grep Mass | awk '{ printf "%8.3f" ,  $4 }' )
   heat=$( cat $fsection | grep Heat | awk '{ printf "%8.3f" ,  $4 }' )
   salt=$( cat $fsection | grep Salt | awk '{ printf "%8.1f" ,  $4 *1. }' )
 
   echo $year $mass $heat $salt  >>  ${CONFCASE}_matrix.mtl
   # JMM remark: previous monitoring added sshmean, tmean smean at the end of the
   #     line but it is not used in the matlab script so I skip it 
 
   expatrie  ${CONFCASE}_matrix.mtl $MONITOR  ${CONFCASE}_matrix.mtl
 
#### cp to web site
   cptoweb ${CONFCASE}_matrix.mtl

   # clean up a little bit
   \rm ${CONFCASE}_matrix.mtl
  fi
 
# Heat and Salt Meridional Transport : Input files : VT, mesh mask, new_maskglo
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  if [ $MHT == 1 ] ; then
# (a) From advection:
#--------------------
   # get VT  files
   rapatrie ${CONFCASE}_y${YEAR}_ANNUAL_VT.nc $MEANY ${CONFCASE}_y${YEAR}_ANNUAL_VT.nc
 
   # get mesh mask files + new_maskglo
   rapatrie  ${MESH_MASK_ID}_byte_mask.nc $IDIR mask.nc
   rapatrie  ${MESH_MASK_ID}_mesh_hgr.nc $IDIR mesh_hgr.nc
   rapatrie  ${MESH_MASK_ID}_mesh_zgr.nc $IDIR mesh_zgr.nc
   rapatrie  new_maskglo.nc $IDIR new_maskglo.nc
 
   # Ascii output file:
   fheat=${CONFCASE}_y${YEAR}_heattrp.dat
   fsalt=${CONFCASE}_y${YEAR}_salttrp.dat
 
   cdfmhst  ${CONFCASE}_y${YEAR}_ANNUAL_VT.nc
 
   expatrie zonal_heat_trp.dat $DIAGS ${CONFCASE}_y${YEAR}_heattrp.dat
   expatrie zonal_salt_trp.dat $DIAGS ${CONFCASE}_y${YEAR}_salttrp.dat
 
   # needed below with the correct name
   cp zonal_heat_trp.dat ${CONFCASE}_y${YEAR}_heattrp.dat
 
# (b) from Surface Heat fluxes
#-----------------------------
    cdfhflx  ${CONFCASE}_y${YEAR}_ANNUAL_gridT.nc
 
    expatrie hflx.out $DIAGS ${CONFCASE}_y${YEAR}_hflx.dat
 
#### Append corresponding lines to matlab file for time series
   # Heat/salt transport
   if [ $(chkfile $MONITOR/${CONFCASE}_heat.mtl ) == present ] ; then
     rapatrie ${CONFCASE}_heat.mtl $MONITOR ${CONFCASE}_heat.mtl
   else 
   # first time 
   # output latitude (North to South) on the first row
     printf "%04d " 0000 > ${CONFCASE}_heat.mtl
     cat zonal_heat_trp.dat  |\
     awk 'BEGIN{s=0} \
        { if (NR > 2) {printf " %8.3f ", $2} \
        }\
        END{ printf "\n" }'   >> ${CONFCASE}_heat.mtl
   fi
 
   # Global Ocean
   printf "%04d "  $YEAR >> ${CONFCASE}_heat.mtl
   cat zonal_heat_trp.dat  |\
   awk 'BEGIN{s=0} \
        { if (NR > 2) {printf " %8.3f ", $3} \
        }\
        END{ printf "\n" }'   >> ${CONFCASE}_heat.mtl
 
   printf "%04d "  $YEAR >> ${CONFCASE}_heat.mtl
   cat hflx.out  |\
   awk 'BEGIN{s=0} \
        { if (NR > 2) {printf " %8.3f ", $3} \
        }\
        END{ printf "\n" }'   >> ${CONFCASE}_heat.mtl
 
 # Atlantic Ocean
   printf "%04d "  $YEAR >> ${CONFCASE}_heat.mtl
   cat zonal_heat_trp.dat  |\
   awk 'BEGIN{s=0} \
        { if (NR > 2) {printf " %8.3f ", $4} \
        }\
        END{ printf "\n" }'   >> ${CONFCASE}_heat.mtl
 
   printf "%04d "  $YEAR >> ${CONFCASE}_heat.mtl
   cat hflx.out  |\
   awk 'BEGIN{s=0} \
        { if (NR > 2) {printf " %8.3f ", $4} \
        }\
        END{ printf "\n" }'   >> ${CONFCASE}_heat.mtl
 
 # Indo-Pacific Ocean
   printf "%04d "  $YEAR >> ${CONFCASE}_heat.mtl
   cat zonal_heat_trp.dat  |\
   awk 'BEGIN{s=0} \
        { if (NR > 2) {printf " %8.3f ", $5 + $6 } \
        }\
        END{ printf "\n" }'   >> ${CONFCASE}_heat.mtl
 
   printf "%04d "  $YEAR >> ${CONFCASE}_heat.mtl
   cat hflx.out  |\
   awk 'BEGIN{s=0} \
        { if (NR > 2) {printf " %8.3f ", $5 } \
        }\
        END{ printf "\n" }'   >> ${CONFCASE}_heat.mtl
 
   expatrie  ${CONFCASE}_heat.mtl $MONITOR  ${CONFCASE}_heat.mtl
 
#### cp to website
   cptoweb ${CONFCASE}_heat.mtl
 
   # clean
   \rm  ${CONFCASE}_heat.mtl
  fi


# MAX and MIN of MOC: requires that MOC files already exists
#^^^^^^^^^^^^^^^^^^^^
  if [ $MAXMOC == 1  ] ; then
   f=moc.nc
   rapatrie ${CONFCASE}_y${YEAR}_MOC.nc $DIAGS $f

   # Ascii output file
   fmaxmoc=${CONFCASE}_y${YEAR}_minmaxmoc.txt
   echo $YEAR > $fmaxmoc
   # GLO
   printf "%s" 'Glo ' >>  $fmaxmoc ; cdfmaxmoc $f glo 20 60 500 2000 | grep Maximum >> $fmaxmoc
   printf "%s" 'Glo ' >>  $fmaxmoc ; cdfmaxmoc $f glo -40 30 2000 5500 | grep Minimum >> $fmaxmoc
   # ATL
   printf "%s" 'Atl ' >>  $fmaxmoc ; cdfmaxmoc $f atl 0 60 500 2000 | grep Maximum >> $fmaxmoc
   printf "%s" 'Atl ' >>  $fmaxmoc ; cdfmaxmoc $f atl -20 40 2000 5500 | grep Minimum  >> $fmaxmoc
   #INP
   printf "%s" 'Inp ' >>  $fmaxmoc ; cdfmaxmoc $f inp 15 50 100 1000 | grep Minimum >> $fmaxmoc
   printf "%s" 'Inp ' >>  $fmaxmoc ; cdfmaxmoc $f inp -30 20 1000 5500  | grep Minimum >> $fmaxmoc
   #AUS
   printf "%s" 'Aus ' >>  $fmaxmoc ; cdfmaxmoc $f glo -70 0 0 2000   | grep Maximum >> $fmaxmoc
   printf "%s" 'Aus ' >>  $fmaxmoc ; cdfmaxmoc $f glo -70 0 2000 5500  | grep Minimum >> $fmaxmoc
   
   expatrie $fmaxmoc $DIAGS $fmaxmoc

   # Max and Min of MOC at some specific latitudes
   fmaxmoc40=${CONFIG}-${CASE}_y${YEAR}_maxmoc40.txt

   echo $YEAR > $fmaxmoc40
   # GLO  MAX at 40 N and 30S
   printf "%s" 'Glo ' >>  $fmaxmoc40 ; cdfmaxmoc $f glo 40 40 500 2000 | grep Maximum >> $fmaxmoc40
   printf "%s" 'Glo ' >>  $fmaxmoc40 ; cdfmaxmoc $f glo -30 -30 500  5500 | grep Maximum >> $fmaxmoc40
   # ATL  MAX at 40N and 30S
   printf "%s" 'Atl ' >>  $fmaxmoc40 ; cdfmaxmoc $f atl 40 40 500 2000 | grep Maximum >> $fmaxmoc40
   printf "%s" 'Atl ' >>  $fmaxmoc40 ; cdfmaxmoc $f atl -30 -30  500 5000 | grep Maximum >> $fmaxmoc40
   #INP  Min at 30 S
   printf "%s" 'Inp ' >>  $fmaxmoc40 ; cdfmaxmoc $f inp -30 -30 1000 5500  | grep Minimum >> $fmaxmoc40
   #AUS  MAX at 50 S
   printf "%s" 'Aus ' >>  $fmaxmoc40 ; cdfmaxmoc $f glo -50 -50 0 2000   | grep Maximum >> $fmaxmoc40

   expatrie $fmaxmoc40 $DIAGS $fmaxmoc40

#### Append to matlab file
   # maxmoc
   if [ $(chkfile $MONITOR/${CONFCASE}_maxmoc.mtl ) == present ] ; then
     rapatrie ${CONFCASE}_maxmoc.mtl $MONITOR ${CONFCASE}_maxmoc.mtl
   fi
 
   year=$( head -1 $fmaxmoc )
   mht=${CONFCASE}_y${YEAR}_heattrp.dat
   # get mht file if necessary
   rapatrie $mht $DIAGS $mht
 
   # GLO
   maxglo=$( cat $fmaxmoc    | grep -e '^Glo' | grep Max | awk '{ printf "%8.3f" ,  $3 }' )
   minglo=$( cat $fmaxmoc    | grep -e '^Glo' | grep Min | awk '{ printf "%8.3f" ,  $3 }' )
 
   # ATL
   maxatl=$( cat $fmaxmoc    | grep -e '^Atl' | grep Max | awk '{ printf "%8.3f" ,  $3 }' )
   minatl=$( cat $fmaxmoc    | grep -e '^Atl' | grep Min | awk '{ printf "%8.3f" ,  $3 }' )
 
   # INP : attention we have 2 Minimum for INP mininp1 and mininp2
   tmp=$( cat $fmaxmoc    | grep -e '^Inp' | grep Min | awk '{ printf "%8.3f" ,  $3 }' )
   mininp1=$( echo $tmp | awk '{print $1}' )
   mininp2=$( echo $tmp | awk '{print $2}' )
 
   # AUS
   maxaus=$( cat $fmaxmoc    | grep -e '^Aus' | grep Max | awk '{ printf "%8.3f" ,  $3 }' )
   minaus=$( cat $fmaxmoc    | grep -e '^Aus' | grep Min | awk '{ printf "%8.3f" ,  $3 }' )
 
   # heattrp at 20 N
   heattrp=$(  cat $mht | awk '{  if ( $2 >= 20 ) { atlmht=$4 ; glomht=$3 } } END { printf "%6.3f %6.3f ", glomht, atlmht }')
   mhtglo=$( echo $heattrp | awk '{print $1}' )
   mhtatl=$( echo $heattrp | awk '{print $2}' )
 

   echo $year $maxglo $minglo $mhtglo $maxatl $minatl $mhtatl  $mininp1 $mininp2 0000  \
              $maxaus $minaus 0000 >>  ${CONFCASE}_maxmoc.mtl

   expatrie ${CONFCASE}_maxmoc.mtl $MONITOR ${CONFCASE}_maxmoc.mtl

#### send it to web site
   cptoweb ${CONFCASE}_maxmoc.mtl
 
   #maxmoc40
   if [ $(chkfile $MONITOR/${CONFCASE}_maxmoc40.mtl ) == present ] ; then
     rapatrie ${CONFCASE}_maxmoc40.mtl $MONITOR ${CONFCASE}_maxmoc40.mtl
   fi

   year=$( head -1 $fmaxmoc40 )
 
   # GLO  max a 40 et -30
   tmp=$( cat $fmaxmoc40    | grep -e '^Glo' | grep Max | awk '{ printf "%8.3f" ,  $3 }' )
   maxglo40n=$( echo $tmp | awk '{print $1}' )
   maxglo30s=$( echo $tmp | awk '{print $2}' )
 
   # ATL  max a 40 et -30
   tmp=$( cat $fmaxmoc40    | grep -e '^Atl' | grep Max | awk '{ printf "%8.3f" ,  $3 }' )
   maxatl40n=$( echo $tmp | awk '{print $1}' )
   maxatl30s=$( echo $tmp | awk '{print $2}' )
 
   # INP : 1 Min at -30 S
   mininp30s=$( cat $fmaxmoc40    | grep -e '^Inp' | grep Min | awk '{ printf "%8.3f" ,  $3 }' )
 
   # AUS 1 max at -50 S
   maxaus50s=$( cat $fmaxmoc40    | grep -e '^Aus' | grep Max | awk '{ printf "%8.3f" ,  $3 }' )

   echo $year $maxglo40n $maxglo30s  $maxatl40n $maxatl30s  $mininp30s $maxaus50s >>  ${CONFCASE}_maxmoc40.mtl

   expatrie ${CONFCASE}_maxmoc40.mtl $MONITOR ${CONFCASE}_maxmoc40.mtl

#### send it to web site
   cptoweb ${CONFCASE}_maxmoc40.mtl
  fi


# DCT :Density Class transport: Input files : gridT, gridU gridV, mesh mask, dens_section.dat
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  if [ $DCT == 1 ] ; then
  # dens_section.dat describe the sections (either zonal or meridional) where the DCT is computed
  cp $P_CTL/dens_section.dat .

  # get mesh mask files
  rapatrie  ${MESH_MASK_ID}_byte_mask.nc $IDIR mask.nc
  rapatrie  ${MESH_MASK_ID}_mesh_hgr.nc $IDIR mesh_hgr.nc
  rapatrie  ${MESH_MASK_ID}_mesh_zgr.nc $IDIR mesh_zgr.nc

  # Required post_processing script: DCT are computed on monthly means, then average is performed
  # for annual values. This process is still done through temporary bimg/dimg files (remnant of the
  # old Clipper times).  By the way, 2 bimgtools are required: bimgmoy4 and bimgcaltrans
  #  In-lining of this script may be confusing. I leave it as an external module.
  cp $CDFTOOLS/JOBS/trpsig_postproc.ksh ./

  # due to the large amount of files that are produced by this diags, we prefer to keep them
  # on a separate directory
  chkdirg ${CONFIG}/${CONFCASE}-TRPSIG/
  chkdirg ${CONFIG}/${CONFCASE}-TRPSIG/$YEAR/
  chkdirg $DIAGS/$YEAR/
  chkdirg $DIAGS/TRPSIG/

  # also need temporary directories in the actual tmpdir:
  chkdir ${CONFIG}
  chkdir ${CONFIG}/${CONFCASE}-TRPSIG
  chkdir ${CONFIG}/${CONFCASE}-TRPSIG/$YEAR/

  TRPSIGY=${CONFIG}/${CONFCASE}-TRPSIG/$YEAR/

  for  m in  1 2 3 4 5 6 7 8 9 10 11 12 ; do
    mm=$(printf "%02d" $m)
    tfich=${CONFCASE}_y${YEAR}m${mm}_gridT.nc 
    ufich=$(echo  $tfich | sed -e 's/gridT/gridU/' )
    vfich=$(echo  $tfich | sed -e 's/gridT/gridV/' )

    #get files on gaya
    rapatrie  $tfich $MEANY  $tfich
    rapatrie  $ufich $MEANY  $ufich
    rapatrie  $vfich $MEANY  $vfich
    
    #retrieve tag time from file name
    tag=$(echo $tfich | sed -e "s/${CONFCASE}_//" -e 's/_gridT.nc//')

    echo $tag > ${CONFCASE}_y${tag}_trpsig_monitor.lst

    cdfsigtrp $tfich $ufich $vfich 21 30 180 -bimg -print  >>  ${CONFCASE}_y${tag}_trpsig_monitor.lst

    # save the monthly log file on gaya for an (improbable) eventual post processing ...
    expatrie ${CONFCASE}_y${tag}_trpsig_monitor.lst $TRPSIGY ${CONFCASE}_y${tag}_trpsig_monitor.lst
    # and create a mirror on the local tmpdir
    mv ${CONFCASE}_y${tag}_trpsig_monitor.lst  $TRPSIGY

    # Idem : save temporary bimg files on gaya and create local mirror
    for  b in *bimg ; do
        mv  $b ${CONFCASE}_y${tag}_$b
        expatrie ${CONFCASE}_y${tag}_$b $TRPSIGY ${CONFCASE}_y${tag}_$b
        mv  ${CONFCASE}_y${tag}_$b  $TRPSIGY
    done
    
    # Idem: for txt files
    mv trpsig.txt ${CONFCASE}_y${tag}_trpsig.txt
    expatrie ${CONFCASE}_y${tag}_trpsig.txt $TRPSIGY ${CONFCASE}_y${tag}_trpsig.txt
    mv  ${CONFCASE}_y${tag}_trpsig.txt $TRPSIGY

    # erase useless files ( monthly averages ) Keep tfich which can be used for MXL
    \rm *.bimg  $ufich $vfich

  # end of month loop
  done

  # Launch post processing   ( by itself a complex script ...)
  # This script retrieve CONFIG name and CASE from the directory name where it runs...
  cd ${CONFIG}/${CONFCASE}-TRPSIG
  $TMPDIR/trpsig_postproc.ksh
  cd $TMPDIR

  # save results on gaya ( as many files as sections in dens_section.dat)
  for f in ${CONFCASE}_y*_trpsig.txt ; do
    expatrie $f $DIAGS/TRPSIG/ $f
  done

#### Append results to mtl file for denmark strait and faroes bank channel
   #trpsig

   # mini and maxi of sigma0
   mini=25.2 ; maxi=28.5   # selected range convenient for DS and FBC
   file=${CONFCASE}_y${YEAR}_01_Denmark_strait_trpsig.txt 
   fil2=$( echo $file | sed -e 's/01_Denmark_strait/02_Faoes_Bank_Channel/' )

   sig=$(   cat $file | grep -v -e '^#' | awk '{ if ( $1 > mini  && $1 < maxi ) {printf "%8.3f"  ,  $1 } }' mini=$mini  maxi=$maxi)
   trp01=$( cat $file | grep -v -e '^#' | awk '{ if ( $1 > mini  && $1 < maxi ) {printf "%13.4e" ,  $2 } }' mini=$mini  maxi=$maxi)
   trp02=$( cat $fil2 | grep -v -e '^#' | awk '{ if ( $1 > mini  && $1 < maxi ) {printf "%13.4e" ,  $2 } }' mini=$mini  maxi=$maxi)

   if [ $(chkfile $MONITOR/${CONFCASE}_TRPSIG.mtl ) == present ] ; then
     rapatrie ${CONFCASE}_TRPSIG.mtl $MONITOR ${CONFCASE}_TRPSIG.mtl
   else
    #first time : add header with sigma classes
    echo 000000  $sig > ${CONFCASE}_TRPSIG.mtl 
   fi
   echo $YEAR $trp01 >> ${CONFCASE}_TRPSIG.mtl
   echo $YEAR $trp02 >> ${CONFCASE}_TRPSIG.mtl

   expatrie ${CONFCASE}_TRPSIG.mtl $MONITOR ${CONFCASE}_TRPSIG.mtl

####  cp to web site
   cptoweb ${CONFCASE}_TRPSIG.mtl
  
   # return to tmpdir
   cd $TMPDIR
   # Erase the TRPSIG tree for this current year
   \rm -r ${CONFIG} \rm ${CONFCASE}_y*_trpsig.txt *.mtl

  fi

# TRACER DIAGS (31/12 of each year) : Input files : ptrcT, mesh mask
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  if [ $TRACER == 1 ] ; then
   # get mesh mask files
   rapatrie  ${MESH_MASK_ID}_byte_mask.nc $IDIR mask.nc
   rapatrie  ${MESH_MASK_ID}_mesh_hgr.nc $IDIR mesh_hgr.nc
   rapatrie  ${MESH_MASK_ID}_mesh_zgr.nc $IDIR mesh_zgr.nc
   rapatrie  new_maskglo.nc $IDIR new_maskglo.nc
 
   # get tracer file from gaya: note that this is from -S dir (5 day average ... to discuss ...)
   rapatrie ${CONFCASE}_y${YEAR}m12d31_ptrcT.nc $SDIRY ${CONFCASE}_y${YEAR}m12d31_ptrcT.nc
 
   # Ascii output file:
   ftrc=${CONFCASE}_y${YEAR}_TRCmean.dat
 
   # Number of mol in the ocean ...
   printf "%04d "  $YEAR   >  $ftrc
 
   # CFC11
   \rm -f tmp1
   cdfmean  ${CONFCASE}_y${YEAR}m12d31_ptrcT.nc  invcfc T > tmp1
   area=$(cat tmp1 |  grep -e 'Mean value at level' | awk ' {print $12}')
   mean=$(cat tmp1 |  grep -e 'Mean value over the ocean' | awk ' {print $6}')
   total=$(echo $mean $area |  awk '{print $1 * $2 }' )
   printf "%s "  $total  >> $ftrc
 
   # B-C14
   \rm -f tmp1
   cdfmean  ${CONFCASE}_y${YEAR}m12d31_ptrcT.nc  invc14 T > tmp1
   area=$(cat tmp1 |  grep -e 'Mean value at level' | awk ' {print $12}')
   mean=$(cat tmp1 |  grep -e 'Mean value over the ocean' | awk ' {print $6}')
   total=$(echo $mean $area |  awk '{print $1 * $2 }' )
   printf "%s \n"  $total  >> $ftrc
 
   expatrie $ftrc $DIAGS $ftrc
 
   # zonal integral of inventories
   cdfzonalsum  ${CONFCASE}_y${YEAR}m12d31_ptrcT.nc  T
 
   # zonal means
   cdfzonalmean  ${CONFCASE}_y${YEAR}m12d31_ptrcT.nc  T
 
   # ncks is required on the prod machine ... !! not standard !!
   # it is used to take only the interesting variables from the results
   ncks -F -d deptht,1,1 -v zocfc11_glo,zobc14_glo,nav_lon,nav_lat zonalmean.nc zonalsurf.nc
 
   # put in ascii format the 1D profiles
   cdfzonalout zonalmean.nc > zonalmean.dat
   cdfzonalout zonalsum.nc >  zonalsum.dat
   cdfzonalout zonalsurf.nc >  zonalsurf.dat
 
   expatrie zonalmean.nc $MEANY ${CONFCASE}_y${YEAR}_TRCzonalmean.nc
   expatrie zonalsum.nc $MEANY ${CONFCASE}_y${YEAR}_TRCzonalsum.nc
 
   expatrie zonalmean.dat $DIAGS ${CONFCASE}_y${YEAR}_TRCzonalmean.dat
   expatrie zonalsum.dat $DIAGS ${CONFCASE}_y${YEAR}_TRCzonalsum.dat
   expatrie zonalsurf.dat $DIAGS ${CONFCASE}_y${YEAR}_TRCzonalsurf.dat
   \rm zonalsurf.nc
 
#### append to matlab file
   # Not done already
 
  fi
