#!/bin/ksh
set -x

CONFIG=CCOONNFF           # set the name of the config
CASE=CCAASSEE                 # set the case of the config

# set the list of years you want to monitor 'at once'  
yinit=YYYY              # initial year 
yend=YYYE               # last year

YEARS=''
while (( $yinit <= $yend )) ; do
 YEARS="$YEARS $yinit "
 yinit=$(( yinit + 1 ))
done

MESH_MASK_ID='ORCA025-G70_noBS_noRS_noPG_noMCBO' # root part of the mesh-mask files
                         # (they must be in the -I directory ( $CONFIG/${CONFIG}-I)
                         #  Standard name is thus : ${MESH_MASK_ID}_byte_mask.nc
                         #                          ${MESH_MASK_ID}_mesh_hgr.nc
                         #                          ${MESH_MASK_ID}_mesh_zgr.nc
#

# define the I-J window for GIB diags and El NINO DIAG
if [ $CONFIG = 'NATL025' ] ; then
  GIBWIN='338 353 239 260'
  # NOT RELEVANT FOR NATL025. Here for compatibility
  NINO12='790 830 459 499'
  NINO3='550 790 479 519 '
  NINO4='350 550 479 519 '
  NINO34='470 670 479 519 '
elif [ $CONFIG = 'ORCA025' ] ; then
  GIBWIN='1094 1109 653 674 '
  NINO12='790 830 459 499'
  NINO3='550 790 479 519 '
  NINO4='350 550 479 519 '
  NINO34='470 670 479 519 '
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
. ./config_def.ksh

cp config_def.ksh $TMPDIR
cp function_def.ksh $TMPDIR

CONFCASE=${CONFIG}-${CASE}
cd $TMPDIR
. ./config_def.ksh

for   YEAR in  $YEARS ; do
   .  $CDFTOOLS/JOBS/monitor_prod.ksh
   # clean the TMPDIR (on WORDIR in fact) from all files for YEAR
   cd $TMPDIR
   find . -name "*${YEAR}*" -exec \rm -rf {} \;
done
