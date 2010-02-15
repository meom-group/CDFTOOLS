#!/bin/ksh
# @ cpu_limit  = 4200
# @ data_limit = 1gb
# Nom du travail LoadLeveler
# @ job_name   = cdfsstconv
# Fichier de sortie standard du travail       
# @ output     = $(job_name).$(jobid)
# Fichier de sortie d'erreur du travail
# @ error      =  $(job_name).$(jobid)
# @ queue                   

set -x
#
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

# chkfile : Usage: chkfile gaya_file
#    check if a file exists on gaya, return present or absent.
chkfile() { rsh gaya " if [ -f $1 ] ; then echo present ;
                       else echo absent ; fi " ; }

# chkdirg  : Usage: chkdirg gaya_directory
#    check the existence of a directory on gaya. Create it if not present
chkdirg() { rsh gaya " if [ ! -d $1 ] ; then mkdir $1 ; fi " ; }

# chkdirw  : Usage: chkdirw web_site_directory
#   check the existence of a dir. on the web site. Create it if not present
chkdirw() { rsh meolipc.hmg.inpg.fr -l apache " if [ ! -d web/DRAKKAR/$1 ] ; 
            then mkdir web/DRAKKAR/$1 ; fi " ; }

# chkdir  : Usage: chkdir local_dir 
#   check the existence of a directory. Create it if not present
chkdir() { if [ ! -d $1 ] ; then mkdir $1 ; fi  ; }


CDFTOOLS=CDFTOOLS-2.1
CONFIG=ATL3
DIRCOO=/cache2/rost011/CLIPPER/GRID
COORD=coordinates.${CONFIG}
IDIR=/cache3/rost005/rcli007/${CONFIG}-I
IDIRNC=${CONFIG}/${CONFIG}-I

chkdirg $CONFIG
chkdirg $IDIRNC

cd $TMPDIR

cp ~/$CDFTOOLS/cdfsstconv ./
# coordinates.diags
rapatrie $COORD $DIRCOO coordinates.diags

year=1992
year2=1995
while (( $year <= $year2 )) ; do
  emp=ECMWF_emp_1d_${year}.${CONFIG}.nc
  sst0=REYNOLDS_sst_1d_${year}.${CONFIG}.nc
  if [ $(chkfile $IDIRNC/$sst0 ) == absent ] ; then
  # get fluxes and STRESS monthly files
  m=1
# while (( $m <= 12 )) ; do
#   mm=$( printf "%02d" $m )
#   flx=ECMWF.Y${year}.M${mm}.FLUX.${CONFIG}.dimg
#   str=ECMWF.Y${year}.M${mm}.STRESS.${CONFIG}.dimg
#   rapatrie $flx $IDIR $flx
#   rapatrie $str $IDIR $str
#   m=$(( m + 1 ))
# done
  # get SST for year -1 year and year+1
    ym1=$(( year - 1 ))
    yp1=$(( year + 1 ))
    y=$ym1
    while (( $y <= $yp1 )) ; do
      sst=REYNOLDS.Y${y}.SST.${CONFIG}.dimg
      if (( $y <= 2000 )) ; then
        rapatrie $sst $IDIR $sst
      fi
      y=$(( y + 1 ))
    done
    ./cdfsstconv $year $CONFIG
    for f in *.nc ; do
       expatrie  $f $IDIRNC $f
       # clean unnecessary files from tmpdir
       \rm $f
    done
       \rm -f  ECMWF*.dimg REYNOLDS.Y${ym1}.SST.${CONFIG}.dimg
    fi
   year=$(( year + 1 ))
done
  


