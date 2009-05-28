#!/bin/ksh
# function_def.ksh file for zahir and gaya. To be used with cdfmoy and cdfvT jobs


#  $Rev: 207 $
#  $Date: 2008-11-21 15:26:12 +0100 (Fri, 21 Nov 2008) $
#  $Id: function_def_zahir.ksh 207 2008-11-21 14:26:12Z rcli002 $

# FROM CDFMOY and CDFVT suite
#############################
# TEST network between service and scratch in batch
chknet() { NET=KO ; ping -c 1 $1 | grep "1 packets transmitted, 1 received, 0% packet loss" && NET=OK ; }
# chkdirg  path : check existence of directory path  on (remote) archiving machine. If it does not exist, create it.
chkdirg() { ssh $REMOTE_USER@${login_node} " if [ ! -d $1 ] ; then mkdir $1 ; fi " ; }
# getmonth  mm type : retrieve all 5 days average files for month mm and grid type 'type', corresponding to current year, current confcase.A
#                     ex: getmonth 04 gridU  : retrieve all april files for gridU
##########################################################################
getmonth() { for f  in $( ssh $REMOTE_USER@${login_node} " ls $SDIR/$SDIRY/${CONFCASE}_y${YEAR}m${1}d*_$2.nc " )  ; do
              file=`basename $f` ; if [ ! -f $file ] ; then if [ ! -f $WORKDIR/${CONFCASE}-S/$YEAR/$file ] ; then scp $REMOTE_USER@${login_node}:$f . ; else if [ -f $WORKDIR/${CONFCASE}-S/$YEAR/$file ] ; then ln -sf $WORKDIR/${CONFCASE}-S/$YEAR/$file . ; else echo 'file not found exit' ; exit ; fi ; fi ; else echo 'file is here'; fi ;
             done  ; }
# get monhtly mean file
getmonthlymean() { for f  in $( ssh $REMOTE_USER@${login_node} " ls $SDIR/$MEANY/${CONFCASE}_y${YEAR}m$2_$1.nc " )  ; do
              file=`basename $f` ; if [ ! -f $file ] ; then if [ ! -f $WORKDIR/${CONFCASE}-MEAN/$YEAR/$file ] ; then scp $REMOTE_USER@${login_node}:$f . ; else ln -sf $WORKDIR/${CONFCASE}-MEAN/$YEAR/$file . ; fi ; else echo 'file is here'; fi ;
             done  ; }
# getannual file
getannualmean() { for f  in $( ssh $REMOTE_USER@${login_node} " ls $SDIR/$MEANY/${CONFCASE}_y${YEAR}_ANNUAL_$1.nc " )  ; do
              file=`basename $f` ; if [ ! -f $file ] ; then if [ ! -f $WORKDIR/${CONFCASE}-MEAN/$YEAR/$file ] ; then scp $REMOTE_USER@${login_node}:$f . ; else ln -sf $WORKDIR/${CONFCASE}-MEAN/$YEAR/$file . ; fi ; else echo 'file is here'; fi ;
             done  ; }
#########################################################################
# getmask file
getmask() { file=$1 ; if [ ! -f $file ] ; then if [ ! -f $WORKDIR/$CONFIG-I/$file ] ; then scp $REMOTE_USER@${login_node}:$SDIR/$IDIR/$file mask.nc ; else ln -sf $WORKDIR/${CONFCASE}-I/$file mask.nc ; fi ; else echo 'file is here'; fi ; }
# getmesh_hgr file
getmesh_hgr() { file=$1 ; if [ ! -f $file ] ; then if [ ! -f $WORKDIR/$CONFIG-I/$file ] ; then scp $REMOTE_USER@${login_node}:$SDIR/$IDIR/$file mesh_hgr.nc ; else ln -sf $WORKDIR/${CONFCASE}-I/$file mesh_hgr.nc ; fi ; else echo 'file is here'; fi ; }
#getmesh_zgr file
getmesh_zgr() { file=$1 ; if [ ! -f $file ] ; then if [ ! -f $WORKDIR/$CONFIG-I/$file ] ; then scp $REMOTE_USER@${login_node}:$SDIR/$IDIR/$file mesh_zgr.nc ; else ln -sf $WORKDIR/${CONFCASE}-I/$file mesh_zgr.nc ; fi ; else echo 'file is here'; fi ; }
#############################################################################
# get levitus file
getlevitus() { file=$1_Levitus-${CONFIG}.nc ; if [ ! -f $file ] ; then if [ ! -f $WORKDIR/$CONFIG-I/$file ] ; then scp $REMOTE_USER@${login_node}:$SDIR/$IDIR/INITIAL/$file levitus.nc ; else cp -f $WORKDIR/${CONFIG}-I/$file levitus.nc ; fi ; else echo 'file is here'; fi; 
if [ $1 = votemper ] ; then flevitus=$T_levitus ; fi ;
if [ $1 = vosaline ] ; then flevitus=$S_levitus ; fi ;
$BIN/cdfmoy levitus.nc ; $BIN/cdfmltmask cdfmoy.nc mask.nc $1 T
mv cdfmoy.nc_masked $flevitus ; rm cdfmoy.nc cdfmoy2.nc ; }
#############################################################################
# rm annual file
rmannualmean() { rm -f ${CONFCASE}_y${YEAR}_ANNUAL_$1.nc ; }
# rm monthly file
rmmonthlymean() { rm -f ${CONFCASE}_y${YEAR}m${2}_$1.nc ; }
# rm mask file
rmmask() { rm -f mask.nc ; }
# rm meshhgr file
rmmesh_hgr() { rm -f mesh_hgr.nc ; }
# rm meshzgr file
rmmesh_zgr() { rm -f mesh_zgr.nc ; }
# chkfile : Usage: chkfile gaya_file
############################################################################
#save file (for big file, used ssh+cp and not scp)
savemeanfile() { ssh $USER@$login_node " cd $R_MONITOR/$YEAR ; cp $1 $SDIR/$MEANY/$1 " ; }
savediagfile() { scp $1 $USER@$login_node:$SDIR/$DIAGS/$1 ; } 
############################################################################
#    check if a file exists on gaya, return present or absent.
chkfile() { ssh $REMOTE_USER@$login_node " if [ -f $1 ] ; then echo present ;\
                       else echo absent ; fi " ; }
# chkdir  : Usage: chkdir local_dir
#   check the existence of a directory. Create it if not present
chkdir() { if [ ! -d $1 ] ; then mkdir $1 ; fi  ; }
############################################################################
#function for TRPSIG
reorganize() {
# for f in *.bimg ; do
#   tmp=${f#${CONFCASE}_y????m??_} ; dir=${tmp%_*.bimg}
#   if [ ! -d $dir ] ; then mkdir $dir ; fi
#   mv $f $dir
# done
 if [ ! -d LST ] ; then  mkdir LST ; fi
 if [ ! -d TRPSIG ] ; then  mkdir TRPSIG ; fi
 mv *.lst LST/. ; mv *trpsig.txt TRPSIG/. ; mv *.bimg TRPSIG/.
 cd TRPSIG
 for f in *.bimg ; do
   tmp=${f#${CONFCASE}_y????m??_} ; dir=${tmp%_*.bimg}
   if [ ! -d $dir ] ; then mkdir $dir ; fi
   mv $f $dir
 done
 cd ..
}
#        we suppose that section name starts with 2 digits 01_ 02_ 03_ etc ...
mean() {
# section name are codes as 2 digit_Capitalized_Name (eg: 01_Denmark_Strait or 07_Bab_el_Mandeb )
cd TRPSIG
for stnam in [0-9][0-9]_[a-zA-Z]* ; do
  cd $stnam
  #printf "%s"  "Working for station $stnam  "
  echo "Working for station $stnam  "
  # note that bimgmoy4 and bimgcaltrans exec are in cdftools-2.0 (extension ...)
  #for d in ???? ; do
    #printf "%4d " $(( YEAR ))
    echo $YEAR
    $BIN/bimgmoy4 ${CONFCASE}_y*trpsig.bimg > /dev/null
    mv moy.bimg ${CONFCASE}_y${YEAR}_${stnam}_trpsig.bimg
   # (2.2) : translate results into txt file foreach section
    $BIN/bimgcaltrans ${CONFCASE}_y${YEAR}_${stnam}_trpsig.bimg > ${CONFCASE}_y${YEAR}_${stnam}_trpsig.txt
    mv ${CONFCASE}_y${YEAR}_${stnam}_trpsig.txt ../.
   cd ../
  done
cd ..
  #printf "\n"
}
#mkvt() {  }
# FROM gaya monitor.ksh
#############################################################################
# cptoweb : Usage: cptoweb  file.mtl
#    rcp the matlab file to the corresponding DATA dir of the website
cptoweb() { rcp $1 \
       apache@meolipc.hmg.inpg.fr:web/DRAKKAR/$CONFIG/$CONFCASE/DATA/ ; }

# chkdirw  : Usage: chkdirw web_site_directory
#   check the existence of a dir. on the web site. Create it if not present
chkdirw() { rsh meolipc.hmg.inpg.fr -l apache " if [ ! -d web/DRAKKAR/$1 ] ; \
            then mkdir web/DRAKKAR/$1 ; fi " ; }





