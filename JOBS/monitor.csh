#!/bin/csh
# This script is intended to be sourced from a main script. Not Stand Alone

# EKE
#-----
  cp $CDFTOOLS/att.txt .
  cp $CDFTOOLS/cdfrmsssh ./
  cp $CDFTOOLS/cdfeke .
  cp $CDFTOOLS/cdfstdevw ./
  chmod 755 cdfeke cdfrmsssh cdfstdevw


  foreach f ( `rsh gaya ls $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_ANNUAL_grid\[UV\]\*nc ` )
   mfget $f ./
  end
   mfget $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_ANNUAL_gridT2.nc  ./

  ./cdfeke ${CONFCASE}_y${YEAR}_ANNUAL_gridU.nc ${CONFCASE}_y${YEAR}_ANNUAL_gridU2.nc \
    ${CONFCASE}_y${YEAR}_ANNUAL_gridV.nc ${CONFCASE}_y${YEAR}_ANNUAL_gridV2.nc ${CONFCASE}_y${YEAR}_ANNUAL_gridT2.nc

  mv eke.nc ${CONFCASE}_y${YEAR}_EKE.nc
  mfput ${CONFCASE}_y${YEAR}_EKE.nc $CONFIG/${CONFCASE}-MEAN/$YEAR/


# RMS SSH and W
#--------------
   mfget $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_ANNUAL_gridT.nc ./

  ./cdfrmsssh  ${CONFCASE}_y${YEAR}_ANNUAL_gridT.nc ${CONFCASE}_y${YEAR}_ANNUAL_gridT2.nc

  mfput rms.nc $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_ANNUAL_RMSSSH.nc
  \rm rms.nc

   mfget $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_ANNUAL_gridW.nc ./
   mfget $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_ANNUAL_gridW2.nc ./

  ./cdfstdevw  ${CONFCASE}_y${YEAR}_ANNUAL_gridW.nc ${CONFCASE}_y${YEAR}_ANNUAL_gridW2.nc

  mfput rmsw.nc $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_ANNUAL_STDEVW.nc
  \rm rmsw.nc

# Global MEANS
#--------------

cp $CDFTOOLS/cdfmean .
chmod 755  cdfmean

mfget ${CONFIG}/${CONFIG}-I/${MESH_MASK_ID}_byte_mask.nc mask.nc
mfget ${CONFIG}/${CONFIG}-I/${MESH_MASK_ID}_mesh_hgr.nc mesh_hgr.nc
mfget ${CONFIG}/${CONFIG}-I/${MESH_MASK_ID}_mesh_zgr.nc mesh_zgr.nc

  echo $YEAR >  ${CONFCASE}_y${YEAR}_SSHMEAN.txt
  echo $YEAR >  ${CONFCASE}_y${YEAR}_TMEAN.txt
  echo $YEAR >  ${CONFCASE}_y${YEAR}_SMEAN.txt

  mfget  ${CONFIG}/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_ANNUAL_gridT.nc ./
  ./cdfmean  ${CONFCASE}_y${YEAR}_ANNUAL_gridT.nc sossheig T >> ${CONFCASE}_y${YEAR}_SSHMEAN.txt
  ./cdfmean  ${CONFCASE}_y${YEAR}_ANNUAL_gridT.nc votemper T >> ${CONFCASE}_y${YEAR}_TMEAN.txt
  ./cdfmean  ${CONFCASE}_y${YEAR}_ANNUAL_gridT.nc vosaline T >> ${CONFCASE}_y${YEAR}_SMEAN.txt

  mfput ${CONFCASE}_y${YEAR}_SSHMEAN.txt  ${CONFIG}/${CONFCASE}-DIAGS
  mfput ${CONFCASE}_y${YEAR}_TMEAN.txt  ${CONFIG}/${CONFCASE}-DIAGS
  mfput ${CONFCASE}_y${YEAR}_SMEAN.txt  ${CONFIG}/${CONFCASE}-DIAGS

# Ice Volume area and extent for m02 m03   m08 m09
#--------------------------------------------------

cp $CDFTOOLS/cdficediags .
chmod 755 cdficediags

mfget $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}m02_icemod.nc ./
mfget $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}m03_icemod.nc ./
mfget $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}m08_icemod.nc ./
mfget $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}m09_icemod.nc ./

echo '###' $YEAR 02 > ${CONFCASE}_y${YEAR}_ice.txt
./cdficediags ${CONFCASE}_y${YEAR}m02_icemod.nc  >> ${CONFCASE}_y${YEAR}_ice.txt
echo '###' $YEAR 03 >> ${CONFCASE}_y${YEAR}_ice.txt
./cdficediags ${CONFCASE}_y${YEAR}m03_icemod.nc  >> ${CONFCASE}_y${YEAR}_ice.txt
echo '###' $YEAR 08 >> ${CONFCASE}_y${YEAR}_ice.txt
./cdficediags ${CONFCASE}_y${YEAR}m08_icemod.nc  >> ${CONFCASE}_y${YEAR}_ice.txt
echo '###' $YEAR 09 >> ${CONFCASE}_y${YEAR}_ice.txt
./cdficediags ${CONFCASE}_y${YEAR}m09_icemod.nc  >> ${CONFCASE}_y${YEAR}_ice.txt

mfput ${CONFCASE}_y${YEAR}_ice.txt ${CONFIG}/${CONFCASE}-DIAGS

# Vertical T-S profiles off the coast of Portugal for Gib monitoring
#-------------------------------------------------------------------
 \rm -f ${CONFCASE}_y${YEAR}_TGIB.txt
     echo $YEAR > ${CONFCASE}_y${YEAR}_TGIB.txt
    ./cdfmean ${CONFCASE}_y${YEAR}_ANNUAL_gridT.nc  votemper T $GIB  0 0  >>  ${CONFCASE}_y${YEAR}_TGIB.txt
 \rm -f ${CONFCASE}_y${YEAR}_SGIB.txt
     echo $YEAR > ${CONFCASE}_y${YEAR}_SGIB.txt
    ./cdfmean ${CONFCASE}_y${YEAR}_ANNUAL_gridT.nc  vosaline T $GIB  0 0 >>  ${CONFCASE}_y${YEAR}_SGIB.txt

mfput ${CONFCASE}_y${YEAR}_TGIB.txt ${CONFIG}/${CONFCASE}-DIAGS
mfput ${CONFCASE}_y${YEAR}_SGIB.txt ${CONFIG}/${CONFCASE}-DIAGS

# El nino indexes
#----------------
 \rm -f ${CONFCASE}_y${YEAR}_NINO.txt

foreach m ( 01 02 03 04 05 06 07 08 09 10 11 12 )
 foreach f ( `rsh gaya ls $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}m${m}_gridT.nc ` )

   mfget $f ./
   set g=`basename $f`

   echo  -n $YEAR $m >>!   ${CONFCASE}_y${YEAR}_NINO.txt

# nino 1+2   [ -90 W -- -80 W, -10 S -- 10 N ]
  ./cdfmean  $g votemper T $NINO12 1 1 | tail -1 | awk '{ printf " %8.5f 0.00", $6 }'  >> ${CONFCASE}_y${YEAR}_NINO.txt
# nino 3     [ -150 W -- -90 W, -5 S -- 5 N ]
  ./cdfmean  $g votemper T $NINO3 1 1  | tail -1 | awk '{ printf " %8.5f 0.00", $6 }'  >> ${CONFCASE}_y${YEAR}_NINO.txt
# nino 4     [ -200 W -- -150 W, -5 S -- 5 N ]
  ./cdfmean  $g votemper T $NINO4 1 1 | tail -1 | awk '{ printf " %8.5f 0.00", $6 }'  >> ${CONFCASE}_y${YEAR}_NINO.txt
# nino 3.4   [ -170 W -- -120 W, -% S -- % N ]
  ./cdfmean  $g votemper T $NINO34 1 1 | tail -1 | awk '{ printf " %8.5f 0.00\n", $6 }'  >> ${CONFCASE}_y${YEAR}_NINO.txt


\rm $g



 end
end

  mfput ${CONFCASE}_y${YEAR}_NINO.txt  ${CONFIG}/${CONFCASE}-DIAGS


# Transport
#----------
set P_CTL=$HOME/RUN_${CONFIG}/${CONFCASE}/CTL

cp $CDFTOOLS/cdftransportiz .
cp $P_CTL/section.dat .

set year=$YEAR

mfget ${CONFIG}/${CONFCASE}-MEAN/$year/${CONFCASE}_y${year}_ANNUAL_VT.nc .

echo $year > ${CONFCASE}_y${year}_section_monitor.txt

./cdftransportiz ${CONFCASE}_y${year}_ANNUAL_VT.nc \
   ${CONFCASE}_y${year}_ANNUAL_gridU.nc \
   ${CONFCASE}_y${year}_ANNUAL_gridV.nc  < section.dat >> ${CONFCASE}_y${year}_section_monitor.txt

 grep -v Give ${CONFCASE}_y${year}_section_monitor.txt | grep -v level | grep -v IMAX | grep -v FROM > tmp
mv -f tmp ${CONFCASE}_y${year}_section_monitor.txt
mfput  ${CONFCASE}_y${year}_section_monitor.txt  ${CONFIG}/${CONFCASE}-DIAGS/

# Heat and Salt Meridional Transport
#------------------------------------

cp $CDFTOOLS/cdfmhst .
mfget ${CONFIG}/${CONFIG}-I/new_maskglo.nc new_maskglo.nc

./cdfmhst  ${CONFCASE}_y${year}_ANNUAL_VT.nc
mv zonal_heat_trp.dat ${CONFCASE}_y${year}_heattrp.dat
mv zonal_salt_trp.dat ${CONFCASE}_y${year}_salttrp.dat

mfput ${CONFCASE}_y${year}_heattrp.dat ${CONFIG}/${CONFCASE}-DIAGS/
mfput ${CONFCASE}_y${year}_salttrp.dat ${CONFIG}/${CONFCASE}-DIAGS/

# heat transport from surface fluxes
#____________________________________
cp  $CDFTOOLS/cdfhflx .
./cdfhflx  ${CONFCASE}_y${year}_ANNUAL_gridT.nc
mv hflx.out ${CONFCASE}_y${year}_hflx.dat

mfput ${CONFCASE}_y${year}_hflx.dat ${CONFIG}/${CONFCASE}-DIAGS/

# MOC
#----

cp $CDFTOOLS/cdfmoc .

./cdfmoc ${CONFCASE}_y${year}_ANNUAL_gridV.nc
mv moc.nc  ${CONFCASE}_y${year}_MOC.nc
mfput  ${CONFCASE}_y${year}_MOC.nc  ${CONFIG}/${CONFCASE}-MEAN/${year}/
rsh rhodes "cd bin ; cat ovtplot.csh | sed -e ""s/YYYY/$year/"" -e""s/CCCC/$CASE/"" -e ""s/FFFF/$CONFIG/""> moctmp.csh ; ./moctmp.csh "


# MAX and MIN of MOC
#-------------------
cp $CDFTOOLS/cdfmaxmoc .
set f=${CONFCASE}_y${year}_MOC.nc
set outfile=${CONFCASE}_y${year}_minmaxmoc.txt
    echo $year > $outfile
# GLO
printf "%s" 'Glo ' >>  $outfile ; ./cdfmaxmoc $f glo 20 60 500 2000 | grep Maximum >> $outfile
printf "%s" 'Glo ' >>  $outfile ; ./cdfmaxmoc $f glo -40 30 2000 5500 | grep Minimum >> $outfile
# ATL
printf "%s" 'Atl ' >>  $outfile ; ./cdfmaxmoc $f atl 0 60 500 2000 | grep Maximum >> $outfile
printf "%s" 'Atl ' >>  $outfile ; ./cdfmaxmoc $f atl -20 40 2000 5500 | grep Minimum  >> $outfile
#INP
printf "%s" 'Inp ' >>  $outfile ; ./cdfmaxmoc $f inp 15 50 100 1000 | grep Minimum >> $outfile
printf "%s" 'Inp ' >>  $outfile ; ./cdfmaxmoc $f inp -30 20 1000 5500  | grep Minimum >> $outfile
#AUS
printf "%s" 'Aus ' >>  $outfile ; ./cdfmaxmoc $f glo -70 0 0 2000   | grep Maximum >> $outfile
printf "%s" 'Aus ' >>  $outfile ; ./cdfmaxmoc $f glo -70 0 2000 5500  | grep Minimum >> $outfile

mfput $outfile ${CONFIG}/${CONFCASE}-DIAGS/

# Max and Min of MOC at some specific latitudes
set f=${CONFCASE}_y${year}_MOC.nc
set outfile=${CONFIG}-${CASE}_y${year}_maxmoc40.txt
  \rm -f $outfile

    echo $year > $outfile
# GLO  MAX at 40 N and 30S
printf "%s" 'Glo ' >>  $outfile ; cdfmaxmoc $f glo 40 40 500 2000 | grep Maximum >> $outfile
printf "%s" 'Glo ' >>  $outfile ; cdfmaxmoc $f glo -30 -30 500  5500 | grep Maximum >> $outfile
# ATL  MAX at 40N and 30S
printf "%s" 'Atl ' >>  $outfile ; cdfmaxmoc $f atl 40 40 500 2000 | grep Maximum >> $outfile
printf "%s" 'Atl ' >>  $outfile ; cdfmaxmoc $f atl -30 -30  500 5000 | grep Maximum >> $outfile
#INP  Min at 30 S
printf "%s" 'Inp ' >>  $outfile ; cdfmaxmoc $f inp -30 -30 1000 5500  | grep Minimum >> $outfile
#AUS  MAX at 50 S
printf "%s" 'Aus ' >>  $outfile ; cdfmaxmoc $f glo -50 -50 0 2000   | grep Maximum >> $outfile

mfput $outfile ${CONFIG}/${CONFIG}-${CASE}-DIAGS/


# Barotropic Transport
#---------------------
cp $CDFTOOLS/cdfpsi .
chmod 755 cdfpsi
  ./cdfpsi ${CONFCASE}_y${YEAR}_ANNUAL_gridU.nc ${CONFCASE}_y${YEAR}_ANNUAL_gridV.nc 

  mv psi.nc ${CONFCASE}_y${YEAR}_PSI.nc
  mfput ${CONFCASE}_y${YEAR}_PSI.nc $CONFIG/${CONFCASE}-MEAN/$YEAR/

DCT:
  if ( $DCT == 1 ) then
# Density Class transport
#-------------------------
cp $CDFTOOLS/cdfsigtrp ./
cp $P_CTL/dens_section.dat .
cp $CDFTOOLS/JOBS/trpsig_postproc.ksh ./

  rsh gaya mkdir ${CONFIG}/${CONFCASE}-TRPSIG/$YEAR/
  rsh gaya mkdir ${CONFIG}/${CONFCASE}-DIAGS/TRPSIG/

  if ( ! -d ${CONFIG} )                           mkdir ${CONFIG}
  if ( ! -d ${CONFIG}/${CONFCASE}-TRPSIG  )       mkdir ${CONFIG}/${CONFCASE}-TRPSIG
  if ( ! -d ${CONFIG}/${CONFCASE}-TRPSIG/$YEAR/ ) mkdir ${CONFIG}/${CONFCASE}-TRPSIG/$YEAR/

  foreach tfich (`rsh gaya ls $CONFIG/${CONFCASE}-MEAN/$YEAR/\*m\?\?_gridT.nc ` )
      set ufich=`echo  $tfich | sed -e 's/gridT/gridU/' `
      set vfich=`echo  $tfich | sed -e 's/gridT/gridV/' `

      mfget $tfich ./
      mfget $ufich ./
      mfget $vfich ./

      set tfich=`basename $tfich`
      set ufich=`basename $ufich`
      set vfich=`basename $vfich`

     set tag=`echo $tfich | sed -e "s/${CONFCASE}_//" -e 's/_gridT.nc//'`

     echo $tag > ${CONFCASE}_y${tag}_trpsig_monitor.lst


    ./cdfsigtrp $tfich $ufich $vfich 21 30 180 -bimg -print  >>  ${CONFCASE}_y${tag}_trpsig_monitor.lst

    mfput  ${CONFCASE}_y${tag}_trpsig_monitor.lst  ${CONFIG}/${CONFCASE}-TRPSIG/$YEAR/
    mv ${CONFCASE}_y${tag}_trpsig_monitor.lst  ${CONFIG}/${CONFCASE}-TRPSIG/$YEAR/

    foreach b (*.bimg)
        mv  $b ${CONFCASE}_y${tag}_$b
        mfput ${CONFCASE}_y${tag}_$b  ${CONFIG}/${CONFCASE}-TRPSIG/$YEAR/
        mv  ${CONFCASE}_y${tag}_$b  ${CONFIG}/${CONFCASE}-TRPSIG/$YEAR/
    end

    \rm *.bimg $tfich  $ufich $vfich

    mv trpsig.txt ${CONFCASE}_y${tag}_trpsig.txt
    mfput ${CONFCASE}_y${tag}_trpsig.txt ${CONFIG}/${CONFCASE}-TRPSIG/$YEAR/
    mv  ${CONFCASE}_y${tag}_trpsig.txt ${CONFIG}/${CONFCASE}-TRPSIG/$YEAR/

  end
# Launch post processing  ( ex range2.ksh sur rhodes)

  cd ${CONFIG}/${CONFCASE}-TRPSIG
  $TMPDIR/trpsig_postproc.ksh
  cd $TMPDIR

# save results on gaya
  mfput ${CONFCASE}_y*_trpsig.txt $CONFIG/${CONFCASE}-DIAGS/TRPSIG/

# Erase the TRPSIG tree for this current year
  \rm -r ${CONFIG} \rm ${CONFCASE}_y*_trpsig.txt

  endif


# MXL Diagnostics
#-----------------
  cp $CDFTOOLS/cdfmxl .
  chmod 755 cdfmxl 

  rsh gaya mkdir $CONFIG/${CONFCASE}-DIAGS/$YEAR

foreach f ( `rsh gaya ls $CONFIG/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}m0\[39\]_gridT.nc ` )
   mfget $f ./
   set g=`basename $f | sed -e 's/gridT/MXL/' `

  ./cdfmxl `basename $f`

  mfput mxl.nc $CONFIG/${CONFCASE}-MEAN/$YEAR/$g

end

 if ( $TRACER == 1 ) then
# TRACER DIAGS (31/12 of each year)
#-------------

cp $CDFTOOLS/cdfzonalmean .
cp $CDFTOOLS/cdfzonalsum .
cp $CDFTOOLS/cdfzonalout .

# Absolute mean of concentration
  echo -n $YEAR ' '  >  ${CONFCASE}_y${YEAR}_TRCmean.dat

  mfget  ${CONFIG}/${CONFCASE}-S/$YEAR/${CONFCASE}_y${YEAR}m12d31_ptrcT.nc ./

  \rm -f tmp1
  ./cdfmean  ${CONFCASE}_y${YEAR}m12d31_ptrcT.nc  invcfc T > tmp1
  set area=`cat tmp1 |  grep -e 'Mean value at level' | awk ' {print $12}'`
  set mean=`cat tmp1 |  grep -e 'Mean value over the ocean' | awk ' {print $6}'`
  set total=` echo $mean $area |  awk '{print $1 * $2 }' `
  echo -n $total ' ' >> ${CONFCASE}_y${YEAR}_TRCmean.dat

  \rm -f tmp1
  ./cdfmean  ${CONFCASE}_y${YEAR}m12d31_ptrcT.nc  invc14 T > tmp1
  set area=`cat tmp1 |  grep -e 'Mean value at level' | awk ' {print $12}'`
  set mean=`cat tmp1 |  grep -e 'Mean value over the ocean' | awk ' {print $6}'`
  set total=` echo $mean $area |  awk '{print $1 * $2 }' `
  echo $total ' ' >> ${CONFCASE}_y${YEAR}_TRCmean.dat

  mfput ${CONFCASE}_y${YEAR}_TRCmean.dat  ${CONFIG}/${CONFCASE}-DIAGS

# zonal integral of inventories
 ./cdfzonalsum  ${CONFCASE}_y${YEAR}m12d31_ptrcT.nc  T

# zonal means
 ./cdfzonalmean  ${CONFCASE}_y${YEAR}m12d31_ptrcT.nc  T

 ncks -F -d deptht,1,1 -v zocfc11_glo,zobc14_glo,nav_lon,nav_lat zonalmean.nc zonalsurf.nc

# put in ascii format the 1D profiles
 ./cdfzonalout zonalmean.nc > zonalmean.dat
 ./cdfzonalout zonalsum.nc >  zonalsum.dat
 ./cdfzonalout zonalsurf.nc >  zonalsurf.dat

mfput zonalmean.nc  ${CONFIG}/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_TRCzonalmean.nc
mfput zonalsum.nc  ${CONFIG}/${CONFCASE}-MEAN/$YEAR/${CONFCASE}_y${YEAR}_TRCzonalsum.nc

mfput zonalmean.dat  ${CONFIG}/${CONFCASE}-DIAGS/${CONFCASE}_y${YEAR}_TRCzonalmean.dat
mfput zonalsum.dat  ${CONFIG}/${CONFCASE}-DIAGS/${CONFCASE}_y${YEAR}_TRCzonalsum.dat
mfput zonalsurf.dat  ${CONFIG}/${CONFCASE}-DIAGS/${CONFCASE}_y${YEAR}_TRCzonalsurf.dat

 endif
