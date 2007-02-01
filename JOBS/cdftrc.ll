#!/bin/csh
# @ cpu_limit  = 3600
# @ data_limit = 1gb
# Nom du travail LoadLeveler
# @ job_name   = cdftrc
# Fichier de sortie standard du travail       
# @ output     = $(job_name).$(jobid)
# Fichier de sortie d'erreur du travail
# @ error      =  $(job_name).$(jobid)
# @ queue                   


set echo

set CONFIG=ORCA025
set CASE=G50

set YEARS=(1950 1951 1952 1953 1954 1955 1956 1957 1958 1959 1960 1961 1962 1963 1964 1965)
#set YEARS=(1963)
set MESH_MASK_ID='ORCA025-G50'  


#
set CONFCASE=${CONFIG}-${CASE}
set CDFTOOLS=~rcli002/CDFTOOLS-2.0

cd $TMPDIR
cp $CDFTOOLS/att.txt .
cp $CDFTOOLS/cdfmean .
cp $CDFTOOLS/cdfzonalmean .
cp $CDFTOOLS/cdfzonalsum .
cp $CDFTOOLS/cdfzonalout .

mfget ${CONFIG}/${CONFIG}-I/${MESH_MASK_ID}_byte_mask.nc mask.nc
mfget ${CONFIG}/${CONFIG}-I/${MESH_MASK_ID}_mesh_hgr.nc mesh_hgr.nc
mfget ${CONFIG}/${CONFIG}-I/${MESH_MASK_ID}_mesh_zgr.nc mesh_zgr.nc

foreach YEAR ( $YEARS )
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

mfput zonalmean.nc  ${CONFIG}/${CONFCASE}-DIAGS/${CONFCASE}_y${YEAR}_TRCzonalmean.nc
mfput zonalsum.nc  ${CONFIG}/${CONFCASE}-DIAGS/${CONFCASE}_y${YEAR}_TRCzonalsum.nc

mfput zonalmean.dat  ${CONFIG}/${CONFCASE}-DIAGS/${CONFCASE}_y${YEAR}_TRCzonalmean.dat
mfput zonalsum.dat  ${CONFIG}/${CONFCASE}-DIAGS/${CONFCASE}_y${YEAR}_TRCzonalsum.dat
mfput zonalsurf.dat  ${CONFIG}/${CONFCASE}-DIAGS/${CONFCASE}_y${YEAR}_TRCzonalsurf.dat


end
