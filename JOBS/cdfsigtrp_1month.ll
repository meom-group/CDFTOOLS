#!/bin/csh
# @ cpu_limit  = 7200
# @ data_limit = 2gb
# Nom du travail LoadLeveler
# @ job_name   = cdfsigtrp1m
# Fichier de sortie standard du travail       
# @ output     = $(job_name).$(jobid)
# Fichier de sortie d'erreur du travail
# @ error      =  $(job_name).$(jobid)
# @ queue                   


set echo


set CONFIG=ORCA05
set CASE=G50

set YEARS=( 1949 \
    1950 1951 1952 1953 1954 1955 1956 1957 1958 1959 \
    1960 1961 1962 1963 1964 1965 1966 1967 1968 1969 \
    1970 1971 1972 1973 1974 1975 1976 1977 1978 1979 \
    1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 \
    1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 \
    2000 2001 2002 2003 2004 )
set MESH_MASK_ID='ORCA05-G50'

set CONFCASE=${CONFIG}-${CASE}

set CDFTOOLS=~rcli002/CDFTOOLS-2.0
set P_CTL=$HOME/RUN_${CONFIG}/${CONFCASE}/CTL


cd $TMPDIR
cp $CDFTOOLS/cdfsigtrp ./
cp $CDFTOOLS/att.txt .
cp $P_CTL/dens_section.dat .

mfget ${CONFIG}/${CONFIG}-I/${MESH_MASK_ID}_mesh_hgr.nc mesh_hgr.nc
mfget ${CONFIG}/${CONFIG}-I/${MESH_MASK_ID}_mesh_zgr.nc mesh_zgr.nc

foreach year ( $YEARS )

  rsh gaya mkdir ${CONFIG}/${CONFCASE}-TRPSIG/$year/

  foreach tfich (`rsh gaya ls $CONFIG/${CONFCASE}-MEAN/$year/\*m\?\?_gridT.nc ` )
      set ufich=`echo  $tfich | sed -e 's/gridT/gridU/' `
      set vfich=`echo  $tfich | sed -e 's/gridT/gridV/' `

      mfget $tfich ./
      mfget $ufich ./
      mfget $vfich ./

      set tfich=`basename $tfich`
      set ufich=`basename $ufich`
      set vfich=`basename $vfich`
  
     set tag=`echo $tfich | sed -e "s/${CONFCASE}_//" -e 's/_gridT.nc//'`

       

     echo $tag > ${CONFCASE}_${tag}_trpsig_monitor.lst


    ./cdfsigtrp $tfich $ufich $vfich 21 30 180 -bimg -print  >>  ${CONFCASE}_y${tag}_trpsig_monitor.lst

    mfput  ${CONFCASE}_y${tag}_trpsig_monitor.lst  ${CONFIG}/${CONFCASE}-TRPSIG/$year/
    foreach b (*.bimg)
        mv  $b ${CONFCASE}_y${tag}_$b
        mfput ${CONFCASE}_y${tag}_$b  ${CONFIG}/${CONFCASE}-TRPSIG/$year/
    end

    \rm *.bimg
    mv trpsig.txt ${CONFCASE}_y${tag}_trpsig.txt
    mfput ${CONFCASE}_y${tag}_trpsig.txt ${CONFIG}/${CONFCASE}-TRPSIG/$year/

  end

end


#
