PROGRAM cdfzonalout
  !!======================================================================
  !!                     ***  PROGRAM  cdfzonalout  ***
  !!=====================================================================
  !!  ** Purpose : Output zonal mean/integral as ascii files
  !!
  !!  ** Method  : Read zonalmean or zonalsum file, determine 1D variable
  !!               and dump them on the standard output.
  !!
  !! History : 2.1  : 02/2006  : J.M. Molines : Original code
  !!           3.0  : 05/2011  : J.M. Molines : Doctor norm + Lic.
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2011, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                               :: jj, jvar, jt  ! dummy loop index
  INTEGER(KIND=4)                               :: ivar          ! variable counter
  INTEGER(KIND=4)                               :: narg, iargc   ! command line 
  INTEGER(KIND=4)                               :: npjglo, npt   ! size of the domain
  INTEGER(KIND=4)                               :: nvarin, nvar  ! variables count
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: ipki          ! input ipk variables
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: id_varin      ! input variables id's

  REAL(KIND=4), DIMENSION (:),      ALLOCATABLE :: tim           ! time counter
  REAL(KIND=4), DIMENSION (:,:),    ALLOCATABLE :: zdumlat       ! latitude for i = north pole
  REAL(KIND=4), DIMENSION (:,:,:),  ALLOCATABLE :: zv            ! data values

  TYPE(variable), DIMENSION(:),     ALLOCATABLE :: stypvar       ! dummy structure

  CHARACTER(LEN=256)                            :: cf_zonal      ! input file name
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_names      ! input variable names
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()

  IF ( narg == 0 ) THEN
     PRINT *,' usage :  cdfzonalout ZONAL-file'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'        This is a formatting program for zonal files, either mean or integral.'
     PRINT *,'        It displays results on the standard output from the input zonal file.'
     PRINT *,'        It only works with 1D zonal variables, skipping 2D variables, that'
     PRINT *,'        cannot be easily displayed !'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'        ZONAL-file : input netcdf zonal file produced by one of the zonal'
     PRINT *,'                     tools.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        none'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'        - Standard output,  structured in columns:'
     PRINT *,'             J  LAT  ( zonal mean, var = 1--> nvar) '
     STOP
  ENDIF

  CALL getarg (1, cf_zonal)
  IF ( chkfile(cf_zonal) ) STOP ! missing file

  nvarin  = getnvar(cf_zonal)
  ALLOCATE ( cv_names(nvarin), ipki(nvarin), id_varin(nvarin), stypvar(nvarin)  )

  cv_names(:) = getvarname(cf_zonal, nvarin, stypvar )
  ipki(:)     = getipk    (cf_zonal, nvarin          )

  ! Open standard output with reclen 2048 for avoid wrapping with ifort
  OPEN(6,FORM='FORMATTED',RECL=2048)
  ! look for 1D var ( f(lat) )
  nvar = 0
  DO jvar = 1,nvarin
     ! skip variables such as nav_lon, nav_lat, time_counter deptht ...
     IF (ipki(jvar) == 0 .OR. ipki(jvar) > 1 ) THEN
        cv_names(jvar)='none'
     ELSE
        nvar = nvar + 1          ! count for elligible  input variables
        id_varin(nvar) = jvar    ! use indirect adressing for those variables
     ENDIF
  END DO

  WRITE(6,*) 'Number of 1D variables :', nvar
  DO jvar=1,nvar
     ivar=id_varin(jvar)
     WRITE(6,*) '     ',TRIM(cv_names(ivar))
  ENDDO

  npjglo = getdim (cf_zonal,cn_y)
  npt    = getdim (cf_zonal,cn_t)

  WRITE(6,*)  'npjglo =', npjglo
  WRITE(6,*)  'npt    =', npt

  ! Allocate arrays
  ALLOCATE ( zv(1,npjglo,nvar), tim(npt)         )
  ALLOCATE ( zdumlat(1,npjglo) )

  zdumlat(:,:) = getvar  (cf_zonal, 'nav_lat', 1, 1, npjglo)
  tim(:)       = getvar1d(cf_zonal, cn_vtimec, npt         )

  DO jt = 1, npt  ! time loop
     ! main elligible variable loop
     DO jvar = 1, nvar
        ivar         = id_varin(jvar)
        zv(:,:,jvar) = getvar(cf_zonal, cv_names(ivar), 1, 1, npjglo, ktime=jt)
     END DO ! next variable

     WRITE(6,*) ' JT = ', jt, ' TIME = ', tim(jt)
     WRITE(6,*) ' J  LAT ', (TRIM(cv_names(id_varin(jvar))),' ',jvar=1,nvar)

     DO jj=npjglo,1,-1
        WRITE(6,*)  jj, zdumlat(1,jj), zv(1,jj,1:nvar)
     ENDDO
  ENDDO


END PROGRAM cdfzonalout
