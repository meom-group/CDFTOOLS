PROGRAM cdfprofile
  !!======================================================================
  !!                     ***  PROGRAM  cdfprofile  ***
  !!=====================================================================
  !!  ** Purpose : extract a vertical profile from a CDFfile
  !!
  !!  ** Method  : read (i,j) position of point to extract
  !!               read varname
  !!               print profile
  !!
  !! History : 2.1  : 06/2005  : J.M. Molines : Original code
  !!           3.0  : 01/2011  : J.M. Molines : Doctor norm + Lic.
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

  INTEGER(KIND=4)                               :: jk, jt, jvar      ! dummy loop index
  INTEGER(KIND=4)                               :: narg, iargc       ! argument numbers
  INTEGER(KIND=4)                               :: ilook, jlook      ! look position
  INTEGER(KIND=4)                               :: npiglo, npjglo    ! size of the domain
  INTEGER(KIND=4)                               :: npk, npt, nvars   ! vetical, time size, number of variables
  INTEGER(KIND=4)                               :: ikx=1, iky=1, ikz ! dims of netcdf output file
  INTEGER(KIND=4)                               :: nboutput=1        ! number of variables to write in cdf output
  INTEGER(KIND=4)                               :: ncout, ierr       ! ncid and error flag for cdfio
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: ipk, id_varout    ! vertical size and id of output variables

  REAL(KIND=4), DIMENSION(:),       ALLOCATABLE :: gdept, rprofile   ! depth and profile values
  REAL(KIND=4), DIMENSION(:),       ALLOCATABLE :: tim               ! time counter
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: v2d, rlon, rlat   ! 2d data array, longitude, latitude
  REAL(KIND=4), DIMENSION(1,1)                  :: rdumlon, rdumlat  ! dummy array for output
  REAL(KIND=4), DIMENSION(1,1)                  :: rdummy            ! dummy array

  CHARACTER(LEN=256)                            :: cldum             ! dummy character variable
  CHARACTER(LEN=256)                            :: cf_in             ! input file 
  CHARACTER(LEN=256)                            :: cf_out='profile.nc'
  CHARACTER(LEN=256)                            :: cv_in, cv_dep     ! variable name and depth name
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_names          ! array of var name

  TYPE(variable), DIMENSION(:),     ALLOCATABLE :: stypvar_input     ! structure of input data
  TYPE(variable), DIMENSION(:),     ALLOCATABLE :: stypvar           ! structure of output data
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg /= 4 ) THEN
     PRINT *,' usage : cdfprofile  I J IN-file IN-var '
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Extract a vertical profile at location I J, for a variable' 
     PRINT *,'       in an input file.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       I   J   : I, J position of the point to extract from file.' 
     PRINT *,'       IN-file : input file to work with.'
     PRINT *,'       IN-var  : variable name whose profile is requested.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        none '
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out)
     PRINT *,'          variable : name given as argument.'
     PRINT *,'       Profile is also written on standard output.'
     STOP
  ENDIF

  CALL getarg (1, cldum ) ; READ(cldum,*) ilook
  CALL getarg (2, cldum ) ; READ(cldum,*) jlook
  CALL getarg (3, cf_in)
  CALL getarg (4, cv_in)

  IF ( chkfile(cf_in) ) STOP ! missing file

  npiglo = getdim (cf_in, cn_x)
  npjglo = getdim (cf_in, cn_y)
  npk    = getdim (cf_in, cn_z, cv_dep)
  nvars  = getnvar(cf_in)
  npt    = getdim (cf_in, cn_t)

  ! Allocate arrays
  ALLOCATE ( v2d (npiglo,npjglo), gdept(npk), rprofile(npk), tim(npt) )
  ALLOCATE ( stypvar_input(nvars) ,cv_names(nvars) )
  ALLOCATE ( stypvar(nboutput), ipk(nboutput), id_varout(nboutput) )
  ALLOCATE ( rlon(npiglo,npjglo), rlat(npiglo,npjglo))

  rlon(:,:)= getvar(cf_in, cn_vlon2d, 1, npiglo, npjglo)
  rlat(:,:)= getvar(cf_in, cn_vlat2d, 1, npiglo, npjglo)

  rdumlon(:,:) = rlon(ilook,jlook)
  rdumlat(:,:) = rlat(ilook,jlook)

  ipk(:) = npk

  cv_names(:) = getvarname(cf_in, nvars, stypvar_input)

  DO jvar = 1, nvars
     IF ( cv_names(jvar) == cv_in ) THEN
        stypvar=stypvar_input(jvar)
        EXIT  ! found cv_in
     ENDIF
  ENDDO

  gdept(:) = getvar1d(cf_in, cv_dep, npk, ierr)
  ikz      = npk

  ! create output fileset
  ncout = create      (cf_out, 'none',  ikx,      iky, npk,     cdep='depth')
  ierr  = createvar   (ncout,  stypvar, nboutput, ipk, id_varout            )
  ierr  = putheadervar(ncout,  cf_in,   ikx,      iky, ikz,     pnavlon=rdumlon, pnavlat=rdumlat, pdep=gdept)

  tim  = getvar1d(cf_in, cn_vtimec, npt     )
  ierr = putvar1d(ncout, tim,       npt, 'T')

  DO jt=1,npt
     DO jk=1,npk
        v2d(:,:)     = getvar(cf_in, cv_in, jk, npiglo, npjglo, ktime=jt)
        rprofile(jk) = v2d(ilook,jlook)
        ! netcdf output 
        rdummy(1,1)  = rprofile(jk)
        ierr         = putvar(ncout, id_varout(1), rdummy, jk, ikx, iky, ktime=jt)
     END DO
     PRINT *, 'FILE : ', TRIM(cf_in), ' TIME = ', jt
     PRINT *, '    ', TRIM(cv_dep),'         ', TRIM(cv_in),'(',ilook,',',jlook,')'

     DO jk=1, npk
        PRINT *, gdept(jk), rprofile(jk)
     END DO
  END DO

  ierr = closeout(ncout)

END PROGRAM cdfprofile
