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
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class file_informations
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                               :: jk, jt, jvar      ! dummy loop index
  INTEGER(KIND=4)                               :: narg, iargc       ! argument numbers
  INTEGER(KIND=4)                               :: ijarg             ! argument counter
  INTEGER(KIND=4)                               :: iilook, ijlook      ! look position
  INTEGER(KIND=4)                               :: npiglo, npjglo    ! size of the domain
  INTEGER(KIND=4)                               :: npk, npt, nvars   ! vetical, time size, number of variables
  INTEGER(KIND=4)                               :: ikx=1, iky=1, ikz ! dims of netcdf output file
  INTEGER(KIND=4)                               :: nboutput=1        ! number of variables to write in cdf output
  INTEGER(KIND=4)                               :: ncout, ierr       ! ncid and error flag for cdfio
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: ipk, id_varout    ! vertical size and id of output variables

  REAL(KIND=4)                                  :: rdep              ! vertical interpolation stuff
  REAL(KIND=4), DIMENSION(:),       ALLOCATABLE :: gdept, rprofile   ! depth and profile values
  REAL(KIND=4), DIMENSION(:),       ALLOCATABLE :: gdepout           ! output depth array
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

  LOGICAL                                       :: l_vert_interp=.FALSE. ! flag for -dep option
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfprofile  -f IN-file -v IN-var -IJ I J [-dep depth ] [-o OUT-file]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Extracts a vertical profile at location I J, for the IN-var variable of'
     PRINT *,'       the IN-file. A depth can optionally be specified.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'      -f  IN-file : input file to work with.'
     PRINT *,'      -v  IN-var  : variable name whose profile is requested.'
     PRINT *,'      -IJ I J     : I, J position of the point to extract from file.' 
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-dep depth] : specify a depth where vertical value will be'
     PRINT *,'                     interpolated.'
     PRINT *,'       [-o OUT-file]: specify output file name instead of ',TRIM(cf_out)
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        none '
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out),' unless -o option is used.'
     PRINT *,'          variable : name given as argument.'
     PRINT *,'       Profile is also written on standard output.'
     STOP
  ENDIF

  ijarg = 1

  DO WHILE ( ijarg < narg )
     CALL getarg (ijarg, cldum ) ; ijarg = ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-f'   ) ; CALL getarg (ijarg, cf_in ) ; ijarg = ijarg+1  
     CASE ( '-v'   ) ; CALL getarg (ijarg, cv_in ) ; ijarg = ijarg+1  
     CASE ( '-IJ'  ) ; CALL getarg (ijarg, cldum ) ; ijarg = ijarg+1  ; READ(cldum,*) iilook
        ;              CALL getarg (ijarg, cldum ) ; ijarg = ijarg+1  ; READ(cldum,*) ijlook
        ! options
     CASE ( '-dep' ) ; CALL getarg (ijarg, cldum ) ; ijarg = ijarg+1  ; READ(cldum,*) rdep ; l_vert_interp = .TRUE.
     CASE ( '-o'   ) ; CALL getarg (ijarg, cf_out) ; ijarg = ijarg+1
     CASE  DEFAULT   ; PRINT *,' ERROR :',TRIM(cldum),' : unknown option.' ; STOP
     END SELECT
  ENDDO

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

  rdumlon(:,:) = rlon(iilook,ijlook)
  rdumlat(:,:) = rlat(iilook,ijlook)

  gdept(:) = getvar1d(cf_in, cv_dep, npk, ierr)
  cv_names(:) = getvarname(cf_in, nvars, stypvar_input)

  CALL CreateOutput

  DO jt=1,npt
     DO jk=1,npk
        v2d(:,:)     = getvar(cf_in, cv_in, jk, npiglo, npjglo, ktime=jt)
        rprofile(jk) = v2d(iilook,ijlook)
        ! netcdf output 
        IF ( .NOT. l_vert_interp ) THEN
           rdummy(1,1)  = rprofile(jk)
           ierr         = putvar(ncout, id_varout(1), rdummy, jk, ikx, iky, ktime=jt)
        ENDIF
     END DO
     PRINT *, 'FILE : ', TRIM(cf_in), ' TIME = ', jt
     PRINT *, '    ', TRIM(cv_dep),'         ', TRIM(cv_in),'(',iilook,',',ijlook,')'

     IF ( l_vert_interp ) THEN
        rdummy(1,1) = vinterp (rprofile, gdept , rdep, npk )
        ierr        = putvar(ncout, id_varout(1), rdummy, 1, ikx, iky, ktime=jt)
        PRINT *, ' Interpolated value is : ', rdummy(1,1)
     ENDIF

     ! Ascii output
     DO jk=1, npk
        PRINT *, gdept(jk), rprofile(jk)
     END DO

  END DO

  ierr = closeout(ncout)

CONTAINS

  REAL(KIND=4) FUNCTION vinterp ( profile, pdept, pdep, kpk)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION vinterp  ***
    !!
    !! ** Purpose :  return the interpolated value at specified depth from
    !!               an  input profile.
    !!
    !! ** Method  :  Linear interpolation
    !!
    !!----------------------------------------------------------------------
    REAL(KIND=4), DIMENSION(kpk), INTENT(in) :: profile  ! vertical profile
    REAL(KIND=4), DIMENSION(kpk), INTENT(in) :: pdept    ! depth array
    REAL(KIND=4),                 INTENT(in) :: pdep     ! required depth
    INTEGER(KIND=4),              INTENT(in) :: kpk      ! number of level in the profile
    !
    INTEGER(KIND=4)                          :: jk       ! loop index
    INTEGER(KIND=4)                          :: ik0, ik1 ! limit index
    REAL(KIND=8)                             :: dalfa    ! weight
    !!----------------------------------------------------------------------
    ! find interpolation limits for required depth
    DO jk = 1, kpk
       IF ( pdept(jk) > pdep ) THEN
          ik0 = jk - 1
          ik1 = ik0 + 1
          EXIT
       ENDIF
    ENDDO
    dalfa   = ( pdep  - pdept(ik0) ) / ( pdept(ik1) - pdept(ik0) ) 
    vinterp =  dalfa * profile (ik1) + ( 1.d0 - dalfa ) * profile (ik0 )

  END FUNCTION vinterp

  SUBROUTINE CreateOutput
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutput  ***
    !!
    !! ** Purpose :  Create netcdf output file(s) 
    !!
    !! ** Method  :  Use stypvar global description of variables
    !!
    !!----------------------------------------------------------------------
    IF ( l_vert_interp ) THEN
       ipk(:) = 1
       ikz    = 1
       ALLOCATE ( gdepout(ikz) )
       gdepout(1) = rdep
    ELSE
       ipk(:) = npk
       ikz    = npk
       ALLOCATE ( gdepout(ikz) )
       gdepout(:) = gdept(:)
    ENDIF


    DO jvar = 1, nvars
       IF ( cv_names(jvar) == cv_in ) THEN
          stypvar=stypvar_input(jvar)
          EXIT  ! found cv_in
       ENDIF
    ENDDO

    ! create output fileset
    ncout = create      (cf_out, 'none',  ikx,      iky, ikz,     cdep='depth')
    ierr  = createvar   (ncout,  stypvar, nboutput, ipk, id_varout            )
    ierr  = putheadervar(ncout,  cf_in,   ikx,      iky, ikz,     pnavlon=rdumlon, pnavlat=rdumlat, pdep=gdepout)

    tim  = getvar1d(cf_in, cn_vtimec, npt     )
    ierr = putvar1d(ncout, tim,       npt, 'T')

  END SUBROUTINE CreateOutput

END PROGRAM cdfprofile
