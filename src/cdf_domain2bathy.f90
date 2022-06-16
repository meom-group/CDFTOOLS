PROGRAM cdf_domain2bathy
  !!======================================================================
  !!                     ***  PROGRAM  cdf_domain2bathy  ***
  !!=====================================================================
  !!  ** Purpose : Create a bathy file from domain_cfg file
  !!
  !!  ** Method  : integrate e3t_0 on the wet points
  !!
  !! History :  4.0  : 10/2020  : J.M. Molines : 
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!----------------------------------------------------------------------

  USE cdfio
  USE modcdfnames

  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2020
  !! $Id$
  !! Copyright (c) 2012, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class domain_file
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                           :: jk, jt              ! dummy loop index
  INTEGER(KIND=4)                           :: ierr,  iko, it   ! working integer
  INTEGER(KIND=4)                           :: narg, iargc, ijarg  ! command line 
  INTEGER(KIND=4)                           :: npiglo, npjglo      ! size of the domain
  INTEGER(KIND=4)                           :: npk, npt            ! size of the domain
  INTEGER(KIND=4)                           :: ncout 
  INTEGER(KIND=4), DIMENSION(1)                :: ipk, id_varout      ! only one output variable
  INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: mbathy           ! working input variable

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e3t                 ! vertical metric
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: rlon, rlat          ! output longitude, latitude

  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dbati               ! bathymetry
  REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: dtim                ! time

  CHARACTER(LEN=256)                        :: cf_in, cf_out       ! input/output file
  CHARACTER(LEN=256)                        :: cv_out              ! output variable name
  CHARACTER(LEN=256)                        :: cldum               ! working variables

  LOGICAL                                   :: lnc4=.FALSE.        ! nc4 flag
  LOGICAL                                   :: lchk                ! file flag
  TYPE(variable), DIMENSION(1)              :: stypvar             ! extension for attributes
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdf_domain2bathy -d DOMAINCFG-file [-nc4] [-o BATI-file]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the bathymetry of a configuration, consistent with the'
     PRINT *,'       domain_cfg file given in argument. This program compute the bathymetry'
     PRINT *,'       from the e3t_0 variable and bottom_level variable.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'         -d DOMAINCFG-file: pass the name of the domain file to work with.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'        [-nc4]  : use netcdf4 output with chunking and deflation'
     PRINT *,'        [-o BATI-file] : use specified output file instead of <IN-var>.nc'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        none'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file :  bati.nc of BATI-file if specified'
     PRINT *,'         variables :  bathy_meter'
     PRINT *,'      '
     STOP 
  ENDIF

  ! browse command line
  ijarg = 1   
  cf_out='none'
  DO WHILE ( ijarg <= narg ) 
     CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1
     SELECT CASE ( cldum)
     CASE ( '-d'    ) ; CALL getarg (ijarg, cf_in  ) ; ijarg = ijarg + 1
        ! options

     CASE ( '-nc4'  ) ; lnc4  = .TRUE. 
     CASE ( '-o'    ) ; CALL getarg (ijarg, cf_out) ; ijarg = ijarg + 1
     CASE DEFAULT     ; PRINT *,' ERROR : ', TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  END DO

  ! Security check
  lchk = chkfile ( cf_in   )
  IF ( lchk ) STOP 99 ! missing files

  IF ( cf_out ==  'none' ) cf_out='bati.nc'

  npiglo = getdim (cf_in, cn_x )
  npjglo = getdim (cf_in, cn_y )
  npk    = getdim (cf_in, 'z' )
  npt    = getdim (cf_in, 't' )

  PRINT *, ' NPIGLO = ', npiglo
  PRINT *, ' NPJGLO = ', npjglo
  PRINT *, ' NPK    = ', npk
  PRINT *, ' NPT    = ', npt

  ! Allocate arrays
  ALLOCATE ( rlon(npiglo,npjglo), rlat(npiglo,npjglo)  )
  ALLOCATE ( mbathy(npiglo,npjglo) )
  ALLOCATE ( dbati(npiglo,npjglo)  )
  ALLOCATE ( e3t(npiglo,npjglo)    )
  ALLOCATE ( dtim(npt)             )

  CALL CreateOutput

  mbathy(:,:) = getvar(cf_in,'bottom_level',1,npiglo, npjglo)
  dbati=0.d0

  DO jk = 1, npk
     e3t(:,:) = getvar(cf_in,'e3t_0',jk, npiglo,npjglo)
     WHERE( mbathy >= jk ) dbati(:,:)=dbati(:,:) + e3t(:,:)
  ENDDO
  ierr = putvar(ncout, id_varout(1) ,REAL(dbati), 1, npiglo, npjglo)
  ierr = closeout(ncout)

CONTAINS

  SUBROUTINE CreateOutput
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutput  ***
    !!
    !! ** Purpose :  Create netcdf output file(s) 
    !!
    !! ** Method  :  Use stypvar global description of variables
    !!
    !!----------------------------------------------------------------------
    ! prepare output variable
    ipk(:)                       = 1
    ! choose chunk size for output ... not easy not used if lnc4=.false. but
    ! anyway ..
    stypvar(1)%ichunk=(/npiglo,MAX(1,npjglo/30),1,1 /)
    cv_out = cn_bathymet

    stypvar(1)%cname             = cv_out
    stypvar(1)%cunits            = 'm'
    stypvar(1)%rmissing_value    = 0.
    stypvar(1)%valid_min         = 0.
    stypvar(1)%valid_max         = 10000.
    stypvar(1)%clong_name        = 'Bathymetry'
    stypvar(1)%cshort_name       =  cv_out
    stypvar(1)%conline_operation = 'N/A'
    stypvar(1)%caxis             = 'TYX'

   rlon =getvar(cf_in, cn_glamt, 1, npiglo, npjglo)  ! nav_lon
   rlat =getvar(cf_in, cn_gphit, 1, npiglo, npjglo)  ! nav_lat

    ncout = create      (cf_out, 'none', npiglo, npjglo, 0 , ld_nc4=lnc4)
    ierr  = createvar   (ncout, stypvar, 1, ipk, id_varout , ld_nc4=lnc4)
    ierr  = putheadervar(ncout, cf_in,   npiglo, npjglo, 0,  pnavlon=rlon, pnavlat=rlat )

    dtim  = getvar1d    (cf_in, cn_vtimec, npt     )
    ierr  = putvar1d    (ncout, dtim,      npt, 'T')

  END SUBROUTINE CreateOutput

END PROGRAM cdf_domain2bathy
