PROGRAM cdfthic
  !!======================================================================
  !!                      ***  PROGRAM  cdfthic   ***
  !!======================================================================
  !!  ** Purpose : Compute water column thickness on either T, U, or V grid.
  !!               Handle time-varying thickness either using the -vvl option
  !!               (i.e. using time-varying e3[tuv]), or using the -ssh option
  !!               (i.e. using constant e3[tuv] and varying ssh).
  !!               Handle full step configuration using the -full option.
  !!
  !!  ** Method  : Compute the vertical sum of e3[tuv]*[tuv]mask 
  !!                 (+ssh if used with the -ssh option)
  !!
  !!               Based on cdfvertmean routine
  !!
  !! History : 0.1  : 09/2017  : N. Jourdain
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  USE modutils
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class integration
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                            :: jk, jvar, jvarin, jt ! dummy loop index
  INTEGER(KIND=4)                            :: it                   ! time index for vvl
  INTEGER(KIND=4)                            :: ierr, iko            ! working integer
  INTEGER(KIND=4)                            :: narg, iargc, ijarg   ! command line 
  INTEGER(KIND=4)                            :: npiglo, npjglo       ! size of the domain
  INTEGER(KIND=4)                            :: npk, npt             ! size of the domain
  INTEGER(KIND=4)                            :: ncout                ! ID for netcdf file
  INTEGER(KIND=4)                            :: nvout=1              ! number of output variables
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: ipk, id_varout       ! levels and varid's of output vars

  REAL(KIND=4), PARAMETER                    :: ppspval= 9999.99     ! missing value
  REAL(KIND=4), PARAMETER                    :: eps=1.e-9            ! epsilon
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE  :: e3_1d                ! vertical metrics in case of full step
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE  :: e3, ssh              ! vertical metric
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE  :: area1, area2, ztmp   ! used to interpolate ssh
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE  :: tmask                ! npiglo x npjglo

  REAL(KIND=8), DIMENSION(:),   ALLOCATABLE  :: dtim                 ! time counter
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE  :: dl_vint1             ! verticall int quantity         

  CHARACTER(LEN=1)                           :: cgrid='T'            ! grod used for thickness computation
  CHARACTER(LEN=256)                         :: cf_in, cn_fe3, cn_ve3! input file and variable
  CHARACTER(LEN=256)                         :: cf_out='thic.nc'     ! output file 
  CHARACTER(LEN=256)                         :: cldum                ! dummy string for command line browsing

  LOGICAL                                    :: lfull  =.FALSE.      ! flag for full step computation
  LOGICAL                                    :: lchk   =.FALSE.      ! flag for missing files
  LOGICAL                                    :: lssh   =.FALSE.      ! flag for using ssh for top layer thickness

  TYPE(variable), DIMENSION(:), ALLOCATABLE  :: stypvar              ! extension for attributes
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfthic [-vvl <file>] [-ssh <file>] [-full]'
     PRINT *,'                 [-U] [-V] [-o <OUT-file>]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the water column thickness at T, U, or V points.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'        (default)        : partial step computation (constant grid metrics)'
     PRINT *,'                           on T points'
     PRINT *,'        -full            : full step computation (constant grid metrics)' 
     PRINT *,'        -vvl INPUT-file  : directly use time-varying vertical metrics in INPUT-file'
     PRINT *,'        -ssh INPUT-file  : use time-varying ssh and initial grid properties' 
     PRINT *,'        -U               : computation on U points instead of default T points'
     PRINT *,'        -V               : computation on V points instead of default T points'
     PRINT *,'        -o OUT-file      : use specified output file instead of ', TRIM(cf_out)
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'      ', TRIM(cn_fzgr),', and ',TRIM(cn_fhgr),', and ',TRIM(cn_fmsk)
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file :  thic.nc (or specified with -o option)'
     PRINT *,'         variables :  thic_T, thic_U, or thic_V (according to grid)'
     PRINT *,'                      ssh_U or ssh_V (if -ssh option with U or V grids)'
     PRINT *,'      '
     STOP 
  ENDIF

  ! browse command line
  ijarg = 1   
  DO WHILE ( ijarg <= narg ) 
     CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1
     SELECT CASE ( cldum)
     CASE ( '-vvl'  ) ; CALL getarg (ijarg, cf_in    ) ; lg_vvl = .TRUE. ; ijarg = ijarg + 1
     CASE ( '-ssh'  ) ; CALL getarg (ijarg, cf_in    ) ; lssh   = .TRUE. ; ijarg = ijarg + 1
     CASE ( '-full' ) ; lfull  = .TRUE. 
     CASE ( '-U'    ) ; cgrid  = 'U'
     CASE ( '-V'    ) ; cgrid  = 'V'
     CASE ( '-o'    ) ; CALL getarg (ijarg, cf_out   ) ; ijarg = ijarg + 1
     CASE DEFAULT     ; PRINT *,' ERROR : ',TRIM(cldum),' : unknown options.' ; STOP 99
     END SELECT
  END DO

  ijarg = 0
  IF ( lg_vvl ) ijarg = ijarg + 1
  IF ( lssh   ) ijarg = ijarg + 1
  IF ( lfull  ) ijarg = ijarg + 1
  IF ( ijarg .GT. 1 ) PRINT *,' ERROR : choose no more than one option among -vvl -ssh -full' ; STOP 99

  IF ( lg_vvl ) THEN
    cn_fe3 = cf_in
    SELECT CASE (cgrid)
      CASE ('T') cn_ve3 = cn_ve3tvvl
      CASE ('U') cn_ve3 = cn_ve3uvvl
      CASE ('V') cn_ve3 = cn_ve3vvvl
    END SELECT
  ELSE
    SELECT CASE (cgrid)
      CASE ('T') cn_fe3 = cn_fe3t ; cn_ve3 = cn_ve3t
      CASE ('U') cn_fe3 = cn_fe3u ; cn_ve3 = cn_ve3u
      CASE ('V') cn_fe3 = cn_fe3v ; cn_ve3 = cn_ve3v
    END SELECT
  ENDIF

  ! Security check
  lchk = chkfile ( cn_fe3  )
  lchk = chkfile ( cn_fzgr ) .OR. lchk
  lchk = chkfile ( cn_fhgr ) .OR. lchk
  lchk = chkfile ( cn_fmsk ) .OR. lchk
  IF ( lchk ) STOP 99 ! missing files

  IF ( lssh .AND. cgrid .NE. 'T' ) nvout=2
  ALLOCATE (stypvar(nvout), ipk(nvout), id_varout(nvout))

  ! log information so far
  PRINT *,' OUTPUT FILE     : ' , TRIM(cf_out)

  npiglo = getdim (cf_in, cn_x )
  npjglo = getdim (cf_in, cn_y )
  npk    = getdim (cf_in, cn_z )
  npt    = getdim (cf_in, cn_t )

  PRINT *, ' NPIGLO = ', npiglo
  PRINT *, ' NPJGLO = ', npjglo
  PRINT *, ' NPK    = ', npk
  PRINT *, ' NPT    = ', npt

  ! Allocate arrays
  ALLOCATE (     dtim( npt                  ) )
  ALLOCATE (    tmask(        npiglo,npjglo ) )
  ALLOCATE (       e3(        npiglo,npjglo ) )
  ALLOCATE ( dl_vint1(        npiglo,npjglo ) )
  IF ( lssh .AND. cgrid .NE. 'T' ) THEN
    ALLOCATE (  area1(        npiglo,npjglo ) )
    ALLOCATE (  area2(        npiglo,npjglo ) )
    ALLOCATE (   ztmp(        npiglo,npjglo ) )
  ENDIF
  IF ( lssh  ) ALLOCATE ( ssh(npiglo,npjglo ) )
  IF ( lfull ) ALLOCATE ( e3_1d( npk ) )

  ! Read vertical axis in full-step-case
  IF ( lfull ) e3_1d(:) = getvare3(cn_fzgr, cn_ve3t1d, npk)

  ! Read area fields to interpolate ssh onto U and V grids
  IF ( lssh .AND. cgrid .NE. 'T' ) THEN
    area1(:,:) = getvar(cn_fhgr, cn_ve1t, 1, npiglo, npjglo)  ! Save memory: e1t is first read into area
    ztmp( :,:) = getvar(cn_fhgr, cn_ve2t, 1, npiglo, npjglo)  ! ..and ztmp is used to read e2t
    area1(:,:) = area1(:,:) * ztmp(:,:)                       ! ..
    SELECT CASE (cgrid)
    CASE ('U')
      area2(:,:) = getvar(cn_fhgr, cn_ve1u, 1, npiglo, npjglo)  ! Save memory: e1u is first read into area
      ztmp( :,:) = getvar(cn_fhgr, cn_ve2u, 1, npiglo, npjglo)  ! ..and ztmp is used to read e2u
      area2(:,:) = area2(:,:) * ztmp(:,:)                       ! ..
    CASE ('V')
      area2(:,:) = getvar(cn_fhgr, cn_ve1v, 1, npiglo, npjglo)  ! Save memory: e1v is first read into area
      ztmp( :,:) = getvar(cn_fhgr, cn_ve2v, 1, npiglo, npjglo)  ! ..and ztmp is used to read e2v
      area2(:,:) = area2(:,:) * ztmp(:,:)                       ! ..
    END SELECT    
  ENDIF

  CALL CreateOutput
  PRINT *, 'Output file initialised ...'

  DO jt=1,npt

     IF ( lg_vvl ) THEN ; it=jt
     ELSE               ; it=1
     ENDIF

     dl_vint1(:,:) = 0.d0

     DO jk = 1, npk
        ! Get values at jk
        tmask(:,:) = getvar(cn_fmsk, cn_tmask, jk, npiglo, npjglo )
        ! Get e3 at level jk ( ps...)
        IF ( lfull ) THEN ; e3(:,:) = e3_1d(jk)
        ELSE              ; e3(:,:) = getvar(cn_fe3, cn_ve3, jk, npiglo, npjglo, ktime=it,  ldiom=.NOT.lg_vvl)
        ENDIF
        ! Vertical integration
        dl_vint1(:,:) = dl_vint1(:,:) + e3(:,:) * tmask(:,:) * 1.d0
     END DO

     IF ( lssh ) THEN
       e3(:,:) = 0.e0
       e3(:,:) = getvar(cf_in, cn_sossheig, 1, npiglo, npjglo, ktime=jt) ! e3 <= ssh (to save memory)
       SELECT CASE (cgrid)
         CASE ('T')
           ssh(:,:) = e3(:,:)
         CASE ('U') ! ssh at U points
           ssh(npiglo,:) = e3(npiglo,:) ! to improve if periodic domain
           ssh(1:npiglo-1,:) = 0.5 * (  e3(1:npiglo-1,:) * area1(1:npiglo-1,:) &
           &                          + e3(2:npiglo  ,:) * area1(2:npiglo  ,:) ) / ( eps + area2(1:npiglo-1,:) )
         CASE ('V') ! ssh at V points
           ssh(:,npjglo) = e3(:,npjglo)
           ssh(:,1:npjglo-1) = 0.5 * (  e3(:,1:npjglo-1) * area1(:,1:npjglo-1) &
           &                          + e3(:,2:npjglo  ) * area1(:,2:npjglo  ) ) / ( eps + area2(:,1:npjglo-1) )
       END SELECT
       dl_vint1(:,:) = dl_vint1(:,:) + ssh(:,:) * 1.d0
     ENDIF

     ! Output to netcdf file 
     ierr = putvar(ncout, id_varout(1), REAL(dl_vint1), jk, 1, 1, ktime=jt)
     IF ( lssh .AND. cgrid .NE. 'T' ) THEN
       ierr = putvar(ncout, id_varout(2), ssh, jk, 1, 1, ktime=jt)
     ENDIF

  END DO  ! loop on time

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
    ipk(:)                       = npk

    stypvar(1)%cname             = 'thic_'//TRIM(cgrid)
    stypvar(1)%cunits            = 'm'
    stypvar(1)%rmissing_value    = ppspval
    stypvar(1)%valid_min         = 0.0
    stypvar(1)%valid_max         = 10000.0
    stypvar(1)%clong_name        = 'Water column thickness'
    stypvar(1)%cshort_name       = 'thickness'
    stypvar(1)%conline_operation = 'N/A'
    stypvar(1)%caxis             = 'TYX'

    IF ( lssh .AND. cgrid .NE. 'T' ) THEN
      stypvar(2)%cname             = 'ssh_'//TRIM(cgrid)
      stypvar(2)%cunits            = 'm'
      stypvar(2)%rmissing_value    = ppspval
      stypvar(2)%valid_min         = 0.0
      stypvar(2)%valid_max         = 100.0
      stypvar(2)%clong_name        = 'Sea Surface Height'
      stypvar(2)%cshort_name       = 'ssh'
      stypvar(2)%conline_operation = 'N/A'
      stypvar(2)%caxis             = 'TYX'
    ENDIF

    ! Initialize output file
    ncout = create      (cf_out, 'none', npiglo, npjglo, 1  )
    ierr  = createvar   (ncout, stypvar, 1, ipk, id_varout  )
    SELECT CASE (cgrid)
      CASE( 'T' ) ierr  = putheadervar(ncout, 'none', npiglo, npjglo, 1, glamt, gphit)
      CASE( 'U' ) ierr  = putheadervar(ncout, 'none', npiglo, npjglo, 1, glamu, gphiu)
      CASE( 'V' ) ierr  = putheadervar(ncout, 'none', npiglo, npjglo, 1, glamv, gphiv)
    END SELECT

    dtim  = getvar1d    (cf_in, cn_vtimec, npt     )
    ierr  = putvar1d    (ncout, dtim,      npt, 'T')

  END SUBROUTINE CreateOutput

END PROGRAM cdfthic
