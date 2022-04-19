PROGRAM cdfsteric_rho1035
  !!======================================================================
  !!                     ***  PROGRAM  cdfsteric  ***
  !!=====================================================================
  !!  ** Purpose : Compute steric height field from gridT file
  !!               at each levels.
  !!               Store the results on a 3D cdf file.
  !!
  !!  ** Method  : the integral of sum [ (1/rho0) * (rho -rhoref) * dz ]
  !!               rho = in situ sea water density
  !!               rhoref = in situ sea water density at
  !!                     T=0°C and S = 35 su for EOS80
  !!              or
  !!                     CT=0°C an dSA =35.16504 g/kg for TEOS10
  !!               rho0 is set to 1035 kg/m3 as recommanded in the NEMO book
  !!
  !! History : 4.2  : 10/2021  :  W. Llovel     : Original code (from cdfdynh_anom.f90)
  !!                  04/2022  : J.M. Molines   : Add some options 
  !!----------------------------------------------------------------------
  USE cdfio
  USE eos, ONLY : sigmai, eos_init
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class derived_fields
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                           :: jk, jt           ! dummy loop index
  INTEGER(KIND=4)                           :: it               ! time index for vvl
  INTEGER(KIND=4)                           :: ierr             ! working integer
  INTEGER(KIND=4)                           :: narg, iargc      ! browse arguments
  INTEGER(KIND=4)                           :: ijarg            ! argument counter
  INTEGER(KIND=4)                           :: npiglo, npjglo   ! size of the domain
  INTEGER(KIND=4)                           :: npk, npt         !  "           "
  INTEGER(KIND=4)                           :: nlev1, nlev2     ! limit of vertical integration
  INTEGER(KIND=4)                           :: ncout            ! ncid of output fileset
  INTEGER(KIND=4), DIMENSION(1)             :: ipk              ! outptut variables : number of levels,
  INTEGER(KIND=4), DIMENSION(1)             :: id_varout        ! ncdf varid's


  REAL(KIND=4)                              :: zsps             ! Missing value for salinity
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: e3t_1d           ! time counter, vertical level spacing
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: temp, zsal       ! Temperature and salinity at current level
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: temp0, zsal0     ! reference temperature and salinity
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: tmask            ! 2D mask at current level
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: rdep, rdepth     ! depth at current level including SSH
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zssh             ! Sea Surface Heigh

  REAL(KIND=8)                              :: drau0 = 1000.d0  ! density of fresh water
  REAL(KIND=8)                              :: drhoref =1035.d0 ! reference density of sea water
  REAL(KIND=8)                              :: dgrav = 9.81d0   ! gravity
  REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: dtim             ! time counter
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dhdy, dterm      ! dynamic height, working array
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dsig0, dsig      ! In situ density (reference, local)

  CHARACTER(LEN=256)                        :: cf_tfil           ! input file name
  CHARACTER(LEN=256)                        :: cf_sfil           ! input file name for Salinity
  CHARACTER(LEN=256)                        :: cf_out='cdfsteric3d.nc' ! output file name
  CHARACTER(LEN=256)                        :: cv_out='vosteric'   ! output file name
  CHARACTER(LEN=256)                        :: cf_out2d='cdfsteric2d.nc' ! output file name
  CHARACTER(LEN=256)                        :: cv_out2d='sosteric'   ! output file name
  CHARACTER(LEN=256)                        :: cldum            ! working char variable

  TYPE(variable) , DIMENSION(1)             :: stypvar          ! structure for attributes

  LOGICAL                                   :: lout  = .FALSE.   ! flag to use non standard output name
  LOGICAL                                   :: lnc4  = .FALSE.   ! Use nc4 with chunking and deflation
  LOGICAL                                   :: limit = .FALSE.   ! flag set if limit are used (2D case)
  LOGICAL                                   :: ll_teos10  = .FALSE. ! teos10 flag
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfsteric_rho1035 -f T-file [-s S-file ] [-limit lev1 lev2] [-vvl] '
     PRINT *,'             [-rhoref RHO-ref] [-teos10] [-o OUT-file] [-nc4]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'        Compute steric height from T-file given as argument.'
     PRINT *,'        In this tool, the cumulated values (from top to bottom) are saved at'
     PRINT *,'        each model level. '
     PRINT *,'        If the -limit option is used, only the 2D integral of the dynamic '
     PRINT *,'        height anomaly between lev1 and lev2 is saved. '
     PRINT *,'        This program replace cdfhdy ( case using -limit option) and cdfhdy3d '
     PRINT *,'        in the standard case.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f T-file : netcdf file with temperature and salinity.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [ -s S-file:] Give salinity file if not the same than T-file.'
     PRINT *,'       [ -rhoref RHO-ref:] Specify the reference density for SeaWater'
     PRINT *,'          Default is ',drhoref, ' kg/m3'
     PRINT *,'       [-limit lev1 lev2] : if specified, the program will only output the'
     PRINT *,'              dynamic height anomaly at lev1 with reference at lev2.'
     PRINT *,'       [-vvl] : use time-varying vertical metrics.'
     PRINT *,'       [-o OUT-file] : Specify output file name instead of ',TRIM(cf_out)
     PRINT *,'       [-nc4] : Use netcdf4 output with chunking and deflation level 1.'
     PRINT *,'                This option is effective only if cdftools are compiled with'
     PRINT *,'                a netcdf library supporting chunking and deflation.'
     PRINT *,'       [-teos10] : use TEOS10 equation of state instead of default EOS80'
     PRINT *,'                 Temperature should be conservative temperature (CT) in deg C.'
     PRINT *,'                 Salinity should be absolute salinity (SA) in g/kg.'
     PRINT *,'      '
     PRINT *,'     OPENMP SUPPORT : yes'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        ',TRIM(cn_fmsk),' and ', TRIM(cn_fzgr)
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) ,' unless -o option is used.'
     PRINT *,'         variables : ', TRIM(cv_out),' ( m )'
     PRINT *,'       If -limit option is used :'
     PRINT *,'       netcdf file : ', TRIM(cf_out2d) ,' unless -o option is used.'
     PRINT *,'         variables : ', TRIM(cv_out2d),' ( m )'
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'       replace old tools cdfhdy and cdfhdy3d.'
     PRINT *,'      '
     STOP
  ENDIF

  cf_sfil = 'none'
  ijarg   = 1
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-f'     ) ; CALL getarg (ijarg, cf_tfil) ; ijarg=ijarg+1
        ! options
     CASE ( '-s'     ) ; CALL getarg (ijarg, cf_sfil) ; ijarg=ijarg+1
     CASE ( '-rhoref') ; CALL getarg (ijarg, cldum  ) ; ijarg=ijarg+1 ; READ(cldum,*) drhoref
     CASE ( '-o'     ) ; CALL getarg (ijarg, cf_out ) ; ijarg=ijarg+1
        ;                lout  = .TRUE.
     CASE ( '-nc4'   ) ; lnc4  = .TRUE.
     CASE ( '-vvl'   ) ; lg_vvl= .TRUE.
     CASE ( '-teos10'  ) ; ll_teos10 = .TRUE.

     CASE ( '-limit' ) ; limit = .TRUE.
        ;                CALL getarg (ijarg, cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) nlev1
        ;                CALL getarg (ijarg, cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) nlev2
     CASE DEFAULT      ; PRINT *,' ERROR : ',TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  ENDDO

  IF ( cf_sfil == 'none') cf_sfil = cf_tfil

  CALL eos_init ( ll_teos10 )

  IF ( lout ) cf_out2d = cf_out  !  name of 2D file if -limit is used

  IF ( chkfile(cf_tfil) .OR. chkfile(cf_sfil) .OR. chkfile(cn_fmsk) .OR. chkfile(cn_fzgr) ) STOP 99 ! missing files
  IF ( lg_vvl ) THEN
    cn_fe3t = cf_tfil
    cn_ve3t = cn_ve3tvvl
  ENDIF

  ! Look for Missing value for salinity
  zsps = getspval(cf_sfil, cn_vosaline)

  npiglo = getdim (cf_tfil,cn_x)
  npjglo = getdim (cf_tfil,cn_y)
  npk    = getdim (cf_tfil,cn_z)
  npt    = getdim (cf_tfil,cn_t)

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  IF ( .NOT. limit ) THEN
     nlev1 = 1
     nlev2 = npk
  ENDIF

  ALLOCATE (temp0(npiglo,npjglo), zsal0(npiglo,npjglo), dsig0(npiglo,npjglo) ,tmask(npiglo,npjglo))
  ALLOCATE (temp(npiglo,npjglo), zsal(npiglo,npjglo), dsig(npiglo,npjglo) , dhdy(npiglo,npjglo), dterm(npiglo,npjglo))
  ALLOCATE (rdep(npiglo,npjglo), rdepth(npiglo,npjglo), zssh(npiglo,npjglo), e3t_1d(npk))
  ALLOCATE (dtim(npt))

  CALL CreateOutput

  ! Temperature and salinity for reference profile
  IF ( .NOT. ll_teos10 ) THEN
     temp0(:,:) =  0.
     zsal0(:,:)  = 35.
  ELSE
     temp0(:,:) =  0.
     zsal0(:,:)  = 35.16504
  END IF

  zssh(:,:)  = getvar(cf_tfil,  cn_sossheig, 1,  npiglo, npjglo)

  !e3t_1d(:)  = getvare3(cn_fzgr, cn_ve3t1d, npk)

  DO jt = 1, npt
     IF ( lg_vvl ) THEN ; it=jt
     ELSE               ; it=1
     ENDIF

     PRINT *,' TIME = ', jt, dtim(jt)/86400.,' days'
     dhdy(:,:)   = 0.
     rdepth(:,:) = 0.

     DO jk = nlev1, nlev2

        IF ( lg_vvl ) THEN ; rdep(:,:) =  getvar(cn_fe3t, cn_ve3t, jk,npiglo,npjglo, ktime=it)
  !      ELSE               ; rdep(:,:)  = e3t_1d(jk)
        ELSE               ; rdep(:,:) = getvar(cn_fe3t, cn_ve3t, jk, npiglo , npjglo)
        ENDIF

        tmask(:,:) = getvar(cn_fmsk, cn_tmask, jk, npiglo, npjglo)

        IF ( jk == 1 .AND. .NOT. lg_vvl ) THEN
           rdep(:,:) = rdep(:,:) + zssh(:,:)
           !rdep(:,:) = rdep(:,:) + 0.
        ENDIF

        ! depth at current level, including ssh (used for computation of rho in situ)
        rdepth(:,:) = rdepth(:,:) + rdep(:,:)
        temp(:,:) = getvar(cf_tfil, cn_votemper,  jk ,npiglo, npjglo, ktime=jt)
        zsal(:,:) = getvar(cf_sfil, cn_vosaline,  jk ,npiglo, npjglo, ktime=jt)

        dsig0 = sigmai(temp0, zsal0, rdepth, npiglo, npjglo)
        dsig  = sigmai(temp , zsal , rdepth, npiglo, npjglo)

        ! we compute the term of the integral : (1/rhoref) * sum [ (rho -
        ! rhoref) * dz ]
        !
!       dterm =  -(1.d0 / drhoref) *  ( ( drau0 + dsig(:,:)) - ( drau0 + dsig0(:,:) ) ) * rdep
        dterm =  -(1.d0 / drhoref) *  (   dsig(:,:) - dsig0(:,:) )  * rdep
        ! in land, it seems appropriate to stop the computation
        WHERE(zsal == zsps ) dterm = 0

        dhdy(:,:) = dhdy(:,:) + dterm(:,:)
        ! masked
        ! remove mask
        !!!!!!!!dhdy(:,:) = dhdy(:,:) * tmask(:,:)

        IF ( .NOT. limit ) ierr = putvar(ncout, id_varout(1) ,REAL(dhdy), jk, npiglo, npjglo, ktime=jt)

     END DO  ! loop to next level
     IF ( limit ) ierr = putvar(ncout, id_varout(1) ,REAL(dhdy), 1 , npiglo, npjglo, ktime=jt)
  END DO  ! next time frame

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
    IF ( limit ) THEN
       ipk(:)                       = 1
       stypvar(1)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
       stypvar(1)%cname             = cv_out2d
       stypvar(1)%cunits            = 'm'
       stypvar(1)%rmissing_value    = 0.
       stypvar(1)%valid_min         = -100.
       stypvar(1)%valid_max         = 100.
       stypvar(1)%clong_name        = 'Dynamical height anomaly'
       stypvar(1)%cshort_name       = cv_out2d
       stypvar(1)%conline_operation = 'N/A'
       stypvar(1)%caxis             = 'TYX'

       ! create output fileset (2D case )
       ncout = create      (cf_out2d, cf_tfil,  npiglo, npjglo, 1       , ld_nc4=lnc4 )
       ierr  = createvar   (ncout,    stypvar, 1,      ipk,    id_varout, ld_nc4=lnc4 )
       ierr  = putheadervar(ncout,    cf_tfil,  npiglo, npjglo, 1                     )
    ELSE
       ipk(:)                       = npk
       stypvar(1)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
       stypvar(1)%cname             = cv_out
       stypvar(1)%cunits            = 'm'
       stypvar(1)%rmissing_value    = 0.
       stypvar(1)%valid_min         = -100.
       stypvar(1)%valid_max         = 100.
       stypvar(1)%clong_name        = 'Dynamical height anomaly'
       stypvar(1)%cshort_name       = cv_out
       stypvar(1)%conline_operation = 'N/A'
       stypvar(1)%caxis             = 'TZYX'

       ! create output fileset (3D case )
       ncout = create      (cf_out, cf_tfil,  npiglo, npjglo, npk     , ld_nc4=lnc4 )
       ierr  = createvar   (ncout,  stypvar, 1,      ipk,    id_varout, ld_nc4=lnc4 )
       ierr  = putheadervar(ncout,  cf_tfil,  npiglo, npjglo, npk       )
    ENDIF

    dtim  = getvar1d(cf_tfil, cn_vtimec, npt)
    ierr  = putvar1d(ncout, dtim,   npt, 'T')

  END SUBROUTINE CreateOutput

END PROGRAM cdfsteric_rho1035


