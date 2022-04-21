PROGRAM cdfsteric
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
  !!                  04/2022  : J.M. Molines   : Generalisation of W. Llovel package
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
  INTEGER(KIND=4)                           :: itt,its          ! time index for halosteric/thermosteric
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
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: rdepth           ! depth at current level including SSH
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e3t              ! thickness of layers (ssh taken into account)
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zssh             ! Sea Surface Heigh

  REAL(KIND=8)                              :: drau0 = 1000.d0  ! density of fresh water
  REAL(KIND=8)                              :: drhoref =1035.d0 ! reference density of sea water
  REAL(KIND=8)                              :: dgrav = 9.81d0   ! gravity
  REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: dtim             ! time counter
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dhdy, dterm      ! dynamic height, working array
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dsig0, dsig      ! In situ density (reference, local)

  CHARACTER(LEN=256)                        :: cf_tfil           ! input file name
  CHARACTER(LEN=256)                        :: cf_sfil           ! input file name for Salinity
  CHARACTER(LEN=256)                        :: cf_sshfil         ! input file name for SSH
  CHARACTER(LEN=256)                        :: cf_out='cdfsteric3d.nc' ! output file name
  CHARACTER(LEN=256)                        :: cv_out='vosteric'   ! output file name
  CHARACTER(LEN=256)                        :: cf_out2d='cdfsteric2d.nc' ! output file name
  CHARACTER(LEN=256)                        :: cv_out2d='sosteric'   ! output file name
  CHARACTER(LEN=256)                        :: cldum            ! working char variable

  TYPE(variable) , DIMENSION(1)             :: stypvar          ! structure for attributes

  LOGICAL                                   :: lout  = .FALSE.   ! flag to use non standard output name
  LOGICAL                                   :: lnc4  = .FALSE.   ! Use nc4 with chunking and deflation
  LOGICAL                                   :: limit = .FALSE.   ! flag set if limit are used (2D case)
  LOGICAL                                   :: lhalo = .FALSE.   ! Halo steric flag
  LOGICAL                                   :: lther = .FALSE.   ! Thermo steric flag
  LOGICAL                                   :: ll_teos10  = .FALSE. ! teos10 flag
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfsteric -t T-file [-s S-file ] [-ssh SSH-file]'
     PRINT *,'             [-sshvar SSH-varname] [-limit lev1 lev2] [-vvl] [-HALO] [-THERMO]'
     PRINT *,'             [-rhoref RHO-ref] [-teos10] [-tvar T-varname] [-svar S-varname]'
     PRINT *,'             [-o OUT-file] [-nc4]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'        Compute steric height from T-file given as argument.'
     PRINT *,'        In this tool, the cumulated values (from top to bottom) are saved at'
     PRINT *,'        each model level. '
     PRINT *,'        If the -limit option is used, only the 2D integral of the dynamic '
     PRINT *,'        height anomaly between lev1 and lev2 is saved. '
     PRINT *,'        This program replace cdfhdy ( case using -limit option) and cdfhdy3d '
     PRINT *,'        in the standard case.'
     PRINT *,'        Options -HALO and -THERMO allows the computation of respectively the'
     PRINT *,'        halosteric and thermosteric part of the dynamical ssh. In these latter'
     PRINT *,'        computations, either temperature (halosteric) or salitnity (thermost.)'
     PRINT *,'        are replaced by climatological values, when computing the density.'
     PRINT *,'        It is from user responsability to provide the correct T or S files.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -t T-file : netcdf file with temperature and salinity.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [ -s S-file:] Give salinity file if not the same than T-file.'
     PRINT *,'       [ -ssh SSH-file:] Give SSH file if not the same than T-file.'
     PRINT *,'       [ -rhoref RHO-ref:] Specify the reference density for SeaWater'
     PRINT *,'          Default is ',drhoref, ' kg/m3'
     PRINT *,'       [-HALO :] compute the HALO steric contribution. In this case, we assume'
     PRINT *,'           that temperature is taken as a reference temperature in T-file,'
     PRINT *,'           fixed in time (1 time frame only in T-file.'
     PRINT *,'       [-THERMO :] compute the THERMO steric contribution. In this case, we'
     PRINT *,'           assume that salinity is taken as a reference salinity in S-file'
     PRINT *,'           fixed in time (1 time frame only in S-file.'
     PRINT *,'       [-limit lev1 lev2] : if specified, the program will only output the'
     PRINT *,'              dynamic height anomaly at lev1 with reference at lev2.'
     PRINT *,'              JMM : to be improved..'
     PRINT *,'       [-vvl] : use time-varying vertical metrics.'
     PRINT *,'       [-o OUT-file] : Specify output file name instead of ',TRIM(cf_out)
     PRINT *,'       [-nc4] : Use netcdf4 output with chunking and deflation level 1.'
     PRINT *,'                This option is effective only if cdftools are compiled with'
     PRINT *,'                a netcdf library supporting chunking and deflation.'
     PRINT *,'       [-teos10] : use TEOS10 equation of state instead of default EOS80'
     PRINT *,'                 Temperature should be conservative temperature (CT) in deg C.'
     PRINT *,'                 Salinity should be absolute salinity (SA) in g/kg.'
     PRINT *,'       [-tvar T-varname :] specify name of temperature variable. [',TRIM(cn_votemper),']'
     PRINT *,'       [-svar S-varname :] specify name of salinity variable. [',TRIM(cn_vosaline),']'
     PRINT *,'       [-sshvar SSH-varname :] specify name of ssh variable. [',TRIM(cn_sossheig),']'
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
  cf_sshfil = 'none'
  ijarg   = 1
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-f','-t') ; CALL getarg (ijarg, cf_tfil   ) ; ijarg=ijarg+1
        ! options
     CASE ( '-s'     ) ; CALL getarg (ijarg, cf_sfil   ) ; ijarg=ijarg+1 
     CASE ( '-ssh'   ) ; CALL getarg (ijarg, cf_sshfil ) ; ijarg=ijarg+1 
     CASE ( '-rhoref') ; CALL getarg (ijarg, cldum     ) ; ijarg=ijarg+1 ; READ(cldum,*) drhoref
     CASE ( '-o'     ) ; CALL getarg (ijarg, cf_out    ) ; ijarg=ijarg+1
        ;                lout  = .TRUE.
     CASE ( '-nc4'   ) ; lnc4  = .TRUE.
     CASE ( '-vvl'   ) ; lg_vvl= .TRUE.
     CASE ( '-teos10') ; ll_teos10 = .TRUE.
     CASE ( '-HALO'  ) ; lhalo     = .TRUE.
     CASE ( '-THERMO') ; lther     = .TRUE.
     CASE ( '-tvar'  ) ; CALL getarg (ijarg, cn_votemper) ; ijarg=ijarg+1
     CASE ( '-svar'  ) ; CALL getarg (ijarg, cn_vosaline) ; ijarg=ijarg+1
     CASE ( '-sshvar') ; CALL getarg (ijarg, cn_sossheig) ; ijarg=ijarg+1

     CASE ( '-limit' ) ; limit = .TRUE.
        ;                CALL getarg (ijarg, cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) nlev1
        ;                CALL getarg (ijarg, cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) nlev2
     CASE DEFAULT      ; PRINT *,' ERROR : ',TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  ENDDO

  IF ( cf_sfil   == 'none') cf_sfil   = cf_tfil
  IF ( cf_sshfil == 'none') cf_sshfil = cf_tfil
  IF ( lhalo ) THEN
     PRINT *,' HALO-STERIC  calculation with :'
  ELSEIF (lther ) THEN
     PRINT *,' THERMO-STERIC  calculation with :'
  ELSE
     PRINT *,' STERIC  calculation with :'
  ENDIF
  PRINT *,'     T-File   : ', TRIM(cf_tfil)
  PRINT *,'        T-var : ', TRIM(cn_votemper)
  PRINT *,'     S-File   : ', TRIM(cf_sfil)
  PRINT *,'        S-var : ', TRIM(cn_vosaline)
  PRINT *,'     Rho-ref  : ', drhoref
  IF ( lg_vvl) THEN
    PRINT *,'     VVL      : YES '
  ELSE
    PRINT *,'     VVL      :  NO '
  ENDIF
  IF ( ll_teos10 ) THEN
    PRINT *,'     EOS      :   TEOS10'
  ELSE
    PRINT *,'     EOS      :   EOS80'
  ENDIF
  PRINT *,'     OUT-file : ', TRIM(cf_out)
  IF ( lnc4) THEN
    PRINT *,'        NC4   :  YES '
  ELSE
    PRINT *,'        NC4   :  NO '
  ENDIF

  CALL eos_init ( ll_teos10 )

  IF ( lout ) cf_out2d = cf_out  !  name of 2D file if -limit is used

  IF ( chkfile(cf_tfil) .OR. chkfile(cf_sfil) .OR. chkfile(cn_fmsk) .OR. chkfile(cn_fzgr) .OR. chkfile(cf_sshfil) ) STOP 99 ! missing files
  IF ( lg_vvl ) THEN
    cn_ve3t = cn_ve3tvvl
    IF (  .NOT. chkvar ( cf_tfil, cn_ve3t )) THEN
       cn_fe3t = cf_tfil
    ELSE
       IF ( .NOT.  chkvar ( cf_sfil, cn_ve3t )) THEN
         cn_fe3t = cf_sfil
       ELSE
         PRINT *,' E R R O R : with vvl option, unable to find e3t(time) in input files...'
         STOP 999
       ENDIF
    ENDIF
    PRINT *, '  Using VVL with ',TRIM(cn_ve3t),' taken in ', TRIM(cn_fe3t)
  ENDIF

  ! Look for Missing value for salinity
  !zsps = getspval(cf_sfil, cn_vosaline)
  !PRINT *,   'SALINITY missing value is :', zsps

  npiglo = getdim (cf_tfil,cn_x)
  npjglo = getdim (cf_tfil,cn_y)
  npk    = getdim (cf_tfil,cn_z)
  IF ( lhalo ) THEN
     npt    = getdim (cf_sfil,cn_t)
  ELSE
     ! valid for thermo steric or steric (T file as many time frames)
     npt    = getdim (cf_tfil,cn_t)
  ENDIF

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  IF ( .NOT. limit ) THEN
     nlev1 = 1
     nlev2 = npk
  ENDIF

  PRINT *,'     NLEV1    : ', nlev1 
  PRINT *,'     NLEV2    : ', nlev2 

  ALLOCATE (temp0(npiglo,npjglo), zsal0(npiglo,npjglo), dsig0(npiglo,npjglo) ,tmask(npiglo,npjglo))
  ALLOCATE (temp(npiglo,npjglo), zsal(npiglo,npjglo), dsig(npiglo,npjglo) , dhdy(npiglo,npjglo), dterm(npiglo,npjglo))
  ALLOCATE (e3t(npiglo,npjglo), rdepth(npiglo,npjglo), zssh(npiglo,npjglo), e3t_1d(npk))
  ALLOCATE (dtim(npt))

  CALL CreateOutput

  ! Temperature and salinity for reference profile
  IF ( .NOT. ll_teos10 ) THEN
     temp0(:,:) =  0.   ! Deg Celsius
     zsal0(:,:)  = 35.  ! PSU
  ELSE
     temp0(:,:) =  0.        ! Deg Celsius (Conservative temperature)
     zsal0(:,:)  = 35.16504  ! g/kg        (Absolute Salinity)
  END IF

  !e3t_1d(:)  = getvare3(cn_fzgr, cn_ve3t1d, npk)

  DO jt = 1, npt
     IF ( lg_vvl ) THEN ; it=jt
     ELSE               ; it=1
     ENDIF
     IF ( .NOT. lg_vvl ) THEN
       zssh(:,:)  = getvar(cf_sshfil,  cn_sossheig, 1,  npiglo, npjglo, ktime=jt)
     ENDIF
     itt=jt ; its=jt
     IF ( lhalo ) itt=1   ! one time frame in temperature file
     IF ( lther ) its=1   ! one time frame in salinity file

     PRINT *,' TIME = ', jt, dtim(jt)/86400.,' days'
     dhdy(:,:)   = 0.d0
     rdepth(:,:) = 0.d0

     DO jk = nlev1, nlev2
        !WRITE(6,'(1H+,i2.2)')  jk

        IF ( lg_vvl ) THEN ; e3t(:,:) = getvar(cn_fe3t, cn_ve3t, jk, npiglo, npjglo, ktime=it )
        ELSE               ; e3t(:,:) = getvar(cn_fe3t, cn_ve3t, jk, npiglo, npjglo )
        ENDIF

        tmask(:,:) = getvar(cn_fmsk, cn_tmask, jk, npiglo, npjglo)

        IF ( jk == 1 .AND. .NOT. lg_vvl ) THEN
           ! add ssh to the first level thickness
           e3t(:,:) = e3t(:,:) + zssh(:,:)
        ENDIF

        ! depth at current level, including ssh (used for computation of rho in situ)
        rdepth(:,:) = rdepth(:,:) + e3t(:,:)

        temp(:,:) = getvar(cf_tfil, cn_votemper,  jk ,npiglo, npjglo, ktime=itt)
        zsal(:,:) = getvar(cf_sfil, cn_vosaline,  jk ,npiglo, npjglo, ktime=its)

        dsig0 = sigmai(temp0, zsal0, rdepth, npiglo, npjglo)
        dsig  = sigmai(temp , zsal , rdepth, npiglo, npjglo)

        ! we compute the term of the integral : (1/rhoref) * sum [ (rho -
        ! rhoref) * dz ]
        !
        dterm(:,:) =  -(1.d0 / drhoref) *  ( dsig(:,:) - dsig0(:,:) )  * e3t(:,:)

        ! mask land values
        dterm(:,:)=dterm(:,:)*tmask(:,:)

        dhdy(:,:) = dhdy(:,:) + dterm(:,:)

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
    CHARACTER(LEN=255) :: cl_longname
    !!----------------------------------------------------------------------
    cl_longname="Steric SSH anomaly"
    IF (lhalo ) THEN 
         cv_out2d='sosshhast'
         cv_out='vosshhast'
         cl_longname="Halosteric SSH anomaly"
    ENDIF
    IF (lther ) THEN 
         cv_out2d='sosshthst'
         cv_out='vosshthst'
         cl_longname="Thermosteric SSH anomaly"
    ENDIF
    IF ( limit ) THEN
       ipk(:)                       = 1
       stypvar(1)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
       stypvar(1)%cname             = cv_out2d
       stypvar(1)%cunits            = 'm'
       stypvar(1)%rmissing_value    = 0.
       stypvar(1)%valid_min         = -100.
       stypvar(1)%valid_max         = 100.
       stypvar(1)%clong_name        = TRIM(cl_longname)
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
       stypvar(1)%clong_name        =  TRIM(cl_longname)
       stypvar(1)%cshort_name       = cv_out
       stypvar(1)%conline_operation = 'N/A'
       stypvar(1)%caxis             = 'TZYX'

       ! create output fileset (3D case )
       ncout = create      (cf_out, cf_tfil,  npiglo, npjglo, npk     , ld_nc4=lnc4 )
       ierr  = createvar   (ncout,  stypvar, 1,      ipk,    id_varout, ld_nc4=lnc4 )
       ierr  = putheadervar(ncout,  cf_tfil,  npiglo, npjglo, npk       )
    ENDIF

    ! read time in the model file, not the reference file !
    IF ( lhalo ) THEN
      dtim  = getvar1d(cf_sfil, cn_vtimec, npt)
    ELSE
      dtim  = getvar1d(cf_tfil, cn_vtimec, npt)
    ENDIF
    ierr  = putvar1d(ncout, dtim,   npt, 'T')

  END SUBROUTINE CreateOutput

END PROGRAM cdfsteric


