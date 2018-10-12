PROGRAM cdfdynh_anom
  !!======================================================================
  !!                     ***  PROGRAM  cdfdynh_anom  ***
  !!=====================================================================
  !!  ** Purpose : Compute dynamical height anomaly field from gridT file
  !!               at each levels.
  !!               Store the results on a 3D cdf file.
  !!
  !!  ** Method  : the integral of (1/g) *10e4 * sum [ delta * dz ]
  !!               with delta = (1/rho - 1/rho0)
  !!               10e4 factor is conversion decibar/pascal
  !!
  !! History : 2.1  : 05/2010  : R. Dussin    : Original code
  !!           3.0  : 12/2010  : J.M. Molines : Doctor norm + Lic.
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  USE cdfio
  USE eos, ONLY : sigmai
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
  REAL(KIND=8)                              :: dgrav = 9.81d0   ! gravity
  REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: dtim             ! time counter
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dhdy, dterm      ! dynamic height, working array
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dsig0, dsig      ! In situ density (reference, local)

  CHARACTER(LEN=256)                        :: cf_tfil           ! input file name
  CHARACTER(LEN=256)                        :: cf_sfil           ! input salinity file (option)
  CHARACTER(LEN=255)                        :: cf_sshfil         ! input SSH file (option)
  CHARACTER(LEN=256)                        :: cf_out='cdfhdy3d.nc' ! output file name
  CHARACTER(LEN=256)                        :: cv_out='vohdy'   ! output file name
  CHARACTER(LEN=256)                        :: cf_out2d='cdfhdy2d.nc' ! output file name
  CHARACTER(LEN=256)                        :: cv_out2d='sohdy'   ! output file name
  CHARACTER(LEN=256)                        :: cldum            ! working char variable

  TYPE(variable) , DIMENSION(1)             :: stypvar          ! structure for attributes

  LOGICAL                                   :: lout  = .FALSE.   ! flag to use non standard output name
  LOGICAL                                   :: lnc4  = .FALSE.   ! Use nc4 with chunking and deflation
  LOGICAL                                   :: limit = .FALSE.   ! flag set if limit are used (2D case)
  LOGICAL                                   :: lchk  = .FALSE.   ! flag for file existence checking
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfdynh_anom -t T-file [-s S-file] [-limit lev1 lev2] [-vvl] ...'
     PRINT *,'             ... [--ssh-file SSH-file] [-o OUT-file] [-nc4]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'        Compute dynamic height anomaly from T-file given as argument.'
     PRINT *,'        In this tool, the cumulated values (from top to bottom) are saved at'
     PRINT *,'        each model level. '
     PRINT *,'        If the -limit option is used, only the 2D integral of the dynamic '
     PRINT *,'        height anomaly between lev1 and lev2 is saved. '
     PRINT *,'        This program replace cdfhdy ( case using -limit option) and cdfhdy3d '
     PRINT *,'        in the standard case.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -t T-file : netcdf file with temperature and salinity.' 
     PRINT *,'             If salinity is not in T-file, use -s option.'
     PRINT *,'             If ssh not in T-file use --ssh-file option.'
     PRINT *,'     '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-s S-file ] : Specify salinity file if not T-file.' 
     PRINT *,'       [--ssh-file SSH-file] : specify the ssh file if not in T-file.'
     PRINT *,'       [-limit lev1 lev2] : if specified, the program will only output the'
     PRINT *,'              dynamic height anomaly at lev1 with reference at lev2.'
     PRINT *,'       [-vvl] : use time-varying vertical metrics.'
     PRINT *,'       [-o OUT-file] : Specify output file name instead of ',TRIM(cf_out)
     PRINT *,'       [-nc4] : Use netcdf4 output with chunking and deflation level 1.'
     PRINT *,'                This option is effective only if cdftools are compiled with'
     PRINT *,'                a netcdf library supporting chunking and deflation.'
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

  cf_sfil   = 'none'
  cf_sshfil = 'none'
  ijarg = 1
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-t','-f'  ) ; CALL getarg (ijarg, cf_tfil) ; ijarg=ijarg+1
        ! options
     CASE ( '-s'       ) ; CALL getarg (ijarg, cf_sfil) ; ijarg=ijarg+1
     CASE ('--ssh-file') ; CALL getarg(ijarg, cf_sshfil); ijarg=ijarg+1
     CASE ( '-o'       ) ; CALL getarg (ijarg, cf_out ) ; ijarg=ijarg+1 
        ;                  lout  = .TRUE.
     CASE ( '-nc4'     ) ; lnc4  = .TRUE.
     CASE ( '-vvl'     ) ; lg_vvl= .TRUE.
     CASE ( '-limit'   ) ; limit = .TRUE.
        ;                  CALL getarg (ijarg, cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) nlev1
        ;                  CALL getarg (ijarg, cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) nlev2
     CASE DEFAULT        ; PRINT *,' ERROR : ',TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  ENDDO

  IF ( cf_sfil   == 'none' ) cf_sfil   = cf_tfil
  IF ( cf_sshfil == 'none' ) cf_sshfil = cf_tfil
  IF ( lout ) cf_out2d = cf_out  !  name of 2D file if -limit is used

  lchk = chkfile(cf_tfil)
  lchk = lchk .OR. chkfile(cf_sfil  )
  lchk = lchk .OR. chkfile(cf_sshfil)
  lchk = lchk .OR. chkfile(cn_fmsk  )
  lchk = lchk .OR. chkfile(cn_fzgr  )

  IF ( lchk ) STOP 99 ! missing files

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
  temp0(:,:) =  0.
  zsal0(:,:)  = 35.

  zssh(:,:)  = getvar(cf_sshfil,  cn_sossheig, 1,  npiglo, npjglo)
  e3t_1d(:)  = getvare3(cn_fzgr, cn_ve3t1d, npk)

  DO jt = 1, npt
     IF ( lg_vvl ) THEN ; it=jt
     ELSE               ; it=1
     ENDIF

     PRINT *,' TIME = ', jt, dtim(jt)/86400.,' days'
     dhdy(:,:)   = 0.
     rdepth(:,:) = 0.

     DO jk = nlev1, nlev2

        IF ( lg_vvl ) THEN ; rdep(:,:) =  getvar(cn_fe3t, cn_ve3t, jk,npiglo,npjglo, ktime=it)
        ELSE               ; rdep(:,:)  = e3t_1d(jk)
        ENDIF

        tmask(:,:) = getvar(cn_fmsk, cn_tmask, jk, npiglo, npjglo)

        IF ( jk == 1 .AND. .NOT. lg_vvl ) THEN
           rdep(:,:) = rdep(:,:) + zssh(:,:)
        ENDIF

        ! depth at current level, including ssh (used for computation of rho in situ)
        rdepth(:,:) = rdepth(:,:) + rdep(:,:)

        temp(:,:) = getvar(cf_tfil, cn_votemper,  jk ,npiglo, npjglo, ktime=jt)
        zsal(:,:) = getvar(cf_sfil, cn_vosaline,  jk ,npiglo, npjglo, ktime=jt)

        dsig0 = sigmai(temp0, zsal0, rdepth, npiglo, npjglo) 
        dsig  = sigmai(temp , zsal , rdepth, npiglo, npjglo) 

        ! we compute the term of the integral : (1/g) *10e4 * sum [ delta * dz ]
        ! with delta = (1/rho - 1/rho0)
        ! 10e4 factor is conversion decibar/pascal
        !
        dterm = ( ( 1.d0 / ( drau0 + dsig(:,:) ) ) - ( 1.d0 / ( drau0 + dsig0(:,:) ) ) ) * 10000.d0 * rdep / dgrav
        ! in land, it seems appropriate to stop the computation
        WHERE(zsal == zsps ) dterm = 0

        dhdy(:,:) = dhdy(:,:) + dterm(:,:)
        ! masked
        dhdy(:,:) = dhdy(:,:) * tmask(:,:)

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

END PROGRAM cdfdynh_anom
