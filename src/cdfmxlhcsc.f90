PROGRAM cdfmxlhcsc
  !!======================================================================
  !!                     ***  PROGRAM  cdfmxlhcsc  ***
  !!=====================================================================
  !!  ** Purpose : Compute mixed layer depth and the heat and salt contents
  !!               in the mixed layer. There is an option to limit this 
  !!               computation between hmin and ml depth. For that, hmin is
  !!               given as last argument (>0) with no arguments, hmin is 
  !!               supposed to be 0.
  !!
  !!  ** Method  : This program is a merge of  cdfmxl, cdfmxlheatc and
  !!               cdfmxlsaltc.
  !!               MXL computation:
  !!                - compute surface properties
  !!                - initialize depths and model levels number
  !!                - from bottom to top compute rho and
  !!                check if rho > rho_surf +rho_c, where rho_c is a 
  !!                density criteria given as argument
  !!               Heat Content and Salt Content:
  !!                 HC = sum ( rho cp T * e1 * e2 * e3 * tmask )
  !!                 SC = sum ( rho    S * e1 * e2 * e3 * tmask )
  !!                 where the sum is limited to the MXL, between hmin and
  !!                 MLD
  !!
  !! History : 2.1  : 04/2007  : M. Juza      : Merging of the programs
  !!           3.0  : 01/2011  : J.M. Molines : Doctor norm + Lic.
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  USE cdfio
  USE eos
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class mixed_layer
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                              :: ji, jj, jk, jt  ! dummy loop index
  INTEGER(KIND=4)                              :: it              ! time index for vvl
  INTEGER(KIND=4)                              :: ik              ! level indirect index
  INTEGER(KIND=4)                              :: narg, iargc     ! browse line
  INTEGER(KIND=4)                              :: ijarg           ! argument counter
  INTEGER(KIND=4)                              :: npiglo, npjglo  ! domain size
  INTEGER(KIND=4)                              :: npk, npt        ! domaine size
  INTEGER(KIND=4)                              :: ncout, ierr     ! ncid of output file an error status
  INTEGER(KIND=4), DIMENSION(3)                :: ipk, id_varout  ! levels and varid's of output vars
  INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: mbathy          ! number of w levels in water <= npk
  INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: nmln            ! last level where rho > rho + val_crit
  !              or rtem > rtem + val_crit
  REAL(KIND=4), PARAMETER                      :: rprho0=1020.    ! reference density
  REAL(KIND=4), PARAMETER                      :: rpcp=4000.      ! specific heat of water
  REAL(KIND=4)                                 :: val             ! criteria value
  REAL(KIND=4)                                 :: hmin = 0.       ! minimum depth for vertical integration
  REAL(KIND=4), DIMENSION(1)                   :: rdep            ! dummy depth output
  REAL(KIND=4), DIMENSION(:),      ALLOCATABLE :: gdepw           ! depth of w points
  REAL(KIND=4), DIMENSION(:),      ALLOCATABLE :: e31d            ! vertical metrics (full step)
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: rtem            ! temperature
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: rsal            ! salinity
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: rho             ! density (sigma-0)
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: rho_surf        ! surface density
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: tem_surf        ! surface temperature
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: hmld            ! mixed layer depth
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: tmask_surf      ! surface tmask
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: tmask           ! land sea mask of temperature
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: e3              ! vertical metrics

  REAL(KIND=8), DIMENSION(:),      ALLOCATABLE :: dtim            ! time counter
  REAL(KIND=8), DIMENSION(:,:),    ALLOCATABLE :: dmxlheatc       ! mxl heat content
  REAL(KIND=8), DIMENSION(:,:),    ALLOCATABLE :: dmxlsaltc       ! mxl salt content

  TYPE(variable), DIMENSION(3)                 :: stypvar         ! output attributes

  CHARACTER(LEN=256)                           :: cf_tfil         ! input file
  CHARACTER(LEN=256)                           :: cf_mld          ! mld input file (option)
  CHARACTER(LEN=256)                           :: cf_out='mxlhcsc.nc' ! output file
  CHARACTER(LEN=256)                           :: cv_mld          ! mld imput variable (option)
  CHARACTER(LEN=256)                           :: criteria        ! type of criteria used for mld
  CHARACTER(LEN=256)                           :: cldum           ! dummy string

  LOGICAL                                      :: lchk            ! flag for missing files
  LOGICAL                                      :: lfull=.FALSE.   ! flag for full step
  LOGICAL                                      :: lnc4 =.FALSE.   ! Use nc4 with chunking and deflation
  LOGICAL                                      :: lmld =.TRUE.    ! Flag to compute mld
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfmxlhcsc -f T-file -C criteria -t THRESH-value [-hmin hmin] ...'
     PRINT *,'             ... [-o OUT-file] [-nc4] [-vvl] [-mld MLD-file MLD-var]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the MiXed Layer depth, the Heat Content and Salt Content.' 
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f T-file : netcdf input file for temperature and salinity (gridT).' 
     PRINT *,'       -C criteria : Specify the type of criteria to use for MXL computation'
     PRINT *,'                can be one of temperature, t,  T for temperature criteria.'
     PRINT *,'                    or one of density, d,  D  for density criteria.'
     PRINT *,'       -t THRESH-value : threshold value for the criteria.'
     PRINT *,'                (eg: 0.2 for temp, 0.01 or 0.03 for dens)'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-hmin hmin ] : limit the vertical integral from hmin to mld. By default,' 
     PRINT *,'                 hmin is set to 0 so that the integral is performed on the'
     PRINT *,'                 whole mixed layer.'
     PRINT *,'       [-o OUT-file ] : specify output file name instead of ',TRIM(cf_out)
     PRINT *,'       [ -nc4 ]     : Use netcdf4 output with chunking and deflation level 1.'
     PRINT *,'                 This option is effective only if cdftools are compiled with'
     PRINT *,'                 a netcdf library supporting chunking and deflation.'
     PRINT *,'       [-vvl ] : use time-varying vertical metrics'
     PRINT *,'       [-mld MLD-file MLD-var]: This option replaces both ''-C'' and ''-t'' '
     PRINT *,'                 arguments, indicating that the mixed layer depth is to be read'
     PRINT *,'                 as the variable MLD-var in the MLD-file. With this option, the'
     PRINT *,'                 heat and salt content can be computed from the MLD estimates'
     PRINT *,'                 provided by ''cdfmxl''.' 
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       ',TRIM(cn_fhgr),' ',TRIM(cn_fzgr),' and ',TRIM(cn_fmsk) 
     PRINT *,'        In case of full step configuration, ',TRIM(cn_fbathylev),' is also required.'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) ,' unless -o option is used'
     PRINT *,'         variables : -  somxl010 (mld based on density criterion 0.01)'
     PRINT *,'          (2D)          or somxl030 (mld on density criterion 0.03)'
     PRINT *,'                        or somxlt02 (mld on temperature criterion -0.2)'
     PRINT *,'                        -  somxlheatc (heat content computed in the MLD)'
     PRINT *,'                        -  somxlsaltc (salt content computed in the MLD)'
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'       cdfmxl, cdfmxlheatc and  cdfmxlsaltc.'
     PRINT *,'      '
     STOP 
  ENDIF


  ijarg = 1 
  DO WHILE ( ijarg <= narg )
     CALL getarg (ijarg, cldum ) ; ijarg = ijarg + 1
     SELECT CASE ( cldum )
     CASE ( '-f '    ) ; CALL getarg (ijarg, cf_tfil  ) ; ijarg = ijarg + 1 
     CASE ( '-C '    ) ; CALL getarg (ijarg, criteria ) ; ijarg = ijarg + 1 
     CASE ( '-t '    ) ; CALL getarg (ijarg, cldum    ) ; ijarg = ijarg + 1  ; READ(cldum,*) val
        ! options
     CASE ( '-hmin ' ) ; CALL getarg (ijarg, cldum    ) ; ijarg = ijarg + 1  ; READ(cldum,*) hmin
     CASE ( '-o'     ) ; CALL getarg (ijarg, cf_out   ) ; ijarg = ijarg + 1  
     CASE ( '-nc4'   ) ; lnc4   = .TRUE.
     CASE ( '-vvl'   ) ; lg_vvl = .TRUE.
     CASE ( '-mld'   ) ; CALL getarg (ijarg, cf_mld   ) ; ijarg = ijarg + 1
        ;                CALL getarg (ijarg, cv_mld   ) ; ijarg = ijarg + 1 
        ;                lmld = .FALSE.
     CASE DEFAULT  ; PRINT *, ' ERROR : Option ',TRIM(cldum),' unknown !' ; STOP 99
     END SELECT
  ENDDO

  lchk = chkfile (cn_fhgr)
  lchk = chkfile (cn_fzgr) .OR. lchk
  lchk = chkfile (cn_fmsk) .OR. lchk
  lchk = chkfile (cf_tfil) .OR. lchk
  IF ( .NOT. lmld ) lchk = chkfile (cf_mld) .OR. lchk
  IF ( lchk ) STOP 99 ! missing files
  IF ( lg_vvl ) THEN
     cn_fe3t = cf_tfil 
     cn_ve3t = cn_ve3tvvl
  ENDIF

  ! read dimensions 
  npiglo = getdim (cf_tfil, cn_x)
  npjglo = getdim (cf_tfil, cn_y)
  npk    = getdim (cf_tfil, cn_z)
  npt    = getdim (cf_tfil, cn_t)

  IF ( lmld ) THEN
     ! Variable Mixed Layer Depth
     SELECT CASE ( criteria)
     CASE ( 'Temperature', 'temperature', 't', 'T' ) ; WRITE(cv_mld,'(a,i2.2)' ) 'somxlt', INT(ABS(val)*10)
     CASE ( 'Density',     'density',     'd', 'D' ) ; WRITE(cv_mld,'(a,i3.3)' ) 'somxl' , INT((val)*1000)
     CASE DEFAULT  ; PRINT *,'ERROR: ',TRIM(criteria),' : criteria not understood' ; STOP 99
     END SELECT
  ENDIF

  ! Allocate arrays
  ALLOCATE (rtem(npiglo,npjglo),rsal(npiglo,npjglo)          )
  ALLOCATE (tmask(npiglo,npjglo),tmask_surf(npiglo,npjglo)   )
  ALLOCATE (mbathy(npiglo,npjglo)                            )
  ALLOCATE (nmln(npiglo,npjglo),hmld(npiglo,npjglo)          )
  ALLOCATE (dmxlheatc(npiglo,npjglo),dmxlsaltc(npiglo,npjglo))
  ALLOCATE (e3(npiglo,npjglo)                                )
  ALLOCATE (gdepw(0:npk), dtim(npt)                          )

  CALL CreateOutput

  ! read mbathy and gdepw use real rtem(:,:) as template (getvar is used for real only)
  INQUIRE (FILE=cn_fbathylev, EXIST=lfull)
  IF ( lfull ) THEN ; rtem(:,:) = getvar(cn_fbathylev, cn_bathylev, 1, npiglo, npjglo) ; ALLOCATE ( e31d(npk) )
  ELSE              ; rtem(:,:) = getvar(cn_fzgr,      cn_mbathy,   1, npiglo, npjglo)
  ENDIF

  mbathy(:,:)  = rtem(:,:)
  gdepw(0)     = 999999.  ! dummy values normaly always masked
  gdepw(1:npk) = getvare3(cn_fzgr, cn_gdepw, npk) 
  IF ( lfull ) e31d = getvare3(cn_fzgr, cn_ve3t1d, npk )

  ! read surface tmask
  tmask_surf(:,:) = getvar(cn_fmsk, cn_tmask, 1, npiglo, npjglo)

  DO jt = 1, npt   ! major time loop
     IF (lg_vvl ) THEN ; it=jt
     ELSE              ; it=1
     ENDIF
     IF ( lmld ) THEN
        ! MXL computation
        !---------------
        ! Initialization to the number of w ocean point mbathy
        nmln(:,:) = mbathy(:,:)


        SELECT CASE ( criteria )
        CASE ( 'temperature', 'Temperature', 'T', 't' ) ! Temperature criteria
           ! temp_surf
           IF (.NOT. ALLOCATED ( tem_surf) ) ALLOCATE (tem_surf(npiglo,npjglo))
           tem_surf(:,:) = getvar(cf_tfil, cn_votemper, 1, npiglo, npjglo, ktime=jt )

           ! Last w-level at which ABS(rtem-tem_surf)>=ABS(val) (starting from jpk-1)
           ! (rtem defined at t-point, thus jk-1 for w-level just above)
           DO jk = npk-1, 2, -1
              rtem(:,:) = getvar(cf_tfil, cn_votemper, jk, npiglo, npjglo, ktime=jt)
              WHERE ( ABS(rtem - tem_surf) > ABS(val) ) nmln = jk
           ENDDO
        CASE ( 'density', 'Density', 'D', 'd' ) ! Density criteria
           ! compute rho_surf
           IF ( .NOT. ALLOCATED( rho_surf ) ) ALLOCATE (rho_surf(npiglo,npjglo) )
           IF ( .NOT. ALLOCATED( rho      ) ) ALLOCATE (rho     (npiglo,npjglo) )
           rtem(:,:) = getvar(cf_tfil, cn_votemper, 1, npiglo, npjglo, ktime=jt)
           rsal(:,:) = getvar(cf_tfil, cn_vosaline, 1, npiglo, npjglo, ktime=jt)
           rho_surf(:,:) = sigma0 (rtem, rsal, npiglo, npjglo )* tmask_surf(:,:)

           ! Last w-level at which rhop>=rho surf+rho_c (starting from jpk-1)
           ! (rhop defined at t-point, thus jk-1 for w-level just above)
           DO jk = npk-1, 2, -1
              rtem( :,:) = getvar(cf_tfil, cn_votemper,  jk ,npiglo, npjglo, ktime=jt)
              rsal( :,:) = getvar(cf_tfil, cn_vosaline,  jk ,npiglo, npjglo, ktime=jt)
              tmask(:,:) = getvar(cn_fmsk, cn_tmask,     jk, npiglo, npjglo          )
              rho(  :,:) = sigma0 (rtem, rsal, npiglo, npjglo ) * tmask(:,:)
              WHERE ( rho > rho_surf + val ) nmln = jk
           ENDDO
        END SELECT

        !! Determine mixed layer depth
        DO jj = 1, npjglo
           DO ji = 1, npiglo
              ik = nmln(ji,jj)
              hmld (ji,jj) = gdepw(ik) * tmask_surf(ji,jj)
           ENDDO
        ENDDO
     ELSE  ! read MLD in file
        hmld(:,:) = getvar(cf_mld, cv_mld, jk, npiglo,npjglo, ktime=jt) * tmask_surf(:,:)
     ENDIF

     !!Compute heat and salt contents in the mixed layer depth
     !!-------------------------------------------------------
     !!
     dmxlheatc(:,:) = 0.d0
     dmxlsaltc(:,:) = 0.d0

     DO jk = 1,npk
        ! Get temperature and salinity at jk
        rtem(:,:)  = getvar(cf_tfil,  cn_votemper, jk, npiglo, npjglo, ktime=jt)
        rsal(:,:)  = getvar(cf_tfil,  cn_vosaline, jk, npiglo, npjglo, ktime=jt)
        tmask(:,:) = getvar(cn_fmsk,  cn_tmask ,   jk, npiglo, npjglo          )

        IF ( lfull ) THEN
           e3(:,:) = e31d(jk)
        ELSE
           ! Get e3 at level jk (ps...)
           e3(:,:) = getvar(cn_fe3t, cn_ve3t, jk ,npiglo, npjglo, ktime=it, ldiom=.NOT.lg_vvl )
        ENDIF

        ! e3 is used as a flag for the mixed layer; it is 0 outside the mixed layer
        e3(:,:) = MAX(0., MIN(e3, hmld-gdepw(jk) ) + MIN(e3, gdepw(jk)+ e3-hmin) - e3)

        ! Heat and salt contents
        dmxlheatc(:,:) = dmxlheatc(:,:) + rtem(:,:) * e3(:,:) * tmask(:,:) *1.d0
        dmxlsaltc(:,:) = dmxlsaltc(:,:) + rsal(:,:) * e3(:,:) * tmask(:,:) *1.d0

     END DO

     !! Heat and salt contents (10^9.J/m2 and 10^6.kg/m2)
     dmxlheatc(:,:) = dmxlheatc(:,:) *rprho0 *rpcp * 1.d-9
     dmxlsaltc(:,:) = dmxlsaltc(:,:) *rprho0       * 1.d-6

     ierr = putvar(ncout, id_varout(1), hmld,            1, npiglo, npjglo, ktime=jt)
     ierr = putvar(ncout, id_varout(2), REAL(dmxlheatc), 1, npiglo, npjglo, ktime=jt)
     ierr = putvar(ncout, id_varout(3), REAL(dmxlsaltc), 1, npiglo, npjglo, ktime=jt)

  END DO ! time loop

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

    rdep(1) = 0.
    ipk(:)  = 1 
    stypvar(1)%ichunk      = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(1)%cname       = TRIM(cv_mld)
    stypvar(1)%cshort_name = TRIM(cv_mld)
    stypvar(1)%cunits      = 'm'
    stypvar(1)%clong_name  = 'Mixed Layer Depth'

    ! Variable Heat Content
    stypvar(2)%ichunk      = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(2)%cname       = 'somxlheatc'
    stypvar(2)%cunits      = '10^9 J/m2'
    stypvar(2)%clong_name  = 'Mixed_Layer_Heat_Content'
    stypvar(2)%cshort_name = 'somxlheatc'

    ! Variable Salt Content
    stypvar(3)%ichunk      = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(3)%cname       = 'somxlsaltc'
    stypvar(3)%cunits      = '10^6 kg/m2'
    stypvar(3)%clong_name  = 'Mixed_Layer_Salt_Content'
    stypvar(3)%cshort_name = 'somxlsaltc'

    ncout = create      (cf_out, cf_tfil, npiglo, npjglo, 1        , ld_nc4=lnc4 )
    ierr  = createvar   (ncout,  stypvar, 3,      ipk,    id_varout, ld_nc4=lnc4 )
    ierr  = putheadervar(ncout,  cf_tfil, npiglo, npjglo, 1, pdep=rdep           )
    dtim  = getvar1d(cf_tfil, cn_vtimec, npt     )
    ierr  = putvar1d(ncout, dtim,       npt, 'T')

  END SUBROUTINE CreateOutput


END PROGRAM cdfmxlhcsc
