PROGRAM cdfdiv
  !!======================================================================
  !!                     ***  PROGRAM  cdfdiv  ***
  !!=====================================================================
  !!  ** Purpose : Compute the divergence for given gridU gridV files
  !!               and variables
  !!
  !!  ** Method  : Use the same stencil than in NEMO code for computing
  !!               vertical velocities
  !!
  !! History :  3.0  : 10/2011  : P. Mathiot : first version, based on cdfw.f90
  !!                 : 05/2015  : J.M. Molines : revised version similar to cdfcurl
  !!         :  4.0  : 03/2017  : J.M. Molines  
  !!
  !!----------------------------------------------------------------------

  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class derived_fields
  !!----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER(KIND=4)                           :: ji, jj, jk, jt     ! dummy loop index
  INTEGER(KIND=4)                           :: it                 ! time index for vvl
  INTEGER(KIND=4)                           :: npiglo, npjglo     ! size of the domain
  INTEGER(KIND=4)                           :: npk, npt           ! size of the domain
  INTEGER(KIND=4)                           :: nlev               ! number of output levels
  INTEGER(KIND=4)                           :: narg, iargc        ! browse command line
  INTEGER(KIND=4)                           :: ijarg              !
  INTEGER(KIND=4)                           :: ncout, ierr        ! browse command line
  INTEGER(KIND=4), DIMENSION(1)             :: ipk, id_varout     ! output variable properties
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: nilev             ! level to be processed

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e1v, e2u, e1t, e2t ! horizontql metrics
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e3u, e3v, e3t      ! vertical metrics
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: un, vn             ! velocity field
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zun, zvn           ! working arrays
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: tim, zdep, gdep    ! time counter
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: e31d               !  vertical metrics (full step)
  REAL(KIND=4)                              :: zmask              ! mask at T point for -T option

  REAL(KIND=8)                              :: dl_pi, dl_omega    ! 3.14159... and earth rotation rad/sec
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dhdivn             ! divergence
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dl_ff              ! Coriolis parameter at T point

  CHARACTER(LEN=256)                        :: cf_ufil, cf_vfil   ! file names
  CHARACTER(LEN=256)                        :: cf_tfil            ! T-file for vvl e3t metrics
  CHARACTER(LEN=256)                        :: cf_out = 'curl.nc' ! output file name
  CHARACTER(LEN=256)                        :: cv_u, cv_v         ! variable names
  CHARACTER(LEN=256)                        :: cldum              ! dummy string

  TYPE (variable), DIMENSION(1)             :: stypvar            ! structure for attibutes

  LOGICAL                                   :: lforcing = .FALSE. ! forcing flag
  LOGICAL                                   :: lchk     = .FALSE. ! flag for missing files
  LOGICAL                                   :: lperio   = .FALSE. ! flag for E-W periodicity
  LOGICAL                                   :: ldblpr   = .FALSE. ! flag for dble precision output
  LOGICAL                                   :: lsurf    = .FALSE. ! flag for 1 lev on C grid.
  LOGICAL                                   :: loverf   = .FALSE. ! flag for 1 lev on C grid.
  LOGICAL                                   :: lfull    = .FALSE. ! full step flag
  LOGICAL                                   :: lnc4     = .FALSE. ! Use nc4 with chunking and deflation

  !!----------------------------------------------------------------------
  CALL ReadCdfNames() 
  narg = iargc()
  IF ( narg < 8 ) THEN
     PRINT *,' usage : cdfdiv -u U-file U-var -v V-file V-var -l levlist  [-8]...'
     PRINT *,'          ... [-surf] [-overf] [-full] [-o OUT-file ] [-nc4] '
     PRINT *,'          ... [-vvl T-file]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the divergence of the flow from the U and V velocity components'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -u U-file U-var : file and variable name for zonal component'
     PRINT *,'       -v V-file V-var : file and variable name for meridional component'
     PRINT *,'       -l levlist    : levels to be processed. If set to 0, assume forcing file'
     PRINT *,'                in input. Example of recognized syntax :'
     PRINT *,'                  -l "1,10,30"  or -l "1-20" or even -l "1-3,10-20,30-"'
     PRINT *,'                  -l  1 . Note that -l "3-" set a levlist from 3 to the bottom'
     PRINT * 
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-8 ]: save in double precision instead of standard simple precision.'
     PRINT *,'       [-surf ] : work with single level C-grid (not forcing)'
     PRINT *,'       [-overf ] : store the ratio curl/f where f is the coriolis parameter'
     PRINT *,'       [-full ] : in case of full step configuration. Default is partial step.'
     PRINT *,'       [-o OUT-file] : specify output file name instead of ',TRIM(cf_out) 
     PRINT *,'       [-nc4 ]     : Use netcdf4 output with chunking and deflation level 1'
     PRINT *,'                 This option is effective only if cdftools are compiled with'
     PRINT *,'                 a netcdf library supporting chunking and deflation.'
     PRINT *,'       [-vvl T-file] : use time-varying e3t, specify T-file for e3t.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        ', TRIM(cn_fhgr),' ',TRIM(cn_fzgr)
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) ,' unless -o option is used.'
     PRINT *,'         variables : div units : s^-1'
     PRINT *,'               or divoverf, no units (if -overf option)'
     STOP
  ENDIF

  ijarg=1
  DO WHILE ( ijarg <= narg ) 
     CALL getarg(ijarg, cldum) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ('-u')
        CALL getarg(ijarg, cf_ufil) ; ijarg=ijarg+1
        CALL getarg(ijarg, cv_u   ) ; ijarg=ijarg+1
     CASE ('-v')
        CALL getarg(ijarg, cf_vfil) ; ijarg=ijarg+1
        CALL getarg(ijarg, cv_v   ) ; ijarg=ijarg+1
     CASE ('-l')
        CALL getarg(ijarg, cldum) ; ijarg=ijarg+1 
        CALL ParseLevel(cldum)  ! fills in array nilev(nlev)
     CASE ('-8'    ) ; ldblpr = .TRUE.
     CASE ('-surf' ) ; lsurf  = .TRUE.
     CASE ('-overf') ; loverf = .TRUE.
     CASE ('-full' ) ; lfull  = .TRUE.
     CASE ('-nc4'  ) ; lnc4   = .TRUE.
     CASE ('-o'    ) ; CALL getarg(ijarg, cf_out ) ; ijarg=ijarg+1
     CASE ('-vvl'  ) ; lg_vvl = .TRUE.
        CALL getarg(ijarg, cf_tfil) ; ijarg=ijarg+1 
     CASE DEFAULT
        PRINT *,  TRIM(cldum), ' : unknown option '
     END SELECT
  ENDDO

  lchk = chkfile(cn_fhgr ) .OR. lchk
  lchk = chkfile(cn_fzgr ) .OR. lchk
  lchk = chkfile(cf_ufil ) .OR. lchk
  lchk = chkfile(cf_vfil ) .OR. lchk
  IF ( lchk ) STOP ! missing files

  IF ( lg_vvl ) THEN
     cn_fe3u = cf_ufil
     cn_fe3v = cf_vfil
     cn_fe3t = cf_tfil
  ENDIF

  npiglo = getdim(cf_ufil,cn_x)
  npjglo = getdim(cf_ufil,cn_y)
  npk    = getdim(cf_ufil,cn_z)
  npt    = getdim(cf_ufil,cn_t) 

  PRINT *, 'NPIGLO = ', npiglo
  PRINT *, 'NPJGLO = ', npjglo
  PRINT *, 'NPK    = ', npk
  PRINT *, 'NPT    = ', npt
  PRINT *, 'NLEV   = ', nlev


  !test if lev exists
  IF ( (npk==0) .AND. (nlev > 0) .AND. .NOT. lsurf ) THEN
     PRINT *, 'Problem : npk = 0 and lev > 0 STOP'
     PRINT *, '  Use -surf option is dealing with single level file on C grid '
     STOP
  END IF

  ! case of 1 level on C-grid
  IF ( lsurf ) THEN
     nlev=1
     IF (ALLOCATED (nilev) ) THEN
        DEALLOCATE(nilev) 
     ENDIF
     ALLOCATE(nilev(nlev) )
     npk = 1 ; nilev(1) =1 
  ENDIF

  ! if forcing field 
  IF ( nilev(1) == 0 .AND. npk==0 ) THEN
     lforcing=.TRUE.
     npk = 1 ; nilev(1)=1
     PRINT *, 'npk =0, assume 1'
  END IF

  IF ( npt==0 ) THEN
     PRINT *, 'npt=0, assume 1'
     npt=1
  END IF
  ! 
  DO jk = 1, nlev
     IF (nilev(jk) >= npk ) THEN
        nlev=jk
        EXIT
     ENDIF
  ENDDO
  PRINT *, 'NLEV', nlev

  ! Allocate the memory
  ALLOCATE ( e1v(npiglo,npjglo), e2u(npiglo,npjglo) )
  ALLOCATE ( e1t(npiglo,npjglo), e2t(npiglo,npjglo) )
  ALLOCATE ( e3u(npiglo,npjglo), e3v(npiglo,npjglo), e3t(npiglo,npjglo) )

  ALLOCATE ( un(npiglo,npjglo)  , vn(npiglo,npjglo)  )
  ALLOCATE ( zun(npiglo,npjglo) , zvn(npiglo,npjglo) )
  ALLOCATE ( dhdivn(npiglo,npjglo) )
  ALLOCATE ( tim(npt) )
  ALLOCATE ( gdep(nlev) , zdep(npk))
  IF ( lfull ) ALLOCATE ( e31d (npk) )

  e1v =  getvar(cn_fhgr, cn_ve1v, 1, npiglo, npjglo)
  e1t =  getvar(cn_fhgr, cn_ve1t, 1, npiglo, npjglo)
  e2u =  getvar(cn_fhgr, cn_ve2u, 1, npiglo, npjglo)
  e2t =  getvar(cn_fhgr, cn_ve2t, 1, npiglo, npjglo)

  ! use zun and zvn to store f latitude and longitude for output
  zun = getvar(cn_fhgr, cn_glamt, 1, npiglo, npjglo)
  zvn = getvar(cn_fhgr, cn_gphit, 1, npiglo, npjglo)

  IF ( loverf ) THEN
     ALLOCATE (dl_ff(npiglo,npjglo) )
     dl_pi = acos(-1.d0)
     dl_omega = 2* dl_pi/86400.d0
     dl_ff = 2* dl_omega* sin ( zvn*dl_pi/180.d0 ) 
  ENDIF

  ! fills in gdep
  IF ( lforcing .OR. lsurf ) THEN
     gdep(1)=0.
  ELSE
     zdep(:) = getvar1d(cf_ufil, cn_vdeptht, npk )
     DO jk=1,nlev
        gdep(jk) = zdep( nilev(jk) )
     ENDDO
  ENDIF
  IF ( lfull ) e31d(:) = getvare3(cn_fzgr, cn_ve3t, npk)

  ! look for  E-W periodicity
  IF ( zun(1,1) == zun(npiglo-1,1) ) lperio = .TRUE.

  CALL CreateOutput


  DO jt=1,npt
     IF (MOD(jt,100)==0 ) PRINT *, jt,'/',npt

     IF ( lg_vvl ) THEN  ; it=jt
     ELSE                ; it=1
     ENDIF
     DO jk = 1, nlev

        zun(:,:) =  getvar(cf_ufil, cv_u, nilev(jk) ,npiglo,npjglo, ktime=jt)
        zvn(:,:) =  getvar(cf_vfil, cv_v, nilev(jk) ,npiglo,npjglo, ktime=jt)

        IF ( lfull ) THEN
           e3u(:,:) = e31d(jk)
           e3v(:,:) = e31d(jk)
           e3t(:,:) = e31d(jk)
        ELSE
           ! e3 metrics at level jk ( Partial steps)
           e3u(:,:) = getvar(cn_fe3u, cn_ve3u, jk, npiglo, npjglo, ktime=it, ldiom=.NOT.lg_vvl )
           e3v(:,:) = getvar(cn_fe3v, cn_ve3v, jk, npiglo, npjglo, ktime=it, ldiom=.NOT.lg_vvl )
           e3t(:,:) = getvar(cn_fe3t, cn_ve3t, jk, npiglo, npjglo, ktime=it, ldiom=.NOT.lg_vvl )
        ENDIF

        IF ( lforcing ) THEN ! for forcing file u and v are on the A grid
           DO ji=1, npiglo-1
              un(ji,:) = 0.5*(zun(ji,:) + zun(ji+1,:))
           END DO
           !
           DO jj=1, npjglo-1
              vn(:,jj) = 0.5*(zvn(:,jj) + zvn(:,jj+1))
           END DO
           ! end compute u and v on U and V point
        ELSE
           un(:,:) = zun(:,:)
           vn(:,:) = zvn(:,:)
        END IF

        dhdivn(:,:) = 0.d0
        DO jj = 2, npjglo -1
           DO ji = 2, npiglo -1
              dhdivn(ji,jj) =   &
                   &  (  e2u(ji,  jj  )*e3u(ji  ,jj  ) * un(ji  ,jj  )*1.d0 - &
                   &     e2u(ji-1,jj  )*e3u(ji-1,jj  ) * un(ji-1,jj  )*1.d0   &
                   &   + e1v(ji  ,jj  )*e3v(ji  ,jj  ) * vn(ji  ,jj  )*1.d0 - &
                   &     e1v(ji  ,jj-1)*e3v(ji  ,jj-1) * vn(ji  ,jj-1)*1.d0  )&
                   & / ( e1t(ji,jj)*e2t(ji,jj) * e3t(ji,jj) )
           END DO
        END DO

        IF ( lperio ) dhdivn(npiglo,:) = dhdivn(2, :)

        ! write dhdivn on file at level k and at time jt
        IF ( loverf ) THEN
           WHERE (dl_ff /= 0.d0) dhdivn=dhdivn/dl_ff
        ENDIF

        ierr = putvar(ncout, id_varout(1), dhdivn, nilev(jk), npiglo, npjglo, ktime=jt)
     ENDDO
  END DO
  ierr = closeout(ncout)

CONTAINS

  SUBROUTINE ParseLevel( cdum )
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE ParsLevel  ***
    !!
    !! ** Purpose :  Parse a string representing a level list with no separator
    !!               or , or - as separator 
    !!
    !! ** Method  :   1,3,5 => list = (1,3,5)
    !!                1-4   => liste = ( 1,2,3,4) 
    !!                1-4,6-8   => liste = ( 1,2,3,4,6,7,8) 
    !!                10-      => liste = (10,11,12,...bottom)
    !!                Allocate and fill nilev array with respective levels
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT (in ) :: cdum

    INTEGER(KIND=4)                 :: jc, jk
    INTEGER(KIND=4)                 :: ilength, ik1, ik2
    INTEGER(KIND=4)                 :: icomma, idash, ipos, ipos1
    INTEGER(KIND=4), DIMENSION(350) :: ilev
    CHARACTER(LEN=256) :: cldum
    CHARACTER(LEN=80), DIMENSION(:), ALLOCATABLE :: clblk
    !!----------------------------------------------------------------------
    ilength=LEN(TRIM(cdum) )

    ! look for , and -
    icomma = 0 ; idash = 0
    DO jc=1,ilength
       IF ( cdum(jc:jc) == ',' ) THEN
          icomma = icomma + 1
       ELSE IF ( cdum(jc:jc) == '-' ) THEN
          idash = idash + 1
       ENDIF
    END DO

    ! no dash nor comma
    IF (icomma == 0 .AND. idash == 0 ) THEN
       nlev=1
       ALLOCATE ( nilev(nlev) )
       READ(cdum,*) nilev(1)
       RETURN
    ENDIF
    ! look for the number of blocks (between commas)
    cldum=cdum
    ALLOCATE (clblk(icomma+1))

    ipos1=1
    DO jc=1,icomma
       ipos=INDEX(cldum ,",")
       clblk(jc)=cldum(1:ipos-1)
       ipos1=ipos+1
       cldum=cldum(ipos1:)
    ENDDO
    clblk(icomma+1)=cldum

    ! now parse block
    nlev=0
    DO jc=1,icomma+1
       ipos=INDEX(clblk(jc),"-")
       IF ( ipos == 0 ) THEN
          nlev=nlev+1
          READ(clblk(jc),* ) ilev(nlev) 
       ELSE IF ( ipos == LEN(TRIM((clblk(jc)))) ) THEN
          READ(clblk(jc)(1:ipos-1),*) ik1
          ik2=300  ! use 300 as a maximum mean while we wait for npk
          PRINT *,' BINGO !', ik1, ik2
          DO jk=ik1,ik2
             nlev=nlev+1
             ilev(nlev)=jk
          ENDDO
       ELSE
          READ(clblk(jc)(1:ipos-1),*) ik1
          READ(clblk(jc)(ipos+1: ),*) ik2
          DO jk=ik1,ik2
             nlev=nlev+1
             ilev(nlev)=jk
          ENDDO
       ENDIF
    ENDDO
    ALLOCATE (nilev(nlev) )
    nilev(:)=ilev(1:nlev)

    DEALLOCATE( clblk)


  END SUBROUTINE ParseLevel

  SUBROUTINE CreateOutput
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutput  ***
    !!
    !! ** Purpose :  Create netcdf output file(s) 
    !!
    !! ** Method  :  Use stypvar global description of variables
    !!
    !!----------------------------------------------------------------------
    ! define new variables for output
    stypvar(1)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(1)%cname             = 'div'
    stypvar(1)%cunits            = 's-1'
    IF (loverf)  stypvar(1)%cname             = 'divoverf'
    IF (loverf)  stypvar(1)%cunits            = '-'

    stypvar(1)%cprecision        ='r4'
    IF ( ldblpr )  stypvar(1)%cprecision     ='r8'
    stypvar(1)%rmissing_value    = 0.
    stypvar(1)%valid_min         = -1000.
    stypvar(1)%valid_max         =  1000.
    stypvar(1)%clong_name        = 'Divergence '
    stypvar(1)%cshort_name       = 'div'
    IF (loverf ) stypvar(1)%clong_name        = 'Divergence normalized by f'
    IF (loverf ) stypvar(1)%cshort_name       = 'divoverf'
    stypvar(1)%conline_operation = 'N/A'
    stypvar(1)%caxis             = 'TYX'

    ipk(1) = nlev  !  nlevel so far

    ! create output fileset
    ncout = create      (cf_out, cf_ufil, npiglo, npjglo, nlev, cdep=cn_vdeptht, ld_nc4=lnc4 )
    ierr  = createvar   (ncout ,   stypvar, 1 ,   ipk,    id_varout,             ld_nc4=lnc4 )
    ierr  = putheadervar(ncout,  'dummy', npiglo, npjglo, nlev, pnavlon=zun, pnavlat=zvn, pdep=gdep)

    tim  = getvar1d(cf_ufil, cn_vtimec, npt      )
    ierr = putvar1d(ncout,   tim,       npt,  'T')

  END SUBROUTINE CreateOutput

END PROGRAM cdfdiv


