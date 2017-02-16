PROGRAM cdfcurl
  !!======================================================================
  !!                     ***  PROGRAM  cdfcurl  ***
  !!=====================================================================
  !!  ** Purpose : Compute the curl on F-points for given gridU gridV 
  !!               files and variables
  !!
  !!  ** Method  : Use the same algorithm than NEMO
  !!
  !! History : 2.1  : 05/2005  : J.M. Molines : Original code
  !!         : 2.1  : 06/2007  : P. Mathiot   : for use with forcing fields
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

  INTEGER(KIND=4)                           :: ji, jj, jk, jt     ! dummy loop index
  INTEGER(KIND=4)                           :: npiglo, npjglo     ! size of the domain
  INTEGER(KIND=4)                           :: npk, npt           ! size of the domain
  INTEGER(KIND=4)                           :: nlev               ! number of output levels
  INTEGER(KIND=4)                           :: narg, iargc        ! browse command line
  INTEGER(KIND=4)                           :: ijarg              !
  INTEGER(KIND=4)                           :: ncout, ierr        ! browse command line
  INTEGER(KIND=4), DIMENSION(1)             :: ipk, id_varout     ! output variable properties
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: nilev             ! level to be processed

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e2v, e1u, e1f, e2f ! horizontql metrics
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: un, vn             ! velocity field
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zun, zvn           ! working arrays
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: fmask              ! mask
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: tim, zdep, gdep    ! time counter
  REAL(KIND=4)                              :: zmask              ! mask at T point for -T option

  REAL(KIND=8)                              :: dl_pi, dl_omega    ! 3.14159... and earth rotation rad/sec
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: drotn              ! curl 
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dl_rotn            ! curl at T point 
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dl_ff              ! Coriolis parameter at F point

  CHARACTER(LEN=2048)                        :: cf_ufil, cf_vfil   ! file names
  CHARACTER(LEN=2048)                        :: cf_out = 'curl.nc' ! output file name
  CHARACTER(LEN=2048)                        :: cv_u, cv_v         ! variable names
  CHARACTER(LEN=2048)                        :: cldum              ! dummy string

  TYPE (variable), DIMENSION(1)             :: stypvar            ! structure for attibutes

  LOGICAL                                   :: lforcing = .FALSE. ! forcing flag
  LOGICAL                                   :: lchk     = .FALSE. ! flag for missing files
  LOGICAL                                   :: lperio   = .FALSE. ! flag for E-W periodicity
  LOGICAL                                   :: ltpoint  = .FALSE. ! flag for T-point output
  LOGICAL                                   :: ldblpr   = .FALSE. ! flag for dble precision output
  LOGICAL                                   :: lsurf    = .FALSE. ! flag for 1 lev on C grid.
  LOGICAL                                   :: loverf   = .FALSE. ! flag for 1 lev on C grid.
  LOGICAL                                   :: lnc4=.false.       ! flag for netcdf4 output with chunking and deflation

  !!----------------------------------------------------------------------
  CALL ReadCdfNames() 

  narg = iargc()
  IF ( narg < 5 ) THEN
     PRINT *,' usage : cdfcurl -u U-file U-var -v V-file V-var -l levlist [-T] [-8]...'
     PRINT *,'           ... [-surf] [-overf] [-nc4] [-o OUT-file ]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the curl of a vector field, at a specified level.'  
     PRINT *,'       If level is specified as 0, assume that the input files are'
     PRINT *,'       forcing files, presumably on A-grid. In this latter case, the'
     PRINT *,'       vector field is interpolated on the C-grid. In any case, the'
     PRINT *,'       curl is computed on the F-point (unless -T option is used).'
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
     PRINT *,'       -T : compute curl at T point instead of default F-point'
     PRINT *,'       -8 : save in double precision instead of standard simple precision.'
     PRINT *,'       -surf : work with single level C-grid (not forcing)'
     PRINT *,'       -overf : store the ratio curl/f where f is the coriolis parameter'
     PRINT *,'       -nc4 : use netcdf4 output with chunking and deflation 1'
     PRINT *,'       -o OUT-file : specify output file name instead of ',TRIM(cf_out) 
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        ', TRIM(cn_fhgr)
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : socurl or socurlt (if -T option), units : s^-1'
     PRINT *,'            or socurloverf, no units (if -overf option)'
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
     CASE ( '-nc4' )
        lnc4 = .true.
     CASE ('-T')
        ltpoint = .true.
     CASE ('-8')
        ldblpr = .true.
     CASE ('-surf')
        lsurf = .true.
     CASE ('-overf')
        loverf = .true.
     CASE ('-o')
        CALL getarg(ijarg, cf_out) ; ijarg=ijarg+1
     CASE DEFAULT
        PRINT *,  TRIM(cldum), ' : unknown option '
     END SELECT
  ENDDO

  lchk = chkfile(cn_fhgr ) .OR. lchk
  lchk = chkfile(cf_ufil ) .OR. lchk
  lchk = chkfile(cf_vfil ) .OR. lchk
  IF ( lchk ) STOP ! missing files


  ! define new variables for output
  stypvar(1)%cname             = 'socurl'
  IF (ltpoint) stypvar(1)%cname             = 'socurlt'
  stypvar(1)%cunits            = 's-1'
  IF (loverf)  stypvar(1)%cname             = 'socurloverf'
  IF (loverf)  stypvar(1)%cunits            = '-'

  stypvar(1)%cprecision        ='r4'
  IF ( ldblpr )  stypvar(1)%cprecision     ='r8'
  stypvar(1)%rmissing_value    = 0.
  stypvar(1)%valid_min         = -1000.
  stypvar(1)%valid_max         =  1000.
  stypvar(1)%clong_name        = 'Relative_Vorticity (curl)'
  stypvar(1)%cshort_name       = 'socurl'
  stypvar(1)%conline_operation = 'N/A'
  stypvar(1)%caxis             = 'TYX'

  ipk(1) = nlev  !  nlevel so far

  npiglo = getdim(cf_ufil,cn_x)
  npjglo = getdim(cf_ufil,cn_y)
  npk    = getdim(cf_ufil,cn_z)
  npt    = getdim(cf_ufil,cn_t) 

  PRINT *, 'npiglo = ',npiglo
  PRINT *, 'npjglo = ',npjglo
  PRINT *, 'npk    = ',npk
  PRINT *, 'npt    = ',npt
  PRINT *, 'nlev   = ',nlev

  ! choose chunk size for output ... not easy not used if lnc4=.false. but anyway ..
  stypvar(1)%ichunk=(/npiglo,MAX(1,npjglo/30),1,1 /)
  print *, stypvar(1)%ichunk

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
     lforcing=.true.
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
  ALLOCATE ( e1u(npiglo,npjglo) , e1f(npiglo,npjglo) )
  ALLOCATE ( e2v(npiglo,npjglo) , e2f(npiglo,npjglo) )
  ALLOCATE ( un(npiglo,npjglo)  , vn(npiglo,npjglo)  )
  ALLOCATE ( zun(npiglo,npjglo) , zvn(npiglo,npjglo) )
  ALLOCATE ( drotn(npiglo,npjglo) , fmask(npiglo,npjglo) )
  ALLOCATE ( tim(npt) )
  ALLOCATE ( gdep(nlev) , zdep(npk))

  IF ( ltpoint) ALLOCATE (dl_rotn(npiglo,npjglo) )

  e1u =  getvar(cn_fhgr, cn_ve1u, 1, npiglo, npjglo)
  e1f =  getvar(cn_fhgr, cn_ve1f, 1, npiglo, npjglo)
  e2v =  getvar(cn_fhgr, cn_ve2v, 1, npiglo, npjglo)
  e2f =  getvar(cn_fhgr, cn_ve2f, 1, npiglo, npjglo)

  ! use zun and zvn to store f latitude and longitude for output
  zun = getvar(cn_fhgr, cn_glamf, 1, npiglo, npjglo)
  zvn = getvar(cn_fhgr, cn_gphif, 1, npiglo, npjglo)

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
    zdep(:) = getvar1d(cf_ufil, cn_vdepthu, npk )
    DO jk=1,nlev
      gdep(jk) = zdep( nilev(jk) )
    ENDDO
  ENDIF


  ! look for  E-W periodicity
  IF ( zun(1,1) == zun(npiglo-1,1) ) lperio = .TRUE.

  ! create output fileset
  ncout = create      (cf_out, cf_ufil, npiglo, npjglo, nlev           , ld_nc4=lnc4)
  ierr  = createvar   (ncout , stypvar, 1,      ipk,    id_varout      , ld_nc4=lnc4)
  ierr  = putheadervar(ncout,  cf_ufil, npiglo, npjglo, nlev, pnavlon=zun, pnavlat=zvn, pdep=gdep)

  tim  = getvar1d(cf_ufil, cn_vtimec, npt      )
  ierr = putvar1d(ncout,   tim,       npt,  'T')

  DO jt=1,npt
     IF (MOD(jt,100)==0 ) PRINT *, jt,'/',npt
     DO jk = 1, nlev
        ! if files are forcing fields
        zun(:,:) =  getvar(cf_ufil, cv_u, nilev(jk) ,npiglo,npjglo, ktime=jt)
        zvn(:,:) =  getvar(cf_vfil, cv_v, nilev(jk) ,npiglo,npjglo, ktime=jt)

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

        ! compute the mask ! must be done every time as we have a level loop now
        ! might be replaced by reading fmask from mask file
        !    IF ( jt==1 ) THEN
        DO jj = 1, npjglo - 1
           DO ji = 1, npiglo - 1
              fmask(ji,jj)=0.
              fmask(ji,jj)= un(ji,jj)*un(ji,jj+1) * vn(ji,jj)*vn(ji+1,jj)
              IF (fmask(ji,jj) /= 0.) fmask(ji,jj)=1.
           ENDDO
        ENDDO
        !    END IF

        drotn(:,:) = 0.d0
        DO jj = 1, npjglo -1 
           DO ji = 1, npiglo -1   ! vector opt.
              drotn(ji,jj) = (  e2v(ji+1,jj  ) * vn(ji+1,jj  ) - e2v(ji,jj) * vn(ji,jj)    &
                   &          - e1u(ji  ,jj+1) * un(ji  ,jj+1) + e1u(ji,jj) * un(ji,jj)  ) &
                   &          * fmask(ji,jj) / ( e1f(ji,jj) * e2f(ji,jj) )
           END DO
        END DO

        IF ( lperio ) drotn(npiglo,:) = drotn(2, :)
        IF ( ltpoint ) THEN
           dl_rotn(:,:) = 0.d0
           DO ji = 2, npiglo
              DO jj = 2, npjglo
                 zmask = fmask(ji,jj)*fmask(ji,jj-1)*fmask(ji-1,jj)*fmask(ji-1,jj-1)
                 dl_rotn(ji,jj) = 0.25*( drotn(ji,jj) + drotn(ji,jj-1) + drotn(ji-1,jj) + drotn(ji-1,jj-1) ) * zmask
              ENDDO
           ENDDO
           IF ( lperio ) dl_rotn(1,:) = dl_rotn(npiglo, :)
           drotn(:,:) = dl_rotn(:,:)

        ENDIF
        ! write drotn on file at level k and at time jt
        IF ( loverf ) THEN
           WHERE (dl_ff /= 0.d0) drotn=drotn/dl_ff
        ENDIF

        ierr = putvar(ncout, id_varout(1), drotn, nilev(jk), npiglo, npjglo, ktime=jt)
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
    CHARACTER(LEN=2048) :: cldum
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

END PROGRAM cdfcurl


