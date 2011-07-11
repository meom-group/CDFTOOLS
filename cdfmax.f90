PROGRAM cdfmax
  !!======================================================================
  !!                     ***  PROGRAM  cdfmax  ***
  !!=====================================================================
  !!  ** Purpose : Find the min/max of a variable of an nc file. Give its 
  !!               location. A sub-area can be specified either horizontally
  !!               and/or vertically.
  !!
  !! History : 2.1  : 11/2006  : J.M. Molines : Original code
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

  INTEGER(KIND=4)                               :: ji, jk, jvar, jt
  INTEGER(KIND=4)                               :: idep
  INTEGER(KIND=4)                               :: narg, iargc, ijarg
  INTEGER(KIND=4)                               :: ni, nj, nk, nt       ! size of the global domain
  INTEGER(KIND=4)                               :: ndim                 ! dimension of the variables
  INTEGER(KIND=4)                               :: ntype                ! type of slab (xy, xz, yz ...)
  INTEGER(KIND=4)                               :: ii1, ii2, ij1, ij2   ! index of min max  
  INTEGER(KIND=4)                               :: niz, njz, nkz, nvars ! size of the domain
  INTEGER(KIND=4)                               :: iimin=1, iimax=0     ! i-limit of the domain
  INTEGER(KIND=4)                               :: ijmin=1, ijmax=0     ! j-limit of the domain
  INTEGER(KIND=4)                               :: ikmin=1, ikmax=0     ! k-limit of the domain
  INTEGER(KIND=4)                               :: itmin=1, itmax=0     ! t-limit of the domain
  INTEGER(KIND=4)                               :: istatus              ! working integer
  INTEGER(KIND=4), DIMENSION(2)                 :: ilmin, ilmax         ! working array for minloc, maxloc

  REAL(KIND=4)                                  :: rfact=1.0            ! multiplying factor
  REAL(KIND=4)                                  :: zspval               ! missing value or spval
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE     :: h                    ! depth 
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE     :: v2d, rlon, rlat      ! data array, longitude, latitude

  CHARACTER(LEN=256)                            :: cf_in                ! input file name
  CHARACTER(LEN=256)                            :: cv_in='none'         ! current variable name
  CHARACTER(LEN=256)                            :: cldum                ! dummy char variable
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_names             ! list of variables in file

  TYPE(variable), DIMENSION(:),     ALLOCATABLE :: stypvar              ! dummy dtructure to read var names

  LOGICAL                                       :: lforcexy=.FALSE.     ! flag for forced horizontal slab
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfmax -f file -var cdfvar ...'
     PRINT *,'      ... [-lev kmin kmax ] [-zoom imin imax jmin jmax] ...'
     PRINT *,'      ... [-time tmin tmax ] [-fact multfact]  [-xy ]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'        Find minimum and maximum of a file as well as their '
     PRINT *,'        respective location. Options allow to restrict the '
     PRINT *,'        finding to a sub area in time and space. This program'
     PRINT *,'        also deal with vertical slabs in a domain.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f file  : input file '
     PRINT *,'       -var cdfvar : input variable'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-lev kmin kmax ] : restrict to level between kmin and kmax. '
     PRINT *,'       [-zoom imin imax jmin jmax] : restrict to sub area specified'
     PRINT *,'                       by the given limits. If the zoomed area is '
     PRINT *,'                       degenerated to a single line, then the vertical'
     PRINT *,'                       slab is considered as domain.'
     PRINT *,'       [-time tmin tmax ] : restrict to the indicated time windows.'
     PRINT *,'       [-fact multfact] : use a multiplicative factor for the output'
     PRINT *,'       [-xy ] : force horizontal slab even in the case of a degenerated'
     PRINT *,'                       zoomed area.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       none' 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       output is done on standard output.'
     STOP
  ENDIF

  ijarg=1
  DO  WHILE (ijarg <=  narg)
     CALL getarg(ijarg, cldum) ; ijarg = ijarg + 1
     SELECT CASE (cldum )
     CASE ( '-f'    )
        CALL getarg(ijarg, cf_in) ; ijarg = ijarg + 1
     CASE ( '-lev'  )
        CALL getarg(ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) ikmin
        CALL getarg(ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) ikmax
     CASE ( '-fact'   )
        CALL getarg(ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) rfact
     CASE ( '-zoom' )
        CALL getarg(ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) iimin
        CALL getarg(ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) iimax
        CALL getarg(ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) ijmin
        CALL getarg(ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) ijmax
     CASE ( '-time' )
        CALL getarg(ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) itmin
        CALL getarg(ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) itmax
     CASE ( '-var'  )
        CALL getarg(ijarg, cv_in) ; ijarg = ijarg + 1 
     CASE ( '-xy'   )
        lforcexy = .TRUE.
     CASE DEFAULT
        PRINT *, cldum,' : unknown option '
        STOP
     END SELECT
  END DO

  IF ( chkfile(cf_in) ) STOP ! missing file

  ni=0 ; nj=0; nk=0; nt=0 

  ni = getdim(cf_in, cn_x, cldum, istatus)
  IF ( istatus == 1 ) THEN 
     ni = getdim(cf_in, 'lon', cldum, istatus)
     IF ( istatus == 1 ) THEN
        PRINT *,' No X or lon dim found ' ; STOP
     ENDIF
  ENDIF
  IF ( iimax == 0 ) iimax = ni

  nj = getdim(cf_in, cn_y, cldum, istatus)
  IF ( istatus == 1 ) THEN 
     nj = getdim(cf_in, 'lat', cldum, istatus)
     IF ( istatus == 1 ) THEN
        PRINT *,' No y or lat dim found ' ; STOP
     ENDIF
  ENDIF
  IF ( ijmax == 0 ) ijmax = nj

  nk=getdim(cf_in, cn_z, cldum, istatus)
  IF ( istatus == 1 ) THEN 
     nk = getdim(cf_in, 'z', cldum, istatus)
     IF ( istatus == 1 ) THEN 
        nk = getdim(cf_in, 'lev', cldum, istatus)
        IF ( istatus == 1 ) THEN
           PRINT *,' No dep or z or lev  dim found ' 
        ENDIF
     ENDIF
  ENDIF
  IF ( ikmax == 0 ) ikmax = nk

  nt = getdim(cf_in, cn_t, cldum, istatus)

  IF ( istatus == 1 ) THEN 
     nt = getdim(cf_in, 'step', cldum, istatus)
     IF ( istatus == 1 ) THEN
        PRINT *,' No time or step dim found ' 
     ENDIF
  ENDIF


  ! fix the size of the zoomed area, or the whole domain if no zoom
  niz = iimax - iimin + 1
  njz = ijmax - ijmin + 1
  nkz = ikmax - ikmin + 1

  IF (nt == 0 ) nt = 1  ! assume a 1 time frame file
  IF ( itmax == 0 ) itmax = nt
  ! allocate arrays
  ALLOCATE (h(nk), rlon(niz,njz), rlat(niz,njz))

  ! Look for variable name starting with dep
  nvars     = getnvar(cf_in)
  ALLOCATE (cv_names(nvars), stypvar(nvars))
  cv_names = getvarname(cf_in,nvars,stypvar)
  DO jvar=1,nvars
     idep = INDEX(cv_names(jvar),'dep') + INDEX(cv_names(jvar),'lev')
     IF (idep /= 0 ) EXIT
  END DO

  IF ( jvar == nvars +1 ) THEN
     ! no depth variable found ... we initialize it to levels
     h = (/(ji,ji=1,nk)/)
  ELSE
     h = getvar1d(cf_in, cv_names(jvar), nk)
  ENDIF
  zspval = getatt(cf_in, cv_in, cn_missing_value)

  ! Allocate memory and define ntype : (1) = horizontal i-j slab eventually many layers.
  !                                    (2) = vertical j-k slab, at a given i
  !                                    (3) = vertical i-k slab, at a given j
  IF ( (niz /= 1 .AND. njz /= 1 ) .OR. lforcexy  ) THEN 
     ALLOCATE (v2d(niz,njz) )
     ntype = 1    ! horizontal x-y slabs
  ELSE  
     IF ( niz == 1 ) THEN
        ALLOCATE (v2d(njz,nkz))
        ntype = 2 ! vertical y-z slab
     ELSE
        ALLOCATE(v2d(niz,nkz))
        ntype = 3 ! vertical x-z slab
     ENDIF
  ENDIF

  ! read latitude, longitude from the header
  rlon = getvar(cf_in, cn_vlon2d, 1, niz, njz, iimin, ijmin)
  rlat = getvar(cf_in, cn_vlat2d, 1, niz, njz, iimin, ijmin)

  DO
     ndim = getvdim(cf_in, cv_in) + 1   ! getvdim gives ndim-1 !
     PRINT *,TRIM(cv_in),' with multiplying factor of ', rfact
     ! ndim <=3 corresponds to purely 2D variables (x,y) or (x,y,t)
     IF ( ndim <= 3 ) THEN
        ikmin = 1 ; ikmax = 1 ; nkz = 1
     ENDIF

     SELECT CASE (ntype)
     CASE (1)
        SELECT CASE (ndim)
        CASE( 2,3,4 )  ! assume x,y,z,t variable
           PRINT 9000,'time  level     dep  MAX:   i    long    j    lat   MaxValue   MIN:     i    long    j   lat    MinValue'
           DO jt=itmin, itmax
              DO jk=ikmin,ikmax
                 v2d(:,:) = getvar(cf_in, cv_in, jk, niz, njz, kimin=iimin, kjmin=ijmin, ktime=jt)
                 ilmax = MAXLOC(v2d,(v2d /= zspval) )
                 ilmin = MINLOC(v2d,(v2d /= zspval) )
                 ii1=ilmax(1) ; ij1=ilmax(2)
                 ii2=ilmin(1) ; ij2=ilmin(2)
                 PRINT 9003, jt, jk, h(jk),ii1+iimin -1, rlon(ii1,ij1),ij1+ijmin -1,rlat(ii1,ij1),v2d(ii1,ij1)*rfact, &
                      &             ii2+iimin -1, rlon(ii2,ij2),ij2+ijmin -1,rlat(ii2,ij2),v2d(ii2,ij2)*rfact
              END DO
           END DO
           EXIT

        CASE DEFAULT
           PRINT *,' Non mapable variables x-y :('
           cv_in='none'
        END SELECT

     CASE (2)
        SELECT CASE (ndim)
        CASE( 4 )  ! assume x,y,z,t variable
           PRINT 9000,' time i-slab  MAX:   i    long   j    lat   k     dep    MaxValue    MIN:  i    &
                & long   j     lat   k     dep    MinValue'
           DO jt=itmin, itmax
              v2d(:,:) = getvaryz(cf_in, cv_in, iimin, njz, nkz, ijmin, ikmin, ktime=jt)
              ilmax = MAXLOC(v2d,(v2d/= zspval) )
              ilmin = MINLOC(v2d,(v2d/= zspval) )
              ii1=ilmax(1) ; ij1=ilmax(2)
              ii2=ilmin(1) ; ij2=ilmin(2)
              PRINT 9002, jt, iimin, iimin, rlon(1,ii1),ii1+ijmin -1,rlat(1,ii1),ij1+ikmin-1, h(ij1+ikmin-1), v2d(ii1,ij1)*rfact, &
                   &             iimin, rlon(1,ii2),ii2+ijmin -1,rlat(1,ii2),ij2+ikmin-1, h(ij2+ikmin-1), v2d(ii2,ij2)*rfact
           END DO
           EXIT
        CASE DEFAULT
           PRINT *,' Non mapable variables x-z or y-z :('
           cv_in='none'
        END SELECT

     CASE (3)
        SELECT CASE (ndim)
        CASE( 4 )  ! assume x,y,z,t variable
           PRINT 9000,' time j-slab  MAX:   i    long   j    lat   k     dep    MaxValue    MIN:  i    &
                & long   j     lat   k     dep    MinValue'
           DO jt=itmin, itmax
              v2d(:,:) = getvarxz(cf_in, cv_in, ijmin, niz, nkz, iimin, ikmin, ktime=jt)
              ilmax = MAXLOC(v2d,(v2d /= zspval) )
              ilmin = MINLOC(v2d,(v2d /= zspval) )
              ii1=ilmax(1) ; ij1=ilmax(2)
              ii2=ilmin(1) ; ij2=ilmin(2)
              PRINT 9002, jt, ijmin, ii1, rlon(ii1,1),ijmin,rlat(ii1,1),ij1+ikmin-1, h(ij1+ikmin-1), v2d(ii1,ij1)*rfact, &
                   &             ii2, rlon(ii2,1),ijmin,rlat(ii2,1),ij2+ikmin-1, h(ij2+ikmin-1), v2d(ii2,ij2)*rfact
           END DO
           EXIT
        CASE DEFAULT
           PRINT *,' Non mapable variables x-z or y-z :('
           cv_in='none'
        END SELECT

     CASE DEFAULT
        PRINT *,' ntype = ',ntype, '  is not defined ' ; STOP
     END SELECT ! ntype
  ENDDO

9000 FORMAT(a)
9002 FORMAT(I5, x,i4,9x, i5,f8.2, i5, f7.2, i5, f8.2, e14.5, 6x, i5,f8.2, i5, f7.2, i5, f8.2, e14.5 )
9003 FORMAT(I5, x,i5,1x,f7.2,5x,i5,f8.2, i5, f7.2, e14.5, 5x,i5,f8.2, i5, f7.2, e14.5)

END PROGRAM cdfmax
