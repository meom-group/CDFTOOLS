PROGRAM cdfzoom
  !!======================================================================
  !!                     ***  PROGRAM  cdfzoom  ***
  !!=====================================================================
  !!  ** Purpose : Extract a sub area of a cdf output file and print it 
  !!               on the screen with an easy to read format.
  !!
  !!  ** Method  : specify the variable name and file on the command line
  !!
  !! History : ---  : 1999     : A. de Miranda : Original code in bimgtools
  !! History : 2.1  : 11/2004  : J.M. Molines  : port to CDFTOOLS
  !!           3.0  : 12/2010  : J.M. Molines  : Doctor norm + Lic.
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
  !
  INTEGER(KIND=4)                           :: ji, jj, jt             ! dummy loop index
  INTEGER(KIND=4)                           :: narg, iargc, ijarg     ! browse line
  INTEGER(KIND=4)                           :: ni, nj, nk, nt, ndim   ! domain dimension
  INTEGER(KIND=4)                           :: niz, njz, nkz          ! size of zoom
  INTEGER(KIND=4)                           :: iimin, iimax           ! i-limits
  INTEGER(KIND=4)                           :: ijmin, ijmax           ! j-limits
  INTEGER(KIND=4)                           :: ikmin, ikmax           ! k-limits
  INTEGER(KIND=4)                           :: itmin, itmax           ! t-limit
  INTEGER(KIND=4)                           :: ikext, ierr            ! 
  INTEGER(KIND=4)                           :: iipmin, iipmax         ! 
  INTEGER(KIND=4)                           :: ijpmin, ijpmax         !
  !
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: v2d                    ! data array
  REAL(KIND=4)                              :: fact                   ! scaling factor
  !
  CHARACTER(LEN=256)                        :: cldum                  ! summy character variable
  CHARACTER(LEN=256)                        :: cf_in                  ! input file name
  CHARACTER(LEN=256)                        :: cv_in='none'           ! variable name
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfzoom -f IN-file -zoom imin imax jmin jmax -v IN-var ...'
     PRINT *,'               ... [-lev kmin kmax ] [-time tmin tmax ] [-fact factor]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'      Displays the numerical values of a zoomed area. By default, all times'
     PRINT *,'      and levels are shown. If the zoomed area is degenerated to a single '
     PRINT *,'      line, then the vertical slab is displayed.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f IN-file : name of input file' 
     PRINT *,'       -zoom imin imax jmin jmax : spatial window definition'
     PRINT *,'       -v IN-var : cdf variable name to work with.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-lev kmin kmax ]  : vertical limits for display.' 
     PRINT *,'       [-time tmin tmax ] : time limits for display.' 
     PRINT *,'       [-fact factor ]    : use a scaling factor for display.'
     PRINT *,'                            Values are DIVIDED by factor'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        none'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       display on standard output'
     STOP 
  ENDIF
  !
  ikext = 1 ; ikmin = 1 ; ikmax = 1 ; itmin = 1 ; itmax = 1
  fact  = 1

  ijarg  = 1
  ! Read command line
  DO  WHILE (ijarg <=  narg)
     CALL getarg(ijarg,cldum) ; ijarg = ijarg + 1
     SELECT CASE ( cldum )
     CASE ( '-f'   ) ; CALL getarg(ijarg, cf_in ) ; ijarg = ijarg + 1
     CASE ( '-v'   ) ; CALL getarg(ijarg, cv_in ) ; ijarg = ijarg + 1
     CASE ( '-zoom') ; CALL getarg(ijarg, cldum ) ; ijarg = ijarg + 1 ; READ(cldum,*) iimin
        ;              CALL getarg(ijarg, cldum ) ; ijarg = ijarg + 1 ; READ(cldum,*) iimax
        ;              CALL getarg(ijarg, cldum ) ; ijarg = ijarg + 1 ; READ(cldum,*) ijmin
        ;              CALL getarg(ijarg, cldum ) ; ijarg = ijarg + 1 ; READ(cldum,*) ijmax
     ! options
     CASE ( '-lev' ) ; CALL getarg(ijarg, cldum ) ; ijarg = ijarg + 1 ; READ(cldum,*) ikmin
        ;              CALL getarg(ijarg, cldum ) ; ijarg = ijarg + 1 ; READ(cldum,*) ikmax
     CASE ( '-time') ; CALL getarg(ijarg, cldum ) ; ijarg = ijarg + 1 ; READ(cldum,*) itmin
        ;              CALL getarg(ijarg, cldum ) ; ijarg = ijarg + 1 ; READ(cldum,*) itmax
     CASE ( '-fact') ; CALL getarg(ijarg, cldum ) ; ijarg = ijarg + 1 ; READ(cldum,*) fact
     CASE DEFAULT    ; PRINT *,' ERROR : ', TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  END DO

  IF ( chkfile (cf_in) ) STOP 99 ! missing file
  !
  ni=0 ; nj=0 ; nk=0 ; nt=0 
  niz  = iimax - iimin + 1
  njz  = ijmax - ijmin + 1
  nkz  = ikmax - ikmin + 1
  ikext= ikmin

  IF ( nkz > 1 ) THEN
     !working with vertical slab, either niz or njz must be 1
     IF ( niz == 1  ) THEN ! y/z slab
     ELSE IF ( njz == 1 ) THEN ! x/z slab
     ELSE
        PRINT *, 'Either niz or njz must be  one'
        STOP 99
     ENDIF
  ENDIF

  ni = getdim(cf_in, cn_x, cldum, ierr)
  IF ( ierr == 1 ) THEN 
     ni = getdim(cf_in, 'lon', cldum, ierr)
     IF ( ierr == 1 ) THEN
        PRINT *,' No X or lon dim found ' ; STOP 99
     ENDIF
  ENDIF

  nj = getdim(cf_in, cn_y, cldum, ierr)
  IF ( ierr == 1 ) THEN 
     nj = getdim(cf_in, 'lat', cldum, ierr)
     IF ( ierr == 1 ) THEN
        PRINT *,' No y or lat dim found ' ; STOP 99
     ENDIF
  ENDIF

  nk = getdim(cf_in, cn_z, cldum, ierr)
  IF ( ierr == 1 ) THEN 
     nk = getdim(cf_in, 'z', cldum, ierr)
     IF ( ierr == 1 ) THEN 
        nk = getdim(cf_in, 'lev', cldum, ierr)
        IF ( ierr == 1 ) THEN
           PRINT *,' No dep or z or lev  dim found ' 
        ENDIF
     ENDIF
  ENDIF

  nt = getdim(cf_in, cn_t, cldum, ierr)
  IF ( ierr == 1 ) THEN 
     nt = getdim(cf_in, 'step', cldum, ierr)
     IF ( ierr == 1 ) THEN
        nt = getdim(cf_in, 'time', cldum, ierr)
        IF ( ierr == 1 ) THEN
           nt = getdim(cf_in, 't', cldum, ierr)
           IF ( ierr == 1 ) THEN
             PRINT *,' No time or step dim found ' 
           ENDIF
        ENDIF
     ENDIF
  ENDIF

  IF ( itmax > nt ) THEN 
     PRINT *,' Not enough time steps in this file' 
     STOP 99
  ENDIF

  IF (nk == 0 ) THEN ; nk = 1 ; ikext = 1 ; ENDIF  ! assume a 2D variable
  IF (nt == 0 ) THEN ; nt = 1             ; ENDIF  ! assume a 1 time frame file

  IF ( nkz == 1 ) THEN 
     ALLOCATE ( v2d(niz,njz) )
  ELSE  
     IF ( niz == 1 ) THEN
        ALLOCATE( v2d(njz,nkz))
     ELSE
        ALLOCATE( v2d(niz,nkz))
     ENDIF
  ENDIF

  DO jt = itmin, itmax
     DO   ! for exit statement
        ndim = getvdim(cf_in, cv_in)+1   ! getvdim gives ndim-1 !
        PRINT *,TRIM(cv_in), ndim, ikext
        SELECT CASE (nkz)
        CASE (1)
           iipmin=iimin ; iipmax=iimax; ijpmin=ijmin; ijpmax=ijmax
           SELECT CASE (ndim)
           CASE( 2 )  ! assume x,y variable
              v2d(:,:) = getvar(cf_in, cv_in, 1,     niz, njz, iimin, ijmin, ktime=jt)
              EXIT
           CASE( 3 )  ! assume x,y,t variable
              v2d(:,:) = getvar(cf_in, cv_in, 1,     niz, njz, iimin, ijmin, ktime=jt)
              EXIT
           CASE( 4 )  ! assume x,y,z,t variable
              v2d(:,:) = getvar(cf_in, cv_in, ikext, niz, njz, iimin, ijmin, ktime=jt)
              EXIT
           CASE DEFAULT
              PRINT *,' Non mapable variables x-y :('
              cv_in='none'
           END SELECT

        CASE DEFAULT
           SELECT CASE (ndim)
           CASE( 4 )  ! assume x,y,z,t variable
              IF ( njz == 1 ) THEN
                 iipmin=iimin ; iipmax=iimax; ijpmin=ikmin; ijpmax=ikmax
                 v2d(:,:) = getvarxz(cf_in, cv_in, ijmin, niz, nkz, iimin, ikmin, ktime=jt)
              ELSE
                 iipmin=ijmin ; iipmax=ijmax; ijpmin=ikmin; ijpmax=ikmax
                 v2d(:,:) = getvaryz(cf_in, cv_in, iimin, njz, nkz, ijmin, ikmin, ktime=jt)
              ENDIF
              EXIT
           CASE DEFAULT
              PRINT *,' Non mapable variables x-z or y-z :('
              cv_in='none'
           END SELECT

        END SELECT ! nkz
     ENDDO

     PRINT *,'IMIN IMAX JMIN JMAX KMIN KMAX TIME', iimin,iimax,ijmin,ijmax,ikmin,ikmax, jt
     PRINT 9001,'      ',(ji,ji=iipmin,iipmax)
     IF (nkz == 1 ) THEN
        DO jj=ijpmax,ijpmin,-1
           PRINT 9000,jj,'  ',(v2d(ji-iipmin+1,jj-ijpmin+1)/fact,ji=iipmin,iipmax)
        END DO
     ELSE
        DO jj=ijpmin,ijpmax
           PRINT 9000,jj,'  ',(v2d(ji-iipmin+1,jj-ijpmin+1)/fact,ji=iipmin,iipmax)
        END DO
     ENDIF
  ENDDO
9000 FORMAT(i4,a,20f12.4)
9001 FORMAT(a,20i12)

END PROGRAM cdfzoom
