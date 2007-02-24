PROGRAM cdfzoom
  !!----------------------------------------------------------------------------
  !!                ***   PROGRAM cdfzoom ***
  !!
  !!  ** Purpose:    Extract a sub area of a cdf output file and print it on the screen
  !!                   with an easy to read format
  !!
  !!  ** Method:     Read command line, open the file get the variable and show the sub area
  !!
  !!  ** Usage :  cdfzoom -f file -zoom imin imax jmin jmax -fact factor -lev klev 
  !!
  !!   History: 
  !!       1999 : Anne de Miranda (bimgzoom)
  !!       2001 : J-M Molines for normalization
  !!       2004 : J-M Molines : support for NetCdf IOIPSL files
  !!       2006 : J-M Molines : included as cdfzoom in cdftools
  !!
  !!----------------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  ! * Module used
  USE cdfio

  ! * Local Variable
  IMPLICIT NONE
  INTEGER ,PARAMETER :: jpk=100, jpt=700
  !
  INTEGER :: numin,jk,ji,jj,jt,jl, jd, i
  INTEGER ::  narg, iargc
  INTEGER :: isdirect
  INTEGER :: ni,nj,nk,nt,icod,ndim
  INTEGER :: niz,njz, nkz, itime, nvars
  INTEGER :: imin, imax, jmin, jmax,kext, istatus, kmin, kmax
  INTEGER :: ipmin, ipmax, jpmin, jpmax
  !
  REAL ,DIMENSION(:),ALLOCATABLE     :: h, rtime
  REAL ,DIMENSION (:,:), ALLOCATABLE :: v2d
  REAL                               :: fact
  !
  CHARACTER(LEN=100) ::  cfilein, cline1, cline2
  CHARACTER(LEN=80) :: cvar='none', cdim
  !!
  !! 1. Initializations:
  !! -------------------
  !!
  narg = iargc()
  IF (narg == 0) THEN
     PRINT *,'usage :cdfzoom -f file '// &
          ' -lev kmin kmax -fact facteur' //  &
          ' -zoom imin imax jmin jmax' // &
          ' -var cdfvarname '
     STOP
  END IF
  !
  kext=1
  fact=1
  numin = 10
  i=1
  ! Read command line
  DO  WHILE (i <=  narg)
     CALL getarg(i,cline1)
     i = i + 1
     IF (cline1 == '-f') THEN
        CALL getarg(i,cline2)
        i = i + 1
        cfilein=cline2
     ELSE IF (cline1 == '-lev') THEN
        CALL getarg(i,cline2)
        i = i + 1
        READ(cline2,*) kmin
        CALL getarg(i,cline2)
        i = i + 1
        READ(cline2,*) kmax
     ELSE IF (cline1 == '-fact') THEN
        CALL getarg(i,cline2)
        i = i + 1
        READ(cline2,*) fact
     ELSE IF (cline1 == '-zoom') THEN
        CALL getarg(i,cline2)
        i = i + 1
        READ(cline2,*) imin
        CALL getarg(i,cline2)
        i = i + 1
        READ(cline2,*) imax
        CALL getarg(i,cline2)
        i = i + 1
        READ(cline2,*) jmin
        CALL getarg(i,cline2)
        i = i + 1
        READ(cline2,*) jmax
     ELSE IF ( cline1 == '-var') THEN
        CALL getarg(i,cvar)
        i = i + 1

     ELSE
        PRINT *, cline1,' : unknown option '
        STOP
     END IF
  END DO
  !
  ni=0 ; nj=0; nk=0; nt=0 
  niz=imax-imin+1
  njz=jmax-jmin+1
  nkz=kmax-kmin+1

  IF (nkz > 1 ) THEN
     !working with vertical slab, either niz or njz must be 1
     IF ( niz == 1  ) THEN ! y/z slab
     ELSE IF ( njz == 1 ) THEN ! x/z slab
     ELSE
        PRINT *, 'Either niz or njz must me  one'
        STOP
     ENDIF
  ENDIF

  ni=getdim(cfilein,'x',cdim,istatus)
  IF ( istatus == 1 ) THEN 
     ni=getdim(cfilein,'lon',cdim,istatus)
     IF ( istatus == 1 ) THEN
        PRINT *,' No X or lon dim found ' ; STOP
     ENDIF
  ENDIF

  nj=getdim(cfilein,'y',cdim,istatus)
  IF ( istatus == 1 ) THEN 
     nj=getdim(cfilein,'lat',cdim,istatus)
     IF ( istatus == 1 ) THEN
        PRINT *,' No y or lat dim found ' ; STOP
     ENDIF
  ENDIF

  nk=getdim(cfilein,'dep',cdim,istatus)
  IF ( istatus == 1 ) THEN 
     nk=getdim(cfilein,'z',cdim,istatus)
     IF ( istatus == 1 ) THEN 
        nk=getdim(cfilein,'lev',cdim,istatus)
        IF ( istatus == 1 ) THEN
           PRINT *,' No dep or z or lev  dim found ' 
        ENDIF
     ENDIF
  ENDIF

  nt=getdim(cfilein,'time',cdim,istatus)
  IF ( istatus == 1 ) THEN 
     nt=getdim(cfilein,'step',cdim,istatus)
     IF ( istatus == 1 ) THEN
        PRINT *,' No time or step dim found ' 
     ENDIF
  ENDIF


  IF (nk == 0 ) nk = 1 ; kext=1  ! assume a 2D variable
  IF (nt == 0 ) nt = 1 ; itime=1 ! assume a 1 time frame file
  ALLOCATE (h(nk), rtime(nt))

  IF (nkz == 1 ) THEN 
     ALLOCATE (v2d(niz,njz) )
  ELSE  
     IF ( niz == 1 ) THEN
        ALLOCATE (v2d(njz,nkz))
     ELSE
        ALLOCATE(v2d(niz,nkz))
     ENDIF
  ENDIF

  DO 
     ndim=getvdim(cfilein,cvar)+1   ! getvdim gives ndim-1 !
     PRINT *,TRIM(cvar), ndim
     SELECT CASE (nkz)
     CASE (1)
        ipmin=imin ; ipmax=imax; jpmin=jmin; jpmax=jmax
        SELECT CASE (ndim)
        CASE( 2 )  ! assume x,y variable
           v2d(:,:)=getvar(cfilein,cvar,1,niz,njz,imin,jmin)
           EXIT
        CASE( 3 )  ! assume x,y,t variable
           v2d(:,:)=getvar(cfilein,cvar,1,niz,njz,imin,jmin)
           EXIT
        CASE( 4 )  ! assume x,y,z,t variable
           v2d(:,:)=getvar(cfilein,cvar,kext,niz,njz,imin,jmin)
           EXIT
        CASE DEFAULT
           PRINT *,' Non mapable variables x-y :('
           cvar='none'
        END SELECT

     CASE DEFAULT
       SELECT CASE (ndim)
        CASE( 4 )  ! assume x,y,z,t variable
           IF ( njz == 1 ) THEN
             ipmin=imin ; ipmax=imax; jpmin=kmin; jpmax=kmax
             v2d(:,:)=getvarxz(cfilein,cvar,jmin,niz,nkz,imin,kmin)
           ELSE
             ipmin=jmin ; ipmax=jmax; jpmin=kmin; jpmax=kmax
             v2d(:,:)=getvaryz(cfilein,cvar,imin,njz,nkz,jmin,kmin)
           ENDIF
           EXIT
        CASE DEFAULT
           PRINT *,' Non mapable variables x-z or y-z :('
           cvar='none'
        END SELECT

     END SELECT ! nkz
  ENDDO

  PRINT *,'IMIN IMAX JMIN JMAX KMIN KMAX', imin,imax,jmin,jmax,kmin,kmax
  PRINT 9001,'      ',(ji,ji=ipmin,ipmax)
  IF (nkz == 1 ) THEN
    DO jj=jpmax,jpmin,-1
       PRINT 9000,jj,'  ',(v2d(ji-ipmin+1,jj-jpmin+1)/fact,ji=ipmin,ipmax)
    END DO
  ELSE
    DO jj=jpmin,jpmax
       PRINT 9000,jj,'  ',(v2d(ji-ipmin+1,jj-jpmin+1)/fact,ji=ipmin,ipmax)
    END DO
  ENDIF
9000 FORMAT(i4,a,20f12.4)
9001 FORMAT(a,20i12)

END PROGRAM cdfzoom
