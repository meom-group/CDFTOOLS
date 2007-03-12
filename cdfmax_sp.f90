PROGRAM cdfmax_sp
  !!----------------------------------------------------------------------------
  !!                ***   PROGRAM cdfmax_sp ***
  !!
  !!  ** Purpose:    Find the min/max of a variable of an nc file. Give its location.
  !!                 A sub-area can be specified either horizontally and/or vertically
  !!
  !!  ** Method:     Read command line, open the file get the variable and display the values
  !!                 this version takes spval into account
  !!
  !!  ** Usage :  cdfmax -f file  -var cdfvarname  [-lev kmin kmax -zoom imin imax jmin jmax ]
  !!
  !!   History: 
  !!       2006 : J-M Molines :  from cdfzoom.
  !!
  !!----------------------------------------------------------------------------
  !!  $Rev: 13 $
  !!  $Date: 2007-02-24 22:15:22 +0100 (Sat, 24 Feb 2007) $
  !!  $Id: cdfmax.f90 13 2007-02-24 21:15:22Z molines $
  !!--------------------------------------------------------------
  ! * Module used
  USE cdfio

  ! * Local Variable
  IMPLICIT NONE
  !
  INTEGER :: ji,jk,jvar, idep, jt
  INTEGER :: narg, iargc
  INTEGER :: ni,nj,nk,nt,ndim, ntype
  INTEGER :: i1,i2,j1,j2
  INTEGER :: niz,njz, nkz, itime, nvars
  INTEGER :: imin=1, imax=0, jmin=1, jmax=0,kext, istatus, kmin=1, kmax=0
  INTEGER :: ipmin, ipmax, jpmin, jpmax
  INTEGER, DIMENSION (2) :: ilmin, ilmax   
  !
  REAL ,DIMENSION(:),ALLOCATABLE     :: h, rtime
  REAL ,DIMENSION (:,:), ALLOCATABLE :: v2d, rlon, rlat
  REAL                               :: rfact=1.0, spval
  !
  CHARACTER(LEN=100) ::  cfilein, cline1, cline2
  CHARACTER(LEN=80) :: cvar='none', cdim
  CHARACTER(LEN=80), DIMENSION(:),ALLOCATABLE :: cvarnames
  TYPE(variable), DIMENSION(:),ALLOCATABLE :: typvar
  !
  LOGICAL :: lvar=.false., lfil=.false., lforcexy=.false.
  LOGICAL, DIMENSION(:,:), ALLOCATABLE :: lmask
  !!
  !! Initializations:
  !!-----------------
  !!
  narg = iargc()
  IF (narg == 0) THEN
     PRINT *,'USAGE :cdfmax_sp -f file '// &
          ' -var cdfvarname ' 
     PRINT *, '      [-lev kmin kmax ' //  &
          ' -zoom imin imax jmin jmax  -fact multfact -xy ]'
     PRINT *, '   -lev and -zoom limit the area for min/max computation'
     PRINT *, '    if not specified : the 3D data is taken '
     PRINT *, '    if either imin=imax or jmin=jmax a vertical slab is considered'
     PRINT *, '     UNLESS -xy option is specified !!! '
     
     STOP
  END IF
  !
  kext=1
  ji=1
  ! Read command line
  DO  WHILE (ji <=  narg)
     CALL getarg(ji,cline1)
     ji = ji + 1
     IF (cline1 == '-f') THEN
        lfil=.true.
        CALL getarg(ji,cline2)
        ji = ji + 1
        cfilein=cline2
     ELSE IF (cline1 == '-lev') THEN
        CALL getarg(ji,cline2)
        ji = ji + 1
        READ(cline2,*) kmin
        CALL getarg(ji,cline2)
        ji = ji + 1
        READ(cline2,*) kmax
     ELSE IF (cline1 == '-fact') THEN
        CALL getarg(ji,cline2)
        ji = ji + 1
        READ(cline2,*) rfact
     ELSE IF (cline1 == '-zoom') THEN
        CALL getarg(ji,cline2)
        ji = ji + 1
        READ(cline2,*) imin
        CALL getarg(ji,cline2)
        ji = ji + 1
        READ(cline2,*) imax
        CALL getarg(ji,cline2)
        ji = ji + 1
        READ(cline2,*) jmin
        CALL getarg(ji,cline2)
        ji = ji + 1
        READ(cline2,*) jmax
     ELSE IF ( cline1 == '-var') THEN
        lvar=.true.
        CALL getarg(ji,cvar)
        ji = ji + 1
     ELSE IF ( cline1 == '-xy') THEN
        lforcexy=.true.
     ELSE
        PRINT *, cline1,' : unknown option '
        STOP
     END IF
  END DO
! IF ( .not. lvar .OR. .not. lfil ) THEN
!       PRINT *,' ERROR : you must specify a variable name with -var option AND a filename (-f) '
!       STOP
! ENDIF
  !
  ! Look for dimensions of the variables in the file
  ni=0 ; nj=0; nk=0; nt=0 

  ni=getdim(cfilein,'x',cdim,istatus)
  IF ( istatus == 1 ) THEN 
     ni=getdim(cfilein,'lon',cdim,istatus)
     IF ( istatus == 1 ) THEN
        PRINT *,' No X or lon dim found ' ; STOP
     ENDIF
  ENDIF
  IF ( imax == 0 ) imax =ni

  nj=getdim(cfilein,'y',cdim,istatus)
  IF ( istatus == 1 ) THEN 
     nj=getdim(cfilein,'lat',cdim,istatus)
     IF ( istatus == 1 ) THEN
        PRINT *,' No y or lat dim found ' ; STOP
     ENDIF
  ENDIF
  IF ( jmax == 0 ) jmax =nj

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
  IF ( kmax == 0 ) kmax =nk

  nt=getdim(cfilein,'time',cdim,istatus)
  IF ( istatus == 1 ) THEN 
     nt=getdim(cfilein,'step',cdim,istatus)
     IF ( istatus == 1 ) THEN
        PRINT *,' No time or step dim found ' 
     ENDIF
  ENDIF


  ! fix the size of the zoomed area, or the whole domain if no zoom
  niz=imax-imin+1
  njz=jmax-jmin+1
  nkz=kmax-kmin+1

  IF (nk == 0 ) nk = 1 ; kext=1  ! assume a 2D variable
  IF (nt == 0 ) nt = 1 ; itime=1 ! assume a 1 time frame file
  ! allocate arrays
  ALLOCATE (h(nk), rtime(nt), rlon(niz,njz),rlat(niz,njz))

  ! Look for variable name starting with dep
  nvars=getnvar(cfilein)
  ALLOCATE (cvarnames(nvars), typvar(nvars))
  cvarnames=getvarname(cfilein,nvars,typvar)
  DO jvar=1,nvars
    idep=INDEX(cvarnames(jvar),'dep') + INDEX(cvarnames(jvar),'lev')
    IF (idep /= 0 ) EXIT
  END DO
  IF ( jvar == nvars +1 ) THEN
    ! no depth variable found ... we initialize it to levels
     h=(/(ji,ji=1,nk)/)
  ELSE
    h=getvar1d(cfilein,cvarnames(jvar),nk)
  ENDIF
 
  ! Allocate memory and define ntype : (1) = horizontal i-j slab eventually many layers.
  !                                    (2) = vertical j-k slab, at a given i
  !                                    (3) = vertical i-k slab, at a given j
  IF ( (niz /= 1 .AND. njz /= 1 )  .OR. lforcexy ) THEN 
     ALLOCATE (v2d(niz,njz) ,lmask(niz,njz))
     ntype=1
  ELSE  
     IF ( niz == 1 ) THEN
        ALLOCATE (v2d(njz,nkz),lmask(njz,nkz))
        ntype=2
     ELSE
        ALLOCATE(v2d(niz,nkz),lmask(niz,nkz))
        ntype=3
     ENDIF
  ENDIF

     ! read latitude, longitude from the header
     rlon=getvar(cfilein,'nav_lon',1,niz,njz,imin,jmin)
     rlat=getvar(cfilein,'nav_lat',1,niz,njz,imin,jmin)

DO
     ndim=getvdim(cfilein,cvar)+1   ! getvdim gives ndim-1 !
     PRINT *,TRIM(cvar),' with multiplying factor of ', rfact
     ! ndim <=3 corresponds to purely 2D variables (x,y) or (x,y,t)
     IF ( ndim <= 3 ) THEN
        kmin=1 ; kmax=1 ; nkz=1
     ENDIF
     spval=getatt(cfilein,cvar,'missing_value')

     SELECT CASE (ntype)
     CASE (1)
        ipmin=imin ; ipmax=imax; jpmin=jmin; jpmax=jmax
        SELECT CASE (ndim)
        CASE( 2,3,4 )  ! assume x,y variable
     PRINT 9000,'time  level     dep  MAX:   i    long    j    lat   MaxValue   MIN:     i    long    j   lat    MinValue'
        DO jt=1,nt
          DO jk =kmin,kmax
           v2d(:,:)=getvar(cfilein,cvar,jk,niz,njz,kimin=imin,kjmin=jmin,ktime=jt)
           lmask(:,:)=.true. ; WHERE ( v2d == spval ) lmask=.false.
           ilmax=maxloc(v2d,lmask)
           ilmin=minloc(v2d,lmask)
           i1=ilmax(1) ; j1=ilmax(2)
           i2=ilmin(1) ; j2=ilmin(2)
           PRINT 9003, jt, jk, h(jk),i1+imin -1, rlon(i1,j1),j1+jmin -1,rlat(i1,j1),v2d(i1,j1)*rfact, &
                   &             i2+imin -1, rlon(i2,j2),j2+jmin -1,rlat(i2,j2),v2d(i2,j2)*rfact
          END DO
         END DO
          EXIT

        CASE DEFAULT
           PRINT *,' Non mapable variables x-y :('
           cvar='none'
        END SELECT

     CASE (2)
       SELECT CASE (ndim)
        CASE( 4 )  ! assume x,y,z,t variable
             ipmin=jmin ; ipmax=jmax; jpmin=kmin; jpmax=kmax
             v2d(:,:)=getvaryz(cfilein,cvar,imin,njz,nkz,jmin,kmin)
           lmask(:,:)=.true. ; WHERE ( v2d == 0 ) lmask=.false.
           ilmax=maxloc(v2d,lmask)
           ilmin=minloc(v2d,lmask)
           i1=ilmax(1) ; j1=ilmax(2)
           i2=ilmin(1) ; j2=ilmin(2)
! sorry for nice identation but if not .. rhodes complains
PRINT 9000,' i-slab  MAX:   i    long   j    lat   k     dep    MaxValue    MIN:  i    long   j     lat   k     dep    MinValue'
PRINT 9002, imin, imin, rlon(1,i1),i1+jmin -1,rlat(1,i1),j1+kmin-1, h(j1+kmin-1), v2d(i1,j1)*rfact, &
         &             imin, rlon(1,i2),i2+jmin -1,rlat(1,i2),j2+kmin-1, h(j2+kmin-1), v2d(i2,j2)*rfact
           EXIT
        CASE DEFAULT
           PRINT *,' Non mapable variables x-z or y-z :('
           cvar='none'
        END SELECT

     CASE (3)
       SELECT CASE (ndim)
        CASE( 4 )  ! assume x,y,z,t variable
             ipmin=imin ; ipmax=imax; jpmin=kmin; jpmax=kmax
             v2d(:,:)=getvarxz(cfilein,cvar,jmin,niz,nkz,imin,kmin)
           lmask(:,:)=.true. ; WHERE ( v2d == 0 ) lmask=.false.
           ilmax=maxloc(v2d,lmask)
           ilmin=minloc(v2d,lmask)
           i1=ilmax(1) ; j1=ilmax(2)
           i2=ilmin(1) ; j2=ilmin(2)
PRINT 9000,' j-slab  MAX:   i    long   j    lat   k     dep    MaxValue    MIN:  i    long   j     lat   k     dep    MinValue'
PRINT 9002, jmin, i1, rlon(i1,1),jmin,rlat(i1,1),j1+kmin-1, h(j1+kmin-1), v2d(i1,j1)*rfact, &
         &             i2, rlon(i2,1),jmin,rlat(i2,1),j2+kmin-1, h(j2+kmin-1), v2d(i2,j2)*rfact
           EXIT
        CASE DEFAULT
           PRINT *,' Non mapable variables x-z or y-z :('
           cvar='none'
        END SELECT
     
     CASE DEFAULT
        PRINT *,' ntype = ',ntype, '  is not defined ' ; STOP
     END SELECT ! ntype
 ENDDO

9000 FORMAT(a)
9001 FORMAT(i4,1x,f7.2,5x,i5,f8.2, i5, f7.2, e14.5, 5x,i5,f8.2, i5, f7.2, e14.5)
9002 FORMAT(i4,9x, i5,f8.2, i5, f7.2, i5, f8.2, e14.5, 6x, i5,f8.2, i5, f7.2, i5, f8.2, e14.5 )
9003 FORMAT(I5, x,i5,1x,f7.2,5x,i5,f8.2, i5, f7.2, e14.5, 5x,i5,f8.2, i5, f7.2, e14.5)

END PROGRAM cdfmax_sp
