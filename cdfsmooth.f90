PROGRAM cdfsmooth
  !!----------------------------------------------------------------------------
  !!        ***   PROGRAM cdfsmooth   ***
  !! 
  !!    ** Purpose : perform a spatial filtering on input file.
  !!             - various filters are available :
  !!               1: Lanczos (default)
  !!               2: hanning
  !!               3: shapiro
  !!              ... : to be completed
  !!
  !!    ** Method : read file level by level and perform a x direction filter, then y direction filter
  !!
  !!   * history:
  !!       Original  : J.M. Molines 1995 for SPEM
  !!             In Dr. Form October 2002
  !!             in cdftools J.M. Molines (July 2007)
  !!----------------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !! * Used moules
  USE cdfio
  IMPLICIT NONE
  !
  INTEGER ::  npiglo, npjglo, npk, nt
  !
  INTEGER :: ji,jj, jk, jt, jvar
  INTEGER :: narg, iargc
  INTEGER :: ncoup, nband
  INTEGER :: nfilter=1
  INTEGER, DIMENSION(:,:), ALLOCATABLE ::  iw    !: flag for bad values (or land masked )
  !
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE ::  v2d,w2d !: raw data,  filtered result
  REAL(KIND=4), DIMENSION(:), ALLOCATABLE   ::  h       !: depth
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE   ::  ec,e    !: weigh in r8, starting index 0 :nband
  REAL(KIND=4) ::  fn, spval
  !
  CHARACTER(LEN=80) :: cfile,cnom, cfilout, cdep, ctim
  ! cdf stuff
  INTEGER    :: nvars, ierr
  INTEGER    :: ncout
  INTEGER , DIMENSION(:), ALLOCATABLE :: id_var , &         !: arrays of var id's
       &                             ipk    , &         !: arrays of vertical level for each var
       &                             id_varout
  REAL(KIND=4), DIMENSION(:), ALLOCATABLE   :: tim
  CHARACTER(LEN=80) ,DIMENSION(:), ALLOCATABLE:: cvarname   !: array of var name

  TYPE (variable), DIMENSION(:), ALLOCATABLE :: typvar
  ! ---


  ! * Initializations
  !
  narg=iargc()
  IF (narg == 0) THEN
     PRINT *,' >>>> cdfsmooth usage : cdfsmooth <filename> n(dx) [ filtertype] '
     PRINT *,'   filename = ncdf file with input data'
     PRINT *,'   n        = number of grid step to filter '
     PRINT *,'   Filtertype = optional argument either '
     PRINT *,'               for Lanczos : Lanczos, l , L '
     PRINT *,'               for Hanning : Hanning, H, h '
     PRINT *,'               for Shapiro : Shapiro, S, s '
     PRINT *,'               for Box : Box, B, b '
     PRINT *,'  output is done on ''filename''.smooth'
     PRINT *,'     where smooth is either L H or S.'
     STOP
  END IF
  !
  CALL getarg(1,cfile)
  CALL getarg(2,cnom)
  READ(cnom,*) ncoup   !  remark: for a spatial filter, fn=dx/lambda where dx is spatial step, lamda is cutting wavelength
  fn=1./ncoup
  nband=NINT(2./fn)    ! Bandwidth of filter is twice the filter span
  ALLOCATE ( ec(0:nband) , e(0:nband) )
  WRITE(cfilout,'(a,a,i3.3)') TRIM(cfile),'L',ncoup   ! default name

  IF ( narg  == 3 ) THEN
     CALL getarg(3,  cnom)
     SELECT CASE ( cnom)
     CASE ( 'Lanczos','L','l') 
        nfilter=1
        WRITE(cfilout,'(a,a,i3.3)') TRIM(cfile),'L',ncoup
        PRINT *,' Working with Lanczos filter'
     CASE ( 'Hanning','H','h')
        nfilter=2
        WRITE(cfilout,'(a,a,i3.3)') TRIM(cfile),'H',ncoup
        PRINT *,' Working with Hanning filter'
     CASE ( 'Shapiro','S','s')
        nfilter=3
        WRITE(cfilout,'(a,a,i3.3)') TRIM(cfile),'S',ncoup
        PRINT *,' Working with Shapiro filter'
     CASE ( 'Box','B','b')
        nfilter=4
        WRITE(cfilout,'(a,a,i3.3)') TRIM(cfile),'B',ncoup
        PRINT *,' Working with Box filter'
     CASE DEFAULT
        PRINT *, TRIM(cnom),' : undefined filter ' ; STOP
     END SELECT
  ENDIF

  CALL filterinit (nfilter, fn,nband)
  ! Look for input file and create outputfile
  npiglo= getdim (cfile,'x')
  npjglo= getdim (cfile,'y')
  npk   = getdim (cfile,'depth',cdtrue=cdep, kstatus=ierr)
  nt    = getdim (cfile,'time',cdtrue=ctim)

  ALLOCATE ( h(npk), v2d(npiglo,npjglo),iw(npiglo,npjglo), w2d(npiglo,npjglo), tim(nt) )
  nvars = getnvar(cfile)
  PRINT *,' nvars =', nvars
  ALLOCATE (cvarname(nvars) )
  ALLOCATE (typvar(nvars) )
  ALLOCATE (id_var(nvars),ipk(nvars),id_varout(nvars) )

  ! get list of variable names and collect attributes in typvar (optional)
  cvarname(:)=getvarname(cfile,nvars,typvar)

  ! ipk gives the number of level or 0 if not a T[Z]YX  variable
  ipk(:)     = getipk (cfile,nvars,cdep=cdep)
  WHERE( ipk == 0 ) cvarname='none'
  typvar(:)%name=cvarname

  ! create output file taking the sizes in cfile
  ncout =create(cfilout, cfile,npiglo,npjglo,npk,cdep=cdep)
   print *,ncout, trim(cfilout)
  ierr= createvar(ncout , typvar,  nvars, ipk, id_varout )
  ierr= putheadervar(ncout , cfile, npiglo, npjglo, npk,cdep=cdep)
  tim=getvar1d(cfile,ctim,nt)
  !
  DO jvar =1,nvars
     IF (cvarname(jvar) == 'nav_lon' .OR. &
          cvarname(jvar) == 'nav_lat' ) THEN
        ! skip these variables
     ELSE
        spval=typvar(jvar)%missing_value
        DO jt=1,nt
           DO jk=1,ipk(jvar)
              v2d(:,:) = getvar(cfile,cvarname(jvar),jk,npiglo,npjglo,ktime=jt)
              iw(:,:) = 1
              WHERE ( v2d == spval ) iw =0
              CALL filter(nfilter,v2d,iw,w2d)
!             CALL lislanczos2d(v2d,iw,w2d,npiglo,npjglo,fn,nband,npiglo,npjglo)
              ierr = putvar(ncout, id_varout(jvar) ,w2d, jk, npiglo, npjglo)
              !
           END DO
        END DO
     ENDIF
  END DO
  ierr=putvar1d(ncout, tim,nt,'T')
  ierr=closeout(ncout)

CONTAINS

  SUBROUTINE filterinit(kfilter, pfn, kband)
    INTEGER, INTENT(in) :: kfilter, kband
    REAL(KIND=4),INTENT(in) :: pfn
    SELECT CASE ( kfilter)
    CASE ( 1 )
       CALL initlanc(pfn,kband)
    CASE ( 2 )
       CALL inithann(pfn,kband)
    CASE ( 3 )
       CALL initshap(pfn, kband)
    CASE ( 4 )
       CALL initbox(pfn, kband)
    END SELECT
  END SUBROUTINE filterinit

  SUBROUTINE filter (kfilter,px,kpx,py)
     INTEGER, INTENT(in) :: kfilter
     INTEGER     , DIMENSION(:,:), INTENT(in) :: kpx
     REAL(KIND=4), DIMENSION(:,:), INTENT(in) :: px
     REAL(KIND=4), DIMENSION(:,:), INTENT(out) :: py

    SELECT CASE ( kfilter)
    CASE ( 1 )
       CALL lislanczos2d(px,kpx,py,npiglo,npjglo,fn,nband,npiglo,npjglo)
    CASE ( 2 )
       print *,' not available'
    CASE ( 3 )
       print *,' not available'
    CASE ( 4 )
       CALL lisbox(px,kpx,py,npiglo,npjglo,fn,nband,npiglo,npjglo)
    END SELECT
  END SUBROUTINE filter

  !!
  SUBROUTINE initlanc(pfn,knj)
    INTEGER, INTENT(in)     :: knj  !: bandwidth
    REAL(KIND=4),INTENT(in) ::  pfn

    ! Local variable
    INTEGER :: ji
    REAL(KIND=8) ::  zpi,zey, zcoef
    !
    !
    zpi=ACOS(-1.)
    zcoef=2*zpi*pfn

    e(0)= 2.*pfn
    DO  ji=1,knj
       e(ji)=SIN(zcoef*ji)/(zpi*ji)
    END DO
    !
    ec(0) = 2*pfn
    DO ji=1,knj
       zey=zpi*ji/knj
       ec(ji)=e(ji)*SIN(zey)/zey
    END DO
    !
  END SUBROUTINE initlanc

  SUBROUTINE inithann(pfn,knj)
    INTEGER, INTENT(in)     :: knj  !: bandwidth
    REAL(KIND=4),INTENT(in) ::  pfn
    PRINT *,' Init hann not done already' ; STOP
  END SUBROUTINE inithann

  SUBROUTINE initshap(pfn,knj)
    INTEGER, INTENT(in)     :: knj  !: bandwidth
    REAL(KIND=4),INTENT(in) ::  pfn
    PRINT *,' Init shap not done already' ; STOP
  END SUBROUTINE initshap

  SUBROUTINE initbox(pfn,knj)
    INTEGER, INTENT(in)     :: knj  !: bandwidth
    REAL(KIND=4),INTENT(in) ::  pfn
    ! dummy init
    ec(:) = 1.
  END SUBROUTINE initbox


  SUBROUTINE lislanczos2d(px,kiw,py,knx,kny,pfn,knj,kpi,kpj)
    !----------------------------------------------
    !	px=input data
    !	kiw = validity of input data
    !	py=output filter
    !	n=number of input/output data
    !	knj= bandwith of the filter
    !	pfn= cutoff frequency
    ! Eric Blayo d'apres une source CLS fournie par F. BLANC. et grosse
    !            optimization.
    !--------------------------------------------
    ! * Arguments
    INTEGER, INTENT(in) :: knx, kny, knj, kpi, kpj
    INTEGER,DIMENSION(:,:),INTENT(in) :: kiw
    REAL(KIND=4), DIMENSION(:,:), INTENT(in) :: px
    REAL(KIND=4), DIMENSION(:,:), INTENT(out) :: py
    REAL(KIND=4), INTENT(in) :: pfn

    ! * local variables
    !
    REAL(KIND=8), DIMENSION(kpi,kpj) ::  ztmpx, ztmpy
    REAL(KIND=8) ::  zyy, zden
    INTEGER :: ji,jj,jmx,jkx
    INTEGER :: ik1x, ik2x,ikkx
    INTEGER :: ifrst=0
    INTEGER :: inxmin,inxmaxi,inymin,inymaxi
    !!
    ! filtering
    inxmin   =  knj
    inxmaxi  =  knx-knj+1
    inymin   =  knj
    inymaxi  =  kny-knj+1
    PRINT *,' filtering parameters'
    PRINT *,'    nx=',knx
    PRINT *,'    nband=',knj
    PRINT *,'    fn=',pfn
    DO jj=1,kny
       DO  jmx=1,knx
          ik1x = -knj
          ik2x =  knj
          !
          IF (jmx <= inxmin)  ik1x = 1-jmx
          IF (jmx >= inxmaxi) ik2x = knx-jmx
          !
          zyy  = 0.d0
          zden = 0.d0
          !
          DO jkx=ik1x,ik2x
             ikkx=ABS(jkx)
             IF (kiw(jkx+jmx,jj)  ==  1) THEN
                zden=zden+ec(ikkx)
                zyy=zyy+ec(ikkx)*px(jkx+jmx,jj)
             END IF
          END DO
          !
          ztmpx(jmx,jj)=zyy/zden
       END DO
    END DO

    DO ji=1,knx
       DO  jmx=1,kny
          ik1x = -knj
          ik2x =  knj
          !
          IF (jmx <= inymin) ik1x=1-jmx
          IF (jmx >= inymaxi) ik2x=kny-jmx
          !
          zyy=0.d0
          zden=0.d0
          !
          DO jkx=ik1x,ik2x
             ikkx=ABS(jkx)
             IF (kiw(ji,jkx+jmx)  ==  1) THEN
                zden=zden+ec(ikkx)
                zyy=zyy+ec(ikkx)*ztmpx(ji,jkx+jmx)
!               zyy=zyy+ec(ikkx)*px(ji,jkx+jmx)
             END IF
          END DO
!         ztmpy(ji,jmx)=zyy/zden
          py(ji,jmx)=zyy/zden
       END DO
    END DO
!     py=0.5*(ztmpx + ztmpy )
    !
  END SUBROUTINE lislanczos2d

  SUBROUTINE lisbox(px,kiw,py,knx,kny,pfn,knj,kpi,kpj)
    !----------------------------------------------
    ! perform a box car 2d filtering, of span knj
    !----------------------------------------------
    ! * Arguments
    INTEGER, INTENT(in) :: knx, kny, knj, kpi, kpj
    INTEGER,DIMENSION(:,:),INTENT(in) :: kiw
    REAL(KIND=4), DIMENSION(:,:), INTENT(in) :: px
    REAL(KIND=4), DIMENSION(:,:), INTENT(out) :: py
    REAL(KIND=4), INTENT(in) :: pfn

    ! Local vaariables
    INTEGER :: ji,jj, ik1x, ik2x, ik1y, ik2y
    REAL(KIND=8) :: den
    LOGICAL, DIMENSION(kpi,kpj) :: lmask

    lmask=.true.
    WHERE (kiw == 0 ) lmask=.false.
    DO ji=1,knx
        ik1x=ji-knj ;  ik2x=ji+knj
        ik1x=MAX(1,ik1x)  ; ik2x=MIN(knx,ik2x)
      DO jj=1,kny
        ik1y=jj-knj ;  ik2y=jj+knj
        ik1y=MAX(1,ik1y)  ; ik2y=MIN(kny,ik2y)
        den=SUM(kiw(ik1x:ik2x,ik1y:ik2y) )
        IF ( den /= 0 ) THEN
          py(ji,jj)= SUM(px(ik1x:ik2x,ik1y:ik2y), mask=lmask(ik1x:ik2x,ik1y:ik2y) )/den
        ELSE
          py(ji,jj) = spval
        ENDIF
      END DO
    END DO

   END SUBROUTINE lisbox

END PROGRAM cdfsmooth
