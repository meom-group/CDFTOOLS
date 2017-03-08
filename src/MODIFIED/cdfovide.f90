PROGRAM cdfovide
  !!======================================================================
  !!                     ***  PROGRAM  cdfovide  ***
  !!=====================================================================
  !!  ** Purpose : Easy tool to perform Temperature, Salinity and velocity
  !!               plots along OVIDE section
  !!               PARTIAL STEPS version
  !!
  !!  ** Method  : Works as a standalone file once compiled
  !!               Inspired by cdffindij, cdftransportiz
  !!
  !! History : 2.1  : 12/2009  : R. Dussin    : Original code
  !!         :  4.0  : 03/2017  : J.M. Molines  
  !!           3.0  : 01/2011  : J.M. Molines : Doctor norm + Lic.
  !!           3.1  : 05/2013  : T. Penduff & R. Dussin  : Saving new variables
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class obsolete
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4) :: narg, iargc         ! command line
  INTEGER(KIND=4) :: npiglo, npjglo, npk ! size of the domain
  INTEGER(KIND=4) :: niter
  INTEGER(KIND=4) :: imin, imax, jmin, jmax, k, ik, jk, jclass
  INTEGER(KIND=4) :: iloc, jloc
  INTEGER(KIND=4) :: iloop, jloop

  INTEGER(KIND=4) :: nsec=0 ! nb total de points le long de la section
  INTEGER(KIND=4), DIMENSION (:), ALLOCATABLE :: isec, jsec ! indices des points a recuperer
  
!	R. Dussin's initial choice 
!  INTEGER(KIND=4), PARAMETER :: nsta=4
!	IFREMER (D. Desbruyeres') choice
  INTEGER(KIND=4), PARAMETER :: nsta=5

  INTEGER(KIND=4), DIMENSION(nsta) :: ista, jsta
  INTEGER(KIND=4), DIMENSION(nsta-1) :: ikeepn

  ! broken line stuff
  INTEGER(KIND=4), PARAMETER :: jpseg=10000
  INTEGER(KIND=4) :: i0,j0,i1,j1, i, j
  INTEGER(KIND=4) :: n,nn, jseg, kk
  INTEGER(KIND=4) :: norm_u, norm_v, ist, jst
  INTEGER(KIND=4)    ::  nxtarg
  INTEGER(KIND=4), DIMENSION(nsta-1,jpseg) :: legs1=0, legs2=0

  REAL(KIND=8), DIMENSION(nsta)             :: rlonsta, rlatsta
  REAL(KIND=8)                              :: xmin, xmax, ymin, ymax, rdis
  REAL(KIND=4)                              :: glamfound, glamin, glamax
  REAL(KIND=8)                              :: glam0, emax
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: glam, gphi, e1, e2
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: e1t, e2t, e1u, e2v, e3t
  
  REAL(KIND=4) :: rxi0, ryj0, rxi1, ryj1
  REAL(KIND=4) :: ai, bi, aj,bj,d
  REAL(KIND=4) :: rxx(jpseg), ryy(jpseg)
  REAL(KIND=4), DIMENSION(jpseg) :: gla !, gphi

  REAL(KIND=8), DIMENSION(jpseg) :: voltrp, heatrp, saltrp
  REAL(KIND=8)                   :: voltrpsum, heatrpsum, saltrpsum
  COMPLEX yypt(jpseg), yypti

  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: e1v, e3v ,gphiv, zv, zvt, zvs !: mask, metrics
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: e2u, e3u ,gphiu, zu, zut, zus !: mask, metrics
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: temper, saline, zonalu, meridv, navlon, navlat
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: ovidetemper, ovidesaline, ovidezonalu, ovidemeridv
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: lonsec, latsec, dumisec, dumjsec
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: e1tsec, e1usec, e1vsec, e2tsec, e2usec, e2vsec
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: e3tsec, e3usec, e3vsec
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: e2join, e3join, dummyvmask
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: glamu, glamv
  REAL(KIND=4), DIMENSION(:),     ALLOCATABLE :: gdepw

  REAL(KIND=8)                              :: tmp,barot

! constants
  REAL(KIND=4)   ::  rau0=1000.,  rcp=4000.

  CHARACTER(LEN=80) :: ctype='F'
  CHARACTER(LEN=80) :: cfilet , cfileu, cfilev, csection 
  LOGICAL  :: lagain, lbord
  LOGICAL  :: ltest=.FALSE.
  LOGICAL  :: lchk
  
! cdf output stuff
  CHARACTER(LEN=80) :: cfileoutnc='ovide.nc'
  TYPE (variable), DIMENSION(:), ALLOCATABLE   :: stypvar
  INTEGER(KIND=4)    :: ierr, ncout
  REAL(KIND=4), DIMENSION(1)                    ::  tim
  INTEGER(KIND=4)    :: nfield=14
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: ipk, id_varout
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg < 3  ) THEN
     PRINT *,'usage : cdfovide gridTfile gridUfile gridVfile '
     PRINT *,'     Files ',TRIM(cn_fhgr),' and ',TRIM(cn_fzgr),' must be in te current directory '
     PRINT *,'     Output on netcdf file ',TRIM(cfileoutnc)
     STOP
  ENDIF

  CALL getarg (1, cfilet)
  CALL getarg (2, cfileu)
  CALL getarg (3, cfilev)

  lchk = chkfile(cn_fhgr)
  lchk = chkfile(cn_fzgr) .OR. lchk
  lchk = chkfile(cfilet ) .OR. lchk
  lchk = chkfile(cfileu ) .OR. lchk
  lchk = chkfile(cfilev ) .OR. lchk
  IF ( lchk ) STOP ! missing files

  ! R. Dussin : Location of leg points that define the 3 legs of OVIDE section
  !rlonsta(1) = -43.00 ; rlatsta(1) = 60.60    ! Greenland
  !rlonsta(2) = -31.30 ; rlatsta(2) = 58.90    ! Reykjanes Ridge
  !rlonsta(3) = -12.65 ; rlatsta(3) = 40.33    ! Off Portugal
  !rlonsta(4) =  -8.70 ; rlatsta(4) = 40.33    ! Lisboa

  ! D. Desbruyeres : Location of leg points that define the 4 legs of the OVIDE section
  rlonsta(1) = -43.70 ; rlatsta(1) = 59.90    ! 
  rlonsta(2) = -30.30 ; rlatsta(2) = 58.90    ! 
  rlonsta(3) = -19.40 ; rlatsta(3) = 44.90    ! 
  rlonsta(4) = -12.65 ; rlatsta(4) = 40.33    ! 
  rlonsta(5) = -08.70 ; rlatsta(5) = 40.33    ! 



  PRINT *, '###########################################################'
  PRINT *, '#                                                          '
  PRINT *, '#   CDF ovide                                              '
  PRINT *, '#                                                          '
  PRINT *, '#                                      \___________________'
  PRINT *, '#                                       \                  '
  PRINT *, '#                                        \                 '
  PRINT *, '#                                         \________________'
  PRINT *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

  !!---------------------------------------------------------------------
  !!  Find the indexes of the legs (from cdffindij) 
  !!---------------------------------------------------------------------

  npiglo = getdim (cn_fhgr,cn_x)
  npjglo = getdim (cn_fhgr,cn_y)
  
  ALLOCATE (glam(npiglo,npjglo), gphi(npiglo,npjglo) )
  ALLOCATE (e1(npiglo,npjglo), e2(npiglo,npjglo) )

  SELECT CASE ( ctype )
  CASE ('T' , 't' )
     glam(:,:) = getvar(cn_fhgr, 'glamt',1,npiglo,npjglo)
     gphi(:,:) = getvar(cn_fhgr, 'gphit',1,npiglo,npjglo)
     e1  (:,:) = getvar(cn_fhgr, 'e1t'  ,1,npiglo,npjglo)
     e2  (:,:) = getvar(cn_fhgr, 'e2t'  ,1,npiglo,npjglo)
  CASE ('U','u' )
     glam(:,:) = getvar(cn_fhgr, 'glamu',1,npiglo,npjglo)
     gphi(:,:) = getvar(cn_fhgr, 'gphiu',1,npiglo,npjglo)
     e1  (:,:) = getvar(cn_fhgr, 'e1u'  ,1,npiglo,npjglo)
     e2  (:,:) = getvar(cn_fhgr, 'e2u'  ,1,npiglo,npjglo)
  CASE ('V','v' )
     glam(:,:) = getvar(cn_fhgr, 'glamv',1,npiglo,npjglo)
     gphi(:,:) = getvar(cn_fhgr, 'gphiv',1,npiglo,npjglo)
     e1  (:,:) = getvar(cn_fhgr, 'e1v'  ,1,npiglo,npjglo)
     e2  (:,:) = getvar(cn_fhgr, 'e2v'  ,1,npiglo,npjglo)
  CASE ('F','f' )
     glam(:,:) = getvar(cn_fhgr, 'glamf',1,npiglo,npjglo)
     gphi(:,:) = getvar(cn_fhgr, 'gphif',1,npiglo,npjglo)
     e1  (:,:) = getvar(cn_fhgr, 'e1f'  ,1,npiglo,npjglo)
     e2  (:,:) = getvar(cn_fhgr, 'e2f'  ,1,npiglo,npjglo)
  CASE DEFAULT
     PRINT *,' ERROR : type of point not known: ', TRIM(ctype)
  END SELECT

  ! work with longitude between 0 and 360 to avoid  the date line.
    WHERE( glam < 0 ) glam=glam+360.
  ! For Orca grid, the longitude of ji=1 is about 70 E
  glam0=glam(1, npjglo/2)
  WHERE( glam < glam0 ) glam=glam+360.

  !! loop on the legs
  DO k = 1,nsta-1

   xmin=rlonsta(k)
   ymin=rlatsta(k)
   xmax=rlonsta(k+1)
   ymax=rlatsta(k+1)

   IF (xmin < 0.) xmin = xmin +360.
   IF (xmax < 0.) xmax = xmax +360.
    
   IF (xmin < glam0) xmin = xmin +360.
   IF (xmax < glam0) xmax = xmax +360.

   lagain = .TRUE.
   niter = 0

  !! --- while loop -----------------------------------------------------------
  DO WHILE (lagain)
     CALL Nearestpoint(xmin,ymin,npiglo,npjglo,glam,gphi,iloc,jloc,lbord)
     ! distance between the target point and the nearest point
     rdis=dist(xmin,glam(iloc,jloc),ymin,gphi(iloc,jloc) ) ! in km
     ! typical grid size (diagonal) in the vicinity of nearest point
     emax= MAX(e1(iloc,jloc),e2(iloc,jloc))/1000.*SQRT(2.) ! in km

!!    rdis = (xmin - glam(iloc,jloc))**2 + (ymin - gphi(iloc,jloc))**2
!!    rdis = SQRT(rdis)
     IF (rdis  > emax ) THEN
         glamfound=glam(iloc,jloc) ; IF (glamfound > 180.)  glamfound=glamfound -360.
        PRINT 9000, 'Long= ',glamfound,' Lat = ',gphi(iloc,jloc)&
             &               , iloc, jloc 
        PRINT *,' Algorithme ne converge pas ', rdis 
        IF ( niter >=  1 ) STOP ' pas de convergence apres iteration'
        lagain = .TRUE.
        jloc = npjglo
        niter = niter +1
     ELSE
        !PRINT '("#  rdis= ",f8.3," km")', rdis
        lagain = .FALSE.
     END IF
  END DO
  !!--------------------------------------------------------------------------

  IF (lbord) THEN
     WRITE (*,*)'Point  Out of domain or on boundary'
  ELSE
     imin=iloc
     jmin=jloc
     !      PRINT 9000, 'Long= ',glam(iloc,jloc),' lat = ',gphi(iloc,jloc), iloc, jloc 
  ENDIF
  
  lagain = .TRUE.
  niter = 0
  !! --- while loop ----------------------------------------------------------------
  DO WHILE (lagain)
     CALL Nearestpoint(xmax,ymax,npiglo,npjglo,glam,gphi,iloc,jloc,lbord)
     ! distance between the target point and the nearest point
     rdis=dist(xmax,glam(iloc,jloc),ymax,gphi(iloc,jloc) ) ! in km
     ! typical grid size (diagonal) in the vicinity of nearest point
     emax= MAX(e1(iloc,jloc),e2(iloc,jloc))/1000.*SQRT(2.) ! in km
!!    rdis = (xmax - glam(iloc,jloc))**2 + (ymax - gphi(iloc,jloc))**2
!!    rdis = SQRT(rdis)
     IF (rdis >  emax ) THEN
         glamfound=glam(iloc,jloc) ; IF (glamfound > 180.)  glamfound=glamfound -360.
        PRINT 9000, 'Long= ',glamfound,' Lat = ',gphi(iloc,jloc) &
             &               , iloc, jloc
        PRINT *,' Algorithme ne converge pas ', rdis
        IF ( niter >= 1 ) STOP ' pas de convergence avres iteration'
        lagain = .TRUE.
        jloc  = npjglo
        niter = niter +1
     ELSE
        !PRINT '("#  rdis= ",f8.3," km")', rdis
        lagain = .FALSE.
     END IF
  END DO !! ---------------------------------------------------------------------

  IF (lbord) THEN
     WRITE (*,*) 'Point  Out of domain or on boundary'
  ELSE
     imax=iloc
     jmax=jloc
     !      PRINT 9000, 'Long= ',glam(iloc,jloc),' lat = ',gphi(iloc,jloc), iloc, jloc
  ENDIF
  
  ista(k)=imin
  jsta(k)=jmin
  ista(k+1)=imax
  jsta(k+1)=jmax

!  PRINT 9001, imin,imax, jmin, jmax
!  glamin=glam(imin,jmin) ;glamax=glam(imax,jmax)
!  IF ( glamin > 180 ) glamin=glamin-360.
!  IF ( glamax > 180 ) glamax=glamax-360.
!  PRINT 9002, glamin, glamax, gphi(imin,jmin),gphi(imax,jmax)
9000 FORMAT(a,f8.2,a,f8.2,2i5)
9001 FORMAT(4i10)
9002 FORMAT(4f10.4)


     !! Find the broken line between P1 (imin,jmin) and P2 (imax, jmax)
     !! ---------------------------------------------------------------
     ! ... Initialization
     i0=imin; j0=jmin; i1=imax;  j1=jmax
     rxi1=i1;  ryj1=j1; rxi0=i0; ryj0=j0

     ! .. Compute equation:  ryj = aj rxi + bj
     IF ( (rxi1 -rxi0) /=  0 ) THEN
        aj = (ryj1 - ryj0 ) / (rxi1 -rxi0)
        bj = ryj0 - aj * rxi0
     ELSE
        aj=10000.
        bj=0.
     END IF


     ! .. Compute equation:  rxi = ai ryj + bi
     IF ( (ryj1 -ryj0) /=  0 ) THEN
        ai = (rxi1 - rxi0 ) / ( ryj1 -ryj0 )
        bi = rxi0 - ai * ryj0
     ELSE
        ai=10000.
        bi=0.
     END IF


     ! ..  Compute the integer pathway:
     n=0
     ! .. Chose the strait line with the smallest slope
     IF (ABS(aj) <=  1 ) THEN
        ! ... Here, the best line is y(x)
        ! ... If i1 < i0 swap points and remember it has been swapped
        IF (i1 <  i0 ) THEN
           i  = i0 ; j  = j0
           i0 = i1 ; j0 = j1
           i1 = i  ; j1 = j
        END IF

        IF ( j1 >= j0 ) THEN
           ist = 1     ; jst = 1
           norm_u =  1 ;  norm_v = -1
        ELSE
           ist = 1     ; jst = 0
           norm_u = -1 ; norm_v = -1
        END IF

        ! ... compute the nearest j point on the line crossing at i
        DO i=i0,i1
           n=n+1
           IF (n > jpseg) STOP 'n > jpseg !'
           j=NINT(aj*i + bj )
           yypt(n) = CMPLX(i,j)
        END DO
     ELSE
        ! ... Here, the best line is x(y)
        ! ... If j1 < j0 swap points and remember it has been swapped
        IF (j1 <  j0 ) THEN
           i  = i0 ; j  = j0
           i0 = i1 ; j0 = j1
           i1 = i  ; j1 = j
        END IF
        IF ( i1 >=  i0 ) THEN
           ist = 1    ;  jst = 1
           norm_u = 1 ;  norm_v = -1
        ELSE
           ist = 0
           jst = 1
           norm_u = 1
           norm_v = 1
        END IF

        ! ... compute the nearest i point on the line crossing at j
        DO j=j0,j1
           n=n+1
           IF (n > jpseg) STOP 'n>jpseg !'
           i=NINT(ai*j + bi)
           yypt(n) = CMPLX(i,j)
        END DO
     END IF

     !!
     !! Look for intermediate points to be added.
     !  ..  The final positions are saved in rxx,ryy
     rxx(1)=REAL(yypt(1))
     ryy(1)=IMAG(yypt(1))
     nn=1

     DO kk=2,n
        ! .. distance between 2 neighbour points
        d=ABS(yypt(kk)-yypt(kk-1))
        ! .. intermediate points required if d > 1
        IF ( d > 1 ) THEN
           CALL interm_pt(yypt,kk,ai,bi,aj,bj,yypti)
           nn=nn+1
           IF (nn > jpseg) STOP 'nn>jpseg !'
           rxx(nn)=REAL(yypti)
           ryy(nn)=IMAG(yypti)
        END IF
        nn=nn+1
        IF (nn > jpseg) STOP 'nn>jpseg !'
        rxx(nn)=REAL(yypt(kk))
        ryy(nn)=IMAG(yypt(kk))
     END DO

  IF (rxx(1) < rxx(nn) ) THEN ;
     legs1(k,:)=rxx
     legs2(k,:)=ryy
  ELSE

    DO iloop=1,nn
      legs1(k,iloop)=rxx(nn-iloop+1)
      legs2(k,iloop)=ryy(nn-iloop+1)
   END DO
  END IF

  ! compute the number of total points
  ikeepn(k)=nn
  nsec = nsec + nn
  END DO !! loop on the legs

! fancy control print
WRITE(*,*) '------------------------------------------------------------'
WRITE(*,9100) 'leg 1 start at ', rlonsta(1) ,'°N ', rlatsta(1), '°W and ends at ', rlonsta(2) ,'°N ', rlatsta(2), '°W'
WRITE(*,9101) 'corresponding to F-gridpoints(', ista(1),',',jsta(1),') and (', ista(2),',',jsta(2),')' 
WRITE(*,*) '------------------------------------------------------------'
WRITE(*,*) '------------------------------------------------------------'
WRITE(*,9100) 'leg 2 start at ', rlonsta(2) ,'°N ', rlatsta(2), '°W and ends at ', rlonsta(3) ,'°N ', rlatsta(3), '°W'
WRITE(*,9101) 'corresponding to F-gridpoints(', ista(2),',',jsta(2),') and (', ista(3),',',jsta(3),')' 
WRITE(*,*) '------------------------------------------------------------'
WRITE(*,*) '------------------------------------------------------------'
WRITE(*,9100) 'leg 3 start at ', rlonsta(3) ,'°N ', rlatsta(3), '°W and ends at ', rlonsta(4) ,'°N ', rlatsta(4), '°W'
WRITE(*,9101) 'corresponding to F-gridpoints(', ista(3),',',jsta(3),') and (', ista(4),',',jsta(4),')' 
WRITE(*,*) '------------------------------------------------------------'
WRITE(*,*) '------------------------------------------------------------'
WRITE(*,9100) 'leg 4 start at ', rlonsta(4) ,'°N ', rlatsta(4), '°W and ends at ', rlonsta(5) ,'°N ', rlatsta(5), '°W'
WRITE(*,9101) 'corresponding to F-gridpoints(', ista(4),',',jsta(4),') and (', ista(5),',',jsta(5),')' 
WRITE(*,*) '------------------------------------------------------------'

9100 FORMAT(a,f6.2,a,f6.2,a,f6.2,a,f6.2,a)
9101 FORMAT(a,i4,a,i4,a,i4,a,i4,a)

  ALLOCATE (isec(nsec), jsec(nsec)) 

  DO k=1, nsta-1
    DO iloop=1, ikeepn(k) 
     jloop=iloop + SUM(ikeepn(1:k)) -ikeepn(k)
     isec(jloop)=legs1(k,iloop)
     jsec(jloop)=legs2(k,iloop)
    END DO
  END DO

  npk = getdim (cfilet,'deptht')

  ! input fields
  ALLOCATE(navlon(npiglo,npjglo), navlat(npiglo,npjglo))
  ALLOCATE(temper(npiglo,npjglo), saline(npiglo,npjglo))
  ALLOCATE(zonalu(npiglo,npjglo), meridv(npiglo,npjglo))
  ALLOCATE(e1v(npiglo,npjglo))
  ALLOCATE(e2u(npiglo,npjglo))
  ALLOCATE(e3u(npiglo,npjglo), e3v(npiglo,npjglo))

  ! output fields
  ALLOCATE(lonsec(1,nsec), latsec(1,nsec) )
  ALLOCATE(dumisec(1,nsec), dumjsec(1,nsec) )
  ALLOCATE(e2usec(1,nsec-1), e3usec(nsec-1,npk) )
  ALLOCATE(e1vsec(1,nsec-1), e3vsec(nsec-1,npk) )
  ALLOCATE(ovidetemper(nsec-1,npk), ovidesaline(nsec-1,npk) )
  ALLOCATE(ovidezonalu(nsec-1,npk), ovidemeridv(nsec-1,npk) )
  ALLOCATE(dummyvmask(nsec-1,npk))
    
  
  e1vsec=-9999.
  e2usec=-9999.
  
  
  dumisec(:,:)=0
  dumjsec(:,:)=0
  

  navlon(:,:) = getvar(cfilet, 'nav_lon' ,1,npiglo,npjglo)
  navlat(:,:) = getvar(cfilet, 'nav_lat' ,1,npiglo,npjglo)
  e1v(:,:)    = getvar(cn_fhgr, 'e1v',1,npiglo,npjglo)
  e2u(:,:)    = getvar(cn_fhgr, 'e2u',1,npiglo,npjglo)

  ! il faut faire un test sur la continuité des segements
  ! on va prendre T et S comme etant la moyenne du point
  ! en dessous et au-dessus du segment pour pouvoir calculer
  ! les fluxs de maniere optimales...

  ! loop on 2d arrays
  DO iloop=1,nsec

    lonsec(1,iloop) = navlon(isec(iloop),jsec(iloop))
    latsec(1,iloop) = navlat(isec(iloop),jsec(iloop))
    dumisec(1,iloop)= isec(iloop)
    dumjsec(1,iloop)= jsec(iloop)

  END DO

	jk=1

  DO iloop=1,nsec-1
	 !PRINT*, 'iloop=', iloop
	
  IF ( jsec(iloop+1) == jsec(iloop) ) THEN ! horizontal segment
    IF ( isec(iloop+1) > isec(iloop) ) THEN ! eastward

       e2usec(jk,iloop) = 0.
       e1vsec(jk,iloop) = e1v(isec(iloop)+1,jsec(iloop))

    ELSE

       e2usec(jk,iloop) = 0.
       e1vsec(jk,iloop) = e1v(isec(iloop),jsec(iloop))

    ENDIF
  ELSEIF ( isec(iloop+1) == isec(iloop) ) THEN ! vertical segment
    IF ( jsec(iloop+1) < jsec(iloop) ) THEN ! southward

       e2usec(jk,iloop) = e2u(isec(iloop),jsec(iloop))
       e1vsec(jk,iloop) = 0.

    ELSE

       e2usec(jk,iloop) = e2u(isec(iloop),jsec(iloop)+1)
       e1vsec(jk,iloop) = 0.

    ENDIF
  ELSE
       PRINT *, 'problem'
       exit 
  ENDIF
  END DO

!	PRINT*,nsec
!	PRINT*, MINVAL(e1v),MAXVAL(e1v),MINVAL(e2u),MAXVAL(e2u)  
!	PRINT*, MINVAL(e1vsec),MAXVAL(e1vsec),MINVAL(e2usec),MAXVAL(e2usec)  
!	PAUSE		
	
  ! loop on 3d arrays
  DO jk=1,npk
  temper(:,:) = getvar(cfilet, 'votemper',jk,npiglo,npjglo)
  saline(:,:) = getvar(cfilet, 'vosaline',jk,npiglo,npjglo)
  zonalu(:,:) = getvar(cfileu, 'vozocrtx',jk,npiglo,npjglo)
  meridv(:,:) = getvar(cfilev, 'vomecrty',jk,npiglo,npjglo)
  e3u(:,:) = getvar(cn_fzgr, 'e3u_ps',jk,npiglo,npjglo, ldiom=.true.)
  e3v(:,:) = getvar(cn_fzgr, 'e3v_ps',jk,npiglo,npjglo, ldiom=.true.)

  DO iloop=1,nsec-1
    IF ( jsec(iloop+1) == jsec(iloop) ) THEN ! horizontal segment
       IF ( isec(iloop+1) > isec(iloop) ) THEN ! eastward

       IF ( min( temper(isec(iloop)+1,jsec(iloop)) , temper(isec(iloop)+1,jsec(iloop)+1) ) == 0. ) THEN
       ovidetemper(iloop,jk) = 0. ; ovidesaline(iloop,jk) = 0.
       ELSE
       ovidetemper(iloop,jk) = 0.5 * ( temper(isec(iloop)+1,jsec(iloop)) + temper(isec(iloop)+1,jsec(iloop)+1) )
       ovidesaline(iloop,jk) = 0.5 * ( saline(isec(iloop)+1,jsec(iloop)) + saline(isec(iloop)+1,jsec(iloop)+1) )
       ENDIF
       ovidezonalu(iloop,jk) = 0.
       ovidemeridv(iloop,jk) = meridv(isec(iloop)+1,jsec(iloop))
       e3usec(iloop,jk) = 0.
       e3vsec(iloop,jk) = e3v(isec(iloop)+1,jsec(iloop))

       ELSE ! westward

       IF ( min( temper(isec(iloop),jsec(iloop)) , temper(isec(iloop),jsec(iloop)+1) ) == 0. ) THEN
       ovidetemper(iloop,jk) = 0. ; ovidesaline(iloop,jk) = 0.
       ELSE
       ovidetemper(iloop,jk) = 0.5 * ( temper(isec(iloop),jsec(iloop)) + temper(isec(iloop),jsec(iloop)+1) )
       ovidesaline(iloop,jk) = 0.5 * ( saline(isec(iloop),jsec(iloop)) + saline(isec(iloop),jsec(iloop)+1) )
       ENDIF
       ovidezonalu(iloop,jk) = 0.
       ovidemeridv(iloop,jk) = meridv(isec(iloop),jsec(iloop))
       e3usec(iloop,jk) = 0.
       e3vsec(iloop,jk) = e3v(isec(iloop),jsec(iloop))

       ENDIF
    ELSEIF ( isec(iloop+1) == isec(iloop) ) THEN ! vertical segment
       IF ( jsec(iloop+1) < jsec(iloop) ) THEN ! southward

       IF ( min( temper(isec(iloop),jsec(iloop)) , temper(isec(iloop)+1,jsec(iloop)) ) == 0. ) THEN
       ovidetemper(iloop,jk) = 0. ; ovidesaline(iloop,jk) = 0.
       ELSE
       ovidetemper(iloop,jk) = 0.5 * ( temper(isec(iloop),jsec(iloop)) + temper(isec(iloop)+1,jsec(iloop)) )
       ovidesaline(iloop,jk) = 0.5 * ( saline(isec(iloop),jsec(iloop)) + saline(isec(iloop)+1,jsec(iloop)) )
       ENDIF
       ovidezonalu(iloop,jk) = zonalu(isec(iloop),jsec(iloop))
       ovidemeridv(iloop,jk) = 0.
       e3usec(iloop,jk) = e3u(isec(iloop),jsec(iloop))
       e3vsec(iloop,jk) = 0.

       ELSE ! northward

       IF ( min( temper(isec(iloop),jsec(iloop)+1) , temper(isec(iloop)+1,jsec(iloop)+1) ) == 0. ) THEN
       ovidetemper(iloop,jk) = 0. ; ovidesaline(iloop,jk) = 0.
       ELSE
       ovidetemper(iloop,jk) = 0.5 * ( temper(isec(iloop),jsec(iloop)+1) + temper(isec(iloop)+1,jsec(iloop)+1) )
       ovidesaline(iloop,jk) = 0.5 * ( saline(isec(iloop),jsec(iloop)+1) + saline(isec(iloop)+1,jsec(iloop)+1) )
       ENDIF
       ovidezonalu(iloop,jk) = zonalu(isec(iloop),jsec(iloop)+1)
       ovidemeridv(iloop,jk) = 0.
       e3usec(iloop,jk) = e3u(isec(iloop),jsec(iloop)+1)
       e3vsec(iloop,jk) = 0.

       ENDIF

    ELSE
       PRINT *, 'problem'
       exit 
    ENDIF

  END DO
  END DO


  ALLOCATE ( stypvar(nfield), ipk(nfield), id_varout(nfield) )

  DO iloop=1,nfield
     ipk(iloop) = npk
  END DO

 ! define new variables for output 
  stypvar(1)%cname= 'votemper'
  stypvar(1)%cunits='deg C'
  stypvar%rmissing_value=0.
  stypvar(1)%valid_min= -2.
  stypvar(1)%valid_max= 40.
  stypvar%scale_factor= 1.
  stypvar%add_offset= 0.
  stypvar%savelog10= 0.
  stypvar(1)%clong_name='Temperature along OVIDE section'
  stypvar(1)%cshort_name='votemper'
  stypvar%conline_operation='N/A'
  stypvar%caxis='TXZ'

  stypvar(2)%cname= 'vosaline'
  stypvar(2)%cunits='PSU'
  stypvar(2)%valid_min= 0.
  stypvar(2)%valid_max= 50.
  stypvar(2)%clong_name='Salinity along OVIDE section'
  stypvar(2)%cshort_name='vosaline'

  stypvar(3)%cname= 'vozocrtx_native'
  stypvar(3)%cunits='m.s-1'
  stypvar(3)%valid_min= -20.
  stypvar(3)%valid_max= 20.
  stypvar(3)%clong_name='Zonal velocity along OVIDE section'
  stypvar(3)%cshort_name='vozocrtx'

  stypvar(4)%cname= 'vomecrty_native'
  stypvar(4)%cunits='m.s-1'
  stypvar(4)%valid_min= -20.
  stypvar(4)%valid_max= 20.
  stypvar(4)%clong_name='Meridionnal velocity along OVIDE section'
  stypvar(4)%cshort_name='vomecrty'

  stypvar(5)%cname= 'isec'
  stypvar(5)%valid_min= 0.
  stypvar(5)%valid_max= npiglo 
  stypvar(6)%cname= 'jsec'
  stypvar(6)%valid_min= 0.
  stypvar(6)%valid_max= npjglo 
  stypvar(7)%cname= 'e2u_native'
  stypvar(7)%valid_min= MINVAL(e2usec(1,:))
  stypvar(7)%valid_max= MAXVAL(e2usec(1,:)) 
  stypvar(8)%cname= 'e1v_native'
  stypvar(8)%valid_min= MINVAL(e1vsec(1,:))
  stypvar(8)%valid_max= MAXVAL(e1vsec(1,:))
  stypvar(9)%cname= 'e3u_native'
  stypvar(9)%valid_min= MINVAL(e3usec(:,:))
  stypvar(9)%valid_max= MAXVAL(e3usec(:,:)) 
  stypvar(10)%cname= 'e3v_native'
  stypvar(10)%valid_min= MINVAL(e3vsec(:,:))
  stypvar(10)%valid_max= MAXVAL(e3vsec(:,:)) 

  stypvar(11)%cname= 'vomecrty'
  stypvar(11)%cunits='m.s-1'
  stypvar(11)%valid_min= -20.
  stypvar(11)%valid_max= 20.
  stypvar(11)%clong_name='Normal velocity along OVIDE section'
  stypvar(11)%cshort_name='vomecrty'

  stypvar(12)%cname= 'e1v'
  stypvar(12)%cunits='m'
  stypvar(12)%valid_min= 0.
  stypvar(12)%valid_max= 1000000.
  stypvar(12)%clong_name='Local horiz. resolution along OVIDE section'
  stypvar(12)%cshort_name='e1v'

  stypvar(13)%cname= 'e3v_ps'
  stypvar(13)%cunits='m'
  stypvar(13)%valid_min= 0.
  stypvar(13)%valid_max= 100000000.
  stypvar(13)%clong_name='Local vert. resolution along OVIDE section'
  stypvar(13)%cshort_name='e3v_ps'

  stypvar(14)%cname= 'vmask'
  stypvar(12)%cunits=''
  stypvar(12)%valid_min= 0.
  stypvar(12)%valid_max= 1.
  stypvar(12)%clong_name='Mask along OVIDE section'
  stypvar(12)%cshort_name='vmask'

  ! create output fileset
   ncout =create(cfileoutnc, 'none', nsec,1, npk,cdep='depthw')
   ierr= createvar(ncout ,stypvar,nfield, ipk,id_varout )
   ierr= putheadervar(ncout, cfilet,nsec,1, npk,pnavlon=lonsec,pnavlat=latsec,pdep=gdepw)
   tim=getvar1d(cfilet,'time_counter',1)
   ierr=putvar1d(ncout,tim,1,'T')


  dummyvmask(:,:) = 1.
  WHERE( ovidesaline(:,:) == 0. ) dummyvmask(:,:) = 0.


	!PRINT*, MINVAL(e1v),MAXVAL(e1v),MINVAL(e2u),MAXVAL(e2u)  
	!PRINT*, MINVAL(e1vsec),MAXVAL(e1vsec),MINVAL(e2usec),MAXVAL(e2usec)  
	!PAUSE

!------------------- BAROTROPIC TRANSPORT
    barot = 0.
    
  DO iloop=1,nsec-1
  DO jk=1,npk
	tmp=(ovidezonalu(iloop,jk)+ovidemeridv(iloop,jk))*&
	    (e2usec(1,iloop)+e1vsec(1,iloop))*&
	    (e3usec(iloop,jk)+e3vsec(iloop,jk))*&
	    dummyvmask(iloop,jk)
    barot=barot+tmp
  ENDDO
	    !jk=1
		!PRINT*,iloop,(ovidezonalu(iloop,jk)+ovidemeridv(iloop,jk)),(e2usec(1,iloop)+e1vsec(1,iloop)),&
		!(e3usec(iloop,jk)+e3vsec(iloop,jk)),dummyvmask(iloop,jk),barot
  ENDDO
	PRINT*, 'BAROTROPIC TRANSPORT = ', barot/1.e6, ' Sv.'
!--------------------------------------------


  ! netcdf output 
    DO jk =1, npk
    ierr = putvar (ncout, id_varout(1), REAL(ovidetemper(:,jk)), jk,nsec-1,1)
    ierr = putvar (ncout, id_varout(2), REAL(ovidesaline(:,jk)), jk,nsec-1,1)
    ierr = putvar (ncout, id_varout(3), REAL(ovidezonalu(:,jk)), jk,nsec-1,1)
    ierr = putvar (ncout, id_varout(4), REAL(ovidemeridv(:,jk)), jk,nsec-1,1)
    ierr = putvar (ncout, id_varout(5), REAL(dumisec(1,:)), jk,1,nsec)
    ierr = putvar (ncout, id_varout(6), REAL(dumjsec(1,:)), jk,1,nsec)
    ierr = putvar (ncout, id_varout(7), REAL(e2usec(1,:)), jk,nsec-1,1)
    ierr = putvar (ncout, id_varout(8), REAL(e1vsec(1,:)), jk,nsec-1,1)
    ierr = putvar (ncout, id_varout(9), REAL(e3usec(:,jk)), jk,nsec-1,1)
    ierr = putvar (ncout, id_varout(10),REAL(e3vsec(:,jk)), jk,nsec-1,1)
  ! along-track normal velocity, horiz. and vert. resolution, and mask
    ierr = putvar (ncout, id_varout(11),REAL(ovidezonalu(:,jk) + ovidemeridv(:,jk)), jk,nsec-1,1)
    ierr = putvar (ncout, id_varout(12),REAL(e2usec(1,:) + e1vsec(1,:)), jk,nsec-1,1)
    ierr = putvar (ncout, id_varout(13),REAL(e3usec(:,jk) + e3vsec(:,jk)), jk,nsec-1,1)
    ierr = putvar (ncout, id_varout(14),REAL(dummyvmask(:,jk)), jk,nsec-1,1)
    END DO

  ierr = closeout(ncout)

    !!--------------------------------------------------------------------
    !!
    !!   SUBROUTINES USED
    !!
    !!--------------------------------------------------------------------
CONTAINS

  SUBROUTINE Nearestpoint(pplon,pplat,kpi,kpj,plam,pphi,kpiloc,kpjloc,ldbord)
    !!----------------------------------------------------------------------------
    !!            ***  SUBROUTINE NEARESTPOINT  ***
    !!
    !!   ** Purpose:  Computes the positions of the nearest i,j in the grid
    !!                from the given longitudes and latitudes
    !!
    !!   ** Method :  Starts on the middle of the grid, search in a 20x20 box, and move
    !!     the box in the direction where the distance between the box and the
    !!     point is minimum
    !!     Iterates ...
    !!     Stops when the point is outside the grid.
    !!     This algorithm does not work on the Mediteranean grid !
    !!
    !!   * history:
    !!        Anne de Miranda et Pierre-Antoine Darbon Jul. 2000 (CLIPPER)
    !!        Jean-Marc Molines : In NEMO form
    !!----------------------------------------------------------------------------
    IMPLICIT NONE
    !* arguments
    REAL(KIND=8),INTENT(in)   ::  pplon,pplat   !: lon and lat of target point
    INTEGER(KIND=4),INTENT (in)       ::  kpi,kpj       !: grid size
    INTEGER(KIND=4),INTENT (inout)    :: kpiloc,kpjloc  !: nearest point location
    REAL(KIND=8),DIMENSION(kpi,kpj),INTENT(in) ::  pphi,plam  !: model grid layout
    LOGICAL                   :: ldbord         !: reach boundary flag

    ! * local variables
    INTEGER(KIND=4) :: ji,jj,i0,j0,i1,j1
    INTEGER(KIND=4) :: itbl
    REAL(KIND=4) ::  zdist,zdistmin,zdistmin0
    LOGICAL, SAVE ::  lbordcell, lfirst=.TRUE.
    !!
    ! Initial values
    kpiloc = kpi/2 ; kpjloc = kpj/2    ! seek from the middle of domain
    itbl = 10                          ! block size for search
    zdistmin=1000000. ; zdistmin0=1000000.
    i0=kpiloc ;  j0=kpjloc
    lbordcell=.TRUE.;   ldbord=.FALSE.

    ! loop until found or boundary reach
    DO  WHILE ( lbordcell .AND. .NOT. ldbord)
       i0=kpiloc-itbl ;  i1=kpiloc+itbl
       j0=kpjloc-itbl ;  j1=kpjloc+itbl

       ! search only the inner domain
       IF (i0 <= 0) i0=2
       IF (i1 > kpi) i1=kpi-1
       IF (j0 <= 0) j0=2
       IF( j1 > kpj) j1=kpj-1

       ! within a block itbl+1 x itbl+1:
       DO jj=j0,j1
          DO ji=i0,i1
             ! compute true distance (orthodromy) between target point and grid point
             zdist=dist(pplon,plam(ji,jj),pplat,pphi(ji,jj) )
             zdistmin=MIN(zdistmin,zdist)
             ! update kpiloc, kpjloc if distance decreases
             IF (zdistmin .NE. zdistmin0 ) THEN
                kpiloc=ji
                kpjloc=jj
             ENDIF
             zdistmin0=zdistmin
          END DO
       END DO
       lbordcell=.FALSE.
       ! if kpiloc, kpjloc belong to block boundary proceed to next block, centered on kpiloc, kpjloc
       IF (kpiloc == i0 .OR. kpiloc == i1) lbordcell=.TRUE.
       IF (kpjloc == j0 .OR. kpjloc == j1) lbordcell=.TRUE.
       ! boundary reach ---> not found
       IF (kpiloc == 2  .OR. kpiloc ==kpi-1) ldbord=.TRUE.
       IF (kpjloc == 2  .OR. kpjloc ==kpj-1) ldbord=.TRUE.
    END DO
  END SUBROUTINE  NEARESTPOINT

SUBROUTINE interm_pt (ydpt,k,pai,pbi,paj,pbj,ydpti)
    !! -----------------------------------------------------
    !!           SUBROUTINE INTERM_PT
    !!           ********************
    !!
    !!   PURPOSE:
    !!   --------
    !!     Find the best intermediate points on a pathway.
    !!
    !!    ARGUMENTS:
    !!    ----------
    !!      ydpt : complex vector of the positions of the nearest points
    !!         k : current working index
    !!       pai ,pbi : slope and original ordinate of x(y)
    !!       paj ,pbj : slope and original ordinate of y(x)
    !!      ydpti : Complex holding the position of intermediate point
    !!
    !!    AUTHOR:
    !!    -------
    !!      19/07/1999 : Jean-Marc MOLINES
    !!      14/01/2005 : J M M in F90
    !!
    !!--------------------------------------------------------------
    !!
    !! 0. Declarations:
    !! ----------------
    IMPLICIT NONE
    COMPLEX, INTENT(in) :: ydpt(*)
    COMPLEX, INTENT(out) :: ydpti
    REAL(KIND=4), INTENT(IN) ::  pai,pbi,paj,pbj
    INTEGER(KIND=4) ,INTENT(in) :: k
    ! ... local
    COMPLEX :: ylptmp1, ylptmp2
    REAL(KIND=4) ::  za0,zb0,za1,zb1,zd1,zd2
    REAL(KIND=4) ::  zxm,zym
    REAL(KIND=4) ::  zxp,zyp
    !!
    !! 1. Compute intermediate points
    !! ------------------------------
    !
    ! ... Determines whether we use y(x) or x(y):
    IF (ABS(paj) <=  1) THEN
       ! ..... y(x)
       ! ... possible intermediate points:
       ylptmp1=ydpt(k-1)+(1.,0.)
       ylptmp2=ydpt(k-1)+CMPLX(0.,SIGN(1.,paj))
       !
       ! ... M is the candidate point:
       zxm=REAL(ylptmp1)
       zym=IMAG(ylptmp1)
       za0=paj
       zb0=pbj
       !
       za1=-1./za0
       zb1=zym - za1*zxm
       ! ... P is the projection of M on the strait line
       zxp=-(zb1-zb0)/(za1-za0)
       zyp=za0*zxp + zb0
       ! ... zd1 is the distance MP
       zd1=(zxm-zxp)*(zxm-zxp) + (zym-zyp)*(zym-zyp)
       !
       ! ... M is the candidate point:
       zxm=REAL(ylptmp2)
       zym=IMAG(ylptmp2)
       za1=-1./za0
       zb1=zym - za1*zxm
       ! ... P is the projection of M on the strait line
       zxp=-(zb1-zb0)/(za1-za0)
       zyp=za0*zxp + zb0
       ! ... zd2 is the distance MP
       zd2=(zxm-zxp)*(zxm-zxp) + (zym-zyp)*(zym-zyp)
       ! ... chose the smallest (zd1,zd2)
       IF (zd2 <=  zd1) THEN
          ydpti=ylptmp2
       ELSE
          ydpti=ylptmp1
       END IF
       !
    ELSE
       !
       ! ... x(y)
       ylptmp1=ydpt(k-1)+CMPLX(SIGN(1.,pai),0.)
       ylptmp2=ydpt(k-1)+(0.,1.)
       zxm=REAL(ylptmp1)
       zym=IMAG(ylptmp1)
       za0=pai
       zb0=pbi
       !
       za1=-1./za0
       zb1=zxm - za1*zym
       zyp=-(zb1-zb0)/(za1-za0)
       zxp=za0*zyp + zb0
       zd1=(zxm-zxp)*(zxm-zxp) + (zym-zyp)*(zym-zyp)
       !
       zxm=REAL(ylptmp2)
       zym=IMAG(ylptmp2)
       za1=-1./za0
       zb1=zxm - za1*zym
       zyp=-(zb1-zb0)/(za1-za0)
       zxp=za0*zyp + zb0
       zd2=(zxm-zxp)*(zxm-zxp) + (zym-zyp)*(zym-zyp)
       IF (zd2 <=  zd1) THEN
          ydpti=ylptmp2
       ELSE
          ydpti=ylptmp1
       END IF
    END IF
  END SUBROUTINE interm_pt

  FUNCTION dist(plona,plonb,plata,platb)
    !!----------------------------------------------------------
    !!           ***  FUNCTION  DIST  ***
    !!
    !!  ** Purpose : Compute the distance (km) between
    !!          point A (lona, lata) and B(lonb,latb)
    !!
    !!  ** Method : Compute the distance along the orthodromy
    !!
    !! * history : J.M. Molines in CHART, f90, may 2007
    !!----------------------------------------------------------
    IMPLICIT NONE
    ! Argument
    REAL(KIND=8), INTENT(in) :: plata, plona, platb, plonb
    REAL(KIND=8) :: dist
    ! Local variables
    REAL(KIND=8),SAVE ::  zlatar, zlatbr, zlonar, zlonbr
    REAL(KIND=8) ::  zpds
    REAL(KIND=8),SAVE :: zux, zuy, zuz
    REAL(KIND=8) :: zvx, zvy, zvz

    REAL(KIND=8), SAVE :: prevlat=-1000., prevlon=-1000, zr, zpi, zconv
    LOGICAL :: lfirst=.TRUE.

    ! initialise some values at first call
    IF ( lfirst ) THEN
       lfirst=.FALSE.
       ! constants
       zpi=ACOS(-1.)
       zconv=zpi/180.  ! for degree to radian conversion
       ! Earth radius
       zr=(6378.137+6356.7523)/2.0 ! km
    ENDIF

    ! compute these term only if they differ from previous call
    IF ( plata /= prevlat .OR. plona /= prevlon) THEN
       zlatar=plata*zconv
       zlonar=plona*zconv
       zux=COS(zlonar)*COS(zlatar)
       zuy=SIN(zlonar)*COS(zlatar)
       zuz=SIN(zlatar)
       prevlat=plata
       prevlon=plona
    ENDIF

    zlatbr=platb*zconv
    zlonbr=plonb*zconv
    zvx=COS(zlonbr)*COS(zlatbr)
    zvy=SIN(zlonbr)*COS(zlatbr)
    zvz=SIN(zlatbr)

    zpds=zux*zvx+zuy*zvy+zuz*zvz

    IF (zpds >= 1.) THEN
       dist=0.
    ELSE
       dist=zr*ACOS(zpds)
    ENDIF
  END FUNCTION dist

END PROGRAM cdfovide
