PROGRAM cdfovide
  !!---------------------------------------------------------------------
  !!               ***  PROGRAM cdfovide  ***
  !!
  !!  **  Purpose: Easy tool to perform Temperature, Salinity and velocity 
  !!               plots along OVIDE section 
  !!               PARTIAL STEPS version
  !!  
  !!  **  Method:  Works as a standalone file once compiled
  !!               Inspired by cdffindij, cdftransportiz 
  !! history :
  !!   Original :  R. Dussin (dec. 2008)
  !!
  !!---------------------------------------------------------------------               
  !! * Modules used
  USE cdfio
  !! * Local variables
  IMPLICIT NONE
  INTEGER :: narg, iargc !: command line
  INTEGER :: npiglo, npjglo, npk !: size of the domain
  INTEGER :: niter, nclass
  INTEGER :: imin, imax, jmin, jmax, k, ik, jk, jclass
  INTEGER :: iloc, jloc
  INTEGER :: iloop, jloop

  INTEGER :: nsec=0 ! nb total de points le long de la section
  INTEGER, DIMENSION (:), ALLOCATABLE :: isec, jsec ! indices des points a recuperer
  
  INTEGER, PARAMETER :: nsta=4
  INTEGER, DIMENSION(nsta) :: ista, jsta
  INTEGER, DIMENSION(nsta-1) :: keepn

  INTEGER ,DIMENSION (:),ALLOCATABLE ::  imeter  !: limit beetween depth level, in m (nclass -1)
  INTEGER ,DIMENSION (:),ALLOCATABLE :: ilev0,ilev1 !: limit in levels  ! nclass
  INTEGER   :: numout = 10, numvtrp=11, numhtrp=12, numstrp=14
  ! broken line stuff
  INTEGER, PARAMETER :: jpseg=10000
  INTEGER :: i0,j0,i1,j1, i, j
  INTEGER :: n,nn, jseg, kk
  INTEGER :: norm_u, norm_v, ist, jst
  INTEGER    ::  nxtarg
  INTEGER, DIMENSION(nsta-1,jpseg) :: legs1=0, legs2=0

  REAL(KIND=8),DIMENSION(nsta)              :: lonsta, latsta
  REAL(KIND=8)                              :: xmin, xmax, ymin, ymax, rdis
  REAL(KIND=4)                              :: glamfound, glamin, glamax
  REAL(KIND=8)                              :: glam0, emax
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: glam, gphi, e1, e2
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: e1t, e2t, e1u, e2v, e3t
  
  REAL(KIND=4) ::  rxi0,ryj0, rxi1, ryj1
  REAL(KIND=4) ::   ai,bi, aj,bj,d
  REAL(KIND=4) ::    rxx(jpseg),ryy(jpseg)
  REAL(KIND=4), DIMENSION(jpseg) :: gla !, gphi

  REAL(KIND=8), DIMENSION(jpseg) :: voltrp, heatrp, saltrp
  REAL(KIND=8)                   :: voltrpsum, heatrpsum, saltrpsum
  COMPLEX yypt(jpseg), yypti

  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::         e1v, e3v ,gphiv, zv, zvt, zvs !: mask, metrics
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::         e2u, e3u ,gphiu, zu, zut, zus !: mask, metrics
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::         temper, saline, zonalu, meridv, navlon, navlat
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::         ovidetemper, ovidesaline, ovidezonalu, ovidemeridv
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::         lonsec, latsec, dumisec, dumjsec
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::         e1tsec, e1usec, e1vsec, e2tsec, e2usec, e2vsec
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::         e3tsec, e3usec, e3vsec
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::         glamu, glamv
  REAL(KIND=4), DIMENSION (:),     ALLOCATABLE ::         gdepw
  REAL(KIND=4)                                 ::   rd1, rd2
  REAL(KIND=4)                                 ::  udum, vdum

  REAL(KIND=8),   DIMENSION (:,:), ALLOCATABLE :: zwku,zwkv,    zwkut,zwkvt,   zwkus,zwkvs
  REAL(KIND=8),   DIMENSION (:,:,:), ALLOCATABLE :: ztrpu, ztrpv, ztrput,ztrpvt, ztrpus,ztrpvs

! constants
  REAL(KIND=4)   ::  rau0=1000.,  rcp=4000.

  CHARACTER(LEN=80) :: coord='coordinates.nc', ctype='F'
  CHARACTER(LEN=80) :: cfilet , cfileu, cfilev, csection 
  CHARACTER(LEN=80) :: coordhgr='mesh_hgr.nc',  coordzgr='mesh_zgr.nc', cdum
  CHARACTER(LEN=80) ,DIMENSION(4)   :: cvarname   !: array of var name for output
  LOGICAL  :: lagain, lbord
  LOGICAL    :: ltest=.FALSE.
  
! cdf output stuff
  CHARACTER(LEN=80) :: cfileoutnc='ovide.nc'
  TYPE(variable)    ,DIMENSION(:), ALLOCATABLE   :: typvar
  INTEGER    :: ierr, ncout
  REAL(KIND=4), DIMENSION (1)                    ::  tim
  INTEGER    :: nfield=10
  INTEGER ,DIMENSION (:),ALLOCATABLE :: ipk, id_varout

  !!---------------------------------------------------------------------
  !!  End of declarations, code begins here
  !!---------------------------------------------------------------------

  !!---------------------------------------------------------------------
  !!  Command line 
  !!---------------------------------------------------------------------

  narg= iargc()
  IF ( narg < 3  ) THEN
     PRINT *,' Usage : cdfovide gridTfile gridUfile gridVfile '
     PRINT *,' Files cordinates.nc, mesh_hgr.nc and mesh_zgr.nc must be in te current directory '
     PRINT *,' Output on netcdf '
     STOP
  ENDIF

  CALL getarg (1, cfilet)
  CALL getarg (2, cfileu)
  CALL getarg (3, cfilev)

  !! We define what are the 3 segments of OVIDE section
  !! so that the user don't have to worry about it
  !! sec1 : (lonsta1,latsta1) -> (lonsta2,latsta2)
  !! and so on

  lonsta(1)=-43.0
  lonsta(2)=-31.3
  lonsta(3)=-12.65
  lonsta(4)=-8.7

  latsta(1)=60.6
  latsta(2)=58.9
  latsta(3)=40.33
  latsta(4)=40.33

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
  !!  Find the indexes of the 3 legs (from cdffindij) 
  !!---------------------------------------------------------------------

  npiglo= getdim (coord,'x')
  npjglo= getdim (coord,'y')
  
  ALLOCATE (glam(npiglo,npjglo), gphi(npiglo,npjglo) )
  ALLOCATE (e1(npiglo,npjglo), e2(npiglo,npjglo) )

  SELECT CASE ( ctype )
  CASE ('T' , 't' )
     glam(:,:) = getvar(coord, 'glamt',1,npiglo,npjglo)
     gphi(:,:) = getvar(coord, 'gphit',1,npiglo,npjglo)
     e1  (:,:) = getvar(coord, 'e1t'  ,1,npiglo,npjglo)
     e2  (:,:) = getvar(coord, 'e2t'  ,1,npiglo,npjglo)
  CASE ('U','u' )
     glam(:,:) = getvar(coord, 'glamu',1,npiglo,npjglo)
     gphi(:,:) = getvar(coord, 'gphiu',1,npiglo,npjglo)
     e1  (:,:) = getvar(coord, 'e1u'  ,1,npiglo,npjglo)
     e2  (:,:) = getvar(coord, 'e2u'  ,1,npiglo,npjglo)
  CASE ('V','v' )
     glam(:,:) = getvar(coord, 'glamv',1,npiglo,npjglo)
     gphi(:,:) = getvar(coord, 'gphiv',1,npiglo,npjglo)
     e1  (:,:) = getvar(coord, 'e1v'  ,1,npiglo,npjglo)
     e2  (:,:) = getvar(coord, 'e2v'  ,1,npiglo,npjglo)
  CASE ('F','f' )
     glam(:,:) = getvar(coord, 'glamf',1,npiglo,npjglo)
     gphi(:,:) = getvar(coord, 'gphif',1,npiglo,npjglo)
     e1  (:,:) = getvar(coord, 'e1f'  ,1,npiglo,npjglo)
     e2  (:,:) = getvar(coord, 'e2f'  ,1,npiglo,npjglo)
  CASE DEFAULT
     PRINT *,' ERROR : type of point not known: ', TRIM(ctype)
  END SELECT

  ! work with longitude between 0 and 360 to avoid  the date line.
    WHERE( glam < 0 ) glam=glam+360.
  ! For Orca grid, the longitude of ji=1 is about 70 E
  glam0=glam(1, npjglo/2)
  WHERE( glam < glam0 ) glam=glam+360.

  !! loop on the 3 legs
  DO k = 1,nsta-1

   xmin=lonsta(k)
   ymin=latsta(k)
   xmax=lonsta(k+1)
   ymax=latsta(k+1)

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
  keepn(k)=nn
  nsec = nsec + nn
  END DO !! loop on the 3 legs

! fancy control print
WRITE(*,*) '------------------------------------------------------------'
WRITE(*,9100) 'leg 1 start at ', lonsta(1) ,'°N ', latsta(1), '°W and ends at ', lonsta(2) ,'°N ', latsta(2), '°W'
WRITE(*,9101) 'corresponding to F-gridpoints(', ista(1),',',jsta(1),') and (', ista(2),',',jsta(2),')' 
WRITE(*,*) '------------------------------------------------------------'
WRITE(*,*) '------------------------------------------------------------'
WRITE(*,9100) 'leg 2 start at ', lonsta(2) ,'°N ', latsta(2), '°W and ends at ', lonsta(3) ,'°N ', latsta(3), '°W'
WRITE(*,9101) 'corresponding to F-gridpoints(', ista(2),',',jsta(2),') and (', ista(3),',',jsta(3),')' 
WRITE(*,*) '------------------------------------------------------------'
WRITE(*,*) '------------------------------------------------------------'
WRITE(*,9100) 'leg 3 start at ', lonsta(3) ,'°N ', latsta(3), '°W and ends at ', lonsta(4) ,'°N ', latsta(4), '°W'
WRITE(*,9101) 'corresponding to F-gridpoints(', ista(3),',',jsta(3),') and (', ista(4),',',jsta(4),')' 
WRITE(*,*) '------------------------------------------------------------'

9100 FORMAT(a,f6.2,a,f6.2,a,f6.2,a,f6.2,a)
9101 FORMAT(a,i4,a,i4,a,i4,a,i4,a)

  ALLOCATE (isec(nsec), jsec(nsec)) 

  DO k=1, nsta-1
    DO iloop=1, keepn(k) 
     jloop=iloop + SUM(keepn(1:k)) -keepn(k)
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
  
  dumisec(:,:)=0
  dumjsec(:,:)=0

  navlon(:,:) = getvar(cfilet, 'nav_lon' ,1,npiglo,npjglo)
  navlat(:,:) = getvar(cfilet, 'nav_lat' ,1,npiglo,npjglo)
  e1v(:,:)    = getvar(coordhgr, 'e1v',1,npiglo,npjglo)
  e2u(:,:)    = getvar(coordhgr, 'e2u',1,npiglo,npjglo)

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

  DO iloop=1,nsec-1

  IF ( jsec(iloop+1) == jsec(iloop) ) THEN ! horizontal segment
    IF ( isec(iloop+1) > isec(iloop) ) THEN ! eastward

       e2usec(iloop,jk) = 0.
       e1vsec(iloop,jk) = e1v(isec(iloop)+1,jsec(iloop))

    ELSE

       e2usec(iloop,jk) = 0.
       e1vsec(iloop,jk) = e1v(isec(iloop),jsec(iloop))

    ENDIF
  ELSEIF ( isec(iloop+1) == isec(iloop) ) THEN ! vertical segment
    IF ( jsec(iloop+1) < jsec(iloop) ) THEN ! southward

       e2usec(iloop,jk) = e2u(isec(iloop),jsec(iloop))
       e1vsec(iloop,jk) = 0.

    ELSE

       e2usec(iloop,jk) = e2u(isec(iloop),jsec(iloop)+1)
       e1vsec(iloop,jk) = 0.
    ENDIF
  ELSE
       PRINT *, 'problem'
       exit 
  ENDIF
  END DO

  ! loop on 3d arrays
  DO jk=1,npk
  temper(:,:) = getvar(cfilet, 'votemper',jk,npiglo,npjglo)
  saline(:,:) = getvar(cfilet, 'vosaline',jk,npiglo,npjglo)
  zonalu(:,:) = getvar(cfileu, 'vozocrtx',jk,npiglo,npjglo)
  meridv(:,:) = getvar(cfilev, 'vomecrty',jk,npiglo,npjglo)
  e3u(:,:) = getvar(coordzgr, 'e3u_ps',jk,npiglo,npjglo, ldiom=.true.)
  e3v(:,:) = getvar(coordzgr, 'e3v_ps',jk,npiglo,npjglo, ldiom=.true.)

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

  ALLOCATE ( typvar(nfield), ipk(nfield), id_varout(nfield) )

  DO iloop=1,nfield
     ipk(iloop) = npk
  END DO

 ! define new variables for output 
  typvar(1)%name= 'votemper'
  typvar(1)%units='deg C'
  typvar%missing_value=0.
  typvar(1)%valid_min= -2.
  typvar(1)%valid_max= 40.
  typvar%scale_factor= 1.
  typvar%add_offset= 0.
  typvar%savelog10= 0.
  typvar(1)%long_name='Temperature along OVIDE section'
  typvar(1)%short_name='votemper'
  typvar%online_operation='N/A'
  typvar%axis='TYZ'

  typvar(2)%name= 'vosaline'
  typvar(2)%units='PSU'
  typvar(2)%valid_min= 0.
  typvar(2)%valid_max= 50.
  typvar(2)%long_name='Salinity along OVIDE section'
  typvar(2)%short_name='vosaline'

  typvar(3)%name= 'vozocrtx'
  typvar(3)%units='m.s-1'
  typvar(3)%valid_min= -20.
  typvar(3)%valid_max= 20.
  typvar(3)%long_name='Zonal velocity along OVIDE section'
  typvar(3)%short_name='vozocrtx'

  typvar(4)%name= 'vomecrty'
  typvar(4)%units='m.s-1'
  typvar(4)%valid_min= -20.
  typvar(4)%valid_max= 20.
  typvar(4)%long_name='Meridionnal velocity along OVIDE section'
  typvar(4)%short_name='vomecrty'

  typvar(5)%name= 'isec'
  typvar(5)%valid_min= 0.
  typvar(5)%valid_max= npiglo 
  typvar(6)%name= 'jsec'
  typvar(6)%valid_min= 0.
  typvar(6)%valid_max= npjglo 
  typvar(7)%name= 'e2u'
  typvar(7)%valid_min= MINVAL(e2usec(1,:))
  typvar(7)%valid_max= MAXVAL(e2usec(1,:)) 
  typvar(8)%name= 'e1v'
  typvar(8)%valid_min= MINVAL(e1vsec(1,:))
  typvar(8)%valid_max= MAXVAL(e1vsec(1,:))
  typvar(9)%name= 'e3u'
  typvar(9)%valid_min= MINVAL(e3usec(:,:))
  typvar(9)%valid_max= MAXVAL(e3usec(:,:)) 
  typvar(10)%name= 'e3v'
  typvar(10)%valid_min= MINVAL(e3vsec(:,:))
  typvar(10)%valid_max= MAXVAL(e3vsec(:,:)) 

  ! create output fileset
   ncout =create(cfileoutnc, 'none', 1,nsec,npk,cdep='depthw')
   ierr= createvar(ncout ,typvar,nfield, ipk,id_varout )
   ierr= putheadervar(ncout, cfilet,1, nsec,npk,pnavlon=lonsec,pnavlat=latsec,pdep=gdepw)
   tim=getvar1d(cfilet,'time_counter',1)
   ierr=putvar1d(ncout,tim,1,'T')

  ! netcdf output 
    DO jk =1, npk
    ierr = putvar (ncout, id_varout(1), REAL(ovidetemper(:,jk)), jk,1,nsec-1)
    ierr = putvar (ncout, id_varout(2), REAL(ovidesaline(:,jk)), jk,1,nsec-1)
    ierr = putvar (ncout, id_varout(3), REAL(ovidezonalu(:,jk)), jk,1,nsec-1)
    ierr = putvar (ncout, id_varout(4), REAL(ovidemeridv(:,jk)), jk,1,nsec-1)
    ierr = putvar (ncout, id_varout(5), REAL(dumisec(1,:)), jk,1,nsec)
    ierr = putvar (ncout, id_varout(6), REAL(dumjsec(1,:)), jk,1,nsec)
    ierr = putvar (ncout, id_varout(7),REAL(e2usec(1,:)), jk,1,nsec-1)
    ierr = putvar (ncout, id_varout(8),REAL(e1vsec(1,:)), jk,1,nsec-1)
    ierr = putvar (ncout, id_varout(9),REAL(e3usec(:,jk)), jk,1,nsec-1)
    ierr = putvar (ncout, id_varout(10),REAL(e3vsec(:,jk)), jk,1,nsec-1)
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
    INTEGER,INTENT (in)       ::  kpi,kpj       !: grid size
    INTEGER,INTENT (inout)    :: kpiloc,kpjloc  !: nearest point location
    REAL(KIND=8),DIMENSION(kpi,kpj),INTENT(in) ::  pphi,plam  !: model grid layout
    LOGICAL                   :: ldbord         !: reach boundary flag

    ! * local variables
    INTEGER :: ji,jj,i0,j0,i1,j1
    INTEGER :: itbl
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
    INTEGER ,INTENT(in) :: k
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
