PROGRAM cdffindij
  !!-------------------------------------------------------------------
  !!               ***  PROGRAM cdffindij  ***
  !!
  !!  **  Purpose  :  return the window index (imin imax jmin jmax )
  !!          for the geographical windows given on input (longmin longmax latmin matmax)
  !!  
  !!  **  Method   :  Read the coordinate/mesh_hgr file and look
  !!                  for the glam, gphi variables
  !!                  Then use a seach algorithm to find the corresponding I J
  !!                 The point type ( T U V F ) is specified on the command line
  !!                 as well as the name of the coordinate/mesh hgr file.
  !!
  !! history ;
  !!  Original :  J.M. Molines (November 2005 )
  !!-------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !! * Modules used
  USE cdfio

  !! * Local variables
  IMPLICIT NONE
  INTEGER :: narg, iargc, niter
  INTEGER :: imin, imax, jmin, jmax
  INTEGER :: iloc, jloc
  INTEGER :: npiglo, npjglo

  REAL(KIND=8)                              :: xmin, xmax, ymin, ymax, rdis
  REAL(KIND=4)                              :: glamfound, glamin, glamax
  REAL(KIND=8)                              :: glam0, emax
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: glam, gphi, e1, e2
 
  CHARACTER(LEN=80) :: cdum, coord='coordinates.nc', ctype='F'

  LOGICAL  :: lagain, lbord
  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg < 4 ) THEN
     PRINT *,' Usage : cdffindij  xmin xmax ymin ymax  [coord_file] [point_type]'
     PRINT *,' return the i,j  position for the zoomed area (nearest point ) '
     PRINT *,' as read in coord_file for the point type specified by point_type'
     PRINT *,' Example : cdffindij  -70 15 -20 25  coordinate_ORCA025.nc F '
     STOP
  ENDIF

  CALL getarg (1, cdum )
  READ(cdum,*) xmin
  CALL getarg (2, cdum )
  READ(cdum,*) xmax
  CALL getarg (3, cdum )
  READ(cdum,*) ymin
  CALL getarg (4, cdum )
  READ(cdum,*) ymax
  ! if 5th argument not given coordinates.nc is assumed
  IF ( narg > 4 ) THEN
     CALL getarg (5, coord )
  ENDIF
  ! if 6th argument not given, assume F point
  IF ( narg == 6 ) THEN
     CALL getarg (6, ctype )
  ENDIF

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

    IF (xmin < 0.) xmin = xmin +360.
    IF (xmax < 0.) xmax = xmax +360.
    
    IF (xmin < glam0) xmin = xmin +360.
    IF (xmax < glam0) xmax = xmax +360.

  lagain = .TRUE.
  niter = 0
  DO WHILE (lagain)
     CALL Nearestpoint(xmin,ymin,npiglo,npjglo,gphi,glam,iloc,jloc,lbord)
     ! distance between the target point and the nearest point
     rdis=dist(xmin,glam(iloc,jloc),ymin,gphi(iloc,jloc) ) ! in km
     ! typical grid size (diagonal) in the vicinity of nearest point
     emax= MAX(e1(iloc,jloc),e2(iloc,jloc))/1000.*SQRT(2.) ! in km

!    rdis = (xmin - glam(iloc,jloc))**2 + (ymin - gphi(iloc,jloc))**2
!    rdis = SQRT(rdis)
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
        PRINT '("#  rdis= ",f8.3," km")', rdis
        lagain = .FALSE.
     END IF
  END DO
  IF (lbord) THEN
     WRITE (*,*)'Point  Out of domain or on boundary'
  ELSE
     imin=iloc
     jmin=jloc
     !      PRINT 9000, 'Long= ',glam(iloc,jloc),' lat = ',gphi(iloc,jloc), iloc, jloc 
  ENDIF
  !
  lagain = .TRUE.
  niter = 0
  DO WHILE (lagain)
     CALL Nearestpoint(xmax,ymax,npiglo,npjglo,gphi,glam,iloc,jloc,lbord)
     ! distance between the target point and the nearest point
     rdis=dist(xmax,glam(iloc,jloc),ymax,gphi(iloc,jloc) ) ! in km
     ! typical grid size (diagonal) in the vicinity of nearest point
     emax= MAX(e1(iloc,jloc),e2(iloc,jloc))/1000.*SQRT(2.) ! in km
!    rdis = (xmax - glam(iloc,jloc))**2 + (ymax - gphi(iloc,jloc))**2
!    rdis = SQRT(rdis)
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
        PRINT '("#  rdis= ",f8.3," km")', rdis
        lagain = .FALSE.
     END IF
  END DO
  IF (lbord) THEN
     WRITE (*,*) 'Point  Out of domain or on boundary'
  ELSE
     imax=iloc
     jmax=jloc
     !      PRINT 9000, 'Long= ',glam(iloc,jloc),' lat = ',gphi(iloc,jloc), iloc, jloc
  ENDIF
  PRINT 9001, imin,imax, jmin, jmax
  glamin=glam(imin,jmin) ;glamax=glam(imax,jmax)
  IF ( glamin > 180 ) glamin=glamin-360.
  IF ( glamax > 180 ) glamax=glamax-360.
  PRINT 9002, glamin, glamax, gphi(imin,jmin),gphi(imax,jmax)
9000 FORMAT(a,f8.2,a,f8.2,2i5)
9001 FORMAT(4i10)
9002 FORMAT(4f10.4)

CONTAINS
  SUBROUTINE Nearestpoint(pplon,pplat,kpi,kpj,pphi,plam,kpiloc,kpjloc,ldbord)
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
END PROGRAM cdffindij
