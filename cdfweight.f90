PROGRAM cdfweight
  !!-------------------------------------------------------------------
  !!               ***  PROGRAM cdfweight  ***
  !!
  !!  **  Purpose  :  return the i,j of the nearest point
  !!  
  !!  **  Method   :  Read the coordinate/mesh_hgr file and look
  !!                  for the glam, gphi variables
  !!                  Then use a seach algorithm to find the corresponding I J
  !!                 The point type ( T U V F ) is specified on the command line
  !!                 as well as the name of the coordinate/mesh hgr file.
  !!
  !! history ;
  !!  Original :  J.M. Molines (November 2005 )
  !!              J.M. Molines (May 2007 for weight)
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
  INTEGER :: numgreg=10, ios, id, idep

  REAL(KIND=8)                       :: xmin, xmax, ymin, ymax, rdis, emax
  REAL(KIND=4)                              :: glam0
  REAL(KIND=8)                              :: glamfound, glamin, gphimin
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: glam, gphi, e1, e2

  CHARACTER(LEN=80) :: cdum, coord='coordinates.nc', ctype='F', cfile

  LOGICAL  :: lagain, lbord
  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg < 1 ) THEN
     PRINT *,' Usage : cdfweight  Greg_File  [coord_file] [point_type]'
     PRINT *,' return the i,j  position for the x,y point  (nearest point ) '
     PRINT *,' as read in coord_file for the point type specified by point_type'
     PRINT *,' Example : cdfweight  -70 15  coordinate_ORCA025.nc F '
     STOP
  ENDIF

  CALL getarg (1, cfile )
  ! if 3rd argument not given coordinates.nc is assumed
  IF ( narg > 1 ) THEN
     CALL getarg (2, coord )
  ENDIF
  ! if 4th argument not given, assume F point
  IF ( narg == 3 ) THEN
     CALL getarg (3, ctype )
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
  WHERE( glam < 0 ) glam(:,:)=glam(:,:)+360.
  glam0=glam(1, npjglo/2)  
  WHERE( glam < glam0 ) glam=glam+360.
  !   glam(:,:)=glam(:,:)-glam0
  ! print *, glam(1,:)
  OPEN(numgreg,FILE=cfile)
  ios=0
  DO WHILE (ios == 0 )
     READ(numgreg,*,iostat=ios) id,ymin,xmin,idep
     IF( ios == 0 ) THEN
        IF (xmin < 0.) xmin = xmin +360.
        IF ( xmin < glam0 ) xmin=xmin+360.

        lagain = .TRUE.
        niter = 0
        DO WHILE (lagain)
           CALL Nearestpoint(xmin,ymin,npiglo,npjglo,gphi,glam,iloc,jloc,lbord)
           !          rdis = (xmin - glam(iloc,jloc))**2 + (ymin - gphi(iloc,jloc))**2
           !          rdis = SQRT(rdis)
           rdis=dist(xmin,glam(iloc,jloc),ymin,gphi(iloc,jloc) ) ! in m
           emax= MAX(e1(iloc,jloc),e2(iloc,jloc))/1000.*SQRT(2.)
           glamin=glam(iloc,jloc) ; gphimin=gphi(iloc,jloc)
           IF (glamin >= 360. ) glamin=glamin - 360.
           IF (rdis  > emax ) THEN
              glamfound=glam(iloc,jloc) ; IF (glamfound > 180.)  glamfound=glamfound -360.
              !       PRINT 9000, 'Long= ',glamfound,' Lat = ',gphi(iloc,jloc)&
              !            &               , iloc, jloc 
              !       PRINT *,' Algorithme ne converge pas ', rdis 
              !       IF ( niter >=  1 ) STOP ' pas de convergence apres iteration'
              IF ( niter <  2 ) THEN
                 lagain = .TRUE.
                 jloc = npjglo-2
                 niter = niter +1
              ELSE
                 lagain = .FALSE.
                 iloc=-1000 ; jloc=-1000
              ENDIF
           ELSE
              !       PRINT *,' rdis = ', rdis
              lagain = .FALSE.
           END IF
        END DO
        IF (lbord) THEN
           !    WRITE (*,*)'Point  Out of domain or on boundary', iloc, jloc
           imin=iloc ; jmin=jloc
        ELSE
           imin=iloc
           jmin=jloc
           !      PRINT 9000, 'Long= ',glam(iloc,jloc),' lat = ',gphi(iloc,jloc), iloc, jloc 
        ENDIF
        !
        IF ( xmin > 360.) xmin=xmin - 360.
        PRINT 9001, id, ymin, xmin, idep ,imin, jmin, rdis, emax, heading(xmin,glamin,ymin,gphimin)
        ! glamin=glam(imin,jmin) 
        ! IF ( glamin > 180 ) glamin=glamin-360.
     ENDIF
  ENDDO
  ! PRINT 9002, glamin, gphi(imin,jmin)
9000 FORMAT(a,f8.2,a,f8.2,2i5)
9001 FORMAT(i10, 2f10.4,3i6,4f10.4)
9002 FORMAT(2f10.2)

CONTAINS
  !!--------------------------------------------------------------------------
  SUBROUTINE Nearestpoint(pplon,pplat,jpi,jpj,gphit,glamt,jpiloc,jpjloc,lbord)
!!!----------------------------------------------------------------------------
!!!           SUBROUTINE NEARESTPOINT
!!!           ***********************
!!!    PURPOSE:
!!!    --------
!!!      Computes the positions of the nearest i,j in the grid 
!!!    from the given longitudes and latitudes
    !!
    !!     METHOD:
    !!     -------
    !!        Starts on the middle of the grid, search in a 20x20 box, and move
    !!     the box in the direction where the distance between the box and the 
    !!     point is minimum
    !!     Iterates ...
    !!     Stops when the point is outside the grid.
    !!     This algorithm does not work on the Mediteranean grid !
    !!
    !!     AUTHORS:
    !!     --------
    !!      Anne de Miranda et Pierre-Antoine Darbon Jul. 2000
    !!
!!!----------------------------------------------------------------------------
    !! 0. Declarations:
    !! ----------------
    !!
    IMPLICIT NONE
    ! ... arguments
    INTEGER :: jpiloc,jpjloc,jpi,jpj
    REAL(KIND=8) ::  pplon,pplat
    REAL(KIND=8) ::  gphit(jpi,jpj),glamt(jpi,jpj)
    LOGICAL lbord
    ! ... local
    INTEGER :: ji,jj,i0,j0,i1,j1
    INTEGER :: itbl
    !
    REAL(KIND=4) ::  zdist,zdistmin,zdistmin0
    !
    LOGICAL, SAVE ::  lbordcell, lfirst=.TRUE.
    !!
    !! 1. Initialisations:
    !! -------------------
    !!
    ! Remneant of the Clipper old good time ...
    ! at all latitudes of the Mediteranean grid set jpiloc to 1 !
    !   IF ((pplat > 32.).AND.(pplat < 47).AND.( pplon > 0.) .AND. (pplon < 30) ) jpiloc=1

    jpiloc = jpi/2
    jpjloc = jpj/2
    itbl = 10
    zdistmin=1000000.
    zdistmin0=1000000.
    i0=jpiloc
    j0=jpjloc
    lbordcell=.TRUE.
    lbord=.FALSE.
    !!
    !! 2. Main loop
    !! ------------
    !!
    DO  WHILE ( lbordcell .AND. .NOT. lbord)
       i0=jpiloc-itbl
       i1=jpiloc+itbl
       j0=jpjloc-itbl
       j1=jpjloc+itbl
       IF (i0 <= 0) i0=1
       IF (i1 > jpi) i1=jpi
       IF (j0 <= 0) j0=1
       IF( j1 > jpj) j1=jpj
       DO jj=j0,j1
          DO ji=i0,i1
             !            IF (pplat < 70 ) THEN
             !            zdist=(pplat- gphit(ji,jj))*(pplat- gphit(ji,jj)) + &
             !                 &        (pplon-glamt(ji,jj))*(pplon-glamt(ji,jj))
             !            ELSE
             zdist=dist(pplon,glamt(ji,jj),pplat,gphit(ji,jj) )
             !            ENDIF
             zdistmin=MIN(zdistmin,zdist)
             IF (zdistmin .NE. zdistmin0 ) THEN
                jpiloc=ji
                jpjloc=jj
             ENDIF
             zdistmin0=zdistmin
          END DO
       END DO
       lbordcell=.FALSE.
       IF (jpiloc == i0 .OR. jpiloc == i1) lbordcell=.TRUE.
       IF (jpjloc == j0 .OR. jpjloc == j1) lbordcell=.TRUE.
       IF (jpiloc == 1  .OR. jpiloc ==jpi) lbord=.TRUE.
       IF (jpjloc == 1  .OR. jpjloc ==jpj) lbord=.TRUE.
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

  FUNCTION heading(plona, plonb, plata, platb)
    IMPLICIT NONE
    ! Argument 
    REAL(KIND=8), INTENT(in) :: plata, plona, platb, plonb
    REAL(KIND=8) :: heading
    REAL(KIND=8) :: zpi, zconv
    REAL(KIND=8) ::  angled, pi,cut_dist
    REAL(KIND=8) ::  xa,xb,ya,yb

    zpi=ACOS(-1.)
    zconv=zpi/180.  ! for degree to radian conversion

    xa=plona*zconv
    xb=plonb*zconv

    ya=-log(tand(45.-plata/2.))
    yb=-log(tand(45.-platb/2.))

    angled=ATAN2((xb-xa),(yb-ya))
    heading=angled*180./zpi
    IF (heading < 0) heading=heading+360.
  END FUNCTION heading

END PROGRAM cdfweight
