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

  REAL(KIND=4)                              :: xmin, xmax, ymin, ymax, rdis
  REAL(KIND=4)                              :: glamfound, glamin, glamax
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: glam, gphi
 
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

  SELECT CASE ( ctype )
  CASE ('T' , 't' )
     glam(:,:) = getvar(coord, 'glamt',1,npiglo,npjglo)
     gphi(:,:) = getvar(coord, 'gphit',1,npiglo,npjglo)
  CASE ('U','u' )
     glam(:,:) = getvar(coord, 'glamu',1,npiglo,npjglo)
     gphi(:,:) = getvar(coord, 'gphiu',1,npiglo,npjglo)
  CASE ('V','v' )
     glam(:,:) = getvar(coord, 'glamv',1,npiglo,npjglo)
     gphi(:,:) = getvar(coord, 'gphiv',1,npiglo,npjglo)
  CASE ('F','f' )
     glam(:,:) = getvar(coord, 'glamf',1,npiglo,npjglo)
     gphi(:,:) = getvar(coord, 'gphif',1,npiglo,npjglo)
  CASE DEFAULT
     PRINT *,' ERROR : type of point not known: ', TRIM(ctype)
  END SELECT
  ! work with longitude between 0 and 360 to avoid  the date line.
    WHERE( glam < 0 ) glam=glam+360.
    IF (xmin < 0.) xmin = xmin +360.
    IF (xmax < 0.) xmax = xmax +360.
    

  lagain = .TRUE.
  niter = 0
  DO WHILE (lagain)
     CALL Nearestpoint(xmin,ymin,npiglo,npjglo,gphi,glam,iloc,jloc,lbord)
     rdis = (xmin - glam(iloc,jloc))**2 + (ymin - gphi(iloc,jloc))**2
     rdis = SQRT(rdis)
     IF (rdis  > 1 ) THEN
         glamfound=glam(iloc,jloc) ; IF (glamfound > 180.)  glamfound=glamfound -360.
        PRINT 9000, 'Long= ',glamfound,' Lat = ',gphi(iloc,jloc)&
             &               , iloc, jloc 
        PRINT *,' Algorithme ne converge pas ', rdis 
        IF ( niter >=  1 ) STOP ' pas de convergence apres iteration'
        lagain = .TRUE.
        jloc = npjglo
        niter = niter +1
     ELSE
        PRINT *,' rdis = ', rdis
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
     rdis = (xmax - glam(iloc,jloc))**2 + (ymax - gphi(iloc,jloc))**2
     rdis = SQRT(rdis)
     IF (rdis >  1 ) THEN
         glamfound=glam(iloc,jloc) ; IF (glamfound > 180.)  glamfound=glamfound -360.
        PRINT 9000, 'Long= ',glamfound,' Lat = ',gphi(iloc,jloc) &
             &               , iloc, jloc
        PRINT *,' Algorithme ne converge pas ', rdis
        IF ( niter >= 1 ) STOP ' pas de convergence avres iteration'
        lagain = .TRUE.
        jloc  = npjglo
        niter = niter +1
     ELSE
        PRINT *,' rdis = ', rdis
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
9001 FORMAT(4i6)
9002 FORMAT(4f10.2)

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
    INTEGER jpiloc,jpjloc,jpi,jpj
    !
    REAL pplon,pplat
    !
    REAL*8 gphit(jpi,jpj),glamt(jpi,jpj)
    !
    LOGICAL lbord
    ! ... local
    INTEGER ji,jj,i0,j0,i1,j1
    INTEGER itbl
    !
    REAL zdist,zdistmin,zdistmin0
    !
    LOGICAL lbordcell
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
             zdist=(pplat- gphit(ji,jj))*(pplat- gphit(ji,jj)) + &
                  &        (pplon-glamt(ji,jj))*(pplon-glamt(ji,jj))
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
END PROGRAM cdffindij
