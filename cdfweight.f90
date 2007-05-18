PROGRAM cdfweight
  !!-------------------------------------------------------------------
  !!               ***  PROGRAM cdfweight  ***
  !!
  !!  **  Purpose  :  return a binary weight file to be used by cdfcoloc
  !!  
  !!  **  Method   : Use Greg Holloway iyxz.txt file type as input, to specify
  !!                 the points to search in the model grid. 
  !!                 Read the coordinate/mesh_hgr file and look
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
  INTEGER :: jk
  INTEGER :: imin,  jmin
  INTEGER :: iloc, jloc, kloc
  INTEGER :: npiglo, npjglo, iquadran, npk
  INTEGER :: numgreg=10, numbin=20, ios
  
  ! Greg Holloway input data
  INTEGER       :: id, idep
  REAL(KIND=8)  :: xmin, ymin

  REAL(KIND=8)                              :: emax, hPp          !: local maximum metrics
  REAL(KIND=8)                              :: glam0              !: longitude of grid point ji=1
  REAL(KIND=8)                              :: glamin, gphimin    !: coordinates of the nearest point  (NP)
  REAL(KIND=8)                              :: glamN, gphiN, hN   !: grid point North of NP, true heding from NP
  REAL(KIND=8)                              :: glamE, gphiE, hE   !: grid point East of NP, true heding from NP
  REAL(KIND=8)                              :: glamS, gphiS, hS   !: grid point South of NP, true heding from NP
  REAL(KIND=8)                              :: glamW, gphiW, hW   !: grid point West of NP, true heding from NP
  REAL(KIND=8), DIMENSION(0:4)              :: glami, gphii       !: the 4 grid points around target (1-4) + the target (0)
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: glam, gphi, e1, e2 !: grid layout and metrics
  REAL(KIND=8), DIMENSION(:)  , ALLOCATABLE :: gdept              !: vertical depth
  REAL(KIND=8)                              :: hP, rdis           !: true heading and distance of target point from NP
  REAL(KIND=8)                              :: alpha, beta, gamma !: reduced coordinates (0-1) in the NP gridcell
                                                                  !: vertical weight

  CHARACTER(LEN=80) :: coord='coordinates.nc', ctype='F', cfile, czgr='mesh_zgr.nc'
  CHARACTER(LEN=80) :: cweight             !: weight file name

  LOGICAL  :: lagain, lbord, ldebug=.false.     !: additional debug print if ldebug=true

  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg < 1 ) THEN
     PRINT *,' Usage : cdfweight  Greg_File  [coord_file] [point_type]'
     PRINT *,' return the i,j  position for the x,y point  (nearest point ) '
     PRINT *,' as read in coord_file for the point type specified by point_type'
     PRINT *, TRIM(czgr),' files must be present in the current dir'
     PRINT *, 'if not given as argument, ',TRIM(coord),' must also be present in current dir'
     PRINT *, 'produce a weight file called weight_point_type.bin'
     PRINT *,' Example : cdfweight  iyxz7904.txt  coordinate_ORCA025.nc F '
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
  npk=    getdim (czgr,'z')

  ALLOCATE (glam(npiglo,npjglo), gphi(npiglo,npjglo) )
  ALLOCATE (e1(npiglo,npjglo), e2(npiglo,npjglo),gdept(npk) )

  ! set name and open  output weight file
  WRITE(cweight,'("weight_",a,".bin")') TRIM(ctype)
  OPEN(numbin, FILE=cweight,FORM='unformatted')

  ! read depth of model T points (hence U and V)
  gdept(:)=getvare3(czgr,'gdept',npk)

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

  ! For Orca grid, the longitude of ji=1 is about 70 E
  glam0=glam(1, npjglo/2)  
  WHERE( glam < glam0 ) glam=glam+360.

  OPEN(numgreg,FILE=cfile)
  ! Greg (Holloway) files are iyxz.txt file
  ios=0
  ! loop for each line of Greg File
  DO WHILE (ios == 0 )
     READ(numgreg,*,iostat=ios) id,ymin,xmin,idep
     IF( ios == 0 ) THEN  ! EOF not reached
        ! look for kloc = k index of point above idep
        kloc=1
        DO jk=1, npk-1
           IF ( idep >= gdept(jk) ) THEN
             kloc=jk 
           ELSE
             EXIT
           ENDIF
        ENDDO
        ! compute gamma such that Vint= (1-gamma) x V(kloc) + gamma x V(kloc +1)
        gamma=(idep - gdept(kloc))/(gdept(kloc+1)-gdept(kloc) )
        IF (kloc == npk -1 ) gamma=0
        IF ( ldebug) print '("DEP", f8.1,i8,f8.0,f8.4)', gdept(kloc), idep, gdept(kloc+1), gamma
        IF ( gamma < 0 ) THEN
           kloc=1
           gamma = 0.
        ENDIF
        IF ( gamma > 1 ) THEN
           kloc=npk -1
           gamma = 0.
        ENDIF

        ! Now deal with horizontal interpolation
        ! set longitude of input point in accordance with glam ( [glam0, 360+glam0 [ )
        IF ( xmin < 0.    ) xmin = xmin + 360.
        IF ( xmin < glam0 ) xmin = xmin + 360.

        lagain = .TRUE. ;   niter = 0
        DO WHILE (lagain)
           CALL Nearestpoint(xmin,ymin,npiglo,npjglo,gphi,glam,iloc,jloc,lbord)
           ! distance between the target point and the nearest point
           rdis=dist(xmin,glam(iloc,jloc),ymin,gphi(iloc,jloc) ) ! in km

           ! typical grid size (diagonal) in the vicinity of nearest point
           emax= MAX(e1(iloc,jloc),e2(iloc,jloc))/1000.*SQRT(2.) ! in km

           ! Latitude and longitude of the neighbours on the grid
           ! define longitudes between 0 and 360 deg
           glamin=MOD(glam(iloc,jloc),360.d0)  ; gphimin=gphi(iloc,jloc)  ! nearest point
           glamN=MOD(glam(iloc,jloc+1),360.d0) ; gphiN=gphi(iloc,jloc+1)  ! N (grid)
           glamE=MOD(glam(iloc+1,jloc),360.d0) ; gphiE=gphi(iloc+1,jloc)  ! E (grid)
           glamS=MOD(glam(iloc,jloc-1),360.d0) ; gphiS=gphi(iloc,jloc-1)  ! S (grid)
           glamW=MOD(glam(iloc-1,jloc),360.d0) ; gphiW=gphi(iloc-1,jloc)  ! W (grid)

           IF (rdis  > emax ) THEN
              ! The nearest point was not found, try one iteration  (jmm ???)
              IF ( niter <  2 ) THEN
                 lagain = .TRUE.
                 jloc = npjglo-2  ! change initial point
                 niter = niter +1
              ELSE
                 ! set iloc, jloc to -1000 -1000 ( flag value)
                 lagain = .FALSE.
                 iloc=-1000 ; jloc=-1000
              ENDIF
           ELSE
              ! The nearest point is found
              lagain = .FALSE.
           END IF
        END DO  ! iteration loop

        ! transfert Nearest point to imin, jmin
        imin=iloc
        jmin=jloc

        ! Restore target point longitude between 0 and 360
        xmin=MOD(xmin,360.d0)

        ! Compute heading of target point and neighbours from the nearest point
        hP=heading(glamin,xmin,gphimin,ymin)    ! target point
        hN=heading(glamin,glamN,gphimin,gphiN)  ! 'north' on the grid
        hE=heading(glamin,glamE,gphimin,gphiE)  ! 'east' on the grid
        hS=heading(glamin,glamS,gphimin,gphiS)  ! 'south' on the grid
        hW=heading(glamin,glamW,gphimin,gphiW)  ! 'west' on the grid

        ! determine the sector in wich the target point is located: ( from 1, to 4 resp. NE, SE, SW, NW  of the grid)
        iquadran=4
        ! to avoid problem with the GW meridian, pass to -180, 180 when working around GW
        IF ( hP > 180 ) THEN 
          hPp=hP-360
        ELSE 
          hPp=hP
        ENDIF

        IF ( hN > hE ) hN=hN -360.
        IF ( hPp > hN .AND. hPp <= hE ) iquadran=1
        IF ( hP > hE .AND. hP <= hS )   iquadran=2
        IF ( hP > hS .AND. hP <= hW )   iquadran=3
        IF ( hP > hW .AND. hPp <= hN)   iquadran=4

        glami(0) = xmin   ; gphii(0) = ymin          ! fill glami, gphii for 0 = target point
        glami(1) = glamin ; gphii(1) = gphimin       !                       1 = nearest point
        IF ( iloc /= -1000 ) THEN
        SELECT CASE ( iquadran )         ! point 2 3 4 are counter clockwise in the respective sector
        CASE ( 1 ) 
          glami(2) = glamE ; gphii(2) = gphiE
          glami(3) = MOD(glam(imin+1,jmin+1), 360.) ; gphii(3) = gphi(imin+1,jmin+1)
          glami(4) = glamN ; gphii(4) = gphiN
        CASE ( 2 )
          glami(2) = glamS ; gphii(2) = gphiS
          glami(3) = MOD(glam(imin+1,jmin-1), 360.) ; gphii(3) = gphi(imin+1,jmin-1)
          glami(4) = glamE ; gphii(4) = gphiE
        CASE ( 3 )
          glami(2) = glamW ; gphii(2) = gphiW
          glami(3) = MOD(glam(imin-1,jmin-1), 360.) ; gphii(3) = gphi(imin-1,jmin-1)
          glami(4) = glamS ; gphii(4) = gphiS
        CASE ( 4 )
          glami(2) = glamN ; gphii(2) = gphiN
          glami(3) = MOD(glam(imin-1,jmin+1), 360.) ; gphii(3) = gphi(imin-1,jmin+1)
          glami(4) = glamW ; gphii(4) = gphiW
        END SELECT

        ! resolve a non linear system of equation for alpha and beta ( the non dimensional coordinates of target point)
        CALL localcoord( alpha, beta, glami, gphii)
        ELSE   ! point is outside the domaine, put dummy values
          alpha=-1000. ; beta=-1000.
        ENDIF
        
        IF (ldebug) PRINT 9001, id, ymin, xmin, idep ,imin, jmin, rdis,  hP, hPp, hN, hE, hS, hW, iquadran, alpha, beta
        ! output both on std output and binary weight file (same info).
        PRINT 9002, id, ymin, xmin, idep ,imin, jmin, kloc, iquadran, alpha, beta, gamma
        WRITE(numbin) id, ymin, xmin, idep ,imin, jmin, kloc, iquadran, hN, alpha, beta, gamma
     ENDIF
  ENDDO
9001 FORMAT(i10, 2f10.4,3i6,7f10.4,I4,2f8.4)
9002 FORMAT(i10, 2f10.4,4i6,I4,3f11.4)
     CLOSE(numbin)

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

  SUBROUTINE localcoord( palpha, pbeta, plam, pphi)
    !!----------------------------------------------------------
    !!           ***  SUBROUTINE  localcoord    ***
    !!
    !!  ** Purpose : Compute the local coordinate in a grid cell
    !!
    !!  ** Method : from N. Daget Web page :
    !!       http://aton.cerfacs.fr/~daget/TECHREPORT/TR_CMGC_06_18_html/node8.html
    !!
    !! * history:
    !!      Original : J.M. Molines ( May 2007)
    !!----------------------------------------------------------
  IMPLICIT NONE
    ! * Arguments
    REAL(KIND=8), DIMENSION(0:4), INTENT(in)  :: plam, pphi
    REAL(KIND=8)                , INTENT(out) :: palpha, pbeta

    ! * Local variables
    REAL(KIND=8) :: zalpha=0.d0 , zbeta=0.d0, zresmax=0.001, zres
    REAL(KIND=8) :: zdeta, zdalp, zdbet
    REAL(KIND=8) :: zdlam, zdphi, z1, z2, z3, z4
    REAL(KIND=8), DIMENSION(2,2):: za
    REAL(KIND=8), DIMENSION(0:4):: zplam
    INTEGER :: itermax=200, niter=0  !: maximum of iteration and iteration counter

    zplam=plam       !: save input longitude in workinh array
    IF ( ldebug ) THEN
      print *,plam(0), pphi(0)
      print *,9999,9999
      print *,plam(1), pphi(1)
      print *,plam(2), pphi(2)
      print *,plam(3), pphi(3)
      print *,plam(4), pphi(4)
      print *,plam(1), pphi(1)
      print *,9999,9999
    ENDIF
   IF ( ABS( zplam(1) -zplam(4) ) >= 180. .OR. ABS( zplam(1) -zplam(2) ) >=180.) THEN
      ! then we are near the 0 deg line and we must work in the frame -180 180 
      WHERE ( zplam >= 180. ) zplam=zplam -360.
   ENDIF

    zres=1000.; zdlam=0.5; zdphi=0.5 ;  zalpha=0.d0 ; zbeta=0.d0; niter=0

    DO WHILE (zres > zresmax .AND. niter < itermax)
     z1=(zplam(2)- zplam(1) )
     z2=(zplam(1) -zplam(4) )
     z3=(zplam(3) -zplam(2) )

     za(1,1) =  z1 + (z2 + z3 )* zbeta 
     za(1,2) = -z2 + (z2 + z3 )* zalpha

     za(2,1) = pphi(2)-pphi(1) +  (pphi(1) -pphi(4) +pphi(3) -pphi(2))* zbeta 
     za(2,2) = pphi(4)-pphi(1) +  (pphi(1) -pphi(4) +pphi(3) -pphi(2))* zalpha

     ! determinant 
     zdeta=det(za(1,1), za(1,2), za(2,1), za(2,2) )

     ! solution of 
     ! | zdlam |        | zdalp |
     ! |       | =  za .|       |
     ! | zdphi |        | zdbet |
     zdalp=det(zdlam,  za(1,2) , zdphi, za(2,2)  )/zdeta
     zdbet=det(za(1,1)  , zdlam, za(2,1)   ,zdphi)/zdeta

     ! compute residual ( loop criteria)
     zres=sqrt(zdalp*zdalp + zdbet*zdbet )
     
     ! Compute alpha and beta from 1rst guess :
     zalpha = zalpha + zdalp
     zbeta  = zbeta  + zdbet

     ! compute corresponding lon/lat for this alpha, beta
     zdlam=zplam(0) - ((1.-zalpha)*(1-zbeta)*zplam(1) + zalpha*(1-zbeta)*zplam(2) + & 
        &                       zalpha*zbeta*zplam(3) + (1-zalpha)*zbeta*zplam(4))
     zdphi=pphi(0) - ((1.-zalpha)*(1-zbeta)*pphi(1) + zalpha*(1-zbeta)*pphi(2) + &
        &                      zalpha*zbeta*pphi(3) + (1-zalpha)*zbeta*pphi(4))
     
     niter=niter + 1  ! increment iteration counter
    END DO   ! loop until zres small enough (or itermax reach )

     palpha = zalpha
     pbeta  = zbeta
  END SUBROUTINE localcoord

  FUNCTION det(p1,p2,p3,p4)
    !!----------------------------------------------------------
    !!          ***  FUNCTION DET   ***
    !!
    !!    ** Purpose : compute determinant 
    !!
    !! * history:
    !!     J.M. Molines may 2007
    !!----------------------------------------------------------
    IMPLICIT NONE
    ! * Arguments
    REAL(KIND=8),INTENT(in) :: p1, p2, p3, p4
    REAL(KIND=8) :: det

    det = p1*p4 - p2*p3
  END FUNCTION det

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

  FUNCTION heading(plona, plonb, plata, platb)
    !!--------------------------------------------------------------
    !!           ***   FUNCTION HEADING   ***
    !!
    !!   ** Purpose: Compute true heading between point a and b
    !!
    !!   ** Method : suppose that the 2 points are not too far away from each other
    !!            so that heading can be computed with loxodromy
    !!
    !!  * history
    !!         J.M. Molines, may 2007
    !!--------------------------------------------------------------
    IMPLICIT NONE
    !*  Arguments
    REAL(KIND=8), INTENT(in) :: plata, plona, platb, plonb
    REAL(KIND=8) :: heading

    ! * Local variables
    REAL(KIND=8) :: zpi, zconv
    REAL(KIND=8) ::  angled, pi,cut_dist
    REAL(KIND=8) ::  xa,xb,ya,yb, xb_xa

    zpi=ACOS(-1.d0)
    zconv=zpi/180.d0  ! for degree to radian conversion

    ! there is a problem if the Greenwich meridian pass between a and b
    IF ( ldebug) print *,' Plonb  Plona ' , plonb, plona
    xa=plona*zconv
    xb=plonb*zconv

    ya=-LOG(tand(45.d0-plata/2.d0))
    yb=-LOG(tand(45.d0-platb/2.d0))

    IF (ldebug) PRINT *,' xa_xb , modulo 2pi', xb-xa, MOD((xb-xa),2*zpi)
    xb_xa=MOD((xb-xa),2*zpi)

    IF ( xb_xa >= zpi ) xb_xa = xb_xa -2*zpi
    IF ( xb_xa <= - zpi ) xb_xa = xb_xa +2*zpi
    IF (ldebug)  print *, 'yb -ya, xb_xa ',yb -ya , xb_xa

    angled=ATAN2(xb_xa,(yb-ya))
    heading=angled*180.d0/zpi
    IF (heading < 0) heading=heading+360.d0
  END FUNCTION heading

END PROGRAM cdfweight
