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
  INTEGER :: jk
  INTEGER :: imin,  jmin
  INTEGER :: iloc, jloc, kloc
  INTEGER :: npiglo, npjglo, iquadran, npk
  INTEGER :: numgreg=10, ios, id, idep, numbin=20

  REAL(KIND=8)                              :: glam0, alpha, beta, gamma
  REAL(KIND=8)                              :: xmin, ymin, hP,   rdis, emax, hPp
  REAL(KIND=8)                              :: glamin, gphimin
  REAL(KIND=8)                              :: glamN, gphiN, hN
  REAL(KIND=8)                              :: glamE, gphiE, hE
  REAL(KIND=8)                              :: glamS, gphiS, hS
  REAL(KIND=8)                              :: glamW, gphiW, hW
  REAL(KIND=8), DIMENSION(0:4)              :: glami, gphii
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: glam, gphi, e1, e2
  REAL(KIND=8), DIMENSION(:)  , ALLOCATABLE :: gdept

  CHARACTER(LEN=80) :: cdum, coord='coordinates.nc', ctype='F', cfile, czgr='mesh_zgr.nc'
  CHARACTER(LEN=80) :: cbin

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
  !  get dim (z) does not work with mesh_zgr because dim z_a exists ...
  !  npk= getdim (coord,'z')
  npk=46

  ALLOCATE (glam(npiglo,npjglo), gphi(npiglo,npjglo) )
  ALLOCATE (e1(npiglo,npjglo), e2(npiglo,npjglo),gdept(npk) )
  WRITE(cbin,'("weight_",a,".bin")') TRIM(ctype)
  OPEN(numbin, FILE=cbin,FORM='unformatted')

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
     IF( ios == 0 ) THEN
        ! look for kloc = k index of point above idep
        kloc=1
        DO jk=1, npk-1
           IF ( idep >= gdept(jk) ) THEN
             kloc=jk 
           ELSE
             EXIT
           ENDIF
        ENDDO
        gamma=(idep - gdept(kloc))/(gdept(kloc+1)-gdept(kloc) )
        IF ( gamma < 0 ) THEN
           kloc=1
           gamma = 0.
        ENDIF
        IF ( gamma > 1 ) THEN
           kloc=45
           gamma = 0.
        ENDIF

        IF ( xmin < 0.    ) xmin = xmin + 360.
        IF ( xmin < glam0 ) xmin = xmin + 360.

        lagain = .TRUE. ;   niter = 0
        DO WHILE (lagain)
           CALL Nearestpoint(xmin,ymin,npiglo,npjglo,gphi,glam,iloc,jloc,lbord)
           ! distance between the target point and the nearest point
           rdis=dist(xmin,glam(iloc,jloc),ymin,gphi(iloc,jloc) ) ! in m

           ! typical grid size (diagonal) in the vicinity of nearest point
           emax= MAX(e1(iloc,jloc),e2(iloc,jloc))/1000.*SQRT(2.)

           ! Latitude and longitude of the neighbours on the grid
           ! define longitudes between 0 and 360 deg
           glamin=MOD(glam(iloc,jloc),360.d0)  ; gphimin=gphi(iloc,jloc)  ! nearest point
           glamN=MOD(glam(iloc,jloc+1),360.d0) ; gphiN=gphi(iloc,jloc+1)
           glamE=MOD(glam(iloc+1,jloc),360.d0) ; gphiE=gphi(iloc+1,jloc)
           glamS=MOD(glam(iloc,jloc-1),360.d0) ; gphiS=gphi(iloc,jloc-1)
           glamW=MOD(glam(iloc-1,jloc),360.d0) ; gphiW=gphi(iloc-1,jloc)

           IF (rdis  > emax ) THEN
              ! The nearest point was not found, try one iteration 
              IF ( niter <  2 ) THEN
                 lagain = .TRUE.
                 jloc = npjglo-2
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

        ! determine the quadran in wich the target point is located: ( from 1, to 4 resp. NE, SE, SW, NW  of the grid)
        iquadran=4
        IF ( hP > 180 ) THEN 
          hPp=hP-360
        ELSE 
          hPp=hP
        ENDIF

        IF ( hN > hE ) hN=hN -360.
!       print *, ( hPp > hN .AND. hP < hE ), (hPp > hN ), ( hPp < hE )
        IF ( hPp > hN .AND. hPp <= hE ) iquadran=1
        IF ( hP > hE .AND. hP <= hS ) iquadran=2
        IF ( hP > hS .AND. hP <= hW ) iquadran=3
        IF ( hP > hW .AND. hPp <= hN) iquadran=4

        glami(0) = xmin   ; gphii(0) = ymin
        glami(1) = glamin ; gphii(1) = gphimin
        IF ( iloc /= -1000 ) THEN
        SELECT CASE ( iquadran )
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
        CALL localcoord( alpha, beta, glami, gphii)
        ELSE
          alpha=-1000. ; beta=-1000.
        ENDIF
        
!       PRINT 9001, id, ymin, xmin, idep ,imin, jmin, rdis,  hP, hPp, hN, hE, hS, hW, iquadran, alpha, beta
        PRINT 9002, id, ymin, xmin, idep ,imin, jmin, kloc, iquadran, alpha, beta, gamma
        WRITE(numbin) id, ymin, xmin, idep ,imin, jmin, kloc, iquadran, alpha, beta, gamma
     ENDIF
  ENDDO
9001 FORMAT(i10, 2f10.4,3i6,7f10.4,I4,2f8.4)
9002 FORMAT(i10, 2f10.4,4i6,I4,3f11.4)
     CLOSE(numbin)

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
    ! Remnant of the Clipper old good time ...
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
!      IF (i0 <= 0) i0=1
!      IF (i1 > jpi) i1=jpi
!      IF (j0 <= 0) j0=1
!      IF( j1 > jpj) j1=jpj
       IF (i0 <= 0) i0=2
       IF (i1 > jpi) i1=jpi-1
       IF (j0 <= 0) j0=2
       IF( j1 > jpj) j1=jpj-1
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
       IF (jpiloc == 2  .OR. jpiloc ==jpi-1) lbord=.TRUE.
       IF (jpjloc == 2  .OR. jpjloc ==jpj-1) lbord=.TRUE.
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
    ! * Argument
    REAL(KIND=8), DIMENSION(0:4), INTENT(in)  :: plam, pphi
    REAL(KIND=8)                , INTENT(out) :: palpha, pbeta

    ! * Local variables
    REAL(KIND=8) :: zalpha=0.d0 , zbeta=0.d0, zresmax=0.001, zres
    REAL(KIND=8) :: zdeta, zdalp, zdbet
    REAL(KIND=8) :: zdlam, zdphi, z1, z2, z3, z4
    REAL(KIND=8), DIMENSION(2,2):: za
    INTEGER :: itermax=200, niter=0
    LOGICAL :: ldebug=.false.

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

    zres=1000.; zdlam=0.5; zdphi=0.5 ;  zalpha=0.d0 ; zbeta=0.d0; niter=0
    DO WHILE (zres > zresmax .AND. niter < itermax)
     za(1,1) = plam(2)-plam(1) +  (plam(1) -plam(4) +plam(3) -plam(2))* zbeta 
     za(1,2) = plam(4)-plam(1) +  (plam(1) -plam(4) +plam(3) -plam(2))* zalpha
     za(2,1) = pphi(2)-pphi(1) +  (pphi(1) -pphi(4) +pphi(3) -pphi(2))* zbeta 
     za(2,2) = pphi(4)-pphi(1) +  (pphi(1) -pphi(4) +pphi(3) -pphi(2))* zalpha

     zdeta=det(za(1,1), za(1,2), za(2,1), za(2,2) )
     z1=(plam(4) -plam(1) ) + (plam(1) -plam(4) +plam(3) -plam(2))*zalpha
     z2=(pphi(4) -pphi(1) ) + (pphi(1) -pphi(4) +pphi(3) -pphi(2))*zalpha
     z3=(plam(2) -plam(1) ) + (plam(1) -plam(4) +plam(3) -plam(2))*zbeta
     z4=(pphi(2) -pphi(1) ) + (pphi(1) -pphi(4) +pphi(3) -pphi(2))*zbeta

     zdalp=det(zdlam,  za(1,2) , zdphi, za(2,2)  )/zdeta
     zdbet=det(za(1,1)  , zdlam, za(2,1)   ,zdphi)/zdeta

     zres=sqrt(zdalp*zdalp + zdbet*zdbet )
     zalpha=zalpha + zdalp
     zbeta=zbeta + zdbet
     zdlam=plam(0) - ((1.-zalpha)*(1-zbeta)*plam(1) + zalpha*(1-zbeta)*plam(2) + zalpha*zbeta*plam(3) + (1-zalpha)*zbeta*plam(4))
     zdphi=pphi(0) - ((1.-zalpha)*(1-zbeta)*pphi(1) + zalpha*(1-zbeta)*pphi(2) + zalpha*zbeta*pphi(3) + (1-zalpha)*zbeta*pphi(4))
     
!    print *, plam, pphi
!    print 9040, niter, zres, zdalp, zdbet, zdlam, zdphi
     niter=niter + 1
    END DO
     palpha=zalpha
     pbeta=zbeta
 9040 FORMAT (I5, 10e10.1)
  END SUBROUTINE localcoord

  FUNCTION det(p1,p2,p3,p4)
    IMPLICIT NONE
    REAL(KIND=8) :: p1, p2, p3, p4
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

    zpi=ACOS(-1.d0)
    zconv=zpi/180.d0  ! for degree to radian conversion

    xa=plona*zconv
    IF (ABS(plonb -plona) > 5 ) xa=(plona-360.)*zconv
!   IF (plona > 350 ) xa=(plona-360.)*zconv
    xb=plonb*zconv

    ya=-LOG(tand(45.d0-plata/2.d0))
    yb=-LOG(tand(45.d0-platb/2.d0))
!   print *, 'yb -ya, xb_xa ',yb -ya , xb -xa

    angled=ATAN2((xb-xa),(yb-ya))
    heading=angled*180.d0/zpi
    IF (heading < 0) heading=heading+360.d0
!   print *, plona, plata, plonb, platb, heading
  END FUNCTION heading

END PROGRAM cdfweight
