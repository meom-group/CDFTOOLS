MODULE modpoly
  !!======================================================================
  !!                     ***  MODULE  modpoly  ***
  !! Determine if a given point is within a polygon or not. This module is
  !! inherited from de finite element mesh generator program (TRIGRID) 
  !!=====================================================================
  !! History :  2.1  : 03/2006  : J.M. Molines : Port from trigrid
  !!            3.0  : 12/2010  : J.M. Molines : Doctor norm + Licence
  !!            4.0  : 03/2017  : J.M. Molines  
  !!
  !!          Use algorithms developped in the late 80's for a finite element
  !!          mesh generator (TRIGRID) by R. Walters, C. Werner et Al.
  !!          Some original comments are maintained for references.
  !!----------------------------------------------------------------------

  !!----------------------------------------------------------------------
  !!   routines      : description
  !!   ReadPoly      : Read polygon file 
  !!   PrepPoly      : Compute elements of the sides of polygon
  !!   InPoly        : Check if a point in in or out of a polygon
  !!   PointSlope    : Internal routine which computes slopes and coeff, of the side
  !!----------------------------------------------------------------------

  IMPLICIT NONE

  PRIVATE

  INTEGER(KIND=4), PUBLIC , PARAMETER :: jpvert = 50  !: Number of vertex per polygon
  INTEGER(KIND=4), PUBLIC , PARAMETER :: jpolys = 20  !: Number of polygons.

  ! - Storage for polygon definitions
  INTEGER(KIND=4)                    :: numpolys            ! number of of polygons currently defined
  INTEGER(KIND=4), DIMENSION(jpolys) :: nvertcnt            ! number of vertices of a given polygon

  REAL(KIND=4), DIMENSION(jpolys,jpvert+1) :: vertx, verty  ! 2dim. array of polygons and their X,Y  coordinates
  REAL(KIND=4)                             :: rmaxx, rmaxy  ! max x,y of polygon coordinates
  REAL(KIND=4)                             :: rminx, rminy  ! min x,y of polygon coordinates
  REAL(KIND=8), DIMENSION(jpvert)          :: slope         ! slope of the sides of polygone
  REAL(KIND=8), DIMENSION(jpvert)          :: ra, rb, rc    ! equation of side of polygon

  PUBLIC  :: ReadPoly
  PUBLIC  :: PrepPoly
  PUBLIC  :: InPoly
  PRIVATE :: PointSlope

  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class system
  !!----------------------------------------------------------------------

CONTAINS 

  SUBROUTINE ReadPoly(cdfront, kpoly, cdarea)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE ReadPoly  ***
    !!
    !! ** Purpose : read an ASCII file with names of polygon area
    !!        and vertices.  
    !!
    !! References :  late 80's trigrid (Walters et Al.)
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*),                 INTENT(in) :: cdfront ! Name of input file
    INTEGER(KIND=4),                 INTENT(out) :: kpoly   ! number of poylgons
    CHARACTER(LEN=*), DIMENSION (:), INTENT(out) :: cdarea  ! Name of the polygonal area

    INTEGER(KIND=4) :: jj                          ! dummy loop index
    INTEGER(KIND=4),DIMENSION(jpolys) :: ipac      ! flag for Pacific area (across date line )
    INTEGER(KIND=4) :: inum=8                      ! logical unit for input file
    INTEGER(KIND=4) :: ipoly                       ! polygon counter
    INTEGER(KIND=4) :: ivert                       ! number of vertices of a polygon
    !!----------------------------------------------------------------------
    OPEN (inum,FILE=cdfront)
    ipoly=0
    !
    DO WHILE (.TRUE.)
       ipoly=ipoly+1
       READ(inum,'(a)',END=995) cdarea(ipoly)                   ! 1rst line of block : name of polygon
       READ(inum,*)nvertcnt(ipoly), ipac(ipoly)                 ! 2nd : number of vertices, 
       ivert=nvertcnt(ipoly)
       READ(inum,*)(vertx(ipoly,jj),verty(ipoly,jj),jj=1,ivert)  ! 3rd : (x,y) pairs foreach vertex
       ! take care of the date line for pacific zone
       IF (ipac(ipoly) == 1 ) THEN
          DO jj=1,ivert 
             IF (vertx(ipoly,jj) < 0 ) vertx(ipoly,jj) = vertx(ipoly,jj) + 360.
          END DO
       ENDIF

       ! Automatically close the polygon
       vertx(ipoly,ivert+1)=vertx(ipoly,1)               
       verty(ipoly,ivert+1)=verty(ipoly,1)
       ! add dummy 0.001 to integer vertex coordinates... to avoid singular problem 
       DO jj=1, ivert+1
         IF ( (vertx(ipoly, jj) - INT( vertx(ipoly, jj) ) ) == 0 ) vertx(ipoly, jj) = vertx(ipoly, jj)+0.001
         IF ( (verty(ipoly, jj) - INT( verty(ipoly, jj) ) ) == 0 ) verty(ipoly, jj) = verty(ipoly, jj)+0.001
       END DO
    ENDDO

995 kpoly=ipoly-1

    CLOSE(inum)

  END SUBROUTINE ReadPoly


  SUBROUTINE PrepPoly ( kpolyid )
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE PrepPoly  ***
    !!
    !! ** Purpose :  determine  polygon information in preparation for
    !!               a call to InPoly. 
    !!
    !! ** Method  :  returns slope and equation of lines (ra, rc, rb)
    !!               as well as the min/max of polygon coordinates
    !!
    !! References :  Trigrid  
    !!----------------------------------------------------------------------
    INTEGER(KIND=4) ,INTENT(in) :: kpolyid   ! polygon Id

    INTEGER(KIND=4) ji        ! dummy loop index
    INTEGER(KIND=4) inumvert   ! number of vertices for polygon kpolyid
    !!----------------------------------------------------------------------
    !     - get slopes & line equations for each polygon boundary
    inumvert = nvertcnt(kpolyid)
    DO ji = 1, inumvert-1
       CALL PointSlope ( slope(ji), vertx(kpolyid,ji), vertx(kpolyid,ji+1),       &
            &                       verty(kpolyid,ji), verty(kpolyid,ji+1),       &
            &                       ra(ji), rb(ji), rc(ji) )
    END DO
    !       - ( ji = 1, inumvert-1 )
    CALL PointSlope ( slope(inumvert), vertx(kpolyid,inumvert), vertx(kpolyid,1), &
            &                          verty(kpolyid,inumvert), verty(kpolyid,1), &
            &                          ra(inumvert), rb(inumvert), rc(inumvert) )

    !     - calculate the max x,y's of polygon
    rmaxx = vertx(kpolyid,1)
    rmaxy = verty(kpolyid,1)
    DO ji = 1, inumvert
       IF (vertx(kpolyid,ji) >  rmaxx)  rmaxx = vertx(kpolyid,ji)
       IF (verty(kpolyid,ji) >  rmaxy)  rmaxy = verty(kpolyid,ji)
    END DO

    !     - calculate the min x,y's of polygon
    rminx = vertx(kpolyid,1)
    rminy = verty(kpolyid,1)
    DO ji = 1, inumvert
       IF (vertx(kpolyid,ji) <  rminx)  rminx = vertx(kpolyid,ji)
       IF (verty(kpolyid,ji) <  rminy)  rminy = verty(kpolyid,ji)
    END DO

  END SUBROUTINE PrepPoly


  SUBROUTINE InPoly ( kpolyid, pxpoint, pypoint, ld_in )
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE InPoly  ***
    !!
    !! ** Purpose :  To see if a point is inside or outside of the specified
    !!        polygon. 
    !!
    !! ** Method  :  Use the equation of the side of the polygon to determine
    !!               if a point is in (ld_in = true) or out (ld_in = false) 
    !!
    !! References : Trigrid 
    !!----------------------------------------------------------------------
    INTEGER(KIND=4), INTENT(in) ::  kpolyid           ! Polygon ID
    REAL(KIND=4),    INTENT(in) ::  pxpoint, pypoint  ! Position to check
    LOGICAL,        INTENT(out) :: ld_in              ! True if in the polygon

    INTEGER(KIND=4) :: ji
    INTEGER(KIND=4) :: icross, inumvert
    REAL(KIND=8)    :: zxpt, zx, zy, zevenodd
    !!----------------------------------------------------------------------
    inumvert = nvertcnt(kpolyid)
    !     - store coordinates of point to test
    zx = pxpoint
    zy = pypoint
    !     - get the number of cross with the polygon boundary
    icross = 0
    !     - see if point falls in the max and min range of polygon
    IF ( zx <=  rmaxx ) THEN
       IF ( zx >=  rminx ) THEN
          IF ( zy <=  rmaxy ) THEN
             IF ( zy >=  rminy ) THEN
                !               - step through the polygon boundaries
                DO ji = 1, inumvert
                   !                 - see if slope = 9999 and if point is on same y axis
                   IF ( slope(ji) ==  9999 ) THEN
                      IF ( zx >=  vertx(kpolyid,ji) ) THEN
                         IF ( ji ==  inumvert ) THEN
                            IF (        ( (zy <=  verty(kpolyid,inumvert) ) .AND.    &
                                 &        (zy >   verty(kpolyid,1)        ) ) .OR.   &
                                 &      ( (zy >=  verty(kpolyid,inumvert) ) .AND.    &
                                 &        (zy <   verty(kpolyid,1)        ) ) ) THEN
                               !                         - it has crossed the polygon boundary
                               icross = icross + 1
!                              if (zy == 398) print *, zx, zy, icross ,'A', ji
                            ENDIF ! ( zy test )
                         ELSEIF (       ( (zy <=  verty(kpolyid,ji)   ) .AND.        &
                              &           (zy >   verty(kpolyid,ji+1) ) ) .OR.       &
                              &         ( (zy >=  verty(kpolyid,ji)   ) .AND.        &
                              &           (zy <   verty(kpolyid,ji+1) ) ) ) THEN
                            !                       - it has crossed the polygon boundary
                            icross = icross + 1
!                              if (zy == 398) print *, zx, zy, icross,'B', ji
                         ENDIF !    ( ji = inumvert )
                      ENDIF   !    ( zx >= vertx(kpolyid,ji) )
                      !                   - see if normal slope (+ or -), and if point is not
                      !                    - higher or lower than y endpoints of the vertices
                   ELSEIF ( slope(ji) .NE. 0 ) THEN
                      zxpt = ( rc(ji) + zy ) / ra(ji) 
                      IF ( ji ==  inumvert ) THEN
                         IF (            ( (zxpt <=  vertx(kpolyid,inumvert) ) .AND.   &
                              &            (zxpt >   vertx(kpolyid,1)        ) ) .OR.  &
                              &          ( (zxpt >=  vertx(kpolyid,inumvert) ) .AND.   &
                              &            (zxpt <   vertx(kpolyid,1)        ) ) ) THEN
                            IF ( zx >=  zxpt) THEN
                               !                          - it has crossed the polygon boundary
                               icross = icross + 1
!                              if (zy == 398) print *, zx, zy, icross,'C', ji
                            ENDIF ! ( zx >= zxpt )
                         ENDIF !  ( zxpt test )
                      ELSEIF (            ( (zxpt <=  vertx(kpolyid,ji)   ) .AND.     &
                           &                (zxpt >   vertx(kpolyid,ji+1) ) ) .OR.    &
                           &              ( (zxpt >=  vertx(kpolyid,ji)   ) .AND.     &
                           &                (zxpt <   vertx(kpolyid,ji+1) ) ) ) THEN
                         IF ( zx >=   zxpt ) THEN
                            !                       - it has crossed the polygon boundary
                            icross = icross + 1
!                              if (zy == 398) print *, zx, zy, icross,'D', ji, slope(ji), zxpt
                         ENDIF ! ( zx >= zxpt )
                      ENDIF ! ( ji = inumvert )
                   ENDIF !  ( zxpt test )
                END DO !  ( ji = 1, inumvert )
                !          - decide how many times scanline crossed poly bounds
                zevenodd = AMOD ( ( icross * 1.0 ), 2.0 )
                IF ( zevenodd .NE. 0 ) THEN
                   !            - point is in polygon
                   ld_in = .TRUE.
                ELSE
                   ld_in = .FALSE.
                ENDIF
                !            - ( zevenodd ne 0 )
             ELSE
                ld_in = .FALSE.
             ENDIF
             !          - ( zy >= rminy )
          ELSE
             ld_in = .FALSE.
          ENDIF
          !           - ( zy <= rmaxy )
       ELSE
          ld_in = .FALSE.
       ENDIF
       !         - ( zx >= rminx )
    ELSE
       ld_in = .FALSE.
    ENDIF
    !       - ( zx <= rmaxx )

  END SUBROUTINE InPoly

  SUBROUTINE PointSlope ( pslup, pvertxa, pvertxb, pvertya, pvertyb, pax, pby, pcnstnt )
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE PointSlope  ***
    !!
    !! ** Purpose :   To get the slope and general equations of lines. 
    !!
    !! ** Method  :  GIVEN: vertxa, vertxb, vertya, vertyb = endpoints of line section
    !!                            to operate on.
    !!             RETURNS: slup = slope of the line section
    !!              ax, by, cnstnt = general eqation of the line section. 
    !!
    !! References :    trigrid
    !!----------------------------------------------------------------------
    REAL(KIND=8), INTENT(out) :: pslup
    REAL(KIND=4),  INTENT(in) :: pvertxa, pvertxb, pvertya, pvertyb
    REAL(KIND=8), INTENT(out) :: pax, pby, pcnstnt

    REAL(KIND=8) :: zvertxa, zvertxb, zvertya, zvertyb
    REAL(KIND=8) :: zrise, zrun
    !!----------------------------------------------------------------------

    zvertxa = pvertxa  ; zvertxb = pvertxb 
    zvertya = pvertya  ; zvertyb = pvertyb 

    zrise = zvertyb - zvertya
    zrun  = zvertxb - zvertxa

    IF ( zrun ==  0 ) THEN
       pslup = 9999
    ELSE
       pslup = zrise / zrun
    ENDIF

    IF ( ABS(pslup) <=  0.001 ) THEN
       pslup = 0.0
    ENDIF

    IF ( pslup ==  0 ) THEN
       pax = pslup
       pby = 1
       pcnstnt = zvertya
    ELSEIF ( pslup ==  9999 ) THEN
       pax = 1
       pby = 0
       pcnstnt = zvertxa
    ELSE
       pax = pslup
       pby = -1
       pcnstnt = ( pslup * zvertxa - zvertya )
    ENDIF

  END SUBROUTINE PointSlope

END MODULE modpoly
