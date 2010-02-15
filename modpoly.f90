MODULE modpoly
  !------------------------------------------------------------------------
  !     *** MODULE modpoly ***
  !
  !  ** Purpose : The main function of this module is  the function InPoly,
  !            which returns a boolean set to true, if the point given in input
  !            is within a defined polygon, set to false in the contrary.
  !
  !  ** Method : Use algorithms developped in the late 80's for a finite element
  !           mesh generator (TRIGRID) by R. Walters, C. Werner et Al.
  !           Some original comments are maintained for references.
  ! - DEFINITIONS
  !     vertx(,) = 2dim. array of polygons and their X coordinates.
  !     verty(,) = 2dim. array of polygons and their Y coordinates.
  !     nvertcnt() = # vertices in a given polygon.
  !     numpolys = # of polygons currently defined.
  !
  !   PARAMETERS set for polygon storage in MASTER1.PAR (before)
  !     jpvert = max # vertices a polygon may have.
  !     jpolys = max # polygons that may be defined.
  !     maxpolydels = array dimension for pd(), the pending deletion array.
  !
  !   history:
  !      Original code : late 80's trigrid (Walters et Al.)
  !        adaptation to model diagnostics : J.M. Molines (03/2006)
  !-----------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------

  ! Module Variables
  IMPLICIT NONE
  PRIVATE

  INTEGER,PUBLIC , PARAMETER :: jpvert = 50, &  !: Number of vertex per polygon
       &                 jpolys = 20    !: Number of polygons.

  ! - Storage for polygon definitions
  INTEGER   :: numpolys
  INTEGER,DIMENSION(jpolys) ::  nvertcnt

  REAL(KIND=4), DIMENSION(jpolys,jpvert+1) :: vertx, verty
  REAL(KIND=4) ::  pmaxx, pmaxy, pminx, pminy
  REAL(KIND=8), DIMENSION(jpvert) ::  slope, a, b, c

  PUBLIC ReadPoly, PrepPoly, InPoly

CONTAINS 

  SUBROUTINE ReadPoly(cdfront,kpoly,cdarea)
    !!-------------------------------------------------------------------------
    !!    *** SUBROUTINE ReadPoly  ***
    !!
    !!  ** Purpose : read an ASCII file with names of polygon area
    !!        and vertices.
    !!
    !!  ** Method :
    !!
    !! history :
    !!    Original code : late 80's trigrid (Walters et Al.)
    !!       adaptation to model diagnostics : J.M. Molines (03/2006)
    !!-------------------------------------------------------------------------
    ! * Arguments
    INTEGER, INTENT(OUT) :: kpoly
    CHARACTER(LEN=*),INTENT(IN) :: cdfront   !: Name of input file
    CHARACTER(LEN=*),DIMENSION (:),INTENT(OUT) :: cdarea  !: Name of the polygonal area

    ! * Local Variables
    INTEGER,DIMENSION(jpolys) :: ipac        !: flag for Pacific area (across date line )
    INTEGER :: numpol=8                      !: logical unit for input file
    INTEGER :: ipoly                         ! polygon counter
    INTEGER :: jj, jmax

    OPEN (numpol,FILE=cdfront)
    ipoly=0
    !
    DO WHILE (.TRUE.)
       ipoly=ipoly+1
       READ(numpol,'(a)',END=995) cdarea(ipoly)                   ! 1rst line of block : name of polygon
       READ(numpol,*)nvertcnt(ipoly),ipac(ipoly)                  ! 2nd : number of vertices, 
       jmax=nvertcnt(ipoly)
       READ(numpol,*)(vertx(ipoly,jj),verty(ipoly,jj),jj=1,jmax)  ! 3rd : (x,y) pairs foreach vertex
       ! take care of the date line for pacific zone
       IF (ipac(ipoly) == 1 ) THEN
          DO jj=1,jmax 
             IF (vertx(ipoly,jj) < 0 ) vertx(ipoly,jj) = vertx(ipoly,jj) + 360.
          END DO
       ENDIF

       ! Automatically close the polygon
       vertx(ipoly,jmax+1)=vertx(ipoly,1)               
       verty(ipoly,jmax+1)=verty(ipoly,1)
    ENDDO
995 kpoly=ipoly-1
    CLOSE(numpol)

  END SUBROUTINE ReadPoly

  SUBROUTINE PrepPoly ( kpolyid )
    !!---------------------------------------------------------------------
    !!     *** SUBROUTINE PrepPoly  ***
    !!
    !!  ** Purpose : To determine certain polygon information in preparation for
    !!       a call to InPoly.
    !!  
    !!  ** Method : 
    !!   GIVEN: polyid  = the id number of polygon to use in Common POLYDEFS
    !!              in PolyStor.Inc.
    !!   RETURNS: In Common /SLOP/:
    !!          slope = array of slopes of polygon sides
    !!          a, b, c = arrays of line equation components of polygon sides
    !!          pmaxx, pmaxy, pminx, pminy = min/max x,y polygon coordinates
    !! 
    !! history:
    !!   WRITTEN: June 1990 by JDM for NODER.
    !!     for model diagnostics, J.M. Molines (03/2006)
    !-----------------------------------------------------------------------
    IMPLICIT NONE

    !! * Arguments
    INTEGER ,INTENT(IN) :: kpolyid

    !! * Local Variables
    INTEGER ji, numvert

    !     - get slopes & line equations for each polygon boundary
    numvert = nvertcnt(kpolyid)
    DO ji = 1, numvert-1
       CALL PointSlope ( slope(ji), vertx(kpolyid,ji),                   &
            &                   vertx(kpolyid,ji+1), verty(kpolyid,ji),       &
            &                   verty(kpolyid,ji+1), a(ji), b(ji), c(ji) )
    END DO
    !       - ( ji = 1, numvert-1 )
    CALL PointSlope ( slope(numvert), vertx(kpolyid,numvert),            &
         &                 vertx(kpolyid,1), verty(kpolyid,numvert),       &
         &                 verty(kpolyid,1), a(numvert),                   &
         &                 b(numvert), c(numvert) )

    !     - calculate the max x,y's of polygon
    pmaxx = vertx(kpolyid,1)
    pmaxy = verty(kpolyid,1)
    DO ji = 1, numvert
       IF (vertx(kpolyid,ji) >  pmaxx)  pmaxx = vertx(kpolyid,ji)
       IF (verty(kpolyid,ji) >  pmaxy)  pmaxy = verty(kpolyid,ji)
    END DO

    !     - calculate the min x,y's of polygon
    pminx = vertx(kpolyid,1)
    pminy = verty(kpolyid,1)
    DO ji = 1, numvert
       IF (vertx(kpolyid,ji) <  pminx)  pminx = vertx(kpolyid,ji)
       IF (verty(kpolyid,ji) <  pminy)  pminy = verty(kpolyid,ji)
    END DO

  END SUBROUTINE PrepPoly


  SUBROUTINE InPoly ( kpolyid, pxpoint, pypoint, ld_in )
    !!------------------------------------------------------------------------
    !!        ***  SUBROUTINE  InPoly  ***
    !!
    !!  ** Purpose: To see if a point is inside or outside of the specified
    !!        polygon.
    !!
    !!  ** Method:
    !!  GIVEN: Polygon data in Common SLOP, this data must be obtained by a 
    !!       call to PrepPoly prior to any calls to InPoly referencing the
    !!       same polygon. When kpolyid changes between calls to InPoly,
    !!       PrepPoly must be called to obtain the new info for the new polygon.
    !!       Passed Arguments;
    !!       kpolyid  = the id number of polygon to use in Common POLYDEFS
    !!              in PolyStor.Inc.
    !!       pxpoint, pypoint = x,y coordinates of point to test.
    !!  RETURNS: ld_in = TRUE if point is in polygon, else FALSE.
    !!
    !!  history:
    !!     Original : TRIGRID  June 1990 by JDM for NODER, based on InOut in SplitMod.For.
    !!        Model diags : J.M. Molines (03/2006)
    !!-----------------------------------------------------------------------
    IMPLICIT NONE

    !* Arguments
    INTEGER,INTENT(IN)      ::  kpolyid           ! Polygon ID
    REAL(KIND=4),INTENT(IN) ::  pxpoint, pypoint  ! Position to check
    LOGICAL,INTENT(OUT)     :: ld_in              ! True if in the polygon

    !* Local Variables
    INTEGER :: icross, ji, numvert
    REAL(KIND=8) ::  zxpt, zx, zy, zevenodd

    numvert = nvertcnt(kpolyid)
    !     - store coordinates of point to test
    zx = pxpoint
    zy = pypoint
    !     - get the number of cross with the polygon boundary
    icross = 0
    !     - see if point falls in the max and min range of polygon
    IF ( zx <=  pmaxx ) THEN
       IF ( zx >=  pminx ) THEN
          IF ( zy <=  pmaxy ) THEN
             IF ( zy >=  pminy ) THEN
                !               - step through the polygon boundaries
                DO ji = 1, numvert
                   !                 - see if slope = 9999 and if point is on same y axis
                   IF ( slope(ji) ==  9999 ) THEN
                      IF ( zx >=  vertx(kpolyid,ji) ) THEN
                         IF ( ji ==  numvert ) THEN
                            IF (        ( (zy <=  verty(kpolyid,numvert) ) .AND.    &
                                 &        (zy >   verty(kpolyid,1)       ) ) .OR.   &
                                 &      ( (zy >=  verty(kpolyid,numvert) ) .AND.    &
                                 &        (zy <   verty(kpolyid,1)       ) ) ) THEN
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
                         ENDIF !    ( ji = numvert )
                      ENDIF   !    ( zx >= vertx(kpolyid,ji) )
                      !                   - see if normal slope (+ or -), and if point is not
                      !                    - higher or lower than y endpoints of the vertices
                   ELSEIF ( slope(ji) .NE. 0 ) THEN
                      zxpt = ( c(ji) + zy ) / a(ji) 
                      IF ( ji ==  numvert ) THEN
                         IF (            ( (zxpt <=  vertx(kpolyid,numvert) ) .AND.   &
                              &            (zxpt >   vertx(kpolyid,1)       ) ) .OR.  &
                              &          ( (zxpt >=  vertx(kpolyid,numvert) ) .AND.   &
                              &            (zxpt <   vertx(kpolyid,1)       ) ) ) THEN
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
                      ENDIF ! ( ji = numvert )
                   ENDIF !  ( zxpt test )
                END DO !  ( ji = 1, numvert )
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
             !          - ( zy >= pminy )
          ELSE
             ld_in = .FALSE.
          ENDIF
          !           - ( zy <= pmaxy )
       ELSE
          ld_in = .FALSE.
       ENDIF
       !         - ( zx >= pminx )
    ELSE
       ld_in = .FALSE.
    ENDIF
    !       - ( zx <= pmaxx )

  END SUBROUTINE InPoly

  SUBROUTINE PointSlope ( pslup, pvertxa, pvertxb, pvertya, pvertyb, pax, pby, pcnstnt )
    !!-------------------------------------------------------------------------
    !!         *** SUBROUTINE PointSlope  ***
    !!
    !!  ** Purpose: To get the slope and general equations of lines.
    !!
    !!  ** Method:
    !!    GIVEN: vertxa, vertxb, vertya, vertyb = endpoints of line section
    !!                            to operate on.
    !!    RETURNS: slup = slope of the line section
    !!         ax, by, cnstnt = general eqation of the line section.
    !!
    !!  history:
    !!    Original : TRIGRID
    !!      for model diags ; J.M. Molines (03/2006)
    !!-----------------------------------------------------------------------
    IMPLICIT NONE

    !* Arguments
    REAL(KIND=4), INTENT(IN) :: pvertxa, pvertxb, pvertya, pvertyb
    REAL(KIND=8), INTENT(OUT) :: pax, pby, pcnstnt
    REAL(KIND=8), INTENT(OUT) :: pslup

    !* Local Variables
    REAL(KIND=8) :: zvertxa, zvertxb, zvertya, zvertyb
    REAL(KIND=8) ::  zrise, zrun

    zvertxa=pvertxa  ; zvertxb=pvertxb 
    zvertya=pvertya  ; zvertyb=pvertyb 

    zrise = zvertyb - zvertya
    zrun = zvertxb - zvertxa

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
