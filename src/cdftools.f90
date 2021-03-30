MODULE cdftools
  !!======================================================================
  !!                     ***  MODULE  cdftools  ***
  !! This module holds subroutine that corresponds to cdftools.
  !! For example cdf_findij is the subroutine equivalent to cdffindij 
  !!=====================================================================
  !! History : 2.1  !  05/2010  : J.M. Molines, A. Melet : Original
  !!           3.0  !  12/2010  : J.M. Molines : Doctor + Lic.
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!   cdf_findij    : the routine version of cdffindij
  !!   NearestPoint : determine the nearest point from a lon lat location
  !!   dist          : compute the distance along othodromic route
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames

  IMPLICIT NONE

  PRIVATE     
  ! list of public subroutines that can be called
  PUBLIC  :: cdf_findij 
  PUBLIC  :: broken_line
  PUBLIC  :: dist

  PRIVATE :: NearestPoint
  PRIVATE :: interm_pt

  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class system
  !!----------------------------------------------------------------------

CONTAINS

  SUBROUTINE cdf_findij ( pxmin, pxmax, pymin, pymax,                        &
       &         kimin, kimax, kjmin, kjmax, cd_coord, cd_point, cd_verbose ,&
       &         plonmin, plonmax, platmin, platmax )
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE cdf_findij  ***
    !!
    !! ** Purpose :  the routine equivalent of cdffindij 
    !!
    !!----------------------------------------------------------------------
    REAL(KIND=4),                  INTENT(in) :: pxmin, pxmax, pymin, pymax !: geographical window in lon-lat
    INTEGER(KIND=4),              INTENT(out) :: kimin, kimax, kjmin, kjmax !: equivalent in model coordinates
    CHARACTER(*), OPTIONAL,        INTENT(in) :: cd_coord                   !: coordinate file name (D: cn_fcoo)
    CHARACTER(*), OPTIONAL,        INTENT(in) :: cd_point                   !: point type           (D: F )
    CHARACTER(*), OPTIONAL,        INTENT(in) :: cd_verbose                 !: verbose flag         (D: N ) Y
    REAL(KIND=4), OPTIONAL,       INTENT(out) :: plonmin, plonmax           !: model lonmn, lonmax
    REAL(KIND=4), OPTIONAL,       INTENT(out) :: platmin, platmax           !: model  latmin, latmax

    INTEGER(KIND=4)                           :: initer
    INTEGER(KIND=4)                           :: imin, imax, jmin, jmax
    INTEGER(KIND=4), SAVE                     :: iloc, jloc
    INTEGER(KIND=4), SAVE                     :: ipiglo, ipjglo
    INTEGER(KIND=4), PARAMETER                :: jp_itermax=15

    REAL(KIND=8)                              :: dl_xmin, dl_xmax, dl_ymin, dl_ymax
    REAL(KIND=8)                              :: dl_dis
    REAL(KIND=8)                              :: dl_glam0, dl_emax
    REAL(KIND=8), SAVE, DIMENSION(:,:), ALLOCATABLE :: dl_glam, dl_gphi, dl_e1, dl_e2

    REAL(KIND=4)                              :: zglamfound, zglamin, zglamax, zgphmin, zgphmax

    CHARACTER(LEN=256)                        :: cl_type='F'
    CHARACTER(LEN=256)                        :: clcoo

    LOGICAL                                   :: ll_again, ll_bnd, ll_verbose=.FALSE.
    !!--------------------------------------------------------------------------
    CALL ReadCdfNames()
    imin = 0 ; jmin = 0

    dl_xmin = pxmin
    dl_xmax = pxmax
    dl_ymin = pymin
    dl_ymax = pymax

    clcoo = cn_fhgr

    IF ( PRESENT( cd_coord)  ) clcoo  = cd_coord
    IF ( PRESENT( cd_point)  ) cl_type= cd_point
    IF ( PRESENT( cd_verbose))   THEN
       IF ( cd_verbose(1:1) == 'Y' .OR. cd_verbose(1:1) == 'y' ) ll_verbose=.TRUE.
    ENDIF

    IF ( .NOT. ALLOCATED (dl_glam) ) THEN 

       IF (chkfile (clcoo) ) STOP 96 ! missing file

       ipiglo= getdim (clcoo, cn_x)
       ipjglo= getdim (clcoo, cn_y)

       ALLOCATE (dl_glam(ipiglo,ipjglo), dl_gphi(ipiglo,ipjglo) )
       ALLOCATE (dl_e1  (ipiglo,ipjglo), dl_e2  (ipiglo,ipjglo) )

       SELECT CASE ( cl_type )
       CASE ('T' , 't' )
          dl_glam(:,:) = getvar(clcoo, cn_glamt, 1, ipiglo, ipjglo)
          dl_gphi(:,:) = getvar(clcoo, cn_gphit, 1, ipiglo, ipjglo)
          dl_e1  (:,:) = getvar(clcoo, cn_ve1t,  1, ipiglo, ipjglo)
          dl_e2  (:,:) = getvar(clcoo, cn_ve2t,  1, ipiglo, ipjglo)
       CASE ('U','u' )
          dl_glam(:,:) = getvar(clcoo, cn_glamu, 1, ipiglo, ipjglo)
          dl_gphi(:,:) = getvar(clcoo, cn_gphiu, 1, ipiglo, ipjglo)
          dl_e1  (:,:) = getvar(clcoo, cn_ve1u,  1, ipiglo, ipjglo)
          dl_e2  (:,:) = getvar(clcoo, cn_ve2u,  1, ipiglo, ipjglo)
       CASE ('V','v' )
          dl_glam(:,:) = getvar(clcoo, cn_glamv, 1, ipiglo, ipjglo)
          dl_gphi(:,:) = getvar(clcoo, cn_gphiv, 1, ipiglo, ipjglo)
          dl_e1  (:,:) = getvar(clcoo, cn_ve1v,  1, ipiglo, ipjglo)
          dl_e2  (:,:) = getvar(clcoo, cn_ve2v,  1, ipiglo, ipjglo)
       CASE ('F','f' )
          dl_glam(:,:) = getvar(clcoo, cn_glamf, 1, ipiglo, ipjglo)
          dl_gphi(:,:) = getvar(clcoo, cn_gphif, 1, ipiglo, ipjglo)
          dl_e1  (:,:) = getvar(clcoo, cn_ve1f,  1, ipiglo, ipjglo)
          dl_e2  (:,:) = getvar(clcoo, cn_ve2f,  1, ipiglo, ipjglo)
       CASE DEFAULT
          PRINT *,' ERROR : type of point not known: ', TRIM(cl_type)
       END SELECT

       ! work with longitude between 0 and 360 to avoid  the date line.
       WHERE( dl_glam < 0 ) dl_glam = dl_glam + 360.d0

       ! For Orca grid, the longitude of ji=1 is about 70 E
       dl_glam0 = dl_glam(1, ipjglo/2)
       WHERE( dl_glam < dl_glam0 ) dl_glam =dl_glam + 360.d0

       IF (dl_xmin < 0.) dl_xmin = dl_xmin + 360.d0
       IF (dl_xmax < 0.) dl_xmax = dl_xmax + 360.d0

       IF (dl_xmin < dl_glam0) dl_xmin = dl_xmin + 360.d0
       IF (dl_xmax < dl_glam0) dl_xmax = dl_xmax + 360.d0
    ENDIF  ! 


    ! deal with xmin, ymin
    ll_again = .TRUE.
    initer = 1
    iloc=ipiglo/2 ; jloc=ipjglo/2

    DO WHILE (ll_again)
       CALL NearestPoint(dl_xmin, dl_ymin, ipiglo, ipjglo, dl_glam, dl_gphi, iloc, jloc, ll_bnd)

       ! distance between the target point and the nearest point
       dl_dis = dist(dl_xmin, dl_glam(iloc,jloc), dl_ymin, dl_gphi(iloc,jloc) ) ! in km

       ! typical grid size (diagonal) in the vicinity of nearest point
       dl_emax= MAX(dl_e1(iloc,jloc), dl_e2(iloc,jloc))/1000.*SQRT(2.) ! in km

       IF (dl_dis  > dl_emax ) THEN
          zglamfound = dl_glam(iloc,jloc) ; IF (zglamfound > 180.)  zglamfound=zglamfound - 360.

          !         PRINT 9000, 'Long= ',zglamfound,' Lat = ',dl_gphi(iloc,jloc) , iloc, jloc 
          !         PRINT *,' Algorithm does''nt converge ', dl_dis

          IF ( initer >= jp_itermax ) THEN
             !            PRINT *, ' no convergence after ', jp_itermax,' iterations'
             iloc     = -1000
             jloc     = -1000
             ll_again = .FALSE.
          ELSE
             ll_again = .TRUE.
             initer   = initer +1
             jloc     = (initer -1)* ipjglo/initer
             iloc     = (initer -1)* ipiglo/jp_itermax
          ENDIF
       ELSE
          IF ( ll_verbose ) THEN
             PRINT '("#  dl_dis= ",f8.3," km")', dl_dis
          ENDIF
          ll_again = .FALSE.
       END IF
    END DO

    IF (ll_bnd ) THEN
       IF (ll_verbose) WRITE (*,*)'Point  Out of domain or on boundary'
    ELSE
       imin=iloc
       jmin=jloc
    ENDIF

    ! deal with xmax, ymax
    IF (  pxmin == pxmax .AND. pymin == pymax ) THEN
       ! job already done with first point
       imax=imin
       jmax=jmin
    ELSE
       ll_again = .TRUE.
       initer = 1
       iloc=ipiglo/2 ; jloc=ipjglo/2

       DO WHILE (ll_again)
          CALL NearestPoint(dl_xmax, dl_ymax, ipiglo, ipjglo, dl_glam, dl_gphi, iloc, jloc, ll_bnd)

          ! distance between the target point and the nearest point
          dl_dis = dist(dl_xmax, dl_glam(iloc,jloc), dl_ymax, dl_gphi(iloc,jloc) ) ! in km

          ! typical grid size (diagonal) in the vicinity of nearest point
          dl_emax = MAX(dl_e1(iloc,jloc),dl_e2(iloc,jloc))/1000.*SQRT(2.) ! in km

          IF (dl_dis >  dl_emax ) THEN
             zglamfound=dl_glam(iloc,jloc) ; IF (zglamfound > 180.)  zglamfound=zglamfound -360.

             !            PRINT 9000, 'Long= ',zglamfound,' Lat = ',dl_gphi(iloc,jloc), iloc, jloc
             !            PRINT *,' Algorithm does''nt converge ', dl_dis

             IF ( initer >= jp_itermax ) THEN
                !               PRINT *, ' no convergence after ', jp_itermax,' iterations'
                iloc     = -1000
                jloc     = -1000
                ll_again = .FALSE.
             ELSE
                ll_again = .TRUE.
                initer   = initer +1
                jloc     = (initer -1)* ipjglo/initer
                iloc     = (initer -1)* ipiglo/jp_itermax
             ENDIF
          ELSE
             IF ( ll_verbose ) THEN
                PRINT '("#  dl_dis= ",f8.3," km")', dl_dis
             ENDIF
             ll_again = .FALSE.
          END IF
       END DO
       IF (ll_bnd ) THEN
          IF (ll_verbose) WRITE (*,*) 'Point  Out of domain or on boundary'
       ELSE
          imax=iloc
          jmax=jloc
       ENDIF
    ENDIF

    IF (ll_verbose) PRINT 9001, imin, imax, jmin, jmax

    kimin   = imin ; kimax = imax ; kjmin   = jmin ; kjmax = jmax
    IF ( ll_bnd .OR. imin < 0 ) THEN
       zglamin= -9999. ; zglamax = -9999
       zgphmin= -9999. ; zgphmax = -9999
    ELSE
       zglamin = dl_glam(imin,jmin)  ; zglamax = dl_glam(imax,jmax)
       zgphmin = dl_gphi(imin,jmin)  ; zgphmax = dl_gphi(imax,jmax)
    ENDIF

    IF ( zglamin > 180 ) zglamin=zglamin-360.
    IF ( zglamax > 180 ) zglamax=zglamax-360.

    IF ( ll_verbose) PRINT 9002, zglamin, zglamax, zgphmin, zgphmax

    IF ( PRESENT (plonmin) ) plonmin=zglamin
    IF ( PRESENT (plonmax) ) plonmax=zglamax
    IF ( PRESENT (platmin) ) platmin=zgphmin
    IF ( PRESENT (platmax) ) platmax=zgphmax

9000 FORMAT(a,f8.2,a,f8.2,2i5)
9001 FORMAT(4i10)
9002 FORMAT(4f10.4)

  END SUBROUTINE cdf_findij

  SUBROUTINE  broken_line( kimin, kimax, kjmin, kjmax, & 
       &                     pxx, pyy, knn , kpi, kpj, knormu, knormv)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE broken_line  ***
    !!
    !! ** Purpose :  determine the broken line between 2 points 
    !!
    !! ** Method  :  Compute the equation of the strait line joining
    !!               the 2 points, on the model I,J domain. Compute
    !!               both y(x) and x(y) and work with the minimum sloping
    !!               line. 
    !!
    !!----------------------------------------------------------------------
    INTEGER(KIND=4),            INTENT(in ) :: kimin, kimax, kjmin, kjmax
    INTEGER(KIND=4),            INTENT(out) :: knn
    INTEGER(KIND=4),            INTENT(in ) :: kpi, kpj
    INTEGER(KIND=4), OPTIONAL,  INTENT(out) :: knormu, knormv
    REAL(KIND=4), DIMENSION(kpi+kpj), INTENT(out) :: pxx, pyy

    INTEGER(KIND=4) :: ji, jj, jk   ! dummy loop index
    INTEGER(KIND=4) :: ii0, ij0, ii1, ij1, iitmp, ijtmp
    INTEGER(KIND=4) :: iist, ijst
    INTEGER(KIND=4) :: ii, ij
    INTEGER(KIND=4) :: idirx, idiry
    INTEGER(KIND=4) :: ipts

    REAL(KIND=4)    :: zxi0, zyj0, zxi1, zyj1
    REAL(KIND=4)    :: zai, zbi
    REAL(KIND=4)    :: zaj, zbj
    REAL(KIND=4)    :: zd

    COMPLEX, DIMENSION(kpi+kpj) :: ylpt           ! array of points coordinates in a section
    COMPLEX                     :: ylpti          ! working point
    !!----------------------------------------------------------------------

    !! Find the broken line between P1 (kimin,kjmin) and P2 (kimax, kjmax)
    ! ... Initialization
    ii0  = kimin ; ij0  = kjmin ; ii1  = kimax ;  ij1 = kjmax
    zxi0 = ii0   ; zyj0 = ij0   ; zxi1 = ii1   ; zyj1 = ij1

    ! compute direction of integrations and signs
    !The transport across the section is the dot product of
    !integral(line){(Mx,My)*dS}
    !Mx=integral(u*dz)  My=integral(v*dz)) and dS=(dy,-dx)}

    !By defining the direction of the integration as
    idirx = SIGN(1,ii1-ii0) !positive to the east or if ii1=ii0
    idiry = SIGN(1,ij1-ij0) !positive to the north or if ij1=ij0

    !Then dS=(e2u*idiry,-e1v*idirx)
    !This will produce the following sign convention:
    !    West-to-est line (dx>0, dy=0)=> -My*dx (-ve for a northward flow)
    !    South-to-north   (dy>0, dx=0)=>  Mx*dy (+ve for an eastward flow)
    IF ( PRESENT (knormu) ) knormu =  idiry
    IF ( PRESENT (knormv) ) knormv = -idirx

    ! .. Compute equation:  ryj = zaj rxi + zbj [valid in the (i,j) plane]
    IF ( (zxi1 -zxi0) /=  0 ) THEN
       zaj = (zyj1 - zyj0 ) / (zxi1 -zxi0)
       zbj = zyj0 - zaj * zxi0
    ELSE
       zaj = 10000.  ! flag value
       zbj = 0.
    END IF
    ! .. Compute equation:  rxi = zai ryj + zbi [valid in the (i,j) plane]
    IF ( (zyj1 -zyj0) /=  0 ) THEN
       zai = (zxi1 - zxi0 ) / ( zyj1 -zyj0 )
       zbi = zxi0 - zai * zyj0
    ELSE
       zai = 10000. ! flag value
       zbi = 0.
    END IF

    ! ..  Compute the integer pathway: a succession of F points
    ipts=0
    ! .. Chose the strait line with the smallest slope
    IF (ABS(zaj) <=  1 ) THEN
       ! ... Here, the best line is y(x)
       ! ... If ii1 < ii0 swap points [ always describe section from left to right ]
       IF (ii1 <  ii0 ) THEN
          iitmp = ii0   ; ijtmp = ij0
          ii0   = ii1   ; ij0   = ij1
          ii1   = iitmp ; ij1   = ijtmp
       END IF

       ! iist,ijst is the grid offset to pass from F point to either U/V point
       IF ( ij1 >= ij0 ) THEN     ! line heading NE
          iist = 1 ; ijst = 1
       ELSE                       ! line heading SE
          iist = 1 ; ijst = 0
       END IF

       ! ... compute the nearest ji point on the line crossing at ji
       DO ji=ii0, ii1
          ipts=ipts+1
          IF (ipts > kpi+kpj) THEN ; PRINT *,' ... ERROR in broken_line : ipts > kpi+kpj !' ; STOP 96
          ENDIF
          ij=NINT(zaj*ji + zbj )
          ylpt(ipts) = CMPLX(ji,ij)
       END DO
    ELSE
       ! ... Here, the best line is x(y)
       ! ... If ij1 < ij0 swap points [ always describe section from bottom to top ]
       IF (ij1 <  ij0 ) THEN
          iitmp = ii0   ; ijtmp = ij0
          ii0   = ii1   ; ij0   = ij1
          ii1   = iitmp ; ij1   = ijtmp
       END IF

       ! iist,ijst is the grid offset to pass from F point to either U/V point
       IF ( ii1 >=  ii0 ) THEN
          iist = 1 ; ijst = 1
       ELSE
          iist = 0 ; ijst = 1
       END IF

       ! ... compute the nearest ji point on the line crossing at jj
       DO jj=ij0,ij1
          ipts=ipts+1
          IF (ipts > kpi+kpj) THEN ; PRINT *,' ... ERROR in broken_line : ipts > kpi+kpj !' ; STOP 96
          ENDIF
          ii=NINT(zai*jj + zbi)
          ylpt(ipts) = CMPLX(ii,jj)
       END DO
    END IF
    !!
    !! Look for intermediate points to be added.
    !  ..  The final positions are saved in pxx,pyy
    pxx(1) = REAL(ylpt(1))
    pyy(1) = IMAG(ylpt(1))
    knn     = 1

    DO jk=2,ipts
       ! .. distance between 2 neighbour points
       zd=ABS(ylpt(jk)-ylpt(jk-1))
       ! .. intermediate points required if zd > 1
       IF ( zd > 1 ) THEN
          CALL interm_pt(ylpt, jk, zai, zbi, zaj, zbj, ylpti)
          knn=knn+1
          IF (knn > kpi+kpj) THEN ; PRINT *,' ... ERROR in broken_line : knn > kpi+kpj !' ; STOP 96
          ENDIF
          pxx(knn) = REAL(ylpti)
          pyy(knn) = IMAG(ylpti)
       END IF
       knn=knn+1
       IF (knn > kpi+kpj) THEN ; PRINT *,' ... ERROR in broken_line : knn > kpi+kpj !' ; STOP 96
       ENDIF
       pxx(knn) = REAL(ylpt(jk))
       pyy(knn) = IMAG(ylpt(jk))
    END DO

  END SUBROUTINE broken_line

  ! Private subroutines and functions
  SUBROUTINE NearestPoint(ddlon, ddlat, kpi, kpj, ddlam, ddphi, kpiloc, kpjloc, ld_bnd)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE NearestPoint  ***
    !!
    !! ** Purpose : Compute the positions of the nearest i,j in the grid
    !!              from the given longitudes and latitudes
    !!
    !! ** Method  : Starts on the middle of the grid, search in a 20x20 box,
    !!              and move the box in the direction where the distance 
    !!              between the box and the point is minimum.
    !!              Iterates ...
    !!              Stops when the point is outside the grid.
    !!
    !! References : P.A. Darbon and A. de Miranda acknowledged for this 
    !!              clever algorithm developped in CLIPPER.
    !!----------------------------------------------------------------------
    REAL(KIND=8),                     INTENT(in) :: ddlon, ddlat      !: lon and lat of target point
    INTEGER(KIND=4),                 INTENT (in) :: kpi, kpj          !: grid size
    REAL(KIND=8), DIMENSION(kpi,kpj), INTENT(in) :: ddlam, ddphi      !: model grid layout
    INTEGER(KIND=4),              INTENT (inout) :: kpiloc, kpjloc    !: nearest point location
    LOGICAL                                      :: ld_bnd            !: reach boundary flag

    INTEGER(KIND=4)                              :: ji, jj
    INTEGER(KIND=4), PARAMETER                   :: jp_blk=10
    INTEGER(KIND=4)                              :: ii0, ij0
    INTEGER(KIND=4)                              :: ii1, ij1
    REAL(KIND=4)                                 :: zdist
    REAL(KIND=4)                                 :: zdistmin, zdistmin0
    LOGICAL, SAVE                                :: ll_bndcell, ll_first=.TRUE.
    !!----------------------------------------------------------------------
    IF ( ll_first ) THEN
       kpiloc = kpi/2 ; kpjloc = kpj/2    ! seek from the middle of domain
       ll_first=.FALSE.
    ENDIF

    zdistmin=1000000. ; zdistmin0=1000000.
    ii0 = kpiloc      ; ij0 = kpjloc
    ll_bndcell=.TRUE. ; ld_bnd=.FALSE.

    ! loop until found or boundary reach
    DO  WHILE ( ll_bndcell .AND. .NOT. ld_bnd )
       ii0 = kpiloc - jp_blk ;  ii1 = kpiloc + jp_blk
       ij0 = kpjloc - jp_blk ;  ij1 = kpjloc + jp_blk

       ! search only the inner domain
       IF (ii0 <= 0 ) ii0 = 2
       IF (ii1 > kpi) ii1 = kpi - 1
       IF (ij0 <= 0 ) ij0 = 2
       IF( ij1 > kpj) ij1 = kpj - 1

       ! within a block jp_blk+1 x jp_blk+1:
       DO jj=ij0,ij1
          DO ji=ii0,ii1
             ! compute true distance (orthodromy) between target point and grid point
             zdist    = dist(ddlon, ddlam(ji,jj), ddlat, ddphi(ji,jj) )
             zdistmin = MIN(zdistmin, zdist)
             ! update kpiloc, kpjloc if distance decreases
             IF (zdistmin /=  zdistmin0 ) THEN
                kpiloc=ji
                kpjloc=jj
             ENDIF
             zdistmin0=zdistmin
          END DO
       END DO

       ll_bndcell=.FALSE.
       ! if kpiloc, kpjloc belong to block boundary proceed to next block, centered on kpiloc, kpjloc
       IF (kpiloc == ii0 .OR. kpiloc == ii1) ll_bndcell=.TRUE.
       IF (kpjloc == ij0 .OR. kpjloc == ij1) ll_bndcell=.TRUE.

       ! boundary reach ---> not found
       IF (kpiloc == 2  .OR. kpiloc ==kpi-1) ld_bnd=.TRUE.
       IF (kpjloc == 2  .OR. kpjloc ==kpj-1) ld_bnd=.TRUE.
    END DO

  END SUBROUTINE  NearestPoint


  REAL(KIND=8) FUNCTION dist(ddlona, ddlonb, ddlata, ddlatb)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION dist  ***
    !!
    !! ** Purpose : Compute the distance (km) between
    !!              point A (lona, lata) and B (lonb, latb)  
    !!
    !! ** Method  : Use of double precision is important. Compute the 
    !!              distance along the orthodromy
    !!
    !!----------------------------------------------------------------------
    REAL(KIND=8), INTENT(in) :: ddlata, ddlona, ddlatb, ddlonb

    REAL(KIND=8), SAVE :: dl_latar, dl_latbr, dl_lonar, dl_lonbr
    REAL(KIND=8)       :: dl_pds
    REAL(KIND=8), SAVE :: dl_ux, dl_uy, dl_uz
    REAL(KIND=8)       :: dl_vx, dl_vy, dl_vz
    REAL(KIND=8), SAVE :: dl_prevlat=-1000.d0
    REAL(KIND=8), SAVE :: dl_prevlon=-1000.d0
    REAL(KIND=8), SAVE :: dl_r, dl_pi, dl_conv

    LOGICAL :: ll_first=.TRUE.
    !!----------------------------------------------------------------------
    ! initialise some values at first call
    IF ( ll_first ) THEN
       ll_first = .FALSE.
       ! constants
       dl_pi   = ACOS(-1.d0)
       dl_conv = dl_pi/180.d0  ! for degree to radian conversion
       ! Earth radius
       dl_r    = (6378.137d0+6356.7523d0)/2.0d0 ! km
    ENDIF

    ! compute these term only if they differ from previous call
    IF ( ddlata /= dl_prevlat .OR. ddlona /= dl_prevlon) THEN
       dl_latar   = ddlata*dl_conv
       dl_lonar   = ddlona*dl_conv
       dl_ux      = COS(dl_lonar)*COS(dl_latar)
       dl_uy      = SIN(dl_lonar)*COS(dl_latar)
       dl_uz      = SIN(dl_latar)
       dl_prevlat = ddlata
       dl_prevlon = ddlona
    ENDIF

    dl_latbr = ddlatb*dl_conv
    dl_lonbr = ddlonb*dl_conv
    dl_vx    = COS(dl_lonbr)*COS(dl_latbr)
    dl_vy    = SIN(dl_lonbr)*COS(dl_latbr)
    dl_vz    = SIN(dl_latbr)

    dl_pds   = dl_ux*dl_vx + dl_uy*dl_vy + dl_uz*dl_vz

    IF (dl_pds >= 1.) THEN
       dist = 0.
    ELSE
       dist = dl_r*ACOS(dl_pds)
    ENDIF

  END FUNCTION dist


  SUBROUTINE interm_pt (ydpt, kk, pai, pbi, paj, pbj, ydpti)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE nterm_pt  ***
    !!
    !! ** Purpose : Find the best intermediate points on a pathway.
    !!
    !! ** Method  : ydpt : complex vector of the positions of the nearest points
    !!               kk  : current working index
    !!          pai, pbi : slope and original ordinate of x(y)
    !!          paj, pbj : slope and original ordinate of y(x)
    !!             ydpti : Complex holding the position of intermediate point 
    !!
    !! ** Reference : 19/07/1999 : J.M. Molines in Clipper
    !!----------------------------------------------------------------------
    COMPLEX, DIMENSION(:), INTENT(in ) :: ydpt
    COMPLEX,               INTENT(out) :: ydpti
    REAL(KIND=4),          INTENT(in ) :: pai, pbi, paj, pbj
    INTEGER(KIND=4),       INTENT(in ) :: kk
    ! ... local
    COMPLEX                            :: ylptmp1, ylptmp2
    REAL(KIND=4)                       :: za0, zb0
    REAL(KIND=4)                       :: za1, zb1
    REAL(KIND=4)                       :: zd1, zd2
    REAL(KIND=4)                       :: zxm, zym
    REAL(KIND=4)                       :: zxp, zyp
    !!----------------------------------------------------------------------
    ! ... Determines whether we use y(x) or x(y):
    IF (ABS(paj) <=  1) THEN
       ! .....  use y(x)
       ! ... possible intermediate points:
       ylptmp1=ydpt(kk-1)+(1.,0.)                 ! M1 
       ylptmp2=ydpt(kk-1)+CMPLX(0.,SIGN(1.,paj))  ! M2
       !
       ! ... M1 is the candidate point:
       zxm=REAL(ylptmp1)
       zym=IMAG(ylptmp1)
       za0=paj
       zb0=pbj
       !
       za1=-1./za0
       zb1=zym - za1*zxm
       ! ... P1 is the projection of M1 on the strait line
       zxp=-(zb1-zb0)/(za1-za0)
       zyp=za0*zxp + zb0
       ! ... zd1 is the distance M1P1
       zd1=(zxm-zxp)*(zxm-zxp) + (zym-zyp)*(zym-zyp)
       !
       ! ... M2 is the candidate point:
       zxm=REAL(ylptmp2)
       zym=IMAG(ylptmp2)
       za1=-1./za0
       zb1=zym - za1*zxm
       ! ... P2 is the projection of M2 on the strait line
       zxp=-(zb1-zb0)/(za1-za0)
       zyp=za0*zxp + zb0
       ! ... zd2 is the distance M2P2
       zd2=(zxm-zxp)*(zxm-zxp) + (zym-zyp)*(zym-zyp)
       ! ... chose the smallest (zd1,zd2)
       IF (zd2 <=  zd1) THEN
          ydpti=ylptmp2   ! use M2
       ELSE
          ydpti=ylptmp1   ! use M1
       END IF
       !
    ELSE   
       ! ...  use x(y)
       ! ... possible intermediate points:
       ylptmp1=ydpt(kk-1)+CMPLX(SIGN(1.,pai),0.)  ! M1
       ylptmp2=ydpt(kk-1)+(0.,1.)                 ! M2
       ! 
       ! ... M1 is the candidate point:
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
END MODULE cdftools
