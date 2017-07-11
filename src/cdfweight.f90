PROGRAM cdfweight
  !!======================================================================
  !!                     ***  PROGRAM  cdfweight  ***
  !!=====================================================================
  !!  ** Purpose : Compute a wheight file for further bi-linear colocalisation
  !!               done with cdfcoloc.
  !!
  !!  ** Method  : Use Greg Holloway iyxz.txt file type as input, to specify
  !!               the points to search in the model grid.
  !!               Read the coordinate/mesh_hgr file and look
  !!               for the glam, gphi variables
  !!               Then use a search algorithm to find the corresponding I J
  !!               The point type ( T U V F ) is specified on the command line
  !!               as well as the name of the coordinate/mesh hgr file.
  !!               If -2d option is used, only horizontal weight are produced.
  !!
  !! History : 2.0  : 11/2005  : J.M. Molines : Original code
  !!                : 05/2007  : J.M. Molines : for weight
  !!           3.0  : 03/2011  : J.M. Molines : Doctor norm + Lic.
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !! SUBROUTINE localcoord( palpha, pbeta, plam, pphi)
  !! FUNCTION det(p1,p2,p3,p4)
  !! FUNCTION heading(plona, plonb, plata, platb)
  !!----------------------------------------------------------------------
  USE cdfio
  USE cdftools       ! cdf_find_ij
  USE modutils
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class data_transformation
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                           :: jk                      ! dummy loop counter
  INTEGER(KIND=4)                           :: idum                    ! dummy working integer
  INTEGER(KIND=4)                           :: narg, iargc, iarg       ! Argument management
  INTEGER(KIND=4)                           :: iimin,  ijmin           ! i j position of target point
  INTEGER(KIND=4)                           :: ikloc                   ! k position of target point
  INTEGER(KIND=4)                           :: npiglo, npjglo, npk     ! domain size
  INTEGER(KIND=4)                           :: iquadran                ! quadran
  INTEGER(KIND=4)                           :: numgreg=10              ! logical unit of ASCII input file
  INTEGER(KIND=4)                           :: numbin=20               ! logical unit of BINARY weight file
  INTEGER(KIND=4)                           :: ios                     ! iostat variable
  ! Greg Holloway input data  ( 5 variables)
  INTEGER(KIND=4)                           :: id                      ! station ID
  REAL(KIND=4)                              :: xmin, ymin, rdep        ! longitude, latitude, depth

  REAL(KIND=8)                              :: dl_xmin, dl_ymin
  REAL(KIND=8)                              :: dl_hPp                  ! local maximum metrics
  REAL(KIND=8)                              :: dl_lam0                 ! longitude of grid point ji=1
  REAL(KIND=8)                              :: dl_lamin, dl_phimin     ! coordinates of the nearest point  (NP)
  REAL(KIND=8)                              :: dl_lamN, dl_phiN, dl_hN ! grid point North of NP, true heading from NP
  REAL(KIND=8)                              :: dl_lamE, dl_phiE, dl_hE ! grid point East of NP, true heading from NP
  REAL(KIND=8)                              :: dl_lamS, dl_phiS, dl_hS ! grid point South of NP, true heading from NP
  REAL(KIND=8)                              :: dl_lamW, dl_phiW, dl_hW ! grid point West of NP, true heading from NP
  REAL(KIND=8), DIMENSION(0:4)              :: dl_lami, dl_phii        ! the 4 grid points around target (1-4)
  !  + the target (0)
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dl_lam, dl_phi          ! grid layout and metrics
  REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: dl_dept                 ! vertical depth
  REAL(KIND=8)                              :: dl_hP                   ! true heading of target point from NP
  REAL(KIND=8)                              :: dl_alpha, dl_beta       ! reduced coordinates (0-1) in the NP gridcell
  REAL(KIND=8)                              :: dl_gam                  ! vertical weight

  CHARACTER(LEN=256)                        :: cf_coord, cf_in         ! file names (in)
  CHARACTER(LEN=256)                        :: cf_weight               ! weight file name (out)
  CHARACTER(LEN=256)                        :: ctype, cldum            ! C-grid type point, dummy character

  LOGICAL                                   :: lldebug = .FALSE.       ! verbose/debug flag
  LOGICAL                                   :: ll2d    = .FALSE.       ! 2D field flag
  LOGICAL                                   :: llchk   = .FALSE.       ! for checking missing files
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  ! default values
  cf_coord = cn_fcoo
  ctype    = 'F'

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfweight  -f IN-file [-c COORD-file] [-p C-type] [-2d] [-v] '
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Produces a weight file for further bilinear colocalisation with ' 
     PRINT *,'       cdfcoloc program. It takes the position of the points to be '
     PRINT *,'       colocated into a simple ascii file. '
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f  IN-file   : input file is a iyxz ASCII file, 1 line per point.'  
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-c COORD-file] : coordinate file [',TRIM(cf_coord),']'
     PRINT *,'       [-p C-type    ] : point type on C-grid (either T U V or F ) [',TRIM(ctype),']'
     PRINT *,'       [-2d ]          : tell cdfweight that only 2D weights will be computed.'
     PRINT *,'       [-v ]           : Verbose mode for extra information (debug mode).'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        ',TRIM(cf_coord),' file if not passed as argument.'
     PRINT *,'        If working with 3D files, ',TRIM(cn_fzgr),' is required.'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       binary weight file : weight_point_type.bin'
     PRINT *,'       standard output : almost the same info that is saved in the binary file'
     PRINT *,'                   When using -v option, even more informations !'
     PRINT *,'      '
     PRINT *,'     SEE ALSO : '
     PRINT *,'        cdfcoloc'
     PRINT *,'      '
     STOP 
  ENDIF

  iarg=1
  DO WHILE (iarg <= narg ) 
     CALL getarg(iarg, cldum) ; iarg=iarg+1
     SELECT CASE ( cldum )
     CASE ('-f' ) ; CALL getarg(iarg, cf_in   ) ; iarg=iarg+1 
        ! options
     CASE ('-c' ) ; CALL getarg(iarg, cf_coord) ; iarg=iarg+1 
     CASE ('-p' ) ; CALL getarg(iarg, ctype   ) ; iarg=iarg+1 
     CASE ('-2d') ; ll2d    = .TRUE.
     CASE ('-v' ) ; lldebug = .TRUE.
     CASE DEFAULT ; PRINT *,' ERROR : ',TRIM(cldum),' unknown option.' ; STOP 99
     END SELECT
  END DO

  llchk = llchk .OR. chkfile(cf_in)
  llchk = llchk .OR. chkfile(cf_coord)
  IF ( .NOT. ll2d ) llchk = llchk .OR. chkfile(cn_fzgr)
  IF ( llchk ) STOP 99 ! missing files

  npiglo = getdim (cf_coord,cn_x)
  npjglo = getdim (cf_coord,cn_y)

  IF ( .NOT. ll2d ) THEN 
     npk = getdim (cn_fzgr,cn_z )
     ALLOCATE (dl_dept(npk) )
     ! read depth of model T points (hence U and V)
     dl_dept(:)=getvare3(cn_fzgr, cn_gdept, npk)
  ENDIF

  ALLOCATE (dl_lam(npiglo,npjglo), dl_phi(npiglo,npjglo) )

  ! set name and open  output weight file
  WRITE(cf_weight,'("weight_",a,".bin")') TRIM(ctype)
  OPEN(numbin, FILE=cf_weight,FORM='unformatted')

  SELECT CASE ( ctype )
  CASE ('T' , 't' )
     dl_lam(:,:) = getvar(cf_coord, cn_glamt, 1, npiglo, npjglo)
     dl_phi(:,:) = getvar(cf_coord, cn_gphit, 1, npiglo, npjglo)
  CASE ('U','u' )
     dl_lam(:,:) = getvar(cf_coord, cn_glamu, 1, npiglo, npjglo)
     dl_phi(:,:) = getvar(cf_coord, cn_gphiu, 1, npiglo, npjglo)
  CASE ('V','v' )
     dl_lam(:,:) = getvar(cf_coord, cn_glamv, 1, npiglo, npjglo)
     dl_phi(:,:) = getvar(cf_coord, cn_gphiv, 1, npiglo, npjglo)
  CASE ('F','f' )
     dl_lam(:,:) = getvar(cf_coord, cn_glamf, 1, npiglo, npjglo)
     dl_phi(:,:) = getvar(cf_coord, cn_gphif, 1, npiglo, npjglo)
  CASE DEFAULT
     PRINT *,' ERROR : type of point not known: ', TRIM(ctype)
  END SELECT
  ! work with longitude between 0 and 360 to avoid  the date line.
  WHERE( dl_lam < 0 ) dl_lam(:,:)=dl_lam(:,:)+360.d0

  ! For Orca grid, the longitude of ji=1 is about 70 E
  dl_lam0 = dl_lam(1, npjglo/2)  
  WHERE( dl_lam < dl_lam0 ) dl_lam=dl_lam+360.d0

  OPEN(numgreg,FILE=cf_in)
  ! Greg (Holloway) files are iyxz.txt file
  ios=0
  ! loop for each line of Greg File
  DO WHILE (ios == 0 )
     READ(numgreg,*,iostat=ios) id,ymin,xmin,rdep
     dl_xmin=xmin ; dl_ymin=ymin
     IF( ios == 0 ) THEN  ! EOF not reached
        IF ( .NOT. ll2d ) THEN
           ! Look for vertical position 
           !    ikloc = k index of point above rdep
           ikloc=1
           DO WHILE ( dl_dept(ikloc) <= rdep .AND. ikloc < npk )
              ikloc = ikloc+1
           ENDDO
           ikloc = ikloc -1 ! up one level

           ! compute dl_gam such that Vint= (1-dl_gam) x V(ikloc) + dl_gam x V(ikloc +1)
           dl_gam=(rdep - dl_dept(ikloc))/(dl_dept(ikloc+1)-dl_dept(ikloc) )

           IF (ikloc == npk -1 ) dl_gam = 0.d0

           IF ( dl_gam < 0 ) THEN
              ikloc=1
              dl_gam = 0.d0
           ENDIF
           IF ( dl_gam > 1 ) THEN
              ikloc=npk -1
              dl_gam = 0.d0
           ENDIF
        ELSE
           dl_gam = 0.d0
        ENDIF

        IF ( lldebug) PRINT '("DEP", f8.1,i8,f8.0,f8.4)', dl_dept(ikloc), rdep, dl_dept(ikloc+1), dl_gam

        ! Now deal with horizontal interpolation
        CALL cdf_findij ( xmin, xmin, ymin, ymin, iimin, idum, ijmin, idum, cd_coord=cf_coord, cd_point=ctype)

        IF ( iimin /= -1000 .AND. ijmin /= -1000 ) THEN
           ! Latitude and longitude of the neighbours on the grid
           ! define longitudes between 0 and 360 deg
           dl_lamin = MOD(dl_lam(iimin  ,ijmin  ),360.d0) ; dl_phimin = dl_phi(iimin  ,ijmin  )  ! nearest point
           dl_lamN  = MOD(dl_lam(iimin  ,ijmin+1),360.d0) ; dl_phiN   = dl_phi(iimin  ,ijmin+1)  ! N (grid)
           dl_lamE  = MOD(dl_lam(iimin+1,ijmin  ),360.d0) ; dl_phiE   = dl_phi(iimin+1,ijmin  )  ! E (grid)
           dl_lamS  = MOD(dl_lam(iimin  ,ijmin-1),360.d0) ; dl_phiS   = dl_phi(iimin  ,ijmin-1)  ! S (grid)
           dl_lamW  = MOD(dl_lam(iimin-1,ijmin  ),360.d0) ; dl_phiW   = dl_phi(iimin-1,ijmin  )  ! W (grid)

           ! Compute heading of target point and neighbours from the nearest point
           dl_hP = heading(dl_lamin, dl_xmin, dl_phimin, dl_ymin)  ! target point
           dl_hN = heading(dl_lamin, dl_lamN, dl_phimin, dl_phiN)  ! 'north' on the grid
           dl_hE = heading(dl_lamin, dl_lamE, dl_phimin, dl_phiE)  ! 'east' on the grid
           dl_hS = heading(dl_lamin, dl_lamS, dl_phimin, dl_phiS)  ! 'south' on the grid
           dl_hW = heading(dl_lamin, dl_lamW, dl_phimin, dl_phiW)  ! 'west' on the grid

           ! determine the sector in wich the target point is located: 
           !  ( from 1, to 4 resp. NE, SE, SW, NW  of the grid)
           iquadran = 4
           ! to avoid problem with the GW meridian, pass to -180, 180 when working around GW
           IF ( dl_hP > 180.d0 ) THEN 
              dl_hPp = dl_hP - 360.d0 
              dl_hPp = dl_hP
           ENDIF

           IF ( dl_hN  > dl_hE                       ) dl_hN = dl_hN - 360.d0
           IF ( dl_hPp > dl_hN .AND. dl_hPp <= dl_hE ) iquadran = 1
           IF ( dl_hP  > dl_hE .AND. dl_hP  <= dl_hS ) iquadran = 2
           IF ( dl_hP  > dl_hS .AND. dl_hP  <= dl_hW ) iquadran = 3
           IF ( dl_hP  > dl_hW .AND. dl_hPp <= dl_hN ) iquadran = 4

           dl_lami(0) = xmin     ; dl_phii(0) = ymin          ! fill dl_lami, dl_phii for 0 = target point
           dl_lami(1) = dl_lamin ; dl_phii(1) = dl_phimin     !                           1 = nearest point
           SELECT CASE ( iquadran )         ! point 2 3 4 are counter clockwise in the respective sector
           CASE ( 1 ) 
              dl_lami(2) = dl_lamE ; dl_phii(2) = dl_phiE
              dl_lami(3) = MOD(dl_lam(iimin+1,ijmin+1), 360.d0) ; dl_phii(3) = dl_phi(iimin+1,ijmin+1)
              dl_lami(4) = dl_lamN ; dl_phii(4) = dl_phiN
           CASE ( 2 )
              dl_lami(2) = dl_lamS ; dl_phii(2) = dl_phiS
              dl_lami(3) = MOD(dl_lam(iimin+1,ijmin-1), 360.d0) ; dl_phii(3) = dl_phi(iimin+1,ijmin-1)
              dl_lami(4) = dl_lamE ; dl_phii(4) = dl_phiE
           CASE ( 3 )
              dl_lami(2) = dl_lamW ; dl_phii(2) = dl_phiW
              dl_lami(3) = MOD(dl_lam(iimin-1,ijmin-1), 360.d0) ; dl_phii(3) = dl_phi(iimin-1,ijmin-1)
              dl_lami(4) = dl_lamS ; dl_phii(4) = dl_phiS
           CASE ( 4 )
              dl_lami(2) = dl_lamN ; dl_phii(2) = dl_phiN
              dl_lami(3) = MOD(dl_lam(iimin-1,ijmin+1), 360.d0) ; dl_phii(3) = dl_phi(iimin-1,ijmin+1)
              dl_lami(4) = dl_lamW ; dl_phii(4) = dl_phiW
           END SELECT

           ! resolve a non linear system of equation for dl_alpha and dl_beta 
           !( the non dimensional coordinates of target point)
           CALL localcoord( dl_alpha, dl_beta, dl_lami, dl_phii)
        ELSE   ! point is outside the domaine, put dummy values
           dl_alpha=-1000.d0 ; dl_beta=-1000.d0
        ENDIF

        IF (lldebug) THEN 
           PRINT 9001, id, ymin, xmin, rdep ,iimin, ijmin,  dl_hP, dl_hPp, dl_hN, &
                &         dl_hE, dl_hS, dl_hW, iquadran, dl_alpha, dl_beta
        ENDIF
        ! output both on std output and binary weight file (same info).
        PRINT 9002, id, ymin, xmin, rdep ,iimin, ijmin, ikloc, iquadran, dl_alpha, dl_beta, dl_gam
        WRITE(numbin) id, ymin, xmin, rdep ,iimin, ijmin, ikloc, iquadran, dl_hN, dl_alpha, dl_beta, dl_gam
     ENDIF
  ENDDO
9001 FORMAT(i10, 3f10.4,2i6,6f10.4,I4,2f8.4)
9002 FORMAT(i10, 3f10.4,3i6,I4,3f11.4)
  CLOSE(numbin)

CONTAINS

  SUBROUTINE localcoord( dpalpha, dpbeta, dplam, dpphi)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE localcoord  ***
    !!
    !! ** Purpose :  compute the local coordinate in a grid cell
    !!
    !! ** Method  :  See reference 
    !!
    !! References : from N. Daget Web page :
    !!      http://aton.cerfacs.fr/~daget/TECHREPORT/TR_CMGC_06_18_html/node8.html 
    !!----------------------------------------------------------------------
    REAL(KIND=8), DIMENSION(0:4), INTENT(in)  :: dplam, dpphi
    REAL(KIND=8)                , INTENT(out) :: dpalpha, dpbeta

    INTEGER(KIND=4)              :: itermax=200    ! maximum of iteration 
    INTEGER(KIND=4)              :: iter=0         ! iteration counter
    REAL(KIND=8)                 :: dlalpha=0.d0   ! working variable, initialized to 1rst guess
    REAL(KIND=8)                 :: dlbeta=0.d0    !  ""                  ""
    REAL(KIND=8)                 :: dlresmax=0.001 ! Convergence criteria
    REAL(KIND=8)                 :: dlres          ! residual 
    REAL(KIND=8)                 :: dldeta
    REAL(KIND=8)                 :: dldalp
    REAL(KIND=8)                 :: dldbet
    REAL(KIND=8)                 :: dldlam
    REAL(KIND=8)                 :: dldphi
    REAL(KIND=8)                 :: dl1, dl2, dl3, dl4
    REAL(KIND=8), DIMENSION(2,2) :: dla
    REAL(KIND=8), DIMENSION(0:4) :: dlplam
    !!----------------------------------------------------------------------
    dlplam=dplam       !: save input longitude in working array
    IF ( lldebug ) THEN
       PRINT *,dplam(0), dpphi(0)
       PRINT *,9999,9999
       PRINT *,dplam(1), dpphi(1)
       PRINT *,dplam(2), dpphi(2)
       PRINT *,dplam(3), dpphi(3)
       PRINT *,dplam(4), dpphi(4)
       PRINT *,dplam(1), dpphi(1)
       PRINT *,9999,9999
    ENDIF
    IF ( ABS( dlplam(1) -dlplam(4) ) >= 180.d0 .OR. ABS( dlplam(1) -dlplam(2) ) >=180.d0) THEN
       ! then we are near the 0 deg line and we must work in the frame -180 180 
       WHERE ( dlplam >= 180.d0 ) dlplam=dlplam -360.d0
    ENDIF

    dlres=1000.; dldlam=0.5; dldphi=0.5 ;  dlalpha=0.d0 ; dlbeta=0.d0; iter=0

    DO WHILE (dlres > dlresmax .AND. iter < itermax)
       dl1=(dlplam(2)- dlplam(1) )
       dl2=(dlplam(1) -dlplam(4) )
       dl3=(dlplam(3) -dlplam(2) )

       dla(1,1) =  dl1 + (dl2 + dl3 )* dlbeta 
       dla(1,2) = -dl2 + (dl2 + dl3 )* dlalpha

       dla(2,1) = dpphi(2)-dpphi(1) +  (dpphi(1) -dpphi(4) +dpphi(3) -dpphi(2))* dlbeta 
       dla(2,2) = dpphi(4)-dpphi(1) +  (dpphi(1) -dpphi(4) +dpphi(3) -dpphi(2))* dlalpha

       ! determinant 
       dldeta=det(dla(1,1), dla(1,2), dla(2,1), dla(2,2) )

       ! solution of 
       ! | zdlam |        | zdalp |
       ! |       | =  za .|       |
       ! | zdphi |        | zdbet |
       dldalp=det(dldlam,   dla(1,2), dldphi,   dla(2,2))/dldeta
       dldbet=det(dla(1,1), dldlam,   dla(2,1), dldphi  )/dldeta

       ! compute residual ( loop criteria)
       dlres=SQRT(dldalp*dldalp + dldbet*dldbet )

       ! Compute alpha and beta from 1rst guess :
       dlalpha = dlalpha + dldalp
       dlbeta  = dlbeta  + dldbet

       ! compute corresponding lon/lat for this alpha, beta
       dldlam = dlplam(0) - ((1.-dlalpha)*(1-dlbeta)*dlplam(1) + dlalpha*(1-dlbeta)*dlplam(2) + & 
            &                      dlalpha*dlbeta*dlplam(3) + (1-dlalpha)*dlbeta*dlplam(4))
       dldphi = dpphi(0)  - ((1.-dlalpha)*(1-dlbeta)*dpphi(1)  + dlalpha*(1-dlbeta)*dpphi(2)  + &
            &                      dlalpha*dlbeta*dpphi(3)  + (1-dlalpha)*dlbeta*dpphi(4))

       iter=iter + 1  ! increment iteration counter
    END DO   ! loop until dlres small enough (or itermax reach )

    dpalpha = dlalpha
    dpbeta  = dlbeta
  END SUBROUTINE localcoord

  FUNCTION det(dp1,dp2,dp3,dp4)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION det  ***
    !!
    !! ** Purpose : compute determinant
    !!
    !! ** Method  : just multiply and add  !
    !!
    !!----------------------------------------------------------------------
    REAL(KIND=8), INTENT(in) :: dp1, dp2, dp3, dp4 ! matrix elements
    REAL(KIND=8)             :: det                ! return value
    !!----------------------------------------------------------------------
    det = dp1*dp4 - dp2*dp3
  END FUNCTION det

END PROGRAM cdfweight
