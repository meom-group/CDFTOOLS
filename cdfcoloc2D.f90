PROGRAM cdfcoloc2D
  !!-------------------------------------------------------------------
  !!               ***  PROGRAM cdfcoloc2D  ***
  !!
  !!  **  Purpose  :  colocalisation for Greg Holloway
  !!  
  !!  **  Method   :  Use the weight files computed with cdfweight
  !!
  !! history ;
  !!  Original :  J.M. Molines (16 may 2007 )
  !!-------------------------------------------------------------------
  !!  $Rev: 62 $
  !!  $Date: 2007-05-18 16:31:17 +0200 (Fri, 18 May 2007) $
  !!  $Id: cdfcoloc.f90 62 2007-05-18 14:31:17Z molines $
  !!--------------------------------------------------------------
  !! * Modules used
  USE cdfio

  !! * Local variables
  IMPLICIT NONE
  INTEGER, PARAMETER :: jptyp=1                       !: number of type to produce ( look to ctype)
  INTEGER :: narg, iargc
  INTEGER :: ji,jj, jk, jid, jtyp                     !: dummy loop index
  INTEGER :: i1, j1, i2, j2, i3, j3, i4, j4, k1, k2   !: working integers
  INTEGER :: nid = 0                !: mooring counter initialize to 0
  INTEGER :: npiglo, npjglo, npk    !: grid size of the model 
  INTEGER :: npkv                   !: vertical dimension of the target variable (either 1 (2D) or npk (3D)
  INTEGER :: numbin=20, numout=30, numskip=31   !: logical unit for I/O files other than NetCdf
  REAL(KIND=8) :: zmin  

  ! variables in the weight file, 1 record per mooring
  INTEGER :: id
  INTEGER :: imin,  jmin, kmin      !: location of horizontal nearest point, vertical above target.
  INTEGER :: iquadran               !: grid sector from 1 to 4 (clockwise, 1=NE) in wich target point
  !  is located with respect to nearest point
  REAL(KIND=8)  :: xmin, ymin
  REAL(KIND=8)  :: alpha, beta, gamma, hN, scale

  REAL(KIND=4)                      :: dep                 !: Working variables
  REAL(KIND=4)                      :: vup, vdo, wup, wdo  !: Working variables
  REAL(KIND=8), DIMENSION(:,:)  , ALLOCATABLE :: v2d,  e   !: 2D working variable and horizontal metric
  REAL(KIND=8), DIMENSION(:,:)  , ALLOCATABLE :: vinterp   !: result array (nid,jptyp)

  INTEGER(KIND=2), DIMENSION(:,:), ALLOCATABLE :: mask   !: 3D working mask

  ! file name
  CHARACTER(LEN=80) :: coord='coordinates.nc', czgr='mesh_zgr.nc', cmask='mask.nc', cfilout='izb.txt'
  CHARACTER(LEN=80) :: cfilskip='izb_skip.txt'
  CHARACTER(LEN=80) :: cweight,cweight_root, cgridt, cgridu, cgridv, cfil
  ! Variable type and name
  CHARACTER(LEN=80) :: cctyp, cvar, cvmask    !: current mooring
  CHARACTER(LEN=3), DIMENSION(jptyp) :: ctype !:  all jptyp defined there

  !!
  ctype=(/'H'/)
  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg /= 2  ) THEN
     PRINT *,' Usage : cdfcoloc2D  root_weight 2D_T_file '
     PRINT *,' return an ascii file izb.txt (N x 3 )'
     PRINT *,' Example : cdfcoloc2D weight etopo1.nc '
     PRINT *,' coordinates.nc,  required in local directory'
     STOP
  ENDIF

  CALL getarg (1, cweight_root )
  CALL getarg (2, cgridt )

  npiglo= getdim (cgridt,'lon')
  npjglo= getdim (cgridt,'lat')
  PRINT *, 'NPIGLO = ', npiglo
  PRINT *, 'NPJGLO = ', npjglo

  ALLOCATE (mask(npiglo, npjglo) )
  ALLOCATE (v2d(npiglo, npjglo), e(npiglo,npjglo) )

  ! open the weight (T) file for counting the moorings ( assumed that all
  !    weight file ( T U V ) have  the same number of stations.
  WRITE(cweight,'(a,a,".bin")') TRIM(cweight_root), '_T'
  OPEN(numbin, FILE=cweight,FORM='unformatted')
  ! Determine the number of records in the weight file
  DO WHILE (1 /= 2 )   ! not able to use iostat here ... sorry 
     READ(numbin, END=100)
     nid=nid+1
  END DO
100 CONTINUE
  PRINT *, nid ,' stations to process...'
  CLOSE (numbin)

  ! allocate result array
  ALLOCATE ( vinterp(nid,jptyp) )

  ! loop on all variables to collocate
  DO jtyp=1,jptyp
     cctyp=TRIM(ctype(jtyp))

     ! depending upon the type, set the weigth file, variable name, mask variable, data file
     !  vertical dimension of output variable and a scale factor
     SELECT CASE ( cctyp)
     CASE ('H')          ! Bottom topography
        WRITE(cweight,'(a,a,".bin")') TRIM(cweight_root), '_T'
        cvar='z' 
        cfil=cgridt
        npkv=1
        scale=1.
     END SELECT

     ! Now enter the generic processing
     PRINT *,'START coloc2D for ', TRIM(cctyp)
     OPEN(numbin, FILE=cweight,FORM='unformatted')

     IF (cfil /= 'none' ) THEN   ! read data (except for Sx and Sy )
           v2d(:,:)=getvar(cfil,cvar,1, npiglo,npjglo)
     ENDIF

     ! compute corresponding mask
     mask = 1
     WHERE ( v2d >= 0 ) mask = 0

     DO jid=1,nid

        READ(numbin) id, ymin, xmin, dep ,imin, jmin, kmin, iquadran, hN, alpha, beta, gamma
        vinterp(jid,jtyp)=interp()
        ! do not scale dummy values
        IF ( vinterp (jid,jtyp) > -99990 )  vinterp (jid,jtyp) = vinterp (jid,jtyp) * scale
     END DO

     CLOSE(numbin)
  END DO   ! Loop on type

  OPEN(numout ,FILE=cfilout)
  OPEN(numskip,FILE=cfilskip)

  WRITE(cweight,'(a,a,".bin")') TRIM(cweight_root), '_T'
  OPEN(numbin, FILE=cweight,FORM='unformatted')
  DO jid=1, nid   ! loop on all stations
     READ(numbin) id, ymin, xmin, dep, imin, jmin, kmin, iquadran, hN
     print *, id,ymin, xmin, dep
     IF ( xmin > 180.) xmin=xmin -360
     ! output only stations with no problems ( vinterp > -99990 )
     zmin=MINVAL(vinterp(jid,:) )
     IF ( zmin > -99990 ) THEN
        WRITE(numout, '(I5, 6f10.1)') id, dep, (vinterp(jid,jtyp),jtyp=1,jptyp)
     ELSE
        ! save discarted stations for control
        WRITE(numskip, '(I5, 6f10.1)') id, dep, (vinterp(jid,jtyp),jtyp=1,jptyp)
     ENDIF
END DO
CLOSE(numbin)
CLOSE(numout)

  PRINT *,' Done.'

CONTAINS

FUNCTION interp ()
 REAL(KIND=8) :: interp

 ! skip out of domain stations
 IF (imin == -1000 ) THEN
    interp=-99999
    RETURN
 ENDIF

 ! choose the 4 interpolation points, according to sector and nearest point (imin, jmin)
 SELECT CASE (iquadran)
 CASE (1)
    i1=imin    ; j1 = jmin
    i2=imin +1 ; j2 = jmin
    i3=imin +1 ; j3 = jmin + 1
    i4=imin    ; j4 = jmin + 1
 CASE (2)
    i1=imin    ; j1 = jmin
    i2=imin    ; j2 = jmin - 1
    i3=imin +1 ; j3 = jmin - 1
    i4=imin +1 ; j4 = jmin
 CASE (3)
    i1=imin    ; j1 = jmin
    i2=imin -1 ; j2 = jmin
    i3=imin -1 ; j3 = jmin - 1
    i4=imin    ; j4 = jmin - 1
 CASE (4)
    i1=imin    ; j1 = jmin
    i2=imin    ; j2 = jmin + 1
    i3=imin -1 ; j3 = jmin + 1
    i4=imin -1 ; j4 = jmin
 END SELECT

 ! kmin is always above target point
 k1=kmin    ; k2 = kmin + 1

 IF (npkv == 1 ) THEN   ! 2D var, do not take care of vertical interpolation
    k1 = 1  ; k2 = 0 ; wdo = 0.
 ENDIF

 ! compute sum of masked weight above target point
 wup = mask(i1,j1)*(1-alpha)*(1-beta) + mask(i2,j2) * alpha *(1-beta) + &
      &  mask(i3,j3)*alpha*beta         + mask(i4,j4) * (1-alpha)*beta

 ! interpolate with non-masked  values, above target point
 vup = v2d(i1,j1)*(1-alpha)*(1-beta) + v2d(i2,j2) * alpha *(1-beta) + &
      &  v2d(i3,j3)*alpha*beta         + v2d(i4,j4) * (1-alpha)*beta


 IF ( wup == 0 ) THEN       ! all points are masked
    interp=-99999.
 ELSE IF ( wdo == 0 ) THEN  ! all points below are masked, or 2D
    interp= vup/wup
 ELSE                       ! general case
    interp= (1 - gamma) * vup/wup + gamma * vdo/wdo
 ENDIF

END FUNCTION interp

SUBROUTINE rotation (pval, pcourse)
  ! This subroutine returns the input vectors on the geographical reference
  ! for vectors in the domain ( only those are processed in the main program)
 REAL(KIND=8), DIMENSION(jptyp), INTENT(inout) :: pval  !: input vectors ( U V Sx Sy )
 REAL(KIND=8), INTENT(in) :: pcourse        !: local direction of the I=cst lines with respect to N (deg)

 ! local variable
 REAL(KIND=8), DIMENSION(2) :: zcomp 
 REAL(KIND=8) :: zconv, zcourse, zpi

 zpi=ACOS(-1.d0) ; zconv=zpi/180.d0
 zcourse=pcourse*zconv  ! heading in radians
 ! Velocity : u=1, v=2
 zcomp(1:2)=pval(1:2)
 pval(1) =  zcomp(1)*COS(zcourse) +zcomp(2)*SIN(zcourse)
 pval(2) = -zcomp(1)*SIN(zcourse) +zcomp(2)*COS(zcourse)

 ! Bottom slope : Sx=3, Sy=4
 zcomp(1:2)=pval(3:4)
 pval(3) =  zcomp(1)*COS(zcourse) +zcomp(2)*SIN(zcourse)
 pval(4) = -zcomp(1)*SIN(zcourse) +zcomp(2)*COS(zcourse)

END SUBROUTINE rotation

END PROGRAM cdfcoloc2D
