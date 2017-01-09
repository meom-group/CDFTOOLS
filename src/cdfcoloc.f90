PROGRAM cdfcoloc
  !!======================================================================
  !!                     ***  PROGRAM  cdfcoloc  ***
  !!=====================================================================
  !!  ** Purpose : Colocates model values on data points. The 3D or 2D
  !!               position of the points are already in the corresponding
  !!               weight file. (Bilinear interpolation).
  !!
  !!  ** Method  : Use the weight file provided as argument and computed
  !!               with cdfweight
  !!
  !! History : 2.1  : 05/2007  : J.M. Molines : Original code
  !!           3.0  : 03/2011  : J.M. Molines : Doctor norm + Lic.
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !! subroutine  rotation    : perform vector rotation to get geographical
  !!                           vector components
  !! subroutine getfld       : decipher the field list given on the command line
  !! subroutine help_message : list available fields
  !! function interp         : perform bilinear interpolation
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2011, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4), PARAMETER                     :: jptyp=16      ! number of available types
  INTEGER(KIND=4)                                :: ntyp          ! number of type to produce ( look to ctype)
  INTEGER(KIND=4)                                :: ji, jj, jk    ! dummy loop index
  INTEGER(KIND=4)                                :: jid, jtyp     ! dummy loop index
  INTEGER(KIND=4)                                :: idum          ! dummy integer
  INTEGER(KIND=4)                                :: narg, iargc, iarg
  INTEGER(KIND=4)                                :: nid = 0       ! mooring counter initialize to 0
  INTEGER(KIND=4)                                :: npiglo, npjglo ! grid size of the model 
  INTEGER(KIND=4)                                :: npk            ! grid size of the model 
  INTEGER(KIND=4)                                :: npkv          ! vertical dimension of the target variable 
  !                                                            !  (either 1 (2D) or npk (3D)
  INTEGER(KIND=4)                                :: numbin  = 20  ! logical unit for I/O files other than NetCdf
  INTEGER(KIND=4)                                :: numout  = 30  ! logical unit for I/O files other than NetCdf
  INTEGER(KIND=4)                                :: numskip = 31  ! logical unit for I/O files other than NetCdf
  ! variables in the weight file, 1 record per mooring
  INTEGER(KIND=4)                                :: id, idep
  INTEGER(KIND=4)                                :: nimin, njmin  ! location of horizontal nearest point
  INTEGER(KIND=4)                                :: nkmin         ! location vertical above target.
  INTEGER(KIND=4)                                :: nquadran      ! grid sector from 1 to 4 (clockwise, 1=NE) 
  !                                                               ! in which target point is located with respect 
  !                                                               ! to nearest point.
  INTEGER(KIND=4)                                :: nSx, nSy      ! index of the Sx and Sy for rotation
  INTEGER(KIND=4)                                :: nU, nV        ! index of the U and V for rotation
  INTEGER(KIND=2), DIMENSION(:,:,:), ALLOCATABLE :: mask          ! 3D working mask
  
  REAL(KIND=4)                                   :: xmin, ymin, rdep
  REAL(KIND=4)                                   :: vup, vdo, wup, wdo  ! Working variables
  REAL(KIND=4), DIMENSION(:,:,:),    ALLOCATABLE :: v3d           ! 3D  ! working variable (unavoidable)

  REAL(KIND=8)                                   :: dxmin, dymin
  REAL(KIND=8)                                   :: dalpha, dbeta, dgama
  REAL(KIND=8)                                   :: dhN, dscale
  REAL(KIND=8)                                   :: dlmin          
  REAL(KIND=8), DIMENSION(:,:),      ALLOCATABLE :: d2d,  de       ! 2D working variable and horizontal metric
  REAL(KIND=8), DIMENSION(:,:),      ALLOCATABLE :: dinterp        ! result array (nid,jptyp)

  ! file name
  CHARACTER(LEN=256)                             :: cf_out  
  CHARACTER(LEN=256)                             :: cf_skip 
  CHARACTER(LEN=256)                             :: cf_weight 
  CHARACTER(LEN=256)                             :: cf_weight_root
  CHARACTER(LEN=256)                             :: cf_gridt   = 'none'
  CHARACTER(LEN=256)                             :: cf_gridtrc = 'none'
  CHARACTER(LEN=256)                             :: cf_diag    = 'none'
  CHARACTER(LEN=256)                             :: cf_gridu   = 'none'
  CHARACTER(LEN=256)                             :: cf_gridv   = 'none'
  CHARACTER(LEN=256)                             :: cf_bathy   = 'none'
  CHARACTER(LEN=256)                             :: cf_in
  CHARACTER(LEN=256)                             :: cf_weight_t
  CHARACTER(LEN=256)                             :: cf_weight_u
  CHARACTER(LEN=256)                             :: cf_weight_v
  CHARACTER(LEN=256)                             :: cctyp, cvar, cvmask    ! current mooring
  CHARACTER(LEN=256)                             :: cldum   ! dummy char variable for line input
  CHARACTER(LEN=256)                             :: ctmplst0 ! current list of type: separated by ,
  CHARACTER(LEN=256)                             :: cformat  ! ASCII format adapted to ntyp
  CHARACTER(LEN=12), DIMENSION(jptyp)            :: ctype  !  all possible type defined there
  CHARACTER(LEN=12), DIMENSION(:),ALLOCATABLE    :: cltype !  actual type used given as argument

  LOGICAL                                        :: llchk
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  ! exhaustive list of supported field     
  ctype    = (/'T       ','S       ','SSH     ','CFCINV  ','CFCCONC ','PENDEP  ',   &
     &        'MXL     ','MXL01   ','MXLT02  ','ISOTHICK','U       ','V       ', &
     &        'Sx      ','Sy      ','H       ','etopo   '/)
  ctmplst0 = 'U,V,Sx,Sy,H'                 ! default list
  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg <  1  ) THEN
     PRINT *,' usage : cdfcoloc  -w root_weight -t gridT -trc TRC_file ...'
     PRINT *,'          ...  -u gridU -v gridV [-l field list ] [-h]'
     PRINT *,'       -w root_weight  : specify the root name of the weight files'
     PRINT *,'                         _T.bin, _U.bin, or _V.bin will be appended '
     PRINT *,'                         to name if necessary.'
     PRINT *,'       -t gridT file   : name of gridT model file'
     PRINT *,'       -trc TRC file   : name of gridT model file'
     PRINT *,'       -d  diag file   : name of specific diagnostic file '
     PRINT *,'       -u gridU file   : name of gridU model file'
     PRINT *,'       -v gridV file   : name of gridV model file'
     PRINT *,'       -b bathy file   : name of etopo like bathymetric file'
     PRINT *,'       -l field list   : list of fields to be colocated, separated by '',''' 
     PRINT *,'                         Default list is :',TRIM(ctmplst0)
     PRINT *,'       -h              : Give the details of available field to colocate.'
     PRINT *,'     Return a column ascii file id dep fields()'
     PRINT *,  TRIM(cn_fmsk),' is required in local directory'
     PRINT *,  TRIM(cn_fcoo),',',TRIM(cn_fzgr),' are also required for slope computation'
     STOP
  ENDIF

  iarg = 1
  DO WHILE ( iarg <= narg ) 
     CALL getarg ( iarg, cldum ) ; iarg = iarg + 1
     SELECT CASE ( cldum )
     CASE ('-w'   ) ; CALL getarg ( iarg, cf_weight_root ) ; iarg = iarg + 1
     CASE ('-t'   ) ; CALL getarg ( iarg, cf_gridt       ) ; iarg = iarg + 1
     CASE ('-trc' ) ; CALL getarg ( iarg, cf_gridtrc     ) ; iarg = iarg + 1
     CASE ('-d'   ) ; CALL getarg ( iarg, cf_diag        ) ; iarg = iarg + 1
     CASE ('-u'   ) ; CALL getarg ( iarg, cf_gridu       ) ; iarg = iarg + 1
     CASE ('-v'   ) ; CALL getarg ( iarg, cf_gridv       ) ; iarg = iarg + 1
     CASE ('-b'   ) ; CALL getarg ( iarg, cf_bathy       ) ; iarg = iarg + 1
     CASE ('-l'   ) ; CALL getarg ( iarg, ctmplst0       ) ; iarg = iarg + 1
     CASE ('-h'   ) ; CALL help_message 
     CASE DEFAULT   ; PRINT *,TRIM(cldum),' : option not available.' ; STOP
     END SELECT
  ENDDO

  ! intepret ctmplst0 to set up cltype list, ntype and build cf_out file name
  CALL getfld( ) 

  idum = INDEX(TRIM(cf_out),'.') - 1
  IF ( idum == -1 ) THEN
     idum = LEN_TRIM(cf_out)
  ENDIF
  cf_skip = cf_out(1:idum)//'_skip.txt'

  WRITE(cf_weight_t,'(a,a,".bin")') TRIM(cf_weight_root), '_T'
  WRITE(cf_weight_u,'(a,a,".bin")') TRIM(cf_weight_root), '_U'
  WRITE(cf_weight_v,'(a,a,".bin")') TRIM(cf_weight_root), '_V'

  ! Check if required files are available
  llchk = .FALSE.

  IF ( cf_bathy /= 'none' ) THEN ! dealing with special case of etopo file
    llchk = llchk .OR. chkfile(cf_bathy    )
    IF (llchk ) STOP ! missing files
    npiglo = getdim (cf_bathy,'lon')
    npjglo = getdim (cf_bathy,'lat')
    npk    = 1
  ELSE
    llchk = llchk .OR. chkfile(cn_fmsk    )
    IF (llchk ) STOP ! missing files
    npiglo = getdim (cn_fmsk,cn_x)
    npjglo = getdim (cn_fmsk,cn_y)
    npk    = getdim (cn_fmsk,cn_z)
  ENDIF

  ALLOCATE (v3d(npiglo, npjglo, npk), mask(npiglo, npjglo, npk) )
  ALLOCATE (d2d(npiglo, npjglo     ), de(npiglo,npjglo) )


  ! loop on all variables to collocate
  DO jtyp=1,ntyp
     cctyp=TRIM(cltype(jtyp))

     ! depending upon the type, set the weigth file, variable name, mask variable, data file
     !  vertical dimension of output variable and a scale factor
     SELECT CASE ( cctyp)
     CASE ('T')          ! temperature, not used for Greg Holloway output
        cf_weight = cf_weight_t
        cf_in     = cf_gridt
        cvar      = cn_votemper 
        cvmask    = 'tmask'
        npkv      = npk
        dscale    = 1.d0
     CASE ('S')          ! salinity, not used for Greg Holloway output
        cf_weight = cf_weight_t
        cf_in     = cf_gridt
        cvar      = cn_vosaline 
        cvmask    = 'tmask'
        npkv      = npk
        dscale    = 1.d0
     CASE ('SSH')        !  SSH, not used for Greg Holloway output
        cf_weight = cf_weight_t
        cf_in     = cf_gridt
        cvar      = cn_sossheig 
        cvmask    = 'tmask'
        npkv      = 1
        dscale    = 100.d0
     CASE ('CFCINV')     !  CFC inventory, not used for Greg Holloway output
        cf_weight = cf_weight_t
        cf_in     = cf_gridtrc
        cvar      = cn_invcfc
        cvmask    = 'tmask'
        npkv      = 1
        dscale    = 1000000.d0
     CASE ('CFCCONC')     !  CFC inventory, not used for Greg Holloway output
        cf_weight = cf_weight_t
        cf_in     = cf_gridtrc
        cvar      = cn_cfc11
        cvmask    = 'tmask'
        npkv      = 1
        dscale    = 1.d0
     CASE ('PENDEP')     !  CFC penetration depth
        cf_weight = cf_weight_t
        cf_in     = cf_diag
        cvar      = cn_pendep
        cvmask    = 'tmask'
        npkv      = 1
        dscale    = 1.d0
     CASE ('MXL','MXL01' )  !  Mixed layer depth
        cf_weight = cf_weight_t
        cf_in     = cf_gridt
        cvar      = cn_somxl010
        cvmask    = 'tmask'
        npkv      = 1
        dscale    = 1.d0
     CASE ('MXLT02' )  !  Mixed layer depth
        cf_weight = cf_weight_t
        cf_in     = cf_gridt
        cvar      = cn_somxlt02
        cvmask    = 'tmask'
        npkv      = 1
        dscale    = 1.d0
     CASE ('ISOTHICK' )  !  Mixed layer depth
        cf_weight = cf_weight_t
        cf_in     = cf_diag
        cvar      = cn_isothick
        cvmask    = 'tmask'
        npkv      = 1
        dscale    = 1.d0
     CASE ('U')          ! Zonal component of velocity 
        cf_weight = cf_weight_u
        cf_in     = cf_gridu
        cvar      = cn_vozocrtx 
        cvmask    = 'umask'
        npkv      = npk
        dscale    = 100.d0  ! to be cm/s in the output
     CASE ('V')             ! Meridional component of velocity
        cf_weight = cf_weight_v
        cf_in     = cf_gridv
        cvar      = cn_vomecrty 
        cvmask    = 'vmask'
        npkv      = npk
        dscale    = 100.d0  ! to be cm/s in the output
     CASE ('Sx')            ! Zonal component of bottom slope
        cf_weight = cf_weight_u
        cf_in     = 'none' 
        cvar      = 'none'
        cvmask    = 'umask'
        npkv      = 1
        dscale    = 100.d0  ! to be in % in the output
        llchk = llchk .OR. chkfile(cn_fcoo    )
        llchk = llchk .OR. chkfile(cn_fzgr    )
        IF ( llchk ) STOP 
        ! Sx is the i-slope of bottom topog: read adequate metric
        !   and compute it on v3d(:,:,1)
        de(:,:)  = getvar(cn_fcoo, cn_ve1u, 1, npiglo, npjglo)     
        d2d(:,:) = getvar(cn_fzgr, 'hdepw', 1, npiglo, npjglo)     
        DO ji=2, npiglo-1
           v3d(ji,:,1) = (d2d(ji+1,:) - d2d(ji,:)) / de(ji,:)
        END DO
     CASE ('Sy')         ! Meridional component of bottom slope
        cf_weight = cf_weight_v
        cf_in     = 'none' 
        cvar      = 'none'
        cvmask    = 'vmask'
        npkv      = 1
        dscale    = 100.d0  ! to be in % in the output
        llchk = llchk .OR. chkfile(cn_fcoo    )
        llchk = llchk .OR. chkfile(cn_fzgr    )
        IF ( llchk ) STOP 
        ! Sy is the j-slope of bottom topog: read adequate metric
        !   and compute it on v3d(:,:,1)
        de(:,:)  = getvar(cn_fcoo, cn_ve2v, 1, npiglo, npjglo)
        d2d(:,:) = getvar(cn_fzgr, 'hdepw', 1, npiglo, npjglo)
        DO jj=2, npjglo-1
           v3d(:,jj,1) = (d2d(:,jj+1) - d2d(:,jj)) / de(:,jj)
        END DO
     CASE ('H')          ! Bottom topography
        cf_weight = cf_weight_t
        cf_in     = cn_fzgr
        cvar      = 'hdepw' 
        cvmask    = 'tmask'
        npkv      = 1
        dscale    = 1.d0
     CASE ('etopo')       ! Bottom topography from external file
        cf_weight = cf_weight_t
        cf_in     = cf_bathy
        cvar      = 'z' 
        cvmask    = 'none'
        npkv      = 1
        dscale    = 1.d0
     END SELECT

     IF (chkfile (cf_weight) .OR. chkfile( cf_in) )  STOP ! missing file

     ! Now enter the generic processing
     PRINT *,'START coloc for ', TRIM(cctyp)
     IF (jtyp == 1 ) THEN ! count number of station and allocate dinterp
        !  assuming  weight file ( T U V ) have  the same number of stations.
        OPEN(numbin, FILE=cf_weight,FORM='unformatted')
        ! Determine the number of records in the weight file
        DO 
           READ(numbin, END=100)
           nid=nid+1
        END DO
100     CONTINUE
        CLOSE(numbin)

        PRINT *, nid ,' stations to process...'
        ! allocate result array
        ALLOCATE ( dinterp(nid,ntyp) )
     ENDIF

     OPEN(numbin, FILE=cf_weight,FORM='unformatted')

     IF (cf_in /= 'none' ) THEN   ! read data (except for Sx and Sy )
        DO jk=1, npkv
           v3d(:,:,jk)=getvar(cf_in,cvar,jk, npiglo,npjglo)
        END DO
     ENDIF

     ! read corresponding mask
     IF ( cvmask == 'none' ) THEN   ! special case of etopo files ( valid values are < 0 ) 
          mask = 1
          WHERE ( v3d >= 0 ) mask = 0
     ELSE
        DO jk=1, npkv
          mask(:,:,jk)=getvar(cn_fmsk,cvmask,jk, npiglo,npjglo)
        END DO
     ENDIF

     DO jid=1,nid
!       READ(numbin) id, dymin, dxmin, idep ,nimin, njmin, nkmin, nquadran, dhN, dalpha, dbeta, dgama
        READ(numbin) id, ymin, xmin, rdep ,nimin, njmin, nkmin, nquadran, dhN, dalpha, dbeta, dgama
        dinterp(jid,jtyp)=interp()
        ! do not scale dummy values
        IF ( dinterp (jid,jtyp) > -99990.d0 )  dinterp (jid,jtyp) = dinterp (jid,jtyp) * dscale
     END DO

     CLOSE(numbin)
  END DO   ! Loop on type

  OPEN(numout,  FILE=cf_out )
  OPEN(numskip, FILE=cf_skip)

  ! need to re-read some informations from the weight file (idep, dhN)
  cf_weight = cf_weight_t
  OPEN(numbin, FILE=cf_weight, FORM='unformatted')
  DO jid=1, nid   ! loop on all stations
     READ(numbin) id, ymin, xmin, rdep, nimin, njmin, nkmin, nquadran, dhN
     IF ( xmin > 180.0) xmin = xmin - 360.0
     ! output only stations with no problems ( dinterp > -99990 )
     dlmin=MINVAL(dinterp(jid,:) )
     IF ( dlmin > -99990.d0 ) THEN
        ! apply vector rotation to have results on the geographic reference system (N-S, E-W )
        IF ( nSx > 0 ) THEN   ! (Sx, Sy pair )
          CALL rotation( dinterp(jid,nSx), dinterp(jid,nSy), dhN)  
        ENDIF
        IF ( nU > 0 ) THEN    ! (U, V pair)
          CALL rotation( dinterp(jid,nU), dinterp(jid,nV), dhN)  
        ENDIF
        WRITE(numout, cformat) id, rdep, (dinterp(jid,jtyp),jtyp=1,ntyp)
     ELSE
        ! save discarted stations for control
        WRITE(numskip, cformat) id, rdep, (dinterp(jid,jtyp),jtyp=1,ntyp)
     ENDIF
  END DO
  CLOSE(numbin)
  CLOSE(numout)

  PRINT *,' Done.'

CONTAINS

  FUNCTION interp ()
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION interp  ***
    !!
    !! ** Purpose : Perform spatial interpolation
    !!
    !! ** Method  : Use the informations in weigth file to perform
    !!              bilinear interpolation  
    !!
    !!----------------------------------------------------------------------
    REAL(KIND=8)    :: interp               ! return value

    INTEGER(KIND=4) :: ii1, ij1, ii2, ij2   ! working integers
    INTEGER(KIND=4) :: ii3, ij3, ii4, ij4   ! working integers
    INTEGER(KIND=4) :: ik1, ik2             ! working integers
    !!----------------------------------------------------------------------
    ! skip out of domain stations (flagged with nimin = -1000)
    IF (nimin == -1000 ) THEN
       interp=-99999.d0
       RETURN
    ENDIF

    ! choose the 4 interpolation points, according to sector and nearest point (nimin, njmin)
    SELECT CASE (nquadran)
    CASE (1)
       ii1=nimin    ; ij1 = njmin
       ii2=nimin +1 ; ij2 = njmin
       ii3=nimin +1 ; ij3 = njmin + 1
       ii4=nimin    ; ij4 = njmin + 1
    CASE (2)
       ii1=nimin    ; ij1 = njmin
       ii2=nimin    ; ij2 = njmin - 1
       ii3=nimin +1 ; ij3 = njmin - 1
       ii4=nimin +1 ; ij4 = njmin
    CASE (3)
       ii1=nimin    ; ij1 = njmin
       ii2=nimin -1 ; ij2 = njmin
       ii3=nimin -1 ; ij3 = njmin - 1
       ii4=nimin    ; ij4 = njmin - 1
    CASE (4)
       ii1=nimin    ; ij1 = njmin
       ii2=nimin    ; ij2 = njmin + 1
       ii3=nimin -1 ; ij3 = njmin + 1
       ii4=nimin -1 ; ij4 = njmin
    END SELECT

    ! nkmin is always above target point
    ik1 = nkmin    ; ik2 = nkmin + 1

    IF (npkv == 1 ) THEN   ! 2D var, do not take care of vertical interpolation
       ik1 = 1  ; ik2 = 0 ; wdo = 0.
    ENDIF

    ! compute sum of masked weight above target point
    wup = mask(ii1,ij1,ik1)*(1-dalpha)*(1-dbeta) + mask(ii2,ij2,ik1) * dalpha *(1-dbeta) + &
         &  mask(ii3,ij3,ik1)*   dalpha*dbeta      + mask(ii4,ij4,ik1) * (1-dalpha)*dbeta

    ! interpolate with non-masked  values, above target point
    vup = v3d(ii1,ij1,ik1)*(1-dalpha)*(1-dbeta) + v3d(ii2,ij2,ik1) * dalpha *(1-dbeta) +   &
         &  v3d(ii3,ij3,ik1)*    dalpha*dbeta     + v3d(ii4,ij4,ik1) * (1-dalpha)*dbeta

    IF (ik2 /= 0 ) THEN   ! for 3D variables
       ! compute sum of masked weight below target point
       wdo = mask(ii1,ij1,ik2)*(1-dalpha)*(1-dbeta) + mask(ii2,ij2,ik2) * dalpha *(1-dbeta) + &
            &  mask(ii3,ij3,ik2)*   dalpha*dbeta      + mask(ii4,ij4,ik2) * (1-dalpha)*dbeta

       ! interpolate with non-masked  values, below target point
       vdo = v3d(ii1,ij1,ik2)*(1-dalpha)*(1-dbeta) + v3d(ii2,ij2,ik2) * dalpha *(1-dbeta) + &
            &  v3d(ii3,ij3,ik2)*dalpha*dbeta       + v3d(ii4,ij4,ik2) * (1-dalpha)*dbeta
    ENDIF

    IF ( wup == 0 ) THEN       ! all points are masked
       interp=-99999.d0
    ELSE IF ( wdo == 0 ) THEN  ! all points below are masked, or 2D
       interp= vup/wup
    ELSE                       ! general case
       interp= (1 - dgama) * vup/wup + dgama * vdo/wdo
    ENDIF

  END FUNCTION interp

  SUBROUTINE rotation (ddu, ddv, ddcourse)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE rotation  ***
    !!
    !! ** Purpose : This subroutine returns the input vectors on the 
    !!              geographical reference
    !!
    !! ** Method  : Projection acording to ddcourse (heading)
    !!
    !!----------------------------------------------------------------------
    REAL(KIND=8), INTENT(inout) :: ddu      ! input u component (along I)
    REAL(KIND=8), INTENT(inout) :: ddv      ! input v component (along J)
    REAL(KIND=8), INTENT(in   ) :: ddcourse ! local direction of the I=cst lines 
    !                                                         ! with respect to N (deg).
    REAL(KIND=8)                :: dlu   ! Local working variables
    REAL(KIND=8)                :: dlv   ! Local working variables
    REAL(KIND=8)                :: dlconv   !      " 
    REAL(KIND=8)                :: dlcourse !      "
    REAL(KIND=8)                :: dlpi     !      "
    !!----------------------------------------------------------------------

    dlpi = ACOS(-1.d0) ; dlconv = dlpi/180.d0
    dlcourse = ddcourse*dlconv  ! heading in radians
    dlu = ddu  ; dlv = ddv 

    ddu    =  dlu*COS(dlcourse) +dlv*SIN(dlcourse)
    ddv    = -dlu*SIN(dlcourse) +dlv*COS(dlcourse)
  END SUBROUTINE rotation

  SUBROUTINE getfld ()
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE getfld  ***
    !!
    !! ** Purpose : decipher ctmplst : looking for ',' separating field
    !!              count the number of field in ctmplst0 : ntyp  
    !!              Set up pairing for vector components.
    !!              Initialize format output
    !!
    !! ** Method  :   use global variables 
    !!
    !!----------------------------------------------------------------------
    INTEGER(KIND=4)    :: jt
    CHARACTER(LEN=256) :: cltmplst
    !!----------------------------------------------------------------------
    
    cltmplst = ctmplst0
    ntyp = 1
    idum = INDEX(cltmplst,',')

    DO WHILE ( idum > 0 )
       cltmplst=cltmplst(idum+1:)
       idum=INDEX(cltmplst,',')
       ntyp = ntyp + 1
    ENDDO

    ALLOCATE (cltype(ntyp) )
    ! populates cltype with individual field
    cltmplst = ctmplst0
    DO jtyp = 1, ntyp
       idum=INDEX(cltmplst,',')
       IF (idum == 0 ) THEN
          cltype(jtyp) = TRIM(cltmplst)
       ELSE
          cltype(jtyp) = cltmplst(1:idum-1)
       ENDIF
       cltmplst=cltmplst(idum+1:)
    ENDDO

    ! check if all fields are supported:
    DO jtyp=1, ntyp
       DO jt =1 , jptyp
          IF ( cltype(jtyp) == TRIM(ctype(jt)) )  EXIT
       ENDDO
       IF ( jt == jptyp + 1 ) THEN
          PRINT *, 'ERROR in field list :', TRIM(cltype(jtyp) ),' not supported'
          STOP
       ENDIF
    ENDDO
    
    ! locate pairing for vector variables
    nSx = -1 ; nSy = -1
    nU  = -1 ; nV  = -1
    DO jtyp = 1, ntyp
      IF (  cltype(jtyp) == 'Sx' ) nSx = jtyp
      IF (  cltype(jtyp) == 'Sy' ) nSy = jtyp
      IF (  cltype(jtyp) == 'U'  ) nU  = jtyp
      IF (  cltype(jtyp) == 'V'  ) nV  = jtyp
    END DO

    IF ( nSx * nSy < 0 ) THEN 
       PRINT *, ' You must specify both Sx and Sy'
       PRINT *, ' in order to perform rotation'
       STOP
    ENDIF
    IF ( nU  * nV  < 0 ) THEN 
       PRINT *, ' You must specify both U and V'
       PRINT *, ' in order to perform rotation'
       STOP
    ENDIF

    ! build output file name
    cf_out='iz'
    DO jtyp=1, ntyp
       cf_out=TRIM(cf_out)//'_'//TRIM(cltype(jtyp))
    ENDDO
    cf_out=TRIM(cf_out)//'.txt'
    
    ! Build output format
    WRITE(cformat,'(a,i2,a)') '(I5,  I6,',ntyp,'e14.6)'

  END SUBROUTINE getfld

  SUBROUTINE help_message ()
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE help_message  ***
    !!
    !! ** Purpose : Print the list of available fields, and describes the
    !!              corresponding required input files  
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=24), DIMENSION(jptyp) :: comments
    CHARACTER(LEN=10), DIMENSION(jptyp) :: crequired
    !!----------------------------------------------------------------------
    PRINT *,' List of available field to process:' 
    PRINT *,'field name  comments                input files'

    
    ! ctype    = (/'T','S','SSH','CFCINV','CFCCONC','PENDEP','MXL','MXL01',
    !              'MXLT02','ISOTHICK','U ','V ','Sx','Sy','H ','etopo'/)
    comments = (/' Potential temperature  ', &
          &      ' Salinity               ', &
          &      ' Sea Surface height     ', &
          &      ' CFC inventory          ', &
          &      ' CFC concentration      ', &
          &      ' Penetration depth      ', &
          &      ' Mixed layer depth s0.01', &
          &      ' Mixed layer depth s0.01', &
          &      ' Mixed layer depth t0.2 ', &
          &      ' Isopycnal thickness    ', &
          &      ' Zonal velocity         ', &
          &      ' Meridional velocity    ', &
          &      ' Zonal bottom slope     ', &
          &      ' Meridional bottom slope', &
          &      ' Local model bathymetry ', &
          &      ' etopo like bathymetry  ' /)
    crequired = (/' -t gridT ',   &
          &       ' -t gridT ',   &
          &       ' -t gridT ',   &
          &       ' -trc TRC ',   &
          &       ' -trc TRC ',   &
          &       ' -d  diag ',   &
          &       ' -t gridT ',   &
          &       ' -t gridT ',   &
          &       ' -t gridT ',   &
          &       ' -d  diag ',   &
          &       ' -u gridU ',   &
          &       ' -v gridV ',   &
          &       ' zgr coord',   &
          &       ' zgr coord',   &
          &       ' zgr      ',   &
          &       ' -b etopo ' /)

    DO jtyp=1, jptyp
      PRINT '( 12a,x,24a,x,10a)', TRIM(ctype(jtyp)), comments(jtyp), crequired(jtyp)
    ENDDO
    STOP

  END SUBROUTINE help_message

END PROGRAM cdfcoloc
