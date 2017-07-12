PROGRAM cdf2regular
  !!======================================================================
  !!                     ***  PROGRAM  cdf2regular  ***
  !!=====================================================================
  !!  ** Purpose : remaps (bin) 3D high resolution (finer than 1x1 deg)
  !!               fields on Levitus-like (ie lon/lat regular) grid.
  !!
  !!  ** Method  : data surface averaging
  !!               It assumes that Levitus grid SW grid cell center 
  !!               is 0.5W,89.5S 
  !!
  !! History : 3.0  : 06/2012  : N. Ferry  : Original code
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  USE cdftools
  USE modutils
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class data_transformation
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                              :: ji, jj, jk, jt     ! dummy loop index
  INTEGER(KIND=4)                              :: jilev,jjlev        ! dummy loop index
  INTEGER(KIND=4)                              :: jvar               ! dummy loop index
  INTEGER(KIND=4)                              :: ii, ij             ! array index (not loop)
  INTEGER(KIND=4)                              :: iilev, ijlev       ! array index (not loop)
  INTEGER(KIND=4)                              :: icount             ! array index (not loop)
  INTEGER(KIND=4)                              :: npiglo, npjglo     ! size of the domain
  INTEGER(KIND=4)                              :: npk, npt           ! size of the domain
  INTEGER(KIND=4)                              :: npilev, npjlev     ! size of the Levitus domain
  INTEGER(KIND=4)                              :: narg, iargc, ijarg ! browse line
  INTEGER(KIND=4)                              :: ncout              ! ncid of output file
  INTEGER(KIND=4)                              :: ierr               ! error status
  INTEGER(KIND=4)                              :: ires               ! resolution in fraction
  INTEGER(KIND=4)                              :: nvars              ! number of variables in the input file
  INTEGER(KIND=4)                              :: iimin, iimax       ! IJ coordinates of the closest points
  INTEGER(KIND=4)                              :: ijmin, ijmax       ! "      "            "
  INTEGER(KIND=4)                              :: imethod=1          ! interpolation method
  INTEGER(KIND=4)                              :: iter_shap=3        ! number of Shapiro iteration
  INTEGER(KIND=4), DIMENSION(:),   ALLOCATABLE :: ipk, id_var        ! levels and varid's of input vars
  INTEGER(KIND=4), DIMENSION(:),   ALLOCATABLE :: ipkout, id_varout  ! levels and varid's of output vars

  REAL(KIND=4)                                 :: zradius=120.       ! Distance (km) for the search bubble (FHZ)
  REAL(KIND=4)                                 :: rlon1, rlon2, rlat1, rlat2, rpos
  REAL(KIND=4)                                 :: gphitmin
  REAL(KIND=4)                                 :: rlev_resol=0.33333333     ! degree
  REAL(KIND=4)                                 :: rlon0=-180         ! degree
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: e1t, e2t           ! horizontal T metrics
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: z_in               ! input field
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: z_fill             ! output 
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: glamt, gphit       ! T longitude latitude
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: rlonlev, rlatlev   ! Levitus grid longitude latitude
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: zbt
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: tmask              ! input mask
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: tmasklev           ! output mask
  REAL(KIND=4), DIMENSION(:),      ALLOCATABLE :: gdept              ! depth axis

  REAL(KIND=8), DIMENSION(:),      ALLOCATABLE :: dtim               ! time counter
  REAL(KIND=8), DIMENSION(:,:),    ALLOCATABLE :: d_out, d_n         ! output field and weighting field

  CHARACTER(LEN=256)                           :: cldum              ! dummy char variable
  CHARACTER(LEN=256)                           :: cf_in              ! input file name
  CHARACTER(LEN=256)                           :: cf_out='regular.nc'! output file name ( output)
  CHARACTER(LEN=256)                           :: cv_nam             ! variable name
  CHARACTER(LEN=256)                           :: ctcalendar         ! time attributes
  CHARACTER(LEN=256)                           :: cttitle            ! time attributes
  CHARACTER(LEN=256)                           :: ctlong_name        ! time attributes
  CHARACTER(LEN=256)                           :: ctaxis             ! time attributes
  CHARACTER(LEN=256)                           :: ctunits            ! time attributes
  CHARACTER(LEN=256)                           :: cttime_origin      ! time attributes

  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_names          ! array of var name
  CHARACTER(LEN=6)                              :: ctyp              ! 'fill' or 'smooth' for shapiro

  TYPE(variable), DIMENSION(:),   ALLOCATABLE  :: stypvar            ! input attributes
  TYPE(variable), DIMENSION(:),   ALLOCATABLE  :: stypvarout         ! output attributes

  LOGICAL                                      :: lchk               ! missing files flag
  LOGICAL                                      :: ltest
  LOGICAL                                      :: l360 = .FALSE.     ! flag for 0-360 layout
  LOGICAL                                      :: lnc4 = .FALSE.     ! Use nc4 with chunking and deflation
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  WRITE(cldum,'("1/",i2)') NINT(1./rlev_resol)

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdf2regular -f IN-file -v VAR-name [-o OUT-file] [-360] ...'
     PRINT *,'        ... [-r TGT-resolution] [-nc4] '
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Remap (by binnig) high resolution fields on a coarser regular grid,'
     PRINT *,'       keeping the same vertical grid as in the input file.'
     PRINT *,'       This program is not suitable for vector fields, as far as transport'
     PRINT *,'       conservation is concerned.'
     PRINT *,'       Default output grid resolution is ',TRIM(cldum),' deg. It can be changed using'
     PRINT *,'       the -r option.'
     PRINT *,'       Note that output file is not masked.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f IN-file  : netcdf input file ' 
     PRINT *,'       -v VAR-name : input variable name to be remapped.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-o OUT-file]: netcdf output file, instead of ',TRIM(cf_out)
     PRINT *,'       [-nc4 ]: Use netcdf4 output with chunking and deflation level 1.'
     PRINT *,'                This option is effective only if cdftools are compiled with'
     PRINT *,'                a netcdf library supporting chunking and deflation.'
     PRINT *,'       [-360 ]: Output file longitudes span [0 -> 360 deg.], instead of'
     PRINT *,'                the default [-180 -> 180 deg.].'
     PRINT *,'       [-r TGT-resolution ]:  Target resolution (in degrees).'
     PRINT *,'     '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       ',TRIM(cn_fhgr)
     PRINT *,'       ',TRIM(cn_fmsk)
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : name given as second argument'
     PRINT *,'         variables : 3d_var_name'
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'       cdf2levitus2d (a particular case)'
     PRINT *,'      '
     STOP 
  ENDIF

  ijarg = 1
  DO WHILE ( ijarg <= narg ) 
     CALL getarg( ijarg, cldum)  ; ijarg = ijarg + 1
     SELECT CASE (cldum)
     CASE ( '-f'  ) ; CALL getarg(ijarg, cf_in)  ; ijarg = ijarg + 1
     CASE ( '-v'  ) ; CALL getarg(ijarg, cv_nam) ; ijarg = ijarg + 1
     CASE ( '-o'  ) ; CALL getarg(ijarg, cf_out) ; ijarg = ijarg + 1
     CASE ( '-r'  ) ; CALL getarg( ijarg, cldum) ; ijarg = ijarg + 1 ; READ(cldum,*) rlev_resol
     CASE ( '-360') ; l360 = .TRUE. ; rlon0 = 0.
     CASE ( '-nc4') ; lnc4 = .TRUE.
     CASE DEFAULT   ; PRINT *,' ERROR : ', TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  ENDDO

  lchk = chkfile (cn_fhgr)
  lchk = chkfile (cn_fmsk) .OR. lchk
  lchk = chkfile (cf_in  ) .OR. lchk
  IF ( lchk ) STOP 99 ! missing files

  npiglo = getdim(cf_in,cn_x)
  npjglo = getdim(cf_in,cn_y)
  npk    = getdim(cf_in,cn_z)
  npt    = getdim(cf_in,cn_t)

  npilev = NINT(360./rlev_resol)
  npjlev = NINT(180./rlev_resol)

  nvars = getnvar(cf_in)
  ! Allocate the memory
  ALLOCATE (cv_names(nvars) )
  ALLOCATE (stypvar(nvars)  , stypvarout(1))
  ALLOCATE (id_var(nvars), ipk(nvars), id_varout(1), ipkout(1))
  ALLOCATE (zbt(npiglo,npjglo) , z_in(npiglo,npjglo) )

  ALLOCATE ( e1t(npiglo,npjglo), e2t(npiglo,npjglo) )
  ALLOCATE ( glamt(npiglo,npjglo), gphit(npiglo,npjglo)  )

  ALLOCATE ( d_out(npilev,npjlev) , d_n(npilev,npjlev) )
  ALLOCATE ( tmask(npiglo,npjglo) , tmasklev(npilev,npjlev))
  ALLOCATE ( rlonlev(npilev,npjlev), rlatlev(npilev,npjlev) )
  ALLOCATE ( gdept(1), dtim(npt) )

  ! Read the metrics from the mesh_hgr file
  e1t = getvar(cn_fhgr, cn_ve1t, 1, npiglo, npjglo)
  e2t = getvar(cn_fhgr, cn_ve2t, 1, npiglo, npjglo)

  ! and the coordinates   from the mesh_hgr file
  glamt = getvar(cn_fhgr, cn_glamt, 1, npiglo, npjglo)
  gphit = getvar(cn_fhgr, cn_gphit, 1, npiglo, npjglo)
  gphitmin = MINVAL(gphit(:,1))
  IF ( l360 ) THEN
     WHERE ( glamt < 0. ) glamt = glamt + 360. 
  ENDIF

  ! get the longitude,latitude,mask from the input Levitus mask file
  rlonlev(:,1) = (/ (rlon0+rlev_resol*ji , ji=1,npilev ) /)
  DO jj=2, npjlev
     rlonlev(:,jj) = rlonlev(:,1)
  ENDDO
  rlatlev(1,:) = (/ (rlev_resol*jj - 90.  , jj=1,npjlev ) /)
  DO ji=2, npilev
     rlatlev(ji,:) = rlatlev(1,:)
  ENDDO
  tmasklev(:,:) = 1.

  CALL CreateOutput

  zbt(:,:) = e1t(:,:) * e2t(:,:)    ! for surface weighting

  DO jt = 1, npt
     DO jk = 1, npk
        PRINT *,'jt = ', jt,' jk = ', jk
        ! get the tmask from the byte_mask file
        tmask(:,:) = getvar(cn_fmsk, cn_tmask, jk, npiglo, npjglo)
        z_in (:,:) = getvar(cf_in,   cv_nam,   jk, npiglo, npjglo, ktime=jt)
        ! Compute spatial mean by bin
        !-----------------------------
        ! Perform bining of the input file on the Levitus grid.
        ! Input area weighted values are summed up into a Levitus 1x1 bin
        d_out(:,:) = 0.d0
        d_n  (:,:) = 0.d0
        DO jj=1,npjglo
           DO ji=1,npiglo
              iilev = MIN( npilev,  INT( (glamt(ji,jj)-rlon0)/rlev_resol ) + 1)
              ijlev = MIN (npjlev , INT( (gphit(ji,jj)+ 90.)/rlev_resol ) + 1)
              IF ( z_in(ji,jj) /=  stypvarout(1)%rmissing_value ) THEN
                 d_out(iilev,ijlev) = d_out(iilev,ijlev) + (z_in(ji,jj)*tmask(ji,jj))*zbt(ji,jj)*tmasklev(iilev,ijlev)*1.d0
                 d_n  (iilev,ijlev) = d_n (iilev,ijlev)  +              tmask(ji,jj) *zbt(ji,jj)*tmasklev(iilev,ijlev)*1.d0
              ENDIF
           ENDDO
        ENDDO

        WHERE ( d_n > 0. )
           d_out = d_out / d_n 
        ELSEWHERE
           d_out = stypvarout(1)%rmissing_value 
        END WHERE

        ! Check if there are points with missing values on Levitus grid
        icount= COUNT( d_out == stypvarout(1)%rmissing_value .AND. tmasklev == 1. )
        IF ( icount /=  0. ) THEN
           ALLOCATE ( z_fill(npilev,npjlev) )
           z_fill(:,:) = 0.
           !
           imethod = 0  ! hard coded here JMM : check if all method works and pass method as argument
           SELECT CASE (imethod)
           CASE ( 0 ) 
              PRINT *, 'no filling, even if required for ',icount
           CASE ( 1 ) ! Method 1: fill missing data with shapiro
              ctyp='fill'
              iter_shap = 3  ! number of shapiro iteration
              CALL shapiro_fill_smooth ( REAL(d_out), npilev, npjlev, iter_shap, ctyp,         &
                   stypvarout(1)%rmissing_value, INT(tmasklev), z_fill )

              DO jjlev = 1 , npjlev
                 DO jilev = 1 , npilev
                    IF ( z_fill(jilev,jjlev) .NE. stypvarout(1)%rmissing_value  &
                         &                  .AND. tmasklev(jilev,jjlev) == 1    &
                         &                  .AND. d_out(jilev,jjlev) == stypvarout(1)%rmissing_value ) &
                         &  d_out(jilev,jjlev) = z_fill(jilev,jjlev)
                 ENDDO
              ENDDO

           CASE ( 2 ) ! Method 2: compute with influence bubble 
              ! For each point of Levitus grid, a data screening is performed 
              ! in a influence bubble of radius zradius, centered on Levitus point
              ! and the weighted average of the data in the bubble is computed
              z_fill(:,:) = 0.
              DO jjlev = 1 , npjlev-1
                 DO jilev = 1 , npilev
                    ierr = 0
                    IF (  tmasklev(jilev,jjlev) == 1 .AND. d_out(jilev,jjlev) == stypvarout(1)%rmissing_value ) THEN 
                       ! for the South pole, no treatment performed if data too far from southern most orca points
                       CALL btoe(rlonlev(jilev,jjlev),rlatlev(jilev,jjlev),rlon1,rlat1,-1.2*zradius,-1.2*zradius)
                       IF ( rlat1 > gphitmin ) THEN
                          ! Search the closest point of ORCA grid for this Levitus point
                          CALL cdf_findij (rlonlev(jilev,jjlev), rlonlev(jilev,jjlev), rlatlev(jilev,jjlev), rlatlev(jilev,jjlev), &
                               &         iimin, iimax, ijmin, ijmax,cd_coord=cn_fhgr,cd_point='T',  cd_verbose='N')

                          ! Next valid grid point going northward on ORCA grid
                          ltest = .TRUE. ; ij = ijmin ; ii =  iimin
                          DO WHILE ( ij <= npjglo .AND. ltest .AND. &
                               & etobd(rlonlev(jilev,jjlev),rlatlev(jilev,jjlev),glamt(ii,ij), gphit(ii,ij)) <= zradius )
                             IF ( tmask(ii,ij) == 1 ) THEN
                                ltest = .FALSE.
                                z_fill(jilev,jjlev) = z_fill(jilev,jjlev) + &
                                     &(z_in(ii,ij)*tmask(ii,ij))*zbt(ii,ij)*tmasklev(jilev,jjlev)
                                d_n  (jilev,jjlev) = d_n (jilev,jjlev)  + &
                                     & tmask(ii,ij) *zbt(ii,ij)*tmasklev(jilev,jjlev)
                                ierr = ierr + 1
                             ENDIF
                             ij = ij + 1
                          END DO
                          ! Next valid grid point going southward on ORCA grid
                          ltest = .TRUE. ; ij = ijmin-1 ; ii =  iimin
                          DO WHILE ( ij >= 1 .AND. ltest .AND. &
                               & etobd(rlonlev(jilev,jjlev),rlatlev(jilev,jjlev),glamt(ii,ij), gphit(ii,ij)) <= zradius ) 
                             IF ( tmask(ii,ij) == 1 ) THEN
                                ltest = .FALSE.
                                z_fill(jilev,jjlev) = z_fill(jilev,jjlev) + &
                                     & (z_in(ii,ij)*tmask(ii,ij))*zbt(ii,ij)*tmasklev(jilev,jjlev)
                                d_n  (jilev,jjlev) = d_n (jilev,jjlev)  +  &
                                     & tmask(ii,ij) *zbt(ii,ij)*tmasklev(jilev,jjlev)
                                ierr = ierr + 1
                             ENDIF
                             ij = ij - 1
                          END DO
                          ! Next valid grid point going westward on ORCA grid
                          ltest = .TRUE. ; ij = ijmin ; ii =  iimin+1 ; IF ( ii > npiglo ) ii = ii - npiglo
                          DO WHILE ( ltest .AND. &
                               & etobd(rlonlev(jilev,jjlev),rlatlev(jilev,jjlev),glamt(ii,ij), gphit(ii,ij)) <= zradius ) 
                             IF ( tmask(ii,ij) == 1 ) THEN
                                ltest = .FALSE.
                                z_fill(jilev,jjlev) = z_fill(jilev,jjlev) + &
                                     & (z_in(ii,ij)*tmask(ii,ij))*zbt(ii,ij)*tmasklev(jilev,jjlev) 
                                d_n  (jilev,jjlev) = d_n (jilev,jjlev)  + &
                                     & tmask(ii,ij) *zbt(ii,ij)*tmasklev(jilev,jjlev)
                                ierr = ierr + 1
                             ENDIF
                             ii = ii + 1
                             IF ( ii > npiglo ) ii = ii - npiglo
                          END DO
                          ! Next valid grid point going eastward on ORCA grid
                          ltest = .TRUE. ; ij = ijmin ; ii =  iimin-1 ; IF ( ii < 1) ii = ii + npiglo
                          DO WHILE ( ltest .AND. &
                               & etobd(rlonlev(jilev,jjlev),rlatlev(jilev,jjlev),glamt(ii,ij), gphit(ii,ij)) <= zradius ) 
                             IF ( tmask(ii,ij) == 1 ) THEN
                                ltest = .FALSE.
                                z_fill(jilev,jjlev) = z_fill(jilev,jjlev) + &
                                     & (z_in(ii,ij)*tmask(ii,ij))*zbt(ii,ij)*tmasklev(jilev,jjlev) 
                                d_n  (jilev,jjlev) = d_n (jilev,jjlev)  + &
                                     & tmask(ii,ij) *zbt(ii,ij)*tmasklev(jilev,jjlev)
                                ierr = ierr + 1
                             ENDIF
                             ii = ii - 1
                             IF ( ii < 1 ) ii = ii + npiglo
                          END DO

                          ! computing d_out value
                          IF ( z_fill(jilev,jjlev) .NE. stypvarout(1)%rmissing_value &
                               & .AND. d_n(jilev,jjlev) > 0 &
                               & .AND. d_out(jilev,jjlev) == stypvarout(1)%rmissing_value ) &
                               &  d_out(jilev,jjlev) = z_fill(jilev,jjlev) / d_n(jilev,jjlev)

                       ENDIF ! rlat1 > gphitmin

                    ENDIF ! tmasklev(jilev,jjlev) == 1 .AND. d_out(jilev,jjlev) == stypvarout(1)%rmissing_value

                 ENDDO
              ENDDO

              ! Case of the North Pole
              ijlev = npjlev
              DO jilev = 1 , npilev
                 ierr = 0

                 IF (  tmasklev(jilev,ijlev) == 1 .AND. d_out(jilev,ijlev) == stypvarout(1)%rmissing_value ) THEN

                    ! Search the closest point of ORCA grid for this Levitus point
                    CALL cdf_findij (rlonlev(jilev,ijlev), rlonlev(jilev,ijlev), rlatlev(jilev,ijlev), rlatlev(jilev,ijlev), &
                         & iimin, iimax, ijmin, ijmax,cd_coord=cn_fhgr,cd_point='T', cd_verbose='N')

                    ! Next valid grid point going southward on ORCA grid
                    ltest = .TRUE. ; ij = ijmin ; ii =  iimin
                    DO WHILE ( ij >= 1 .AND. ltest .AND. &
                         & etobd(rlonlev(jilev,ijlev),rlatlev(jilev,ijlev),glamt(ii,ij), gphit(ii,ij)) <= zradius ) 
                       IF ( tmask(ii,ij) == 1 ) THEN
                          ltest = .FALSE.
                          z_fill(jilev,ijlev) = z_fill(jilev,ijlev) + &
                               & (z_in(ii,ij)*tmask(ii,ij))*zbt(ii,ij)*tmasklev(jilev,ijlev) 
                          d_n  (jilev,ijlev) = d_n (jilev,ijlev)  +  &
                               & tmask(ii,ij) *zbt(ii,ij)*tmasklev(jilev,ijlev)
                          ierr = ierr + 1
                       ENDIF
                       ij = ij - 1
                    END DO
                    ! Next valid grid point going westward on ORCA grid
                    ltest = .TRUE. ; ij = ijmin ; ii =  iimin+1 ; IF ( ii > npiglo) ii = ii - npiglo
                    DO WHILE ( ltest .AND. &
                         & etobd(rlonlev(jilev,ijlev),rlatlev(jilev,ijlev),glamt(ii,ij), gphit(ii,ij)) <= zradius ) 
                       IF ( tmask(ii,ij) == 1 ) THEN
                          ltest = .FALSE.
                          z_fill(jilev,ijlev) = z_fill(jilev,ijlev) + &
                               & (z_in(ii,ij)*tmask(ii,ij))*zbt(ii,ij)*tmasklev(jilev,ijlev) 
                          d_n  (jilev,ijlev) = d_n (jilev,ijlev)  + &
                               & tmask(ii,ij) *zbt(ii,ij)*tmasklev(jilev,ijlev)
                          ierr = ierr + 1
                       ENDIF
                       ii = ii + 1
                       IF ( ii > npiglo ) ii = ii - npiglo
                    END DO
                    ! Next valid grid point going eastward on ORCA grid
                    ltest = .TRUE. ; ij = ijmin ; ii =  iimin-1 ; IF ( ii < 1 ) ii = ii + npiglo
                    DO WHILE ( ltest .AND. &
                         & etobd(rlonlev(jilev,ijlev),rlatlev(jilev,ijlev),glamt(ii,ij), gphit(ii,ij)) <= zradius ) 
                       IF ( tmask(ii,ij) == 1 ) THEN
                          ltest = .FALSE.
                          z_fill(jilev,ijlev) = z_fill(jilev,ijlev) + &
                               & (z_in(ii,ij)*tmask(ji,ij))*zbt(ii,ij)*tmasklev(jilev,ijlev) 
                          d_n  (jilev,ijlev) = d_n (jilev,ijlev)  + &
                               & tmask(ii,ij) *zbt(ii,ij)*tmasklev(jilev,ijlev)
                          ierr = ierr + 1
                       ENDIF
                       ii = ii - 1 ;  IF ( ii < 1 ) ii = ii + npiglo
                    END DO
                    ! Going "northward" and crossing the pole
                    ltest = .TRUE.
                    rpos = 1.
                    ij = ijmin + 1 * rpos ; ii =  iimin
                    IF ( ij > npjglo ) THEN
                       ii = ii + npiglo/2 ; IF ( ii > npiglo ) ii = ii - npiglo
                       rpos = -1.
                       ij = ij + 1 * rpos
                    ENDIF
                    DO WHILE ( ij >= 1 .AND. ltest .AND. &
                         & etobd(rlonlev(jilev,ijlev),rlatlev(jilev,ijlev),glamt(ii,ij), gphit(ii,ij)) <= zradius )
                       IF ( tmask(ii,ij) == 1 ) THEN
                          ltest = .FALSE.
                          z_fill(jilev,ijlev) = z_fill(jilev,ijlev) + &
                               &(z_in(ii,ij)*tmask(ii,ij))*zbt(ii,ij)*tmasklev(jilev,ijlev) 
                          d_n  (jilev,ijlev) = d_n (jilev,ijlev)  + &
                               & tmask(ii,ij) *zbt(ii,ij)*tmasklev(jilev,ijlev)
                          ierr = ierr + 1
                       ENDIF
                       ij = ij + 1 * rpos
                       IF ( ij > npjglo ) THEN
                          ii = ii + npiglo/2 ; IF ( ii > npiglo ) ii = ii - npiglo
                          rpos = -1.
                          ij = ij + 1 * rpos
                       ENDIF
                    END DO

                    ! computing d_out value
                    IF ( z_fill(jilev,ijlev) .NE. stypvarout(1)%rmissing_value  &
                         & .AND. d_n(jilev,ijlev) > 0  ) &
                         &  d_out(jilev,ijlev) = z_fill(jilev,ijlev) / d_n(jilev,ijlev)

                 ENDIF

              ENDDO
           CASE DEFAULT 
              PRINT *, ' METHOD ', imethod ,'is not recognized in this program'
              STOP 99

           END SELECT  ! imethod
           IF ( ALLOCATED(z_fill) ) DEALLOCATE( z_fill )

        ENDIF !  filling points
        ! ----------------------------------------------------------------------------------------
        ! write 
        ierr = putvar(ncout, id_varout(1), REAL(d_out(:,:)), jk, npilev, npjlev, ktime=jt)
     ENDDO ! jk

  END DO ! loop on time

  ierr = closeout(ncout)

CONTAINS

  REAL(KIND=4) FUNCTION etobd(plonr, platr, plon, plat)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION etobd  ***
    !!
    !! ** Purpose :
    !!           Compute the beta plane distance between two lon/lat values 
    !!
    !! ** Method  : 
    !!           Return the distance (km) between 2 points given in the
    !!           arguments with their latitudes and longitudes
    !!----------------------------------------------------------------------
    REAL(KIND=4), INTENT(in) :: plonr, platr, plon, plat  ! lon/lat of the 2 input points

    REAL(KIND=8), PARAMETER  :: dp_p1 = 1.745329D-02  ! radians per degree
    REAL(KIND=8), PARAMETER  :: dp_p2 = 111.1940D0    ! dp_p1 * earth radius in km (6370.949)

    REAL(KIND=8)             :: dl_rx, dl_ry          ! working double prec variables
    REAL(KIND=8)             :: dl_r0
    REAL(KIND=8)             :: dl_r1, dl_r2
    !!----------------------------------------------------------------------
    dl_r0 = DBLE(plonr) ; IF ( dl_r0 < 0.d0 ) dl_r0 = dl_r0 + 360.d0
    dl_r1 = DBLE(plon ) ; IF ( dl_r1 < 0.d0 ) dl_r1 = dl_r1 + 360.d0

    dl_rx = dl_r1 - dl_r0
    IF ( dl_rx > 180. ) THEN
       dl_rx = dl_rx - 360.d0
    ELSE IF ( dl_rx < -180.d0 ) THEN
       dl_rx = dl_rx + 360.d0
    ELSE IF ( ABS(dl_rx) == 180.d0 ) THEN
       dl_rx = 180.d0
    ENDIF
    dl_r2 = DBLE(plat) - DBLE(platr)

    dl_rx = dp_p2 * dl_rx * COS(dp_p1*platr)
    dl_ry = dp_p2 * dl_r2

    etobd = REAL(SQRT(dl_rx*dl_rx + dl_ry*dl_ry))

  END FUNCTION etobd


  SUBROUTINE btoe(plonr, platr, plon, plat, plx, ply)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE btoe  ***
    !!
    !! ** Purpose :  Return position (lon/lat) of a point located at
    !!               a given distance from the original point
    !!
    !! ** Method  :  Trigo ... 
    !!
    !!----------------------------------------------------------------------
    REAL(KIND=4), INTENT(in)  :: plonr, platr, plx, ply
    REAL(KIND=4), INTENT(out) :: plon, plat

    REAL(KIND=8), PARAMETER   :: dp_p1= 1.745329D-02  ! radians per degree
    REAL(KIND=8), PARAMETER   :: dp_p2= 111.1940D0    ! dp_p1 * earth radius in km (6370.949)
    REAL(KIND=8)              :: dl_rx, dl_ry
    REAL(KIND=8)              :: dl_r0
    REAL(KIND=8)              :: dl_r1, dl_r2, dl_r3
    !!----------------------------------------------------------------------
    dl_r0 = DBLE(plonr) ; IF ( dl_r0 < 0.d0 ) dl_r0 = dl_r0 + 360.d0
    dl_rx = DBLE(plx  )
    dl_r2 = DBLE(platr)
    dl_r1 = dl_r0 + dl_rx / ( dp_p2 * COS(dp_p1*dl_r2) )
    dl_r3 = dl_r2 + ply   / dp_p2
    plon  = REAL(dl_r1)
    IF ( plon <    0. ) plon = plon + 360.
    IF ( plon >= 360. ) plon = plon - 360.
    plat = REAL(dl_r3)
    IF ( plat >  90. ) plat =  90.
    IF ( plat < -90. ) plat = -90.

  END SUBROUTINE btoe

  SUBROUTINE CreateOutput
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutput  ***
    !!
    !! ** Purpose :  Create netcdf output file(s) 
    !!
    !! ** Method  :  Use stypvar global description of variables
    !!
    !!----------------------------------------------------------------------
  ! get list of variable names and collect attributes in stypvar (optional)
  cv_names(:) = getvarname(cf_in, nvars, stypvar)

  id_var(:)   = (/(jvar, jvar=1,nvars)/)
  ! ipk gives the number of level or 0 if not a T[Z]YX  variable
  ipk(:)     = getipk(cf_in, nvars)
  WHERE( ipk == 0 ) cv_names='none'
  stypvar(:)%cname = cv_names

  ! select variables to output:
  ii=1
  DO jk=1,nvars
     IF ( TRIM(cv_names(jk)) == TRIM(cv_nam) ) THEN
        ipkout(ii) = ipk(jk)
        stypvarout(ii) = stypvar(jk)
        stypvar(ii)%ichunk  = (/npilev,MAX(1,npjlev/30),1,1 /)
     ENDIF
  ENDDO

  ! create output fileset
  ncout = create      (cf_out,  cf_in, npilev, npjlev, npk ,cdlonvar='lon', cdlatvar='lat', ld_nc4=lnc4 )
  ierr  = createvar   (ncout ,  stypvarout, 1,      ipkout,    id_varout                  , ld_nc4=lnc4 )
  ierr  = putheadervar(ncout ,  cf_in, npilev, npjlev, npk , pnavlon=rlonlev, pnavlat=rlatlev )

  dtim = getvar1d(cf_in, cn_vtimec, npt     )
  ierr = putvar1d(ncout, dtim,      npt, 'T')

  ierr = gettimeatt(cf_in, cn_vtimec, ctcalendar, cttitle, ctlong_name, ctaxis, ctunits, cttime_origin )
  ierr = puttimeatt(ncout, cn_vtimec, ctcalendar, cttitle, ctlong_name, ctaxis, ctunits, cttime_origin )
  ierr = putvar1d( ncout, rlonlev(:,1), npilev, 'X')
  ierr = putvar1d( ncout, rlatlev(1,:), npjlev, 'Y')

  END SUBROUTINE CreateOutput

END PROGRAM cdf2regular
