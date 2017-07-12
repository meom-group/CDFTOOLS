PROGRAM cdf2levitusgrid2d
  !!======================================================================
  !!                     ***  PROGRAM  cdf2levitusgrid2d  ***
  !!=====================================================================
  !!  ** Purpose : remaps (bin) 2D high resolution (finer than 1x1 deg)
  !!               fields on Levitus 2D 1x1 deg grid
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

  INTEGER(KIND=4)                              :: ji, jj, jk, jvar,jt! dummy loop index
  INTEGER(KIND=4)                              :: jilev,jjlev        ! dummy loop index
  INTEGER(KIND=4)                              :: numvar0            ! dummy loop index
  INTEGER(KIND=4)                              :: ii, ij             ! array index (not loop)
  INTEGER(KIND=4)                              :: iilev, ijlev       ! array index (not loop)
  INTEGER(KIND=4)                              :: npiglo, npjglo     ! size of the domain
  INTEGER(KIND=4)                              :: npk, npt           ! size of the domain
  INTEGER(KIND=4)                              :: npilev, npjlev     ! size of the Levitus domain
  INTEGER(KIND=4)                              :: narg, iargc, ijarg ! browse line
  INTEGER(KIND=4)                              :: ireq              ! ,andatory arguments counter. 
  INTEGER(KIND=4)                              :: ncout              ! ncid of output file
  INTEGER(KIND=4)                              :: ierr               ! error status
  INTEGER(KIND=4)                              :: nvars              ! number of variables in the input file
  INTEGER(KIND=4)                              :: nvarsout           ! number of variables in the output file
  INTEGER(KIND=4)                              :: iter_shap=3        ! number of Shapiro iteration
  INTEGER(KIND=4), DIMENSION(:),   ALLOCATABLE :: ipk, id_var        ! levels and varid's of input vars
  INTEGER(KIND=4), DIMENSION(:),   ALLOCATABLE :: ipkout, id_varout  ! levels and varid's of output vars

  REAL(KIND=4)                                 :: rlon1, rlon2, rlat1, rlat2, rpos
  REAL(KIND=4)                                 :: gphitmin
  REAL(KIND=4), DIMENSION(:),      ALLOCATABLE :: gdept              ! depth axis
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: e1t, e2t           ! horizontal T metrics
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: z_in               ! input field
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: z_fill             ! output 
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: glamt, gphit       ! T longitude latitude
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: rlonlev, rlatlev   ! Levitus grid longitude latitude
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: zbt
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: tmask              ! input mask
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: tmasklev           ! output mask

  REAL(KIND=8), DIMENSION(:),      ALLOCATABLE :: dtim               ! time counter
  REAL(KIND=8), DIMENSION(:,:),    ALLOCATABLE :: d_out, d_n         ! output field and weighting field

  CHARACTER(LEN=256)                           :: cf_in              ! input file name
  CHARACTER(LEN=256)                           :: cf_out             ! output file name ( output)
  CHARACTER(LEN=256)                           :: cf_levitus_mask='levitus_mask.nc'   ! Levitus mask filename
  CHARACTER(LEN=256)                           :: cv_nam             ! variable name
  CHARACTER(LEN=256)                           :: cldum              ! dummy string
  CHARACTER(LEN=256)                           :: ctcalendar         ! time attributes
  CHARACTER(LEN=256)                           :: cttitle            ! time attributes
  CHARACTER(LEN=256)                           :: ctlong_name        ! time attributes
  CHARACTER(LEN=256)                           :: ctaxis             ! time attributes
  CHARACTER(LEN=256)                           :: ctunits            ! time attributes
  CHARACTER(LEN=256)                           :: cttime_origin      ! time attributes

  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE:: cv_names          ! array of var name
  CHARACTER(LEN=6)                             :: ctyp              ! 'fill' or 'smooth' for shapiro

  TYPE(variable), DIMENSION(:),   ALLOCATABLE  :: stypvar            ! input attributes
  TYPE(variable), DIMENSION(:),   ALLOCATABLE  :: stypvarout         ! output attributes

  LOGICAL                                      :: lchk               ! missing files flag
  LOGICAL                                      :: ltest
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdf2levitusgrid2d -f IN-file -o OUT-file -v VAR-name2D'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Remaps (by binning) 2D high resolution (i.e. finer than 1x1 deg) '
     PRINT *,'       fields on Levitus 2D 1x1 deg grid. This program does not work for'
     PRINT *,'       vector fields.'
     PRINT *,'       It assumes that the southwestern-most grid cell of the target grid'
     PRINT *,'       (Levitus 1 deg) is centered at (0.5W,89.5S).' 
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f IN-file  : netcdf input file ' 
     PRINT *,'       -o OUT-file : netcdf output file ' 
     PRINT *,'       -v VAR-name2D : input variable name for interpolation '
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       ',TRIM(cn_fhgr)
     PRINT *,'       ',TRIM(cn_fmsk)
     PRINT *,'       ',TRIM(cf_levitus_mask)
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : name given as second argument'
     PRINT *,'         variables : 2d_var_name'
     STOP 
  ENDIF

  ijarg = 1 ; ireq = 0 
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum)  ; ijarg = ijarg + 1
     SELECT CASE ( cldum )
     CASE ( '-f' ) ; CALL getarg(ijarg, cf_in ) ; ijarg = ijarg + 1 ; ireq=ireq+1
     CASE ( '-o' ) ; CALL getarg(ijarg, cf_out) ; ijarg = ijarg + 1 ; ireq=ireq+1
     CASE ( '-v' ) ; CALL getarg(ijarg, cv_nam) ; ijarg = ijarg + 1 ; ireq=ireq+1
     CASE DEFAULT  ; PRINT *,', ERROR : ', TRIM(cldum), ' :  unknown option.' ; STOP 99
     END SELECT
  ENDDO

  IF ( ireq /= 3 ) THEN ; PRINT *, ' ERROR : missing arguments.'; STOP 99  ; ENDIF

  lchk = chkfile (cn_fhgr)
  lchk = chkfile (cn_fmsk)         .OR. lchk
  lchk = chkfile (cf_levitus_mask) .OR. lchk
  lchk = chkfile (cf_in)           .OR. lchk
  IF ( lchk ) STOP 99 ! missing files

  npiglo = getdim(cf_in,cn_x)
  npjglo = getdim(cf_in,cn_y)
  npk    = getdim(cf_in,cn_z)
  npt    = getdim(cf_in,cn_t)
  npilev = getdim(cf_levitus_mask,cn_x)
  npjlev = getdim(cf_levitus_mask,cn_y)

  nvars = getnvar(cf_in)
  ! Allocate the memory
  ALLOCATE (cv_names(nvars) )
  ALLOCATE (stypvar(nvars)  , stypvarout(1))
  ALLOCATE (id_var(nvars), ipk(nvars), id_varout(1), ipkout(1))
  ALLOCATE ( zbt(npiglo,npjglo) , z_in(npiglo,npjglo) )
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
  WHERE ( glamt < 0. ) glamt = glamt + 360. 

  ! get the tmask from the byte_mask file
  tmask(:,:) = getvar(cn_fmsk, cn_tmask, 1, npiglo, npjglo)

  CALL CreateOutput

  ! Get input 2D field ( numvar0 is the varid of the requested var from all the var in cf_in).
  z_in(:,:) = getvar(cf_in, cv_names(numvar0), 1, npiglo, npjglo)

  zbt(:,:) = e1t(:,:) * e2t(:,:)    ! for surface weighting

  DO jt = 1, npt
     PRINT *,'jt = ', jt
     ! Compute spatial mean by bin
     !-----------------------------
     ! Perform bining of the input file on the Levitus grid.
     ! Input area weighted values are summed up into a Levitus 1x1 bin
     d_out(:,:) = 0.d0
     d_n  (:,:) = 0.d0
     DO jj=1,npjglo
        DO ji=1,npiglo
           iilev = MIN( 360, INT( glamt(ji,jj) ) + 1)
           ijlev = MIN (180 , INT( gphit(ji,jj) + 90. ) + 1)
           IF ( z_in(ji,jj) /=  stypvarout(1)%rmissing_value ) THEN
              d_out(iilev,ijlev) = d_out(iilev,ijlev) + (z_in(ji,jj)*tmask(ji,jj))*zbt(ji,jj)*tmasklev(iilev,ijlev)*1.d0
              d_n  (iilev,ijlev) = d_n (iilev,ijlev)  +              tmask(ji,jj) *zbt(ji,jj)*tmasklev(iilev,ijlev)*1.d0
           ENDIF
        ENDDO
     ENDDO

     WHERE ( d_n > 0. ) ; d_out = d_out / d_n 
     ELSEWHERE          ; d_out = stypvarout(1)%rmissing_value 
     END WHERE

     ! Check if there are points with missing values on Levitus grid
     IF ( COUNT( d_out == stypvarout(1)%rmissing_value .AND. tmasklev == 1. ) /=  0. ) THEN
        ALLOCATE ( z_fill(npilev,npjlev) )
        z_fill(:,:) = 0.
        ctyp='fill'
        iter_shap = 3  ! number of shapiro iteration
        CALL shapiro_fill_smooth ( REAL(d_out), npilev, npjlev, iter_shap, ctyp,       &
           &                       stypvarout(1)%rmissing_value, INT(tmasklev), z_fill )

        DO jjlev = 1 , npjlev
           DO jilev = 1 , npilev
             IF (       z_fill(jilev,jjlev)   /=   stypvarout(1)%rmissing_value  &
                & .AND. tmasklev(jilev,jjlev) == 1                               &
                & .AND. d_out(jilev,jjlev)    == stypvarout(1)%rmissing_value )  THEN
                d_out(jilev,jjlev) = z_fill(jilev,jjlev)
             ENDIF
           ENDDO
        ENDDO
        IF ( ALLOCATED(z_fill) ) DEALLOCATE( z_fill )
     ENDIF !  filling points
     ! ----------------------------------------------------------------------------------------
     ! write 
     ierr = putvar(ncout, id_varout(1), REAL(d_out(:,:)), 1, npilev, npjlev, ktime=jt)

  END DO ! loop on time

  ierr = closeout(ncout)

CONTAINS

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

    ! select variable to output: (only one passed as argument)
    ii=1
    DO jvar=1,nvars
       IF ( TRIM(cv_names(jvar)) == TRIM(cv_nam) ) THEN
          ipkout(ii)     = ipk(jvar)
          stypvarout(ii) = stypvar(jvar)
          stypvarout(ii)%rmissing_value=getspval ( cf_in, TRIM(cv_nam) )
          PRINT*, 'rmissing_value = ', stypvarout(ii)%rmissing_value
          nvarsout = ii
          numvar0 = jvar
          EXIT ! found !
       ENDIF
    ENDDO

    ! JMM this Levitus mask file is not standard ? Must be provided some-how ...
    !      variable names hard coded here : OK 
    ! get the longitude,latitude,mask from the input Levitus mask file
    rlonlev(:,:)  = getvar(cf_levitus_mask, 'nav_lon',  1, npilev, npjlev)
    rlatlev(:,:)  = getvar(cf_levitus_mask, 'nav_lat' , 1, npilev, npjlev)
    tmasklev(:,:) = getvar(cf_levitus_mask, 'mask',     1, npilev, npjlev)

    ! create output fileset
    ncout = create      (cf_out, cf_levitus_mask, npilev, npjlev, 0 ,cdlonvar='lon', cdlatvar='lat'  )
    ierr  = createvar   (ncout ,  stypvarout, 1,      ipkout,    id_varout    )
    ierr  = putheadervar(ncout ,  'dummy', npilev, npjlev, 0 , pnavlon=rlonlev, pnavlat=rlatlev )

    dtim = getvar1d(cf_in, cn_vtimec, npt     )
    ierr = putvar1d(ncout, dtim,      npt, 'T')
    ierr = gettimeatt(cf_in, cn_vtimec, ctcalendar, cttitle, ctlong_name, ctaxis, ctunits, cttime_origin )
    ierr = puttimeatt(ncout, cn_vtimec, ctcalendar, cttitle, ctlong_name, ctaxis, ctunits, cttime_origin )
    ierr = putvar1d( ncout, rlonlev(:,1), npilev, 'X')
    ierr = putvar1d( ncout, rlatlev(1,:), npjlev, 'Y')

  END SUBROUTINE CreateOutput

END PROGRAM cdf2levitusgrid2d
