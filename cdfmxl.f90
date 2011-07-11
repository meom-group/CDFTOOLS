PROGRAM cdfmxl
  !!======================================================================
  !!                     ***  PROGRAM  cdfmxl  ***
  !!=====================================================================
  !!  ** Purpose : Compute mixed layer depth
  !!
  !!  ** Method  : - compute surface properties
  !!               - initialize depths and model levels number
  !!               - from bottom to top compute rho and
  !!                 check if rho > rho_surf +rho_c
  !!                 where rho_c is a density criteria given as argument
  !!
  !! History : 2.1  : 10/2005  : J.M. Molines : Original code
  !!           3.0  : 01/2011  : J.M. Molines : Doctor norm + Lic.
  !!----------------------------------------------------------------------
  USE cdfio
  USE eos
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2011, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                              :: ji, jj, jk, jt ! dummy loop index
  INTEGER(KIND=4)                              :: ik1, ik2, ikt  ! k vertical index of mixed layers 
  INTEGER(KIND=4)                              :: narg, iargc    ! browse line
  INTEGER(KIND=4)                              :: npiglo, npjglo ! domain size
  INTEGER(KIND=4)                              :: npk, npt       ! domain size
  INTEGER(KIND=4)                              :: ncout, ierr    ! ncid of output file, error status
  INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: mbathy         ! number of w levels in water <= npk
  INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: nmln1          ! last level where rho > rho + rho_c1
  INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: nmln2          ! last level where rho > rho + rho_c1
  INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: nmlnt          ! last level where rho > rho + rho_c1
  INTEGER(KIND=4), DIMENSION(3)                :: ipk, id_varout ! levels and varid's of output vars

  REAL(KIND=4)                                 :: rho_c1=0.01    ! 1rst density criterium
  REAL(KIND=4)                                 :: rho_c2=0.03    ! 2nd density criterium
  REAL(KIND=4)                                 :: temp_c=-0.2    ! temperature criterium
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: rtem           ! temperature
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: rsal           ! salinity
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: rho            ! density
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: rho_surf       ! surface density
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: tem_surf       ! surface temperature
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: tmask_surf     ! surface tmask
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: tmask          ! level tmask
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: hmlp1          ! mxl depth based on density criterium 1
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: hmlp2          ! mxl depth based on density criterium 2
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: hmlt           ! mxl depth based on temperature criterium
  REAL(KIND=4), DIMENSION(:),      ALLOCATABLE :: gdepw          ! depth of w levels
  REAL(KIND=4), DIMENSION(:),      ALLOCATABLE :: tim            ! time counter
  REAL(KIND=4), DIMENSION(1)                   :: rdep           ! dummy depth for output

  CHARACTER(LEN=256)                           :: cf_tfil        ! input T file
  CHARACTER(LEN=256)                           :: cf_out='mxl.nc'! output file name

  TYPE(variable), DIMENSION(3)                 :: stypvar        ! structure for attributes

  LOGICAL                                      :: lexist         ! flag for existence of bathy_level file
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfmxl T-file'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute 3 estimates of the mixed layer depth from temperature'
     PRINT *,'       and salinity given in the input file, based on 3 different criteria:'
     PRINT *,'        1- Density criterium (0.01 kg/m3 difference between surface and MLD)' 
     PRINT *,'        2- Density criterium (0.03 kg/m3 difference between surface and MLD)' 
     PRINT *,'        3- Temperature criterium (0.2 C absolute difference between surface '
     PRINT *,'           and MLD)' 
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       T-file : input netcd file (gridT)' 
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        ',TRIM(cn_fzgr)
     PRINT *,'         In case of FULL STEP configuration, ',TRIM(cn_fbathylev),' is also required.'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : somxl010 = mld on density criterium 0.01'
     PRINT *,'                     somxl030 = mld on density criterium 0.03'
     PRINT *,'                     mld on temperature criterium -0.2'
     STOP
  ENDIF

  CALL getarg (1, cf_tfil)

  IF ( chkfile(cf_tfil) .OR. chkfile(cn_fzgr) ) STOP ! missing file

  ! read dimensions 
  npiglo = getdim (cf_tfil,cn_x)
  npjglo = getdim (cf_tfil,cn_y)
  npk    = getdim (cf_tfil,cn_z)
  npt    = getdim (cf_tfil,cn_t)

  rdep(1) = 0.

  ipk(:)                    = 1
  stypvar(1)%cname          = 'somxl010'
  stypvar(2)%cname          = 'somxl030'
  stypvar(3)%cname          = 'somxlt02'
  stypvar%cunits            = 'm'
  stypvar%rmissing_value    = 0.
  stypvar%valid_min         = 0.
  stypvar%valid_max         = 7000.
  stypvar(1)%clong_name     = 'Mixed_Layer_Depth_on_0.01_rho_crit'
  stypvar(2)%clong_name     = 'Mixed_Layer_Depth_on_0.03_rho_crit'
  stypvar(3)%clong_name     = 'Mixed_Layer_Depth_on_-0.2_temp_crit'
  stypvar(1)%cshort_name    = 'somxl010'
  stypvar(2)%cshort_name    = 'somxl030'
  stypvar(3)%cshort_name    = 'somxlt02'
  stypvar%conline_operation = 'N/A'
  stypvar%caxis             = 'TYX'

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  ALLOCATE (rtem(npiglo,npjglo), rsal(npiglo,npjglo), rho(npiglo,npjglo)    )
  ALLOCATE (rho_surf(npiglo,npjglo), tem_surf(npiglo,npjglo)                )
  ALLOCATE (tmask(npiglo,npjglo), tmask_surf(npiglo,npjglo)                 )
  ALLOCATE (hmlp1(npiglo,npjglo), hmlp2(npiglo,npjglo), hmlt(npiglo,npjglo) )
  ALLOCATE (mbathy(npiglo,npjglo)                                           )
  ALLOCATE (nmln1(npiglo,npjglo), nmln2(npiglo,npjglo), nmlnt(npiglo,npjglo))
  ALLOCATE (gdepw(npk), tim(npt)                                            )

  ! read mbathy and gdepw use real rtem(:,:) as template (getvar is used for real only)
  IF ( chkfile( cn_fbathylev)  ) THEN
     PRINT *, 'Read mbathy in ', TRIM(cn_fzgr),' ...'
     rtem(:,:) = getvar(cn_fzgr,      'mbathy',    1, npiglo, npjglo)
  ELSE
     rtem(:,:) = getvar(cn_fbathylev, cn_bathylev, 1, npiglo, npjglo)
  ENDIF

  mbathy(:,:) = rtem(:,:)
  gdepw(:)    = getvare3(cn_fzgr, cn_gdepw, npk)

  ncout = create      (cf_out, cf_tfil, npiglo, npjglo, 1           )
  ierr  = createvar   (ncout,  stypvar, 3,      ipk,    id_varout   )
  ierr  = putheadervar(ncout,  cf_tfil, npiglo, npjglo, 1, pdep=rdep)

  tim  = getvar1d(cf_tfil, cn_vtimec, npt     )
  ierr = putvar1d(ncout,   tim,       npt, 'T')

  DO jt=1,npt
     ! read surface T and S and deduce land-mask from salinity
     rtem( :,:) = getvar(cf_tfil, cn_votemper, 1, npiglo, npjglo, ktime=jt )
     rsal (:,:) = getvar(cf_tfil, cn_vosaline, 1, npiglo, npjglo, ktime=jt )
     IF (jt == 1 ) THEN
        tmask(:,:) = 1.;  WHERE ( rsal == 0. ) tmask = 0.
        tmask_surf(:,:) = tmask(:,:)
     ENDIF

     ! compute rho_surf
     rho_surf(:,:) = sigma0 (rtem, rsal, npiglo, npjglo )* tmask(:,:)
     tem_surf(:,:) = rtem(:,:)

     ! Initialization to the number of w ocean point mbathy
     nmln1(:,:) = mbathy(:,:)
     nmln2(:,:) = mbathy(:,:)
     nmlnt(:,:) = mbathy(:,:)

     ! compute mixed layer depth
     ! Last w-level at which rhop>=rho surf+rho_c (starting from jpk-1)
     ! (rhop defined at t-point, thus jk-1 for w-level just above)
     DO jk = npk-1, 2, -1
        rtem (:,:) = getvar(cf_tfil, cn_votemper, jk ,npiglo, npjglo, ktime=jt)
        rsal (:,:) = getvar(cf_tfil, cn_vosaline, jk ,npiglo, npjglo, ktime=jt)
        tmask(:,:) = 1.  ; WHERE ( rsal == 0. ) tmask = 0.
        rho  (:,:) = sigma0 (rtem, rsal, npiglo, npjglo )* tmask(:,:)

        DO jj = 1, npjglo
           DO ji = 1, npiglo
              IF( rho(ji,jj)  > rho_surf(ji,jj) + rho_c1 )   nmln1(ji,jj) = jk
              IF( rho(ji,jj)  > rho_surf(ji,jj) + rho_c2 )   nmln2(ji,jj) = jk
              IF( ABS(rtem(ji,jj) - tem_surf(ji,jj)) > ABS( temp_c)  )   nmlnt(ji,jj) = jk
           END DO
        END DO
     END DO

     ! Mixed layer depth
     DO jj = 1, npjglo
        DO ji = 1, npiglo
           ik1 = nmln1(ji,jj) ; ik2 = nmln2(ji,jj) ; ikt = nmlnt(ji,jj)
           hmlp1 (ji,jj) = gdepw(ik1) * tmask_surf(ji,jj)
           hmlp2 (ji,jj) = gdepw(ik2) * tmask_surf(ji,jj)
           hmlt (ji,jj)  = gdepw(ikt) * tmask_surf(ji,jj)
        END DO
     END DO

     ierr = putvar(ncout, id_varout(1), hmlp1, 1, npiglo, npjglo, ktime=jt)
     ierr = putvar(ncout, id_varout(2), hmlp2, 1, npiglo, npjglo, ktime=jt)
     ierr = putvar(ncout, id_varout(3), hmlt , 1, npiglo, npjglo, ktime=jt)

  END DO ! time loop

  ierr = closeout(ncout)

END PROGRAM cdfmxl
