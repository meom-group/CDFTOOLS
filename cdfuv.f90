PROGRAM cdfuv
  !!======================================================================
  !!                     ***  PROGRAM  cdfuv  ***
  !!=====================================================================
  !!  ** Purpose : Compute the average values for the products 
  !!               u.v at T-point
  !!
  !!  ** Method  : pass the CONFIG name and a series of tags as arguments.
  !!
  !! History : 2.1  : 11/2004  : J.M. Molines : Original code
  !!           2.1  : 02/2010  : J.M. Molines : handle multiframes input files.
  !!           3.0  : 04/2011  : J.M. Molines : Doctor norm + Lic.
  !!                : 10/2012  : M. Balmaseda : Split T and S file eventually
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  USE modutils       ! SetFileName function
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2011, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                           :: ji, jj, jk, jt, jtt  ! dummy loop index
  INTEGER(KIND=4)                           :: ierr                 ! working integer
  INTEGER(KIND=4)                           :: narg, iargc          ! command line
  INTEGER(KIND=4)                           :: npiglo,npjglo        ! size of the domain
  INTEGER(KIND=4)                           :: npk, npt             ! size of the domain
  INTEGER(KIND=4)                           :: ntframe              ! Cumul of time frame
  INTEGER(KIND=4)                           :: ncout                ! ncid of output file
  INTEGER(KIND=4), DIMENSION(4)             :: ipk, id_varout       ! level and varid's of output vars

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zu, zv               ! Velocity component
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zworku, zworkv       ! working arrays
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zmean                ! temporary mean value for netcdf write
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zlon                 ! longitude of T points to check periodicity
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: tim                  ! time counter of individual files
  REAL(KIND=4), DIMENSION(1)                :: timean               ! mean time

  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dcumuluv             ! Arrays for cumulated values
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dcumulu              ! Arrays for cumulated values
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dcumulv              ! Arrays for cumulated values
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dupvp                ! Arrays for U'.V'
  REAL(KIND=8)                              :: dtotal_time          ! cumulated time

  CHARACTER(LEN=256)                        :: cf_tfil              ! Temperature file for reference only
  CHARACTER(LEN=256)                        :: cf_ufil              ! zonal velocity file
  CHARACTER(LEN=256)                        :: cf_vfil              ! meridional velocity file
  CHARACTER(LEN=256)                        :: cf_out='uv.nc'       ! output file
  CHARACTER(LEN=256)                        :: config               ! configuration name
  CHARACTER(LEN=256)                        :: ctag                 ! current tag to work with               
  CHARACTER(LEN=256)                        :: cl_name              ! temporary variable name

  TYPE (variable), DIMENSION(4)             :: stypvar              ! structure for attributes

  LOGICAL                                   :: lcaltmean            ! flag for mean time computation
  LOGICAL                                   :: lperio=.false.       ! flag for E-W periodicity
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  !!  Read command line
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfuv CONFIG-CASE ''list_of_tags'' '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the time average values for U.V  product, at T point.' 
     PRINT *,'       Mean U and V values at T points, and mean U''.V'' product are '
     PRINT *,'       saved as well.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       CONFIG-CASE is the config name of a given experiment (eg ORCA025-G70)'
     PRINT *,'            The program will look for gridU and gridV files for' 
     PRINT *,'            this config (grid_U and grid_V are also accepted).'
     PRINT *,'       list_of_tags : a list of time tags that will be used for time'
     PRINT *,'            averaging. e.g. y2000m01d05 y2000m01d10 ...'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        none'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'       variables : ',TRIM(cn_vouv), '  : Mean U.V at T point'
     PRINT *,'                   ',TRIM(cn_vozocrtx)//'_t : Mean U at T point'
     PRINT *,'                   ',TRIM(cn_vomecrty)//'_t : Mean V at T point'
     PRINT *,'                   ',TRIM(CN_VOUV)//'_prime : Mean U''.V'' at T point'
     STOP
  ENDIF

  !! Initialisation from 1st file (all file are assume to have the same geometry)
  CALL getarg (1, config)
  CALL getarg (2, ctag  )

  cf_tfil = SetFileName( config, ctag, 'T')
  cf_ufil = SetFileName( config, ctag, 'U')

   

  npiglo = getdim (cf_ufil,cn_x)
  npjglo = getdim (cf_ufil,cn_y)
  npk    = getdim (cf_ufil,cn_z)



  ipk(:)= npk  ! all variables (input and output are 3D)
  ! define output variables
  stypvar%rmissing_value    = 0.
  stypvar%valid_min         = -100.
  stypvar%valid_max         = 100.
  stypvar%conline_operation = 'N/A'
  stypvar%caxis             = 'TZYX'

  stypvar(1)%cname          = cn_vouv                  ; stypvar(1)%cunits        = 'm2/s2'
  stypvar(1)%clong_name     = 'U.V product at T point' ; stypvar(1)%cshort_name   = cn_vouv

  cl_name = TRIM(cn_vozocrtx)//'_t'
  stypvar(2)%cname          = cl_name                  ; stypvar(2)%cunits        = 'm/s'
  stypvar(2)%clong_name     = 'Mean U at T point '     ; stypvar(2)%cshort_name   = cl_name

  cl_name = TRIM(cn_vomecrty)//'_t'
  stypvar(3)%cname          = cl_name                  ; stypvar(3)%cunits        = 'm/s'
  stypvar(3)%clong_name     = 'Mean V at T point '     ; stypvar(3)%cshort_name   = cl_name

  cl_name = TRIM(cn_vouv)//'_prime'
  stypvar(4)%cname          = cl_name                     ; stypvar(3)%cunits        = 'm2/s2'
  stypvar(4)%clong_name     = 'Uprime .Vprime at T point' ; stypvar(3)%cshort_name   = cl_name

  PRINT *, 'npiglo =', npiglo
  PRINT *, 'npjglo =', npjglo
  PRINT *, 'npk    =', npk

  ALLOCATE( dcumuluv(npiglo,npjglo) )
  ALLOCATE( dcumulu(npiglo,npjglo) )
  ALLOCATE( dcumulv(npiglo,npjglo) )
  ALLOCATE( dupvp(npiglo,npjglo) )
  ALLOCATE( zu(npiglo,npjglo),     zv(npiglo,npjglo) )
  ALLOCATE( zworku(npiglo,npjglo), zworkv(npiglo,npjglo) )
  ALLOCATE( zmean(npiglo,npjglo))
  ALLOCATE( zlon(npiglo,npjglo))

   ! check for E_W periodicity

  zlon(:,:) = getvar(cf_tfil, cn_vlon2d, 1, npiglo, npjglo )
  IF ( zlon(1,1) ==  zlon(npiglo-1,1) ) THEN
     lperio = .TRUE.
     PRINT *,' E-W periodicity detected '
  ENDIF

  ! create output fileset
  ncout = create      (cf_out, cf_tfil, npiglo, npjglo, npk, ld_xycoo=.TRUE. )
  ierr  = createvar   (ncout , stypvar, 4,      ipk,    id_varout            )
  ierr  = putheadervar(ncout,  cf_tfil, npiglo, npjglo, npk, ld_xycoo=.TRUE. )
  
  lcaltmean=.TRUE.
  DO jk = 1, npk
     PRINT *,'level ',jk
     dcumuluv(:,:) = 0.d0 ;  dtotal_time  = 0.d0 ; ntframe = 0
     dcumulu(:,:)  = 0.d0 ;  dcumulv(:,:) = 0.d0

     DO jt = 2, narg           ! loop on tags
        CALL getarg (jt, ctag)

        cf_ufil = SetFileName( config, ctag, 'U', ld_stop=.TRUE. )
        npt = getdim (cf_ufil, cn_t)
        IF ( lcaltmean ) THEN
           ALLOCATE ( tim(npt) )
           tim = getvar1d(cf_ufil, cn_vtimec, npt)
           dtotal_time = dtotal_time + SUM(tim(1:npt) )
           DEALLOCATE( tim )
        END IF

        ! assume U and V file have same time span ...
        cf_ufil = SetFileName( config, ctag, 'U' )
        cf_vfil = SetFileName( config, ctag, 'V' )

        DO jtt = 1, npt  ! loop on time frame in a single file
           ntframe = ntframe+1
           zu(:,:)    = getvar(cf_ufil, cn_vozocrtx, jk, npiglo, npjglo, ktime=jtt )
           zv(:,:)    = getvar(cf_vfil, cn_vomecrty, jk, npiglo, npjglo, ktime=jtt )

           ! U and V velocities at T point
           zworku(:,:) = 0. ; zworkv(:,:) = 0.
           DO ji=2, npiglo
              DO jj = 1, npjglo
                 zworku(ji,jj) = 0.5 * ( zu(ji,jj) + zu(ji-1,jj) )  ! U at T point
              END DO
           END DO

           DO ji=1, npiglo
              DO jj = 2, npjglo 
                 zworkv(ji,jj) = 0.5 * ( zv(ji,jj) + zv(ji,jj-1) )  ! V at T point
              END DO
           END DO

           dcumuluv(:,:) = dcumuluv(:,:) + zworku(:,:) * zworkv(:,:)*1.d0
           dcumulu(:,:)  = dcumulu(:,:)  + zworku(:,:)*1.d0
           dcumulv(:,:)  = dcumulv(:,:)  + zworkv(:,:)*1.d0

        END DO  !jtt
     END DO  ! jt
     ! finish with level jk ; compute mean (assume spval is 0 )
     zmean(:,:) = dcumuluv(:,:)/ntframe
     IF ( lperio ) zmean(1, :) = zmean(npiglo-1,:)
     ierr = putvar(ncout, id_varout(1), zmean, jk,npiglo, npjglo, kwght=ntframe )

     zmean(:,:) = dcumulu(:,:)/ntframe
     IF ( lperio ) zmean(1, :) = zmean(npiglo-1,:)
     ierr = putvar(ncout, id_varout(2), zmean, jk,npiglo, npjglo, kwght=ntframe )

     zmean(:,:) = dcumulv(:,:)/ntframe
     IF ( lperio ) zmean(1, :) = zmean(npiglo-1,:)
     ierr = putvar(ncout, id_varout(3), zmean, jk,npiglo, npjglo, kwght=ntframe )

     dupvp(:,:) = dcumuluv(:,:)/ntframe - dcumulu(:,:)*dcumulv(:,:)/ntframe/ntframe
     zmean(:,:) = dupvp(:,:)
     IF ( lperio ) zmean(1, :) = zmean(npiglo-1,:)
     ierr = putvar(ncout, id_varout(4), zmean, jk,npiglo, npjglo, kwght=ntframe )

     IF (lcaltmean )  THEN
        timean(1) = dtotal_time/ntframe
        ierr      = putvar1d(ncout, timean, 1, 'T')
     END IF
     lcaltmean=.FALSE. ! tmean already computed

  END DO  ! loop to next level

  ierr = closeout(ncout)

END PROGRAM cdfuv
