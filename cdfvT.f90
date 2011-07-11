PROGRAM cdfvT
  !!======================================================================
  !!                     ***  PROGRAM  cdfvT  ***
  !!=====================================================================
  !!  ** Purpose : Compute the average values for the products 
  !!               V.T, V.S, U.T and U.S, used afterward for heat and
  !!               salt transport.
  !!
  !!  ** Method  : pass the CONFIG name and a series of tags as arguments.
  !!
  !! History : 2.1  : 11/2004  : J.M. Molines : Original code
  !!           2.1  : 02/2010  : J.M. Molines : handle multiframes input files.
  !!           3.0  : 04/2011  : J.M. Molines : Doctor norm + Lic.
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

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ztemp, zsal          ! Array to read a layer of data
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zu, zv               ! Velocity component
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zworku, zworkv       ! working arrays
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zmean                ! temporary mean value for netcdf write
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: tim                  ! time counter of individual files
  REAL(KIND=4), DIMENSION(1)                :: timean               ! mean time

  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dcumulut, dcumulus   ! Arrays for cumulated values
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dcumulvt, dcumulvs   ! Arrays for cumulated values
  REAL(KIND=8)                              :: dtotal_time          ! cumulated time

  CHARACTER(LEN=256)                        :: cf_tfil              ! TS file name
  CHARACTER(LEN=256)                        :: cf_ufil              ! zonal velocity file
  CHARACTER(LEN=256)                        :: cf_vfil              ! meridional velocity file
  CHARACTER(LEN=256)                        :: cf_out='vt.nc'       ! output file
  CHARACTER(LEN=256)                        :: config               ! configuration name
  CHARACTER(LEN=256)                        :: ctag                 ! current tag to work with               

  TYPE (variable), DIMENSION(4)             :: stypvar              ! structure for attributes

  LOGICAL                                   :: lcaltmean            ! flag for mean time computation
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  !!  Read command line
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfvT CONFIG ''list_of_tags'' '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the time average values for second order products ' 
     PRINT *,'       V.T, V.S, U.T and U.S used in heat and salt transport computation.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       CONFIG is the config name of a given experiment (eg ORCA025-G70)'
     PRINT *,'            The program will look for gridT, gridU and gridV files for' 
     PRINT *,'            this config ( grid_T, grid_U and grid_V are also accepted).'
     PRINT *,'       list_of_tags : a list of time tags that will be used for time'
     PRINT *,'            averaging. e.g. y2000m01d05 y2000m01d10 ...'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        none'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'       variables : ',TRIM(cn_vozout),', ',TRIM(cn_vozous),', ',TRIM(cn_vomevt),' and ',TRIM(cn_vomevs)
     STOP
  ENDIF

  !! Initialisation from 1st file (all file are assume to have the same geometry)
  CALL getarg (1, config)
  CALL getarg (2, ctag  )

  cf_tfil = SetFileName( config, ctag, 'T')

  npiglo = getdim (cf_tfil,cn_x)
  npjglo = getdim (cf_tfil,cn_y)
  npk    = getdim (cf_tfil,cn_z)

  ipk(:)= npk  ! all variables (input and output are 3D)
  ! define output variables
  stypvar%rmissing_value    = 0.
  stypvar%valid_min         = -100.
  stypvar%valid_max         = 100.
  stypvar%conline_operation = 'N/A'
  stypvar%caxis             = 'TZYX'

  stypvar(1)%cname          = cn_vomevt       ; stypvar(1)%cunits        = 'm.DegC.s-1'
  stypvar(2)%cname          = cn_vomevs       ; stypvar(2)%cunits        = 'm.PSU.s-1'
  stypvar(3)%cname          = cn_vozout       ; stypvar(3)%cunits        = 'm.DegC.s-1'
  stypvar(4)%cname          = cn_vozous       ; stypvar(4)%cunits        = 'm.PSU.s-1'

  stypvar(1)%clong_name     = 'Meridional_VT' ; stypvar(1)%cshort_name   = cn_vomevt
  stypvar(2)%clong_name     = 'Meridional_VS' ; stypvar(2)%cshort_name   = cn_vomevs
  stypvar(3)%clong_name     = 'Zonal_UT'      ; stypvar(3)%cshort_name   = cn_vozout
  stypvar(4)%clong_name     = 'Zonal_US'      ; stypvar(4)%cshort_name   = cn_vozous

  PRINT *, 'npiglo =', npiglo
  PRINT *, 'npjglo =', npjglo
  PRINT *, 'npk    =', npk

  ALLOCATE( dcumulut(npiglo,npjglo), dcumulus(npiglo,npjglo) )
  ALLOCATE( dcumulvt(npiglo,npjglo), dcumulvs(npiglo,npjglo) )
  ALLOCATE( zu(npiglo,npjglo),    zv(npiglo,npjglo) )
  ALLOCATE( zworku(npiglo,npjglo),   zworkv(npiglo,npjglo) )
  ALLOCATE( ztemp(npiglo,npjglo),    zsal(npiglo,npjglo) )
  ALLOCATE( zmean(npiglo,npjglo))

  ! create output fileset
  ncout = create      (cf_out, cf_tfil, npiglo, npjglo, npk       )
  ierr  = createvar   (ncout , stypvar, 4,      ipk,    id_varout )
  ierr  = putheadervar(ncout,  cf_tfil, npiglo, npjglo, npk       )

  lcaltmean=.TRUE.
  DO jk = 1, npk
     PRINT *,'level ',jk
     dcumulut(:,:) = 0.d0 ;  dcumulvt(:,:) = 0.d0 ; dtotal_time = 0.d0
     dcumulus(:,:) = 0.d0 ;  dcumulvs(:,:) = 0.d0 ; ntframe = 0

     DO jt = 2, narg           ! loop on tags
        CALL getarg (jt, ctag)

        cf_tfil = SetFileName( config, ctag, 'T' )

        npt = getdim (cf_tfil, cn_t)
        IF ( lcaltmean ) THEN
           ALLOCATE ( tim(npt) )
           tim = getvar1d(cf_tfil, cn_vtimec, npt)
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
           ztemp(:,:) = getvar(cf_tfil, cn_votemper, jk, npiglo, npjglo, ktime=jtt )
           zsal(:,:)  = getvar(cf_tfil, cn_vosaline, jk, npiglo, npjglo, ktime=jtt )

           ! temperature at u point, v points
           zworku(:,:) = 0. ; zworkv(:,:) = 0.
           DO ji=1, npiglo-1
              DO jj = 1, npjglo -1
                 zworku(ji,jj) = 0.5 * ( ztemp(ji,jj) + ztemp(ji+1,jj) )  ! temper at Upoint
                 zworkv(ji,jj) = 0.5 * ( ztemp(ji,jj) + ztemp(ji,jj+1) )  ! temper at Vpoint
              END DO
           END DO

           dcumulut(:,:) = dcumulut(:,:) + zworku(:,:) * zu(:,:)*1.d0
           dcumulvt(:,:) = dcumulvt(:,:) + zworkv(:,:) * zv(:,:)*1.d0

           ! salinity at u points, v points
           zworku(:,:) = 0. ; zworkv(:,:) = 0.
           DO ji=1, npiglo-1
              DO jj = 1, npjglo -1
                 zworku(ji,jj) = 0.5 * ( zsal(ji,jj) + zsal(ji+1,jj) )  ! salinity  at Upoint
                 zworkv(ji,jj) = 0.5 * ( zsal(ji,jj) + zsal(ji,jj+1) )  ! salinity  at Vpoint
              END DO
           END DO

           dcumulus(:,:) = dcumulus(:,:) + zworku(:,:) * zu(:,:)*1.d0
           dcumulvs(:,:) = dcumulvs(:,:) + zworkv(:,:) * zv(:,:)*1.d0

        END DO  !jtt
     END DO  ! jt
     ! finish with level jk ; compute mean (assume spval is 0 )
     zmean(:,:) = dcumulvt(:,:)/ntframe
     ierr = putvar(ncout, id_varout(1), zmean, jk,npiglo, npjglo, kwght=ntframe )

     zmean(:,:) = dcumulvs(:,:)/ntframe
     ierr = putvar(ncout, id_varout(2), zmean, jk,npiglo, npjglo, kwght=ntframe )

     zmean(:,:) = dcumulut(:,:)/ntframe
     ierr = putvar(ncout, id_varout(3), zmean, jk,npiglo, npjglo, kwght=ntframe )

     zmean(:,:) = dcumulus(:,:)/ntframe
     ierr = putvar(ncout, id_varout(4), zmean, jk,npiglo, npjglo, kwght=ntframe )

     IF (lcaltmean )  THEN
        timean(1) = dtotal_time/ntframe
        ierr      = putvar1d(ncout, timean, 1, 'T')
     END IF
     lcaltmean=.FALSE. ! tmean already computed

  END DO  ! loop to next level

  ierr = closeout(ncout)

END PROGRAM cdfvT
