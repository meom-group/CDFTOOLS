PROGRAM cdf_remake_bathy 
  !!======================================================================
  !!                     ***  PROGRAM  cdf_remake_bathy  ***
  !!=====================================================================
  !!  ** Purpose : Recompute bathymetry from e3t vertical intregration
  !!
  !! History : 4.0  : 09/2019  : Q. Jamet, J.-M. Molines
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames   ! for cdf variable names
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class domain_file
  !!----------------------------------------------------------------------
  IMPLICIT NONE


  INTEGER(KIND=4), PARAMETER                   :: wp=4
  INTEGER(KIND=4), PARAMETER                   :: jpvarout = 3             ! number of output variables
  INTEGER(KIND=4)                              :: ncout                    ! ncid of output file
  INTEGER(KIND=4), DIMENSION(jpvarout)         :: id_varout                ! id of output variables
  INTEGER(KIND=4), DIMENSION(jpvarout)         :: ipk                      ! level of output variables
  INTEGER(KIND=4)                              :: ierr                     ! working integer
  INTEGER(KIND=4)                              :: narg, iargc, ijarg       !
  INTEGER(KIND=4)                              :: ji, jj, jk, jt           ! dummy loop indices
  INTEGER(KIND=4)                              :: npi, npj, npk, npt       ! size of the domain

  REAL(wp), DIMENSION(:,:)    , ALLOCATABLE    :: e3t_0, e3u_0, e3v_0      ! vet. metrics at rest (without vvl)
  REAL(wp), DIMENSION(:,:)    , ALLOCATABLE    :: ht_0, hu_0, hv_0         ! Reference ocean depth at T-, U-, V-points
  REAL(wp), DIMENSION(:,:)    , ALLOCATABLE    :: tmask, umask, vmask      ! Mask at T-, U-, V-points

  CHARACTER(LEN=255)                           :: cf_mz                    ! vert. mesh  netcdf file name
  CHARACTER(LEN=255)                           :: cf_mask                  ! mask        netcdf file name
  CHARACTER(LEN=256)                           :: cf_out='bathy_gdepw_0.nc'! output file name 
  CHARACTER(LEN=256)                           :: cldum                    ! dummy character variable
  CHARACTER(LEN=256)                           :: cglobal                  ! global attribute

  TYPE (variable), DIMENSION(jpvarout)         :: stypvar                  ! structure for attibutes

  LOGICAL                                      :: lchk                     ! check missing files
  LOGICAL                                      :: lnc4  =.FALSE.           ! full step flag
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdf_remake_bathy -mz ZGR-file -msk MASK-file [-o OUT-file] [-nc4]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'      Remake bathymetry from vertical integration of e3t, read in mesh_zgr.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -mz ZGR-file   : netcdf file for vertical mesh'
     PRINT *,'       -msk MASK-file : netcdf file for mask'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       -o OUT-file    : netcdf output file'
     PRINT *,'       -nc4           : Use netcdf4/HDF5 with chunking and deflation for output'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : '//TRIM(cf_out)//' unless -o option is used.'
     PRINT *,'       variables   : gdepw_0, hu_0, hv_0 '
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'      cdf_domain2bathy (can probably be merged)'
     PRINT *,'      '
     STOP
  ENDIF

  cglobal = 'Partial step computation'


  ijarg = 1
  DO WHILE ( ijarg <= narg ) 
     CALL getarg(ijarg, cldum) ; ijarg = ijarg + 1
     SELECT CASE (cldum)
     CASE ('-mz'    ) ; CALL getarg( ijarg, cf_mz    ) ; ijarg=ijarg+1
     CASE ('-msk'   ) ; CALL getarg( ijarg, cf_mask  ) ; ijarg=ijarg+1
        ! options
     CASE ( '-o'    ) ; CALL getarg(ijarg, cf_out ) ; ijarg = ijarg + 1
     CASE ( '-nc4'  ) ; lnc4    = .TRUE.
     CASE DEFAULT     ; PRINT *,' ERROR : ', TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  END DO


  !-- get dimensions (assuming all files have the same dimension that U-file) --
  npi = getdim (cf_mz, cn_x)
  npj = getdim (cf_mz, cn_y)
  npk = getdim (cf_mz, 'z')
  npt = 1

  PRINT *, '==============='
  PRINT *, 'npi=      ', npi
  PRINT *, 'npj=      ', npj
  PRINT *, 'npk=      ', npk
  PRINT *, 'npt=      ', npt
  PRINT *, '==============='

  !-- Allocate arrays --
  ALLOCATE( e3t_0(npi, npj), e3u_0(npi, npj), e3v_0(npi, npj) )
  ALLOCATE( ht_0(npi, npj) , hu_0(npi, npj) , hv_0(npi, npj)  )
  ALLOCATE( tmask(npi, npj), umask(npi, npj), vmask(npi, npj) )

  !-- Creat output netcdf files to fill in --
  CALL CreateOutput

  !-- compute (from domain.F90) --
  ht_0(:,:) = 0.0_wp                       ! Reference ocean depth at T-points
  hu_0(:,:) = 0.0_wp                       ! Reference ocean depth at U-points
  hv_0(:,:) = 0.0_wp                       ! Reference ocean depth at V-points
  DO jk = 1, npk
     !- T-points -
     e3t_0(:,:)           = getvar(cf_mz  , cn_ve3t0  , jk, npi, npj )
     tmask(:,:)           = getvar(cf_mask, cn_tmask  , jk, npi, npj ) 
     ht_0(:,:)            = ht_0(:,:) + e3t_0(:,:) * tmask(:,:)
     !- U-points -
     e3u_0(:,:)           = getvar(cf_mz  , cn_ve3u0  , jk, npi, npj )
     umask(:,:)           = getvar(cf_mask, cn_umask  , jk, npi, npj ) 
     hu_0(:,:)            = hu_0(:,:) + e3u_0(:,:) * umask(:,:)
     !- V-points -
     e3v_0(:,:)           = getvar(cf_mz  , cn_ve3v0  , jk, npi, npj )
     vmask(:,:)           = getvar(cf_mask, cn_vmask  , jk, npi, npj ) 
     hv_0(:,:)            = hv_0(:,:) + e3v_0(:,:) * vmask(:,:)
  END DO

  ierr = putvar(ncout, id_varout(1), ht_0, 1, npi, npj , 1 )
  ierr = putvar(ncout, id_varout(2), hu_0, 1, npi, npj , 1 )
  ierr = putvar(ncout, id_varout(3), hv_0, 1, npi, npj , 1 )

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
    REAL(KIND=8), DIMENSION(npt)        :: dltim                     ! time 
    ! define new variables for output
    ipk(:)                       = 1 
    !stypvar(1)%ichunk            = (/npi,MAX(1,npj/30),1,1 /)
    stypvar(1)%cname             = 'gdepw_0'
    stypvar(1)%cunits            = 'm'
    stypvar(1)%rmissing_value    = 0.
    stypvar(1)%valid_min         = 0.
    stypvar(1)%valid_max         = 9000.
    stypvar(1)%clong_name        = 'Bathymetry as vert. integ.of e3t_0'
    stypvar(1)%cshort_name       = 'Bathymetry'
    stypvar(1)%conline_operation = 'On T-grid'
    stypvar(1)%caxis             = 'nav_lon nav_lat'
    !
    stypvar(2)%ichunk            = (/npi,MAX(1,npj/30),1,1 /)
    stypvar(2)%cname             = 'hu_0'
    stypvar(2)%cunits            = 'm'
    stypvar(2)%rmissing_value    = 0.
    stypvar(2)%valid_min         = 0.
    stypvar(2)%valid_max         = 9000.
    stypvar(2)%clong_name        = 'Bathymetry as vert. integ.of e3u_0'
    stypvar(2)%cshort_name       = 'Bathymetry'
    stypvar(2)%conline_operation = 'On U-grid'
    stypvar(2)%caxis             = 'nav_lon nav_lat'
    !
    stypvar(3)%ichunk            = (/npi,MAX(1,npj/30),1,1 /)
    stypvar(3)%cname             = 'hv_0'
    stypvar(3)%cunits            = 'm'
    stypvar(3)%rmissing_value    = 0.
    stypvar(3)%valid_min         = 0.
    stypvar(3)%valid_max         = 9000.
    stypvar(3)%clong_name        = 'Bathymetry as vert. integ.of e3v_0'
    stypvar(3)%cshort_name       = 'Bathymetry'
    stypvar(3)%conline_operation = 'On V-grid'
    stypvar(3)%caxis             = 'nav_lon nav_lat'

    ! create output fileset
    ncout   = create      (cf_out  , cf_mz ,  npi, npj, 1 , ld_nc4=lnc4   )
    ierr    = createvar   (ncout , stypvar ,  jpvarout, ipk , id_varout   , ld_nc4=lnc4 )
    ierr    = putheadervar(ncout , cf_mz   ,  npi, npj, 1              )

    dltim = getvar1d(cf_mz, cn_vtimec,      npt     )
    ierr  = putvar1d(ncout  , dltim,        npt, 'T')
    ierr  = putvar1d(ncout  , dltim,        npt, 'T')


  END SUBROUTINE CreateOutput

END PROGRAM cdf_remake_bathy
