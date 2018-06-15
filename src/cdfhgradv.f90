PROGRAM cdfhgradv
  !!======================================================================
  !!                     ***  PROGRAM  cdfhgradv ***
  !!=====================================================================
  !!  ** Purpose : Compute the norm of horizontal gradient for a variable
  !!
  !!  ** Method  :
  !!
  !! History :   : 4.0  : 06/2018  : J.M. Molines  (from cdfhgradb)
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!   CreateOutput  :  perform all the stuff linked with output file 
  !!----------------------------------------------------------------------
  USE cdfio
  USE eos
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class derived_fields
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4), PARAMETER                  :: jp_varout = 1     ! number of output variables
  INTEGER(KIND=4)                             :: jk, jt            ! dummy loop index
  INTEGER(KIND=4)                             :: ji, jj, jlev      ! dummy loop index
  INTEGER(KIND=4)                             :: narg, ijarg, iargc! command line
  INTEGER(KIND=4)                             :: npiglo, npjglo    ! size of the domain
  INTEGER(KIND=4)                             :: npk, npkk,  npt   ! size of the domain
  INTEGER(KIND=4)                             :: ncout             ! ncid of output variable
  INTEGER(KIND=4)                             :: ierr              ! error status
  INTEGER(KIND=4), DIMENSION(jp_varout)       :: ipk, id_varout    ! output variable

  REAL(KIND=4)                                :: zdep              ! temporary scalar
  REAL(KIND=4), DIMENSION(:),     ALLOCATABLE :: gdept             ! depth of Tlevels
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: umask, vmask, tmask ! relevant mask
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: e1u, e2v          ! metrics
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: zv                ! Variables

  REAL(KIND=8)                                :: dgrav=9.81d0
  REAL(KIND=8), DIMENSION(:),     ALLOCATABLE :: dtim              ! time variable
  REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: dgradv_xu, dgradv_yv! Temperature gradient, native grid
  REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: dgradv_xt, dgradv_yt! Temperature gradient, t-point

  REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: dgradv              ! norm of the buyoyancy gradient, t-point


  CHARACTER(LEN=256)                         :: cf_in               ! input file name for T and S
  CHARACTER(LEN=256)                         :: cf_out = 'grad.nc'  ! output file name
  CHARACTER(LEN=256)                         :: cldum               ! dummy character variable
  CHARACTER(LEN=256)                         :: cv_in               ! input variable

  TYPE(variable), DIMENSION(jp_varout)       :: stypvar             ! output data structure

  LOGICAL                                    :: lchk = .FALSE.      ! flag for missing files
  LOGICAL                                    :: lnc4 = .FALSE.      ! flag for netcdf4 output with chunking and deflation

  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfhgradv -f IN-file  -v VAR-name [-o OUT-file] [-nc4] ...'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the norm of the horizontal gradient of a variable.'
     PRINT *,'       Results are saved at T points.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f IN-file  : input file'
     PRINT *,'       -v VAR-name : input variable'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'        [-o OUT-file]  : specify the name of output file instead of '
     PRINT *,'                   ',TRIM(cf_out)
     PRINT *,'        [-nc4]         : use netcdf4 chunking and deflation on output.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       ', TRIM(cn_fhgr),' ',TRIM(cn_fmsk),' and ',TRIM(cn_fzgr)
     PRINT *,'      '
     PRINT *,'     OPENMP SUPPORTED : yes '
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out),' ( unless specified with -o option)' 
     PRINT *,'         ',jp_varout,' variables : '
     PRINT *,'              vohgradb: norm of the horizontal buoyancy gradient at t-point'
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'        cdfbuoyflx  '
     PRINT *,'      '
     STOP 
  ENDIF

  ijarg =1 
  DO WHILE ( ijarg <= narg )
     CALL getarg (ijarg, cldum ) ; ijarg = ijarg + 1 
     SELECT CASE ( cldum) 
     CASE ( '-f '  ) ; CALL getarg (ijarg, cf_in  ) ; ijarg = ijarg+1
     CASE ( '-v '  ) ; CALL getarg (ijarg, cv_in  ) ; ijarg = ijarg+1
     CASE ( '-nc4 ') ; lnc4 = .TRUE.
     CASE ( '-o'   ) ; CALL getarg (ijarg, cf_out ) ; ijarg = ijarg+1
     CASE DEFAULT    ; PRINT *,' ERROR : ',TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  ENDDO

  lchk = ( lchk .OR. chkfile(cf_in) )
  lchk = ( lchk .OR. chkfile(cn_fhgr) )
  lchk = ( lchk .OR. chkfile(cn_fzgr) )
  lchk = ( lchk .OR. chkfile(cn_fmsk) )
  IF (lchk ) STOP 99 ! missing file

  npiglo = getdim (cf_in, cn_x)
  npjglo = getdim (cf_in, cn_y)
  npk    = getdim (cf_in, cn_z)
  npt    = getdim (cf_in, cn_t)

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  npkk=npk
  IF ( npk == 0 ) npkk = 1

  !!  Allocate arrays
  ALLOCATE (dtim(npt) )
  ALLOCATE (gdept(npkk) )
  ALLOCATE (e1u  (npiglo,npjglo), e2v  (npiglo,npjglo) )
  ALLOCATE (umask(npiglo,npjglo), vmask(npiglo,npjglo), tmask(npiglo,npjglo))
  ALLOCATE (zv   (npiglo,npjglo))
  ALLOCATE (dgradv_xu(npiglo,npjglo), dgradv_yv(npiglo,npjglo))
  ALLOCATE (dgradv_xt(npiglo,npjglo), dgradv_yt(npiglo,npjglo))
  ALLOCATE (dgradv(npiglo,npjglo))

  CALL CreateOutput

  e1u =  getvar(cn_fhgr, cn_ve1u, 1, npiglo, npjglo)
  e2v =  getvar(cn_fhgr, cn_ve2v, 1, npiglo, npjglo)
  gdept(:)    = getvare3(cn_fzgr, cn_gdept, npkk )

  DO jt = 1,npt
     DO jk = 1, npkk

        PRINT *,'level ',jk
        zdep = gdept(jk)

        ! read files
        zv(:,:) = getvar(cf_in, cv_in, jk,   npiglo, npjglo, ktime=jt)

        umask(:,:) = getvar(cn_fmsk, cn_umask , jk, npiglo, npjglo )
        vmask(:,:) = getvar(cn_fmsk, cn_vmask , jk, npiglo, npjglo )
        tmask(:,:) = getvar(cn_fmsk, cn_tmask , jk, npiglo, npjglo )

        ! zonal grad located at U point at current level
        dgradv_xu(:,:) = 0.d0
        dgradv_xu(1:npiglo-1,:) = 1.d0 / e1u(1:npiglo-1,:) * &
             &                 ( zv(2:npiglo,:) - zv(1:npiglo-1,:) ) * umask(1:npiglo-1,:)

        ! meridional grad located at V point at current level
        dgradv_yv(:,:) = 0.d0
        dgradv_yv(:,1:npjglo-1) = 1.d0 / e2v(:,1:npjglo-1) * &
             &                 ( zv(:,2:npjglo) - zv(:,1:npjglo-1) ) * vmask(:,1:npjglo-1)

        ! get temperature and salinity gradients at t-point 
        dgradv_xt(:,:) = 0.d0
        dgradv_yt(:,:) = 0.d0

        DO ji=2,npiglo
           DO jj=2,npjglo
              dgradv_xt(ji,jj) = 0.5*(dgradv_xu(ji-1,jj)+dgradv_xu(ji,jj))
              dgradv_yt(ji,jj) = 0.5*(dgradv_yv(ji,jj-1)+dgradv_yv(ji,jj))
           ENDDO
        ENDDO


        ! get the norm of the buoyancy gradients at t-point
        dgradv(:,:) =  SQRT( dgradv_xt(:,:) * dgradv_xt(:,:) + dgradv_yt(:,:) * dgradv_yt(:,:) )

        ! write
        ierr = putvar(ncout, id_varout(1), REAL(dgradv), jk, npiglo, npjglo, ktime=jt)
     END DO
  END DO

  ierr = closeout(ncout)

CONTAINS
  SUBROUTINE CreateOutput
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutput  ***
    !!
    !! ** Purpose :  Set up all things required for the output file, create
    !!               the file and write the header part. 
    !!
    !! ** Method  :  Use global module variables 
    !!
    !!----------------------------------------------------------------------
    ipk(:) = npkk  !  3D

    stypvar(1)%cname             = 'vohgradb'
    stypvar(1)%cunits            = 's^{-2}'
    stypvar(1)%rmissing_value    = -1000.
    stypvar(1)%valid_min         = -1.
    stypvar(1)%valid_max         = 1.
    stypvar(1)%clong_name        = 'Horizontal gradient of buoyancy.'
    stypvar(1)%cshort_name       = 'vohgradb'
    stypvar(1)%conline_operation = 'N/A'
    stypvar(1)%ichunk            = (/npiglo, MAX(1,npjglo/30), 1, 1 /)

    ! create output fileset
    ncout = create      (cf_out, cf_in, npiglo, npjglo, npk,         ld_nc4=lnc4    )
    ierr  = createvar   (ncout,  stypvar, jp_varout, ipk, id_varout, ld_nc4=lnc4    )
    ierr  = putheadervar(ncout,  cf_in, npiglo, npjglo, npk       )

    dtim = getvar1d(cf_in,  cn_vtimec,  npt     )
    ierr = putvar1d(ncout,  dtim,       npt, 'T')

  END SUBROUTINE CreateOutput

END PROGRAM cdfhgradv
