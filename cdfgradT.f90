PROGRAM cdfgradT
  !!======================================================================
  !!                     ***  PROGRAM  cdfgradT  ***
  !!=====================================================================
  !!  ** Purpose :
  !!
  !!  ** Method  :
  !!
  !! History : 3.0  : 05/2013  : N. Ducousso
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                            :: jk, jt, jvar      ! dummy loop index
  INTEGER(KIND=4)                            :: narg, iargc       ! command line
  INTEGER(KIND=4)                            :: ijarg, ireq       ! command line
  INTEGER(KIND=4)                            :: npiglo, npjglo    ! size of the domain
  INTEGER(KIND=4)                            :: npk, npt          ! size of the domain
  INTEGER(KIND=4)                            :: ncout             ! ncid of output variable
  INTEGER(KIND=4)                            :: ierr              ! error status
  INTEGER(KIND=4)                            :: iup= 1, icurr= 2  ! 
  INTEGER(KIND=4), DIMENSION(6)              :: ipk, id_varout    ! output variable

  REAL(KIND=4), DIMENSION (:),      ALLOCATABLE :: tim              
  REAL(KIND=4), DIMENSION (:,:,:),  ALLOCATABLE :: zt, zs
  REAL(KIND=4), DIMENSION (:,:),    ALLOCATABLE :: umask, vmask, wmask
  REAL(KIND=4), DIMENSION (:,:),    ALLOCATABLE :: e1u, e2v, e3w 

  REAL(KIND=8), DIMENSION (:,:),    ALLOCATABLE :: gradt_x, gradt_y, gradt_z
  REAL(KIND=8), DIMENSION (:,:),    ALLOCATABLE :: grads_x, grads_y, grads_z
 
  CHARACTER(LEN=256)                         :: cf_tfil             ! input file name
  CHARACTER(LEN=256)                         :: cf_out = 'gradT.nc' ! output file name
  CHARACTER(LEN=256), DIMENSION(2)           :: cv_namesi           ! input variable names

  TYPE(variable), DIMENSION(6)               :: stypvar             ! output data structure

  LOGICAL                                    :: lchk = .FALSE.      ! flag for missing files
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  cv_namesi(1) = cn_votemper
  cv_namesi(2) = cn_vosaline

  narg= iargc()
  IF ( narg /= 1 ) THEN
     PRINT *,' usage : cdfgradT T-file'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'      '
     STOP
  ENDIF

  CALL getarg (1, cf_tfil)
  IF (chkfile(cf_tfil) ) STOP ! missing file

  npiglo = getdim (cf_tfil, cn_x)
  npjglo = getdim (cf_tfil, cn_y)
  npk    = getdim (cf_tfil, cn_z)
  npt    = getdim (cf_tfil, cn_t)

  !!  Create output variables
  ipk(:) = npk  !  3D

  stypvar(1)%cname             = 'vozogradt'
  stypvar(1)%cunits            = ''
  stypvar(1)%rmissing_value    = -1000.
  stypvar(1)%valid_min         = -1.
  stypvar(1)%valid_max         = 1.
  stypvar(1)%clong_name        = 'zonal temper gradient'
  stypvar(1)%cshort_name       = 'vozogradt'
  stypvar(1)%conline_operation = 'N/A'
  stypvar(1)%caxis             = 'TZYX'

  stypvar(2)%cname             = 'vomegradt'
  stypvar(2)%cunits            = ''
  stypvar(2)%rmissing_value    = -1000.
  stypvar(2)%valid_min         = -1.
  stypvar(2)%valid_max         = 1.
  stypvar(2)%clong_name        = 'meridional temper gradient'
  stypvar(2)%cshort_name       = 'vomegradt'
  stypvar(2)%conline_operation = 'N/A'
  stypvar(2)%caxis             = 'TZYX'

  stypvar(3)%cname             = 'vovegradt'
  stypvar(3)%cunits            = ''
  stypvar(3)%rmissing_value    = -1000.
  stypvar(3)%valid_min         = -1.
  stypvar(3)%valid_max         = 1.
  stypvar(3)%clong_name        = 'vertical temper gradient'
  stypvar(3)%cshort_name       = 'vovegradt'
  stypvar(3)%conline_operation = 'N/A'
  stypvar(3)%caxis             = 'TZYX'

  stypvar(4)%cname             = 'vozograds'
  stypvar(4)%cunits            = ''
  stypvar(4)%rmissing_value    = -1000.
  stypvar(4)%valid_min         = -1.
  stypvar(4)%valid_max         = 1.
  stypvar(4)%clong_name        = 'zonal saline gradient'
  stypvar(4)%cshort_name       = 'vozograds'
  stypvar(4)%conline_operation = 'N/A'
  stypvar(4)%caxis             = 'TZYX'

  stypvar(5)%cname             = 'vomegrads'
  stypvar(5)%cunits            = ''
  stypvar(5)%rmissing_value    = -1000.
  stypvar(5)%valid_min         = -1.
  stypvar(5)%valid_max         = 1.
  stypvar(5)%clong_name        = 'meridional saline gradient'
  stypvar(5)%cshort_name       = 'vomegrads'
  stypvar(5)%conline_operation = 'N/A'
  stypvar(5)%caxis             = 'TZYX'

  stypvar(6)%cname             = 'vovegrads'
  stypvar(6)%cunits            = ''
  stypvar(6)%rmissing_value    = -1000.
  stypvar(6)%valid_min         = -1.
  stypvar(6)%valid_max         = 1.
  stypvar(6)%clong_name        = 'vertical saline gradient'
  stypvar(6)%cshort_name       = 'vovegrads'
  stypvar(6)%conline_operation = 'N/A'

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  !!  Allocate arrays
  ALLOCATE (tim(npt) )
  ALLOCATE (e1u(npiglo,npjglo), e2v(npiglo,npjglo), e3w(npiglo,npjglo))
  ALLOCATE (umask(npiglo,npjglo), vmask(npiglo,npjglo), wmask(npiglo,npjglo))
  ALLOCATE (zt(npiglo,npjglo,2), zs(npiglo,npjglo,2))
  ALLOCATE (gradt_x(npiglo,npjglo), gradt_y(npiglo,npjglo), gradt_z(npiglo,npjglo))
  ALLOCATE (grads_x(npiglo,npjglo), grads_y(npiglo,npjglo), grads_z(npiglo,npjglo))

  ! create output fileset
  ncout = create      (cf_out, cf_tfil, npiglo, npjglo, npk       )
  ierr  = createvar   (ncout,  stypvar, 6,      ipk,    id_varout )
  ierr  = putheadervar(ncout,  cf_tfil, npiglo, npjglo, npk       )

  tim  = getvar1d(cf_tfil, cn_vtimec, npt     )
  ierr = putvar1d(ncout,  tim,       npt, 'T')

  e1u =  getvar(cn_fhgr, cn_ve1u, 1, npiglo, npjglo)
  e2v =  getvar(cn_fhgr, cn_ve2v, 1, npiglo, npjglo)

  DO jt = 1,npt
  DO jk = npk, 1, -1  !! Main loop : (2 levels of T are required : iup, icurr)
     
     PRINT *,'level ',jk
     
     ! read files
     IF (jk == 1) THEN
        zt(:,:,iup)   = getvar(cf_tfil, cn_votemper, jk, npiglo, npjglo, ktime=jt)
        zt(:,:,icurr) = getvar(cf_tfil, cn_votemper, jk, npiglo, npjglo, ktime=jt)
        zs(:,:,iup)   = getvar(cf_tfil, cn_vosaline, jk, npiglo, npjglo, ktime=jt)
        zs(:,:,icurr) = getvar(cf_tfil, cn_vosaline, jk, npiglo, npjglo, ktime=jt)
     ELSE        
        zt(:,:,iup)   = getvar(cf_tfil, cn_votemper, jk-1, npiglo, npjglo, ktime=jt)
        zt(:,:,icurr) = getvar(cf_tfil, cn_votemper, jk, npiglo, npjglo, ktime=jt)
        zs(:,:,iup)   = getvar(cf_tfil, cn_vosaline, jk-1, npiglo, npjglo, ktime=jt)
        zs(:,:,icurr) = getvar(cf_tfil, cn_vosaline, jk, npiglo, npjglo, ktime=jt)
     END IF

     e3w(:,:) = getvar(cn_fzgr, 'e3w_ps', jk, npiglo, npjglo, ldiom=.true.)
     
     umask(:,:) = getvar(cn_fmsk, 'umask' , jk, npiglo, npjglo )
     vmask(:,:) = getvar(cn_fmsk, 'vmask' , jk, npiglo, npjglo )
     wmask(:,:) = getvar(cn_fmsk, 'tmask' , jk, npiglo, npjglo )

     ! zonal grad located at U point
     gradt_x(:,:) = 0.
     gradt_x(1:npiglo-1,:) = 1. / e1u(1:npiglo-1,:) * &
          &                 ( zt(2:npiglo,:,icurr) - zt(1:npiglo-1,:,icurr) ) * umask(1:npiglo-1,:)
     grads_x(:,:) = 0.
     grads_x(1:npiglo-1,:) = 1. / e1u(1:npiglo-1,:) * &
          &                 ( zs(2:npiglo,:,icurr) - zs(1:npiglo-1,:,icurr) ) * umask(1:npiglo-1,:)

     ! meridional grad located at V point
     gradt_y(:,:) = 0.
     gradt_y(:,1:npjglo-1) = 1. / e2v(:,1:npjglo-1) * &
          &                 ( zt(:,2:npjglo,icurr) - zt(:,1:npjglo-1,icurr) ) * vmask(:,1:npjglo-1)
     grads_y(:,:) = 0.
     grads_y(:,1:npjglo-1) = 1. / e2v(:,1:npjglo-1) * &
          &                 ( zs(:,2:npjglo,icurr) - zs(:,1:npjglo-1,icurr) ) * vmask(:,1:npjglo-1)

     ! vertical grad located at W point
     gradt_z(:,:) = 0.
     gradt_z(:,:) = 1. / e3w(:,:) * ( zt(:,:,iup) - zt(:,:,icurr) ) * wmask(:,:)
     grads_z(:,:) = 0.
     grads_z(:,:) = 1. / e3w(:,:) * ( zs(:,:,iup) - zs(:,:,icurr) ) * wmask(:,:)

     ! write
     ierr = putvar(ncout, id_varout(1), REAL(gradt_x), jk, npiglo, npjglo, ktime=jt)
     ierr = putvar(ncout, id_varout(2), REAL(gradt_y), jk, npiglo, npjglo, ktime=jt)
     ierr = putvar(ncout, id_varout(3), REAL(gradt_z), jk, npiglo, npjglo, ktime=jt)
     ierr = putvar(ncout, id_varout(4), REAL(grads_x), jk, npiglo, npjglo, ktime=jt)
     ierr = putvar(ncout, id_varout(5), REAL(grads_y), jk, npiglo, npjglo, ktime=jt)
     ierr = putvar(ncout, id_varout(6), REAL(grads_z), jk, npiglo, npjglo, ktime=jt)

  END DO
  END DO

  ierr = closeout(ncout)

END PROGRAM cdfgradT
