PROGRAM cdf_dynadv_ubs 
  !!======================================================================
  !!                     ***  PROGRAM  cdf_dynadv_ubs  ***
  !!=====================================================================
  !!  ** Purpose : Compute the now momentum advection trend in flux form.
  !!
  !!  ** Method  : Adapt NEMO dynadv_ubs.F90 to CDFTOOLS (cf below for further details)
  !!               following the parameter used in the configuration eNATL60.
  !!
  !!
  !! History : 4.0  : 09/2019  : Q. Jamet & J.M. Molines : Original code
  !!                : 05/2021  : Q. Jamet & J.M. Molines : Turn the computation layer per layer
  !!                                                       to avoid memory issues.
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames   ! for cdf variable names
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class derived_fields
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4), PARAMETER                   :: wp=4
  INTEGER(KIND=4), PARAMETER                   :: jpnvarout = 2            ! number of output variables

  INTEGER(KIND=4)                              :: ji, jj, jk, jt           ! dummy loop indices
  INTEGER(KIND=4)                              :: it                       ! time index for vvl
  INTEGER(KIND=4)                              :: ierr                     ! working integer
  INTEGER(KIND=4)                              :: narg, iargc, ijarg       ! 
  INTEGER(KIND=4)                              :: npiglo, npjglo, npk, npt ! size of the domain
  INTEGER(KIND=4)                              :: npim1, npjm1, npkm1      ! index for computation of grad/div
  INTEGER(KIND=4)                              :: nkkm1=1, nkk=2, nkkp1=3  ! for swapping the levels
  INTEGER(KIND=4)                              :: npkk=3                   ! number of level to load (nkkm1, nkk, nkkp1)
  INTEGER(KIND=4)                              :: ncout_u, ncout_v, ncout_ke ! ncid of output file
  INTEGER(KIND=4), DIMENSION(jpnvarout)        :: ipk                      ! level of output variables
  INTEGER(KIND=4), DIMENSION(jpnvarout)        :: id_varout_u              ! id of output variables (u-comp)
  INTEGER(KIND=4), DIMENSION(jpnvarout)        :: id_varout_v              ! id of output variables (v-comp)
  INTEGER(KIND=4), DIMENSION(jpnvarout)        :: id_varout_ke             ! id of output variables (ke-comp)

  REAL(wp), PARAMETER                          :: pp_gamma1 = 1._wp/3._wp  ! =1/4 quick      ; =1/3  3rd order UBS
  REAL(wp), PARAMETER                          :: pp_gamma2 = 1._wp/32._wp ! =0   2nd order  ; =1/32 4th order centred

  REAL(wp), DIMENSION(:),     ALLOCATABLE      :: deptht, depthu, depthv   ! z-grid (t, u,v)
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: rlon_t, rlat_t           ! t-grid hor.
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: rlon_u, rlat_u           ! u-grid hor.
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: rlon_v, rlat_v           ! v-grid hor.
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: ht_0                     ! Reference ocean depth at T-points
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: sshn                     ! now sea surface height
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: e1t, e2t                 ! horizontal metric, t-pts
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: e1u, e2u                 ! horizontal metric, u-pts
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: e1v, e2v                 ! horizontal metric, v-pts
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: e12t, r1_e12u, r1_e12v   ! face area at t-pts and inverse at u- v- pts
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: e3t_0, e3u_0, e3v_0      ! vet. metrics at rest (without vvl)
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: e3u, e3v, e3t            ! vet. metrics, u- v- t- pts
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: tmask, fmask             ! mesh masks
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: umask, vmask             ! mesh masks
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: wn                       ! 3D vert. velocity (now)
  !REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: un, vn                   ! 3D hz. velocity (now)
  !REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: ub, vb                   ! 3D hz. velocity (before)
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: adv_h_u, adv_z_u         ! hor. and vert. advection of u-mom. (outputs)
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: adv_h_v, adv_z_v         ! hor. and vert. advection of v-mom. (outputs)
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: adv_h_ke, adv_z_ke       ! hor. and vert. advection of KE     (outputs)

  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET  :: un, vn               ! 3D hz. velocity (now)
  REAL(wp), DIMENSION(:,:,:), POINTER          :: ub, vb                   ! 3D hz. velocity (before)

  CHARACTER(LEN=256)                           :: cf_tfil                  ! temperature netcdf file name (for mesh only)
  CHARACTER(LEN=256)                           :: cf_ufil                  ! zonal vel  netcdf file name
  CHARACTER(LEN=256)                           :: cf_vfil                  ! merid vel  netcdf file name
  CHARACTER(LEN=255)                           :: cf_wfil                  ! vert. vel  netcdf file name
  CHARACTER(LEN=256)                           :: cf_sshfil                ! ssh        netcdf file name (for vvl)
  CHARACTER(LEN=255)                           :: cf_mh                    ! mesh       netcdf file name
  CHARACTER(LEN=255)                           :: cf_mz                    ! mesh       netcdf file name
  CHARACTER(LEN=255)                           :: cf_mask                  ! mask       netcdf file name
  CHARACTER(LEN=255)                           :: cf_bathy                 ! bathymetry netcdf file name
  CHARACTER(LEN=256)                           :: cf_out_u ='adv_u.nc'     ! output file name (u-comp)
  CHARACTER(LEN=256)                           :: cf_out_v ='adv_v.nc'     ! output file name (v-comp)
  CHARACTER(LEN=256)                           :: cf_out_ke='adv_ke.nc'    ! output file name (ke-comp)
  CHARACTER(LEN=256)                           :: cldum                    ! dummy character variable
  CHARACTER(LEN=256)                           :: cglobal                  ! global attribute

  TYPE (variable), DIMENSION(jpnvarout)        :: stypvar                  ! structure for attibutes (u-comp)
  TYPE (variable), DIMENSION(jpnvarout)        :: stypvar2                 ! structure for attibutes (v-comp)
  TYPE (variable), DIMENSION(jpnvarout)        :: stypvar3                 ! structure for attibutes (ke-comp)

  LOGICAL                                      :: l_w   =.FALSE.           ! flag for vertical location of bn2
  LOGICAL                                      :: lchk                     ! check missing files
  LOGICAL                                      :: lnc4  =.FALSE.           ! full step flag
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdf_dynadv_ubs -t T-file -u U-file -v V-file -w W-file -ssh SSH-file ...'
     PRINT *,'          -mh MESH-file -mz MESZ-file -mask MASK-file -bathy BATHY-file ...'
     PRINT *,'          -o_u OUT-file-u -o_v OUT-file-v -o_ke OUT-file-ke'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'      Compute the now momentum advection trend in flux form.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -t T-file          : netcdf file for temperature (for mesh only)'
     PRINT *,'       -u U-file          : netcdf file for zonal velocity'
     PRINT *,'       -v V-file          : netcdf file for meridional velocity'
     PRINT *,'       -w W-file          : netcdf file for vertical velocity'
     PRINT *,'       -ssh SSH-file      : netcdf file for SSH (for vvl recomputation of vert. grid)'
     PRINT *,'       -mh MESH-file      : netcdf file for horizontal mesh'
     PRINT *,'       -mz MESZ-file      : netcdf file for vertical mesh'
     PRINT *,'       -mask MASK-file    : netcdf file for mask'
     PRINT *,'       -bathy BATHY-file  : netcdf file for model bathymetry'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       -o_u OUT-Ufile      : netcdf file for advection term for u-momentum'
     PRINT *,'       -o_v OUT-Vfile      : netcdf file for advection term for v-momentum'
     PRINT *,'       -o_ke OUT-KEfile    : netcdf file for advection term for KE'
     PRINT *,'       -nc4                : use netcdf4/HDF5 with chunking and deflation on output'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : '
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'      '
     STOP
  ENDIF

  cglobal = 'Partial step computation'

  ijarg = 1
  DO WHILE ( ijarg <= narg ) 
     CALL getarg(ijarg, cldum) ; ijarg = ijarg + 1
     SELECT CASE (cldum)
     CASE ('-t'        ) ; CALL getarg( ijarg, cf_tfil ) ; ijarg=ijarg+1
     CASE ('-u'        ) ; CALL getarg( ijarg, cf_ufil ) ; ijarg=ijarg+1
     CASE ('-v'        ) ; CALL getarg( ijarg, cf_vfil ) ; ijarg=ijarg+1
     CASE ('-w'        ) ; CALL getarg( ijarg, cf_wfil ) ; ijarg=ijarg+1
     CASE ('-ssh'      ) ; CALL getarg( ijarg, cf_sshfil); ijarg=ijarg+1
     CASE ('-mh'       ) ; CALL getarg( ijarg, cf_mh   ) ; ijarg=ijarg+1
     CASE ('-mz'       ) ; CALL getarg( ijarg, cf_mz   ) ; ijarg=ijarg+1
     CASE ('-mask'     ) ; CALL getarg( ijarg, cf_mask ) ; ijarg=ijarg+1
     CASE ('-bathy'    ) ; CALL getarg( ijarg, cf_bathy ) ; ijarg=ijarg+1
        ! options
     CASE ( '-o_u'    ) ; CALL getarg(ijarg, cf_out_u ) ; ijarg = ijarg + 1
     CASE ( '-o_v'    ) ; CALL getarg(ijarg, cf_out_v ) ; ijarg = ijarg + 1
     CASE ( '-o_ke'   ) ; CALL getarg(ijarg, cf_out_ke) ; ijarg = ijarg + 1
     CASE ( '-nc4'    ) ; lnc4    = .TRUE.
     CASE DEFAULT     ; PRINT *,' ERROR : ', TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  END DO

  !-- get dimensions (all files must have the same dimension that U-file) --
  npiglo = getdim (cf_ufil, cn_x)
  npjglo = getdim (cf_ufil, cn_y)
  npk    = getdim (cf_ufil, cn_z)
  npt    = getdim (cf_ufil, cn_t)
  npim1 = npiglo-1
  npjm1 = npjglo-1
  npkm1 = npk-1
  
  !-- summary --
  PRINT *, 'npiglo =', npiglo
  PRINT *, 'npjglo =', npjglo
  PRINT *, 'npk    =', npk
  PRINT *, 'npt    =', npt

  !-- Allocate --
  ! mesh
  ALLOCATE( deptht(npk) , depthu(npk) , depthv(npk)                )
  ALLOCATE( rlon_t(npiglo, npjglo)     , rlat_t(npiglo, npjglo)    )
  ALLOCATE( rlon_u(npiglo, npjglo)     , rlat_u(npiglo, npjglo)    )
  ALLOCATE( rlon_v(npiglo, npjglo)     , rlat_v(npiglo, npjglo)    )
  ALLOCATE( ht_0(npiglo, npjglo)                                   )
  ALLOCATE( e1t(npiglo, npjglo)           , e2t(npiglo, npjglo)    )
  ALLOCATE( e1u(npiglo, npjglo)           , e2u(npiglo, npjglo)    )
  ALLOCATE( e1v(npiglo, npjglo)           , e2v(npiglo, npjglo)    )
  ALLOCATE( e12t(npiglo, npjglo)                                   )
  ALLOCATE( r1_e12u(npiglo, npjglo)      , r1_e12v(npiglo, npjglo) )
  ALLOCATE( e3t_0(npiglo, npjglo)    , e3t(npiglo, npjglo)         )
  ALLOCATE( e3u_0( npiglo, npjglo)   , e3v_0(npiglo, npjglo)       )
  ALLOCATE( e3u(npiglo, npjglo)      , e3v(npiglo, npjglo)         )
  ALLOCATE( tmask(npiglo, npjglo)    , fmask(npiglo, npjglo)       )
  ALLOCATE( umask(npiglo, npjglo)    , vmask(npiglo, npjglo)       )
  !! variables
  ALLOCATE( sshn(npiglo, npjglo)                                   )
  ALLOCATE( un(npiglo, npjglo, npkk) , vn(npiglo, npjglo, npkk)    )
  !ALLOCATE( ub(npiglo, npjglo, npkk), vb(npiglo, npjglo, npkk)    )
  ALLOCATE( wn(npiglo, npjglo, npkk)                               )
  !
  ALLOCATE( adv_h_u(npiglo, npjglo)  , adv_z_u(npiglo, npjglo)     )
  ALLOCATE( adv_h_v(npiglo, npjglo)  , adv_z_v(npiglo, npjglo)     )
  ALLOCATE( adv_h_ke(npiglo, npjglo) , adv_z_ke(npiglo, npjglo)    )


  !!-- loading -- 
  PRINT *, '-- LOAD VARIABLES --'
  rlon_t    = getvar(cf_tfil,   cn_vlon2d, 1, npiglo, npjglo)
  rlat_t    = getvar(cf_tfil,   cn_vlat2d, 1, npiglo, npjglo)
  deptht    = getvar1d(cf_tfil, cn_vdeptht , npk)

  rlon_u    = getvar(cf_ufil,   cn_vlon2d, 1, npiglo, npjglo)
  rlat_u    = getvar(cf_ufil,   cn_vlat2d, 1, npiglo, npjglo)
  depthu    = getvar1d(cf_ufil, cn_vdepthu , npk)

  rlon_v    = getvar(cf_vfil,   cn_vlon2d, 1, npiglo, npjglo)
  rlat_v    = getvar(cf_vfil,   cn_vlat2d, 1, npiglo, npjglo)
  depthv    = getvar1d(cf_vfil, cn_vdepthv , npk)

  ht_0(:,:) = getvar(cf_bathy, 'gdepw_0', 1, npiglo, npjglo )
  !- hz. mesh -
  e1t(:,:)     = getvar(cf_mh  , cn_ve1t  , 1, npiglo, npjglo)
  e2t(:,:)     = getvar(cf_mh  , cn_ve2t  , 1, npiglo, npjglo)
  e1u(:,:)     = getvar(cf_mh  , cn_ve1u  , 1, npiglo, npjglo)
  e2u(:,:)     = getvar(cf_mh  , cn_ve2u  , 1, npiglo, npjglo)
  e1v(:,:)     = getvar(cf_mh  , cn_ve1v  , 1, npiglo, npjglo)
  e2v(:,:)     = getvar(cf_mh  , cn_ve2v  , 1, npiglo, npjglo)
  e12t(:,:)    = e1t(:,:) * e2t(:,:)
  r1_e12u(:,:) = 1._wp / (e1u(:,:) * e2u(:,:))
  r1_e12v(:,:) = 1._wp / (e1v(:,:) * e2v(:,:))

  !-- Creat output netcdf files to fill in --
  PRINT *, '-- Creat output --'
  CALL CreateOutput

  PRINT *, 'toto'

  DO jk = 1, npkm1
   PRINT *, '-- klayer: ', jk

   !-- load vert. mesh (at rest) and masks (dommsk.f90) --
   e3t_0(:,:) = getvar(cf_mz  , 'e3t_0' , jk, npiglo, npjglo)
   e3u_0(:,:) = e3t_0(:,:)
   e3v_0(:,:) = e3t_0(:,:)
   !DO jk = 1,npk                         ! Computed as the minimum of neighbooring scale factors
      DO jj = 1, npjm1
         DO ji = 1, npim1   ! vector opt.
            e3u_0 (ji,jj) = MIN( e3t_0(ji,jj), e3t_0(ji+1,jj) )
            e3v_0 (ji,jj) = MIN( e3t_0(ji,jj), e3t_0(ji,jj+1) )
         END DO
      END DO
   !END DO
   tmask(:,:) = getvar(cf_mask, cn_tmask , jk, npiglo, npjglo)
   umask(:,:) = getvar(cf_mask, cn_umask , jk, npiglo, npjglo )
   vmask(:,:) = getvar(cf_mask, cn_vmask , jk, npiglo, npjglo )
   !! fmask = 2 on lateral boundaries for no-slip bdy conditions on vorticity !!
   fmask(:,:) = getvar(cf_mask, cn_fmask , jk, npiglo, npjglo )

   DO jt = 1,npt
     !PRINT *, '======= time-step = ', jt

     !-- recomputed vert. metric due to non-linear free surface (VVL) --
     ! from domvvl.F90 -->> e3t = e3t_0*(1+sshn/ht_0)
     !IF ( jk == 1 ) THEN
       !PRINT *, '-- Recompute vert. mesh --'
     sshn(:,:)     = getvar(cf_sshfil  , cn_sossheig, 1, npiglo, npjglo, ktime=jt )
     !ENDIF
     e3t(:,:) = e3t_0(:,:) * (1 + sshn/ht_0)
     !- at u- and v- pts (domvvl.F90) -
     e3u(:,:) = e3u_0(:,:)
     e3v(:,:) = e3v_0(:,:)
     !DO jk = 1, npk
        DO jj = 1, npjm1
           DO ji = 1, npim1   ! vector opt.
              e3u(ji,jj) = e3u_0(ji,jj) + 0.5_wp * umask(ji,jj) * r1_e12u(ji,jj)                 &
                 &                       * (   e12t(ji  ,jj) * ( e3t(ji  ,jj) - e3t_0(ji  ,jj) )    &
                 &                           + e12t(ji+1,jj) * ( e3t(ji+1,jj) - e3t_0(ji+1,jj) ) )
              e3v(ji,jj) = e3v_0(ji,jj) + 0.5_wp * vmask(ji,jj) * r1_e12v(ji,jj)                 &
                 &                       * (   e12t(ji,jj  ) * ( e3t(ji,jj  ) - e3t_0(ji,jj  ) )    &
                 &                           + e12t(ji,jj+1) * ( e3t(ji,jj+1) - e3t_0(ji,jj+1) ) )
           END DO
        END DO
     !END DO

     !-- Load variables --
     !PRINT *, '-- Load hz. vel --'
     IF ( jk == 1 ) THEN
        un(:,:,nkkm1) = 0._wp
        vn(:,:,nkkm1) = 0._wp
        wn(:,:,nkkm1) = 0._wp
        un(:,:,nkk  ) = getvar(cf_ufil, cn_vozocrtx, jk  , npiglo, npjglo, ktime=jt )
        vn(:,:,nkk  ) = getvar(cf_vfil, cn_vomecrty, jk  , npiglo, npjglo, ktime=jt )
        wn(:,:,nkk  ) = getvar(cf_wfil, cn_vovecrtz, jk  , npiglo, npjglo, ktime=jt )
        un(:,:,nkkp1) = getvar(cf_ufil, cn_vozocrtx, jk+1, npiglo, npjglo, ktime=jt )
        vn(:,:,nkkp1) = getvar(cf_vfil, cn_vomecrty, jk+1, npiglo, npjglo, ktime=jt )
        wn(:,:,nkkp1) = getvar(cf_wfil, cn_vovecrtz, jk+1, npiglo, npjglo, ktime=jt )
     ELSE
        un(:,:,nkkm1) = un(:,:,nkk  )
        vn(:,:,nkkm1) = vn(:,:,nkk  )
        wn(:,:,nkkm1) = wn(:,:,nkk  )
        un(:,:,nkk  ) = un(:,:,nkkp1)
        vn(:,:,nkk  ) = vn(:,:,nkkp1)
        wn(:,:,nkk  ) = wn(:,:,nkkp1)
        un(:,:,nkkp1) = getvar(cf_ufil, cn_vozocrtx, jk+1, npiglo, npjglo, ktime=jt )
        vn(:,:,nkkp1) = getvar(cf_vfil, cn_vomecrty, jk+1, npiglo, npjglo, ktime=jt )
        wn(:,:,nkkp1) = getvar(cf_wfil, cn_vovecrtz, jk+1, npiglo, npjglo, ktime=jt )
     ENDIF
     !un(:,:,:)   = getvar3d(cf_ufil, cn_vozocrtx, npiglo, npjglo, npk, ktime=jt )
     !vn(:,:,:)   = getvar3d(cf_vfil, cn_vomecrty, npiglo, npjglo, npk, ktime=jt )
     !wn(:,:,:)   = getvar3d(cf_wfil, cn_vovecrtz, npiglo, npjglo, npk, ktime=jt )
     ub => un
     vb => vn
     !ub(:,:,:)   = getvar3d(cf_ufil, cn_vozocrtx, npiglo, npjglo, npk, ktime=jt-1 )
     !vb(:,:,:)   = getvar3d(cf_vfil, cn_vomecrty, npiglo, npjglo, npk, ktime=jt-1 )


     !-- Advection trends --
     !PRINT *, '-- Compute advection trends --'
     CALL dyn_adv_ubs( jt )

     !-- Construct KE --
     adv_h_ke(:,:) = 0._wp
     adv_z_ke(:,:) = 0._wp
     CALL trd_ken( adv_h_u, adv_h_v, adv_h_ke ) 
     CALL trd_ken( adv_z_u, adv_z_v, adv_z_ke ) 

     !-- output trends --
     !DO jk = 1, npk
        !- u-mom advection -
        ierr = putvar(ncout_u, id_varout_u(1), adv_h_u(:,:), jk, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout_u, id_varout_u(2), adv_z_u(:,:), jk, npiglo, npjglo, ktime=jt )
        !- v-mom advection -
        ierr = putvar(ncout_v, id_varout_v(1), adv_h_v(:,:), jk, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout_v, id_varout_v(2), adv_z_v(:,:), jk, npiglo, npjglo, ktime=jt )
        !- KE-mom advection -
        ierr = putvar(ncout_ke, id_varout_ke(1), adv_h_ke(:,:), jk, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout_ke, id_varout_ke(2), adv_z_ke(:,:), jk, npiglo, npjglo, ktime=jt )
     !ENDDO

  ENDDO		!jt-loop

  ENDDO		!jk-loop
 
  ierr = closeout(ncout_u)
  ierr = closeout(ncout_v)
  ierr = closeout(ncout_ke)

CONTAINS

   SUBROUTINE dyn_adv_ubs( kt )
!!----------------------------------------------------------------------
!!                  ***  ROUTINE dyn_adv_ubs  ***
!!
!! ** Purpose :   Compute the now momentum advection trend in flux form
!!              and the general trend of the momentum equation.
!!
!! ** Method  :   The scheme is the one implemeted in ROMS. It depends
!!      on two parameter gamma1 and gamma2. The former control the
!!      upstream baised part of the scheme and the later the centred
!!      part:     gamma1 = 0    pure centered  (no diffusive part)
!!                       = 1/4  Quick scheme
!!                       = 1/3  3rd order Upstream biased scheme
!!                gamma2 = 0    2nd order finite differencing
!!                       = 1/32 4th order finite differencing
!!      For stability reasons, the first term of the fluxes which cor-
!!      responds to a second order centered scheme is evaluated using
!!      the now velocity (centered in time) while the second term which
!!      is the diffusive part of the scheme, is evaluated using the
!!      before velocity (forward in time).
!!      Default value (hard coded in the begining of the module) are
!!      gamma1=1/3 and gamma2=1/32.
!!
!! ** Action : - (ua,va) updated with the 3D advective momentum trends
!!
!! Reference : Shchepetkin & McWilliams, 2005, Ocean Modelling.
!!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt     ! ocean time-step index
! JMM REM  jk is an unused argument... As an argument, (integer) the name should start by k ( jk --> kk)
!      BUT, even better to suppress it as it is not used....
!
      INTEGER  ::   ji, jj, jk            ! dummy loop indices
      REAL(wp) ::   zbu, zbv    ! temporary scalars
      REAL(wp) ::   zui, zvj, zfuj, zfvi, zl_u, zl_v   ! temporary scalars
      REAL(wp), POINTER, DIMENSION(:,:    ) ::  zfu, zfv
      REAL(wp), POINTER, DIMENSION(:,:,:  ) ::  zfw
      REAL(wp), POINTER, DIMENSION(:,:    ) ::  zfu_t, zfv_t, zfu_f, zfv_f
      REAL(wp), POINTER, DIMENSION(:,:,:  ) ::  zfu_uw, zfv_vw
      REAL(wp), POINTER, DIMENSION(:,:,:  ) ::  zlu_uu, zlv_vv, zlu_uv, zlv_vu
!!----------------------------------------------------------------------
!
      ALLOCATE( zfu(npiglo, npjglo)       , zfv(npiglo, npjglo)       )
      ALLOCATE( zfw(npiglo, npjglo, npkk)                             )
      ALLOCATE( zfu_t(npiglo, npjglo)     , zfv_t(npiglo, npjglo)     )
      ALLOCATE( zfu_f(npiglo, npjglo)     , zfv_f(npiglo, npjglo)     )
      ALLOCATE( zfu_uw(npiglo, npjglo, npkk)    , zfv_vw(npiglo, npjglo, npkk)    )
      !
      ALLOCATE( zlu_uu(npiglo, npjglo, 2) , zlv_vv(npiglo, npjglo, 2) )
      ALLOCATE( zlu_uv(npiglo, npjglo, 2) , zlv_vu(npiglo, npjglo, 2) )

!
      zfu_t(:,:)    = 0._wp
      zfv_t(:,:)    = 0._wp
      zfu_f(:,:)    = 0._wp
      zfv_f(:,:)    = 0._wp
      zfu_uw(:,:,:) = 0._wp
      zfv_vw(:,:,:) = 0._wp
!
      zlu_uu(:,:,:) = 0._wp
      zlv_vv(:,:,:) = 0._wp
      zlu_uv(:,:,:) = 0._wp
      zlv_vu(:,:,:) = 0._wp

!                                      ! =========================== !
   !   DO jk = 1, npkm1                       !  Laplacian of the velocity  !
!                                   ! =========================== !
!                                         ! horizontal volume fluxes
         zfu(:,:) = e2u(:,:) * e3u(:,:) * un(:,:,nkk)
         zfv(:,:) = e1v(:,:) * e3v(:,:) * vn(:,:,nkk)
!
         DO jj = 2, npjm1                          ! laplacian
            DO ji = 2, npim1   ! vector opt.
!
               zlu_uu(ji,jj,1) = ( ub (ji+1,jj  ,nkk) - 2.*ub (ji,jj,nkk) + ub (ji-1,jj  ,nkk) ) * umask(ji,jj)
               zlv_vv(ji,jj,1) = ( vb (ji  ,jj+1,nkk) - 2.*vb (ji,jj,nkk) + vb (ji  ,jj-1,nkk) ) * vmask(ji,jj)
               zlu_uv(ji,jj,1) = ( ub (ji  ,jj+1,nkk) - ub (ji  ,jj  ,nkk) ) * fmask(ji  ,jj  )   &
                  &               - ( ub (ji  ,jj  ,nkk) - ub (ji  ,jj-1,nkk) ) * fmask(ji  ,jj-1)
               zlv_vu(ji,jj,1) = ( vb (ji+1,jj  ,nkk) - vb (ji  ,jj  ,nkk) ) * fmask(ji  ,jj  )   &
                  &               - ( vb (ji  ,jj  ,nkk) - vb (ji-1,jj  ,nkk) ) * fmask(ji-1,jj  )
!
               zlu_uu(ji,jj,2) = ( zfu(ji+1,jj  ) - 2.*zfu(ji,jj) + zfu(ji-1,jj  ) ) * umask(ji,jj)
               zlv_vv(ji,jj,2) = ( zfv(ji  ,jj+1) - 2.*zfv(ji,jj) + zfv(ji  ,jj-1) ) * vmask(ji,jj)
               zlu_uv(ji,jj,2) = ( zfu(ji  ,jj+1) - zfu(ji  ,jj  ) ) * fmask(ji  ,jj  )   &
                  &               - ( zfu(ji  ,jj  ) - zfu(ji  ,jj-1) ) * fmask(ji  ,jj-1)
               zlv_vu(ji,jj,2) = ( zfv(ji+1,jj  ) - zfv(ji  ,jj  ) ) * fmask(ji  ,jj  )   &
                  &               - ( zfv(ji  ,jj  ) - zfv(ji-1,jj  ) ) * fmask(ji-1,jj  )
            END DO
         END DO
   !   END DO

!                                      ! ====================== !
!                                      !  Horizontal advection  !
   !   DO jk = 1, npkm1                       ! ====================== !
!                                         ! horizontal volume fluxes
         zfu(:,:) = 0.25 * e2u(:,:) * e3u(:,:) * un(:,:,nkk)
         zfv(:,:) = 0.25 * e1v(:,:) * e3v(:,:) * vn(:,:,nkk)
!
         DO jj = 1, npjm1                          ! horizontal momentum fluxes at T- and F-point
            DO ji = 1, npim1   ! vector opt.
               zui = ( un(ji,jj,nkk) + un(ji+1,jj  ,nkk) )
               zvj = ( vn(ji,jj,nkk) + vn(ji  ,jj+1,nkk) )
!
               IF (zui > 0) THEN   ;   zl_u = zlu_uu(ji  ,jj,1)
               ELSE                ;   zl_u = zlu_uu(ji+1,jj,1)
               ENDIF
               IF (zvj > 0) THEN   ;   zl_v = zlv_vv(ji,jj  ,1)
               ELSE                ;   zl_v = zlv_vv(ji,jj+1,1)
               ENDIF
!
               zfu_t(ji+1,jj  ) = ( zfu(ji,jj) + zfu(ji+1,jj  )                               &
                  &                    - pp_gamma2 * ( zlu_uu(ji,jj,2) + zlu_uu(ji+1,jj  ,2) )  )   &
                  &                * ( zui - pp_gamma1 * zl_u)
               zfv_t(ji  ,jj+1) = ( zfv(ji,jj) + zfv(ji  ,jj+1)                               &
                  &                    - pp_gamma2 * ( zlv_vv(ji,jj,2) + zlv_vv(ji  ,jj+1,2) )  )   &
                  &                * ( zvj - pp_gamma1 * zl_v)
!
               zfuj = ( zfu(ji,jj) + zfu(ji  ,jj+1) )
               zfvi = ( zfv(ji,jj) + zfv(ji+1,jj  ) )
               IF (zfuj > 0) THEN   ;    zl_v = zlv_vu( ji  ,jj  ,1)
               ELSE                 ;    zl_v = zlv_vu( ji+1,jj,1)
               ENDIF
               IF (zfvi > 0) THEN   ;    zl_u = zlu_uv( ji,jj  ,1)
               ELSE                 ;    zl_u = zlu_uv( ji,jj+1,1)
               ENDIF
!
               zfv_f(ji  ,jj  ) = ( zfvi - pp_gamma2 * ( zlv_vu(ji,jj,2) + zlv_vu(ji+1,jj  ,2) )  )   &
                  &                * ( un(ji,jj,nkk) + un(ji  ,jj+1,nkk) - pp_gamma1 * zl_u )
               zfu_f(ji  ,jj  ) = ( zfuj - pp_gamma2 * ( zlu_uv(ji,jj,2) + zlu_uv(ji  ,jj+1,2) )  )   &
                  &                * ( vn(ji,jj,nkk) + vn(ji+1,jj  ,nkk) - pp_gamma1 * zl_v )
            END DO
         END DO
         DO jj = 2, npjm1                          ! divergence of horizontal momentum fluxes
            DO ji = 2, npim1   ! vector opt.
               zbu = e1u(ji,jj) * e2u(ji,jj) * e3u(ji,jj)
               zbv = e1v(ji,jj) * e2v(ji,jj) * e3v(ji,jj)
!
               adv_h_u(ji,jj)  = - (  zfu_t(ji+1,jj  ) - zfu_t(ji  ,jj  )    &
                  &                   + zfv_f(ji  ,jj  ) - zfv_f(ji  ,jj-1)  ) / zbu
               adv_h_v(ji,jj) = - (  zfu_f(ji  ,jj  ) - zfu_f(ji-1,jj  )    &
                  &                  + zfv_t(ji  ,jj+1) - zfv_t(ji  ,jj  )  ) / zbv
!
               adv_h_u(ji,jj) = adv_h_u(ji,jj) * umask(ji,jj)
               adv_h_v(ji,jj) = adv_h_v(ji,jj) * vmask(ji,jj)
            END DO
         END DO
   !   END DO

!                                      ! ==================== !
!                                      !  Vertical advection  !
   !   DO jk = 1, npkm1                       ! ==================== !
!                                         ! Vertical volume fluxes
         zfw(:,:,nkk) = 0.25 * e1t(:,:) * e2t(:,:) * wn(:,:,nkk)
         zfw(:,:,nkkp1) = 0.25 * e1t(:,:) * e2t(:,:) * wn(:,:,nkkp1)
!
   !      IF( jk == 1 ) THEN                        ! surface/bottom advective fluxes
   !       ... moved after interior fluxes ...
   !      ELSE                                      ! interior fluxes
            DO jj = 2, npjm1
               DO ji = 2, npim1   ! vector opt.
                  zfu_uw(ji,jj,nkk) = ( zfw(ji,jj,nkk)+ zfw(ji+1,jj  ,nkk) ) * ( un(ji,jj,nkk) + un(ji,jj,nkkm1) )
                  zfv_vw(ji,jj,nkk) = ( zfw(ji,jj,nkk)+ zfw(ji  ,jj+1,nkk) ) * ( vn(ji,jj,nkk) + vn(ji,jj,nkkm1) )
                  zfu_uw(ji,jj,nkkp1) = ( zfw(ji,jj,nkkp1)+ zfw(ji+1,jj  ,nkkp1) ) * ( un(ji,jj,nkkp1) + un(ji,jj,nkk) )
                  zfv_vw(ji,jj,nkkp1) = ( zfw(ji,jj,nkkp1)+ zfw(ji  ,jj+1,nkkp1) ) * ( vn(ji,jj,nkkp1) + vn(ji,jj,nkk) )
               END DO
            END DO
    !     ENDIF
         IF( jk == npkm1 ) THEN
            zfu_uw(:,:,nkkp1) = 0.e0                      ! Bottom  value : flux set to zero
            zfv_vw(:,:,nkkp1) = 0.e0
         ENDIF
         IF ( jk == 1 ) THEN                        ! Surface value :
!            IF( lk_vvl ) THEN                                ! variable volume : flux set to zero
               zfu_uw(:,:, nkk ) = 0.e0
               zfv_vw(:,:, nkk ) = 0.e0
!            ELSE                                             ! constant volume : advection through the surface
!               ...
!            ENDIF
         ENDIF

   !   END DO
   !   DO jk = 1, npkm1                             ! divergence of vertical momentum flux divergence
         DO jj = 2, npjm1
            DO ji = 2, npim1   ! vector opt.
               adv_z_u(ji,jj) =  - ( zfu_uw(ji,jj,nkk) - zfu_uw(ji,jj,nkkp1) )    &
                  &  / ( e1u(ji,jj) * e2u(ji,jj) * e3u(ji,jj) )
               adv_z_v(ji,jj) =  - ( zfv_vw(ji,jj,nkk) - zfv_vw(ji,jj,nkkp1) )    &
                  &  / ( e1v(ji,jj) * e2v(ji,jj) * e3v(ji,jj) )
!
               adv_z_u(ji,jj) = adv_z_u(ji,jj) * umask(ji,jj)
               adv_z_v(ji,jj) = adv_z_v(ji,jj) * vmask(ji,jj)
            END DO
         END DO
   !   END DO
!
      DEALLOCATE( zfu_t , zfv_t  )
      DEALLOCATE( zfu_f , zfv_f  )
      DEALLOCATE( zfu   , zfv    )
      DEALLOCATE( zfw            )
      DEALLOCATE( zlu_uu, zlv_vv )
      DEALLOCATE( zlu_uv, zlv_vu )

   END SUBROUTINE dyn_adv_ubs


   SUBROUTINE trd_ken( putrd, pvtrd, pktrd )
!!---------------------------------------------------------------------
!!                  ***  ROUTINE trd_ken  ***
!!
!! ** Purpose :   output 3D Kinetic Energy trends using IOM
!!
!! ** Method  : - apply lbc to the input masked velocity trends
!!              - compute the associated KE trend:
!!          zke = 0.5 * (  mi-1[ un * putrd * bu ] + mj-1[ vn * pvtrd * bv]  ) / bt
!!      where bu, bv, bt are the volume of u-, v- and t-boxes.
!!              - vertical diffusion case (jpdyn_zdf):
!!          diagnose separately the KE trend associated with wind stress
!!              - bottom friction case (jpdyn_bfr):
!!          explicit case (ln_bfrimp=F): bottom trend put in the 1st level
!!                                       of putrd, pvtrd
!!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(in ) :: putrd, pvtrd   ! U and V masked trends
      REAL(wp), DIMENSION(:,:), INTENT(out) :: pktrd          ! KE trend
!
      INTEGER ::   ji, jj, jk       ! dummy loop indices
!
      REAL(wp)                              :: zrau0 = 1026._wp    ! volumic mass of reference     [kg/m3] (from phycst.F90)
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: zbu, zbv   ! volume of u- and v-boxes
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: zr1_bt    ! inverse of t-box volume
!!----------------------------------------------------------------------
!
!
      ALLOCATE( zbu(npiglo,npjglo) , zbv(npiglo,npjglo) , zr1_bt(npiglo,npjglo) )
!
!      IF ( lk_vvl .AND. kt /= nkstp ) THEN   ! Variable volume: set box volume at the 1st call of kt time step
!         nkstp = kt
   !      DO jk = 1, npkm1
            zbu   (:,:) =  e1u(:,:) * e2u(:,:) * e3u(:,:)
            zbv   (:,:) =  e1v(:,:) * e2v(:,:) * e3v(:,:)
            zr1_bt(:,:) = 1._wp / ( e12t(:,:) * e3t(:,:) ) * tmask(:,:)
   !      END DO
!      ENDIF
!
   !   pktrd(:,:,npk) = 0._wp
      pktrd(1,: ) = 0._wp
      pktrd(:,1 ) = 0._wp
   !   DO jk = 1, npkm1
         DO jj = 2, npjglo
            DO ji = 2, npiglo
               pktrd(ji,jj) = 0.5_wp * zrau0 *( un(ji  ,jj,nkk) * putrd(ji  ,jj) * zbu(ji  ,jj)  &
                  &                           + un(ji-1,jj,nkk) * putrd(ji-1,jj) * zbu(ji-1,jj)  &
                  &                           + vn(ji,jj  ,nkk) * pvtrd(ji,jj  ) * zbv(ji,jj  )  &
                  &                           + vn(ji,jj-1,nkk) * pvtrd(ji,jj-1) * zbv(ji,jj-1)  ) * zr1_bt(ji,jj)
            END DO
         END DO
   !   END DO

   END SUBROUTINE trd_ken



  SUBROUTINE CreateOutput
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutput  ***
    !!
    !! ** Purpose :  Create netcdf output file(s)
    !!
    !! ** Method  :  Use stypvar global description of variables
    !!
    !!----------------------------------------------------------------------
     REAL(KIND=8), DIMENSION(npt) :: dltim
    ! define new variables for output
    ipk(:)                        = npk
    stypvar(1)%ichunk             = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(1)%cname              = 'advh_u'
    stypvar(1)%cunits             = 'm/s^2'
    stypvar(1)%rmissing_value     = 99999.
    stypvar(1)%valid_min          = -1.
    stypvar(1)%valid_max          = 1.
    stypvar(1)%clong_name         = 'Horizontal divergence of u-momentum fluxes'
    stypvar(1)%cshort_name        = 'advh_u'
    stypvar(1)%conline_operation  = 'On u-grid'
    stypvar(1)%caxis              = 'time depthu nav_lon_u nav_lat_u'
    !
    stypvar(2)%ichunk             = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(2)%cname              = 'advz_u'
    stypvar(2)%cunits             = 'm/s^2'
    stypvar(2)%rmissing_value     = 99999.
    stypvar(2)%valid_min          = -1.
    stypvar(2)%valid_max          = 1.
    stypvar(2)%clong_name         = 'Divergence of vertical u-momentum fluxes'
    stypvar(2)%cshort_name        = 'advz_u'
    stypvar(2)%conline_operation  = 'On u-grid'
    stypvar(2)%caxis              = 'time depthu nav_lon_u nav_lat_u'

    stypvar2(1)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar2(1)%cname             = 'advh_v'
    stypvar2(1)%cunits            = 'm/s^2'
    stypvar2(1)%rmissing_value    = 99999.
    stypvar2(1)%valid_min         = -1.
    stypvar2(1)%valid_max         = 1.
    stypvar2(1)%clong_name        = 'Horizontal divergence of v-momentum fluxes'
    stypvar2(1)%cshort_name       = 'advh_v'
    stypvar2(1)%conline_operation = 'On v-grid'
    stypvar2(1)%caxis             = 'time depthv nav_lon_v nav_lat_v'
    !
    stypvar2(2)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar2(2)%cname             = 'advz_v'
    stypvar2(2)%cunits            = 'm/s^2'
    stypvar2(2)%rmissing_value    = 99999.
    stypvar2(2)%valid_min         = -1.
    stypvar2(2)%valid_max         = 1.
    stypvar2(2)%clong_name        = 'Divergence of vertical v-momentum fluxes'
    stypvar2(2)%cshort_name       = 'advz_v'
    stypvar2(2)%conline_operation = 'On v-grid'
    stypvar2(2)%caxis             = 'time depthv nav_lon_v nav_lat_v'

    stypvar3(1)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar3(1)%cname             = 'advh_ke'
    stypvar3(1)%cunits            = 'm^2/s^3'
    stypvar3(1)%rmissing_value    = 99999.
    stypvar3(1)%valid_min         = -1.
    stypvar3(1)%valid_max         = 1.
    stypvar3(1)%clong_name        = 'Horizontal divergence of KE fluxes'
    stypvar3(1)%cshort_name       = 'advh_ke'
    stypvar3(1)%conline_operation = 'On t-grid'
    stypvar3(1)%caxis             = 'time deptht nav_lon_t nav_lat_t'
    !
    stypvar3(2)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar3(2)%cname             = 'advz_ke'
    stypvar3(2)%cunits            = 'm^2/s^3'
    stypvar3(2)%rmissing_value    = 99999.
    stypvar3(2)%valid_min         = -1.
    stypvar3(2)%valid_max         = 1.
    stypvar3(2)%clong_name        = 'Divergence of vertical KE fluxes'
    stypvar3(2)%cshort_name       = 'advz_ke'
    stypvar3(2)%conline_operation = 'On t-grid'
    stypvar3(2)%caxis             = 'time deptht nav_lon_t nav_lat_t'


    

    ! create output fileset
    ncout_u = create      (cf_out_u, cf_ufil ,  npiglo, npjglo, npk, cdep=cn_vdepthu   , ld_nc4=lnc4 )
    ierr    = createvar   (ncout_u , stypvar ,  jpnvarout, ipk , id_varout_u           , ld_nc4=lnc4 )
    ierr    = putheadervar(ncout_u , cf_ufil ,  npiglo, npjglo, npk, rlon_u, rlat_u,  depthu   )

    ncout_v = create      (cf_out_v, cf_vfil ,  npiglo, npjglo, npk, cdep=cn_vdepthv  ,  ld_nc4=lnc4 )
    ierr    = createvar   (ncout_v , stypvar2,  jpnvarout, ipk , id_varout_v           , ld_nc4=lnc4 )
    ierr    = putheadervar(ncout_v , cf_vfil ,  npiglo, npjglo, npk, rlon_v, rlat_v,  depthv   )

    ncout_ke= create      (cf_out_ke, cf_tfil ,  npiglo, npjglo, npk, cdep=cn_vdeptht  , ld_nc4=lnc4 )
    ierr    = createvar   (ncout_ke , stypvar3,  jpnvarout, ipk , id_varout_ke         , ld_nc4=lnc4 )
    ierr    = putheadervar(ncout_ke , cf_tfil ,  npiglo, npjglo, npk, rlon_t, rlat_t, deptht   )


    dltim = getvar1d(cf_ufil , cn_vtimec,   npt     )
    ierr = putvar1d(ncout_u , dltim,        npt, 'T')
    ierr = putvar1d(ncout_v , dltim,        npt, 'T')
    ierr = putvar1d(ncout_ke, dltim,        npt, 'T')


  END SUBROUTINE CreateOutput


END PROGRAM
