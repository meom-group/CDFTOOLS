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
  !! @class Equation_of_state
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4), PARAMETER                   :: wp=4
  INTEGER(KIND=4), PARAMETER                   :: pnvarout = 2             ! number of output variables
  INTEGER(KIND=4)                              :: ji, jj, jk, jt           ! dummy loop indices
  INTEGER(KIND=4)                              :: it                       ! time index for vvl
  INTEGER(KIND=4)                              :: ierr                     ! working integer
  INTEGER(KIND=4)                              :: narg, iargc, ijarg       ! 
  INTEGER(KIND=4)                              :: jpiglo, jpjglo, jpk, jpt ! size of the domain
  INTEGER(KIND=4)                              :: jpim1, jpjm1, jpkm1      ! index for computation of grad/div
  INTEGER(KIND=4)                              :: jkkm1=1, jkk=2, jkkp1=3  ! for swapping the levels
  INTEGER(KIND=4)                              :: jpkk=3                   ! number of level to load (jkkm1, jkk, jkkp1)
  INTEGER(KIND=4)                              :: ncout_u, ncout_v, ncout_ke ! ncid of output file
  INTEGER(KIND=4), DIMENSION(pnvarout)         :: ipk                      ! level of output variables
  INTEGER(KIND=4), DIMENSION(pnvarout)         :: id_varout_u              ! id of output variables (u-comp)
  INTEGER(KIND=4), DIMENSION(pnvarout)         :: id_varout_v              ! id of output variables (v-comp)
  INTEGER(KIND=4), DIMENSION(pnvarout)         :: id_varout_ke             ! id of output variables (ke-comp)

  REAL(wp), PARAMETER                          :: gamma1 = 1._wp/3._wp     ! =1/4 quick      ; =1/3  3rd order UBS
  REAL(wp), PARAMETER                          :: gamma2 = 1._wp/32._wp    ! =0   2nd order  ; =1/32 4th order centred
  REAL(wp), DIMENSION(:)    , ALLOCATABLE      :: dtim                     ! time
  REAL(wp), DIMENSION(:)    , ALLOCATABLE      :: deptht, depthu, depthv   ! z-grid (t, u,v)
  REAL(wp), DIMENSION(:,:)  , ALLOCATABLE      :: nav_lon_t, nav_lat_t     ! t-grid hor.
  REAL(wp), DIMENSION(:,:)  , ALLOCATABLE      :: nav_lon_u, nav_lat_u     ! u-grid hor.
  REAL(wp), DIMENSION(:,:)  , ALLOCATABLE      :: nav_lon_v, nav_lat_v     ! v-grid hor.
  REAL(wp), DIMENSION(:,:)  , ALLOCATABLE      :: ht_0                     ! Reference ocean depth at T-points
  REAL(wp), DIMENSION(:,:)  , ALLOCATABLE      :: sshn                     ! now sea surface height
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: e1t, e2t                 ! horizontal metric, t-pts
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: e1u, e2u                 ! horizontal metric, u-pts
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: e1v, e2v                 ! horizontal metric, v-pts
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: e12t, r1_e12u, r1_e12v   ! face area at t-pts and inverse at u- v- pts
  REAL(wp), DIMENSION(:,:), ALLOCATABLE      :: e3t_0, e3u_0, e3v_0      ! vet. metrics at rest (without vvl)
  REAL(wp), DIMENSION(:,:), ALLOCATABLE      :: e3u, e3v, e3t            ! vet. metrics, u- v- t- pts
  REAL(wp), DIMENSION(:,:), ALLOCATABLE      :: tmask, fmask             ! mesh masks
  REAL(wp), DIMENSION(:,:), ALLOCATABLE      :: umask, vmask             ! mesh masks
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: wn                       ! 3D vert. velocity (now)
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET  :: un, vn               ! 3D hz. velocity (now)
  REAL(wp), DIMENSION(:,:,:), POINTER          :: ub, vb                   ! 3D hz. velocity (before)
  !REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: un, vn                   ! 3D hz. velocity (now)
  !REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: ub, vb                   ! 3D hz. velocity (before)
  REAL(wp), DIMENSION(:,:), ALLOCATABLE      :: adv_h_u, adv_z_u         ! hor. and vert. advection of u-mom. (outputs)
  REAL(wp), DIMENSION(:,:), ALLOCATABLE      :: adv_h_v, adv_z_v         ! hor. and vert. advection of v-mom. (outputs)
  REAL(wp), DIMENSION(:,:), ALLOCATABLE      :: adv_h_ke, adv_z_ke       ! hor. and vert. advection of KE     (outputs)

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

  TYPE (variable), DIMENSION(pnvarout)         :: stypvar                  ! structure for attibutes (u-comp)
  TYPE (variable), DIMENSION(pnvarout)         :: stypvar2                 ! structure for attibutes (v-comp)
  TYPE (variable), DIMENSION(pnvarout)         :: stypvar3                 ! structure for attibutes (ke-comp)

  LOGICAL                                      :: l_w   =.FALSE.           ! flag for vertical location of bn2
  LOGICAL                                      :: lchk                     ! check missing files
  LOGICAL                                      :: lfull =.FALSE.           ! full step flag
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
     PRINT *,'       -o_u OUT-file      : netcdf file for advection term for u-momentum'
     PRINT *,'       -o_v OUT-file      : netcdf file for advection term for v-momentum'
     PRINT *,'       -o_ke OUT-file     : netcdf file for advection term for KE'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
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
     CASE ( '-full' ) ; lfull   = .TRUE. ; cglobal = 'full step computation'
     CASE ( '-o_u'    ) ; CALL getarg(ijarg, cf_out_u ) ; ijarg = ijarg + 1
     CASE ( '-o_v'    ) ; CALL getarg(ijarg, cf_out_v ) ; ijarg = ijarg + 1
     CASE ( '-o_ke'   ) ; CALL getarg(ijarg, cf_out_ke) ; ijarg = ijarg + 1
     CASE ( '-nc4'  ) ; lnc4    = .TRUE.
     CASE DEFAULT     ; PRINT *,' ERROR : ', TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  END DO

  !-- get dimensions (all files must have the same dimension that U-file) --
  jpiglo = getdim (cf_ufil, cn_x)
  jpjglo = getdim (cf_ufil, cn_y)
  jpk    = getdim (cf_ufil, cn_z)
  jpt    = getdim (cf_ufil, cn_t)
  jpim1 = jpiglo-1
  jpjm1 = jpjglo-1
  jpkm1 = jpk-1
  
  !-- summary --
  PRINT *, 'jpiglo =', jpiglo
  PRINT *, 'jpjglo =', jpjglo
  PRINT *, 'jpk    =', jpk
  PRINT *, 'jpt    =', jpt

  !-- Allocate --
  ! mesh
  ALLOCATE( deptht(jpk)                   , depthu(jpk)                   , depthv(jpk)                    )
  ALLOCATE( nav_lon_t(jpiglo, jpjglo)     , nav_lat_t(jpiglo, jpjglo)      )
  ALLOCATE( nav_lon_u(jpiglo, jpjglo)     , nav_lat_u(jpiglo, jpjglo)      )
  ALLOCATE( nav_lon_v(jpiglo, jpjglo)     , nav_lat_v(jpiglo, jpjglo)      )
  ALLOCATE( ht_0(jpiglo, jpjglo)                                           )
  ALLOCATE( e1t(jpiglo, jpjglo)           , e2t(jpiglo, jpjglo)            )
  ALLOCATE( e1u(jpiglo, jpjglo)           , e2u(jpiglo, jpjglo)            )
  ALLOCATE( e1v(jpiglo, jpjglo)           , e2v(jpiglo, jpjglo)            )
  ALLOCATE( e12t(jpiglo, jpjglo)                                           )
  ALLOCATE( r1_e12u(jpiglo, jpjglo)       , r1_e12v(jpiglo, jpjglo)        )
  ALLOCATE( e3t_0(jpiglo, jpjglo)    , e3t(jpiglo, jpjglo)          )
  ALLOCATE( e3u_0( jpiglo, jpjglo)   , e3v_0(jpiglo, jpjglo)     )
  ALLOCATE( e3u(jpiglo, jpjglo)      , e3v(jpiglo, jpjglo)       )
  ALLOCATE( tmask(jpiglo, jpjglo)    , fmask(jpiglo, jpjglo)     )
  ALLOCATE( umask(jpiglo, jpjglo)    , vmask(jpiglo, jpjglo)     )
  !! variables
  ALLOCATE( sshn(jpiglo, jpjglo)                                           )
  ALLOCATE( un(jpiglo, jpjglo, jpkk)       , vn(jpiglo, jpjglo, jpkk)        )
  !ALLOCATE( ub(jpiglo, jpjglo, jpkk)       , vb(jpiglo, jpjglo, jpkk)        )
  ALLOCATE( wn(jpiglo, jpjglo, jpkk)                                        )
  !
  ALLOCATE( adv_h_u(jpiglo, jpjglo)  , adv_z_u(jpiglo, jpjglo)   )
  ALLOCATE( adv_h_v(jpiglo, jpjglo)  , adv_z_v(jpiglo, jpjglo)   )
  ALLOCATE( adv_h_ke(jpiglo, jpjglo) , adv_z_ke(jpiglo, jpjglo)  )


  !!-- loading -- 
  PRINT *, '-- LOAD VARIABLES --'
  nav_lon_t    = getvar(cf_tfil, 'nav_lon', 1, jpiglo, jpjglo)
  nav_lat_t    = getvar(cf_tfil, 'nav_lat', 1, jpiglo, jpjglo)
  deptht       = getvar1d(cf_tfil, cn_vdeptht , jpk)
  nav_lon_u    = getvar(cf_ufil, 'nav_lon', 1, jpiglo, jpjglo)
  nav_lat_u    = getvar(cf_ufil, 'nav_lat', 1, jpiglo, jpjglo)
  depthu       = getvar1d(cf_ufil, cn_vdepthu , jpk)
  nav_lon_v    = getvar(cf_vfil, 'nav_lon', 1, jpiglo, jpjglo)
  nav_lat_v    = getvar(cf_vfil, 'nav_lat', 1, jpiglo, jpjglo)
  depthv       = getvar1d(cf_vfil, cn_vdepthv , jpk)
  ht_0(:,:)    = getvar(cf_bathy, 'gdepw_0', 1, jpiglo, jpjglo )
  !- hz. mesh -
  e1t(:,:)     = getvar(cf_mh  , 'e1t'  , 1, jpiglo, jpjglo)
  e2t(:,:)     = getvar(cf_mh  , 'e2t'  , 1, jpiglo, jpjglo)
  e1u(:,:)     = getvar(cf_mh  , 'e1u'  , 1, jpiglo, jpjglo)
  e2u(:,:)     = getvar(cf_mh  , 'e2u'  , 1, jpiglo, jpjglo)
  e1v(:,:)     = getvar(cf_mh  , 'e1v'  , 1, jpiglo, jpjglo)
  e2v(:,:)     = getvar(cf_mh  , 'e2v'  , 1, jpiglo, jpjglo)
  e12t(:,:)    = e1t(:,:) * e2t(:,:)
  r1_e12u(:,:) = 1._wp / (e1u(:,:) * e2u(:,:))
  r1_e12v(:,:) = 1._wp / (e1v(:,:) * e2v(:,:))

  !-- Creat output netcdf files to fill in --
  PRINT *, '-- Creat output --'
  CALL CreateOutput

  PRINT *, 'toto'

  DO jk = 1, jpkm1
   PRINT *, '-- klayer: ', jk

   !-- load vert. mesh (at rest) and masks (dommsk.f90) --
   e3t_0(:,:) = getvar(cf_mz  , 'e3t_0' , jk, jpiglo, jpjglo)
   e3u_0(:,:) = e3t_0(:,:)
   e3v_0(:,:) = e3t_0(:,:)
   !DO jk = 1,jpk                         ! Computed as the minimum of neighbooring scale factors
      DO jj = 1, jpjm1
         DO ji = 1, jpim1   ! vector opt.
            e3u_0 (ji,jj) = MIN( e3t_0(ji,jj), e3t_0(ji+1,jj) )
            e3v_0 (ji,jj) = MIN( e3t_0(ji,jj), e3t_0(ji,jj+1) )
         END DO
      END DO
   !END DO
   tmask(:,:) = getvar(cf_mask, 'tmask' , jk, jpiglo, jpjglo)
   umask(:,:) = getvar(cf_mask, 'umask' , jk, jpiglo, jpjglo )
   vmask(:,:) = getvar(cf_mask, 'vmask' , jk, jpiglo, jpjglo )
   !! fmask = 2 on lateral boundaries for no-slip bdy conditions on vorticity !!
   fmask(:,:) = getvar(cf_mask, 'fmask' , jk, jpiglo, jpjglo )

   DO jt = 1,jpt
     !PRINT *, '======= time-step = ', jt

     !-- recomputed vert. metric due to non-linear free surface (VVL) --
     ! from domvvl.F90 -->> e3t = e3t_0*(1+sshn/ht_0)
     !IF ( jk == 1 ) THEN
       !PRINT *, '-- Recompute vert. mesh --'
     sshn(:,:)     = getvar(cf_sshfil  , cn_sossheig, 1, jpiglo, jpjglo, ktime=jt )
     !ENDIF
     e3t(:,:) = e3t_0(:,:) * (1 + sshn/ht_0)
     !- at u- and v- pts (domvvl.F90) -
     e3u(:,:) = e3u_0(:,:)
     e3v(:,:) = e3v_0(:,:)
     !DO jk = 1, jpk
        DO jj = 1, jpjm1
           DO ji = 1, jpim1   ! vector opt.
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
        un(:,:,jkkm1) = 0._wp
        vn(:,:,jkkm1) = 0._wp
        wn(:,:,jkkm1) = 0._wp
        un(:,:,jkk  ) = getvar(cf_ufil, cn_vozocrtx, jk  , jpiglo, jpjglo, ktime=jt )
        vn(:,:,jkk  ) = getvar(cf_vfil, cn_vomecrty, jk  , jpiglo, jpjglo, ktime=jt )
        wn(:,:,jkk  ) = getvar(cf_wfil, cn_vovecrtz, jk  , jpiglo, jpjglo, ktime=jt )
        un(:,:,jkkp1) = getvar(cf_ufil, cn_vozocrtx, jk+1, jpiglo, jpjglo, ktime=jt )
        vn(:,:,jkkp1) = getvar(cf_vfil, cn_vomecrty, jk+1, jpiglo, jpjglo, ktime=jt )
        wn(:,:,jkkp1) = getvar(cf_wfil, cn_vovecrtz, jk+1, jpiglo, jpjglo, ktime=jt )
     ELSE
        un(:,:,jkkm1) = un(:,:,jkk  )
        vn(:,:,jkkm1) = vn(:,:,jkk  )
        wn(:,:,jkkm1) = wn(:,:,jkk  )
        un(:,:,jkk  ) = un(:,:,jkkp1)
        vn(:,:,jkk  ) = vn(:,:,jkkp1)
        wn(:,:,jkk  ) = wn(:,:,jkkp1)
        un(:,:,jkkp1) = getvar(cf_ufil, cn_vozocrtx, jk+1, jpiglo, jpjglo, ktime=jt )
        vn(:,:,jkkp1) = getvar(cf_vfil, cn_vomecrty, jk+1, jpiglo, jpjglo, ktime=jt )
        wn(:,:,jkkp1) = getvar(cf_wfil, cn_vovecrtz, jk+1, jpiglo, jpjglo, ktime=jt )
     ENDIF
     !un(:,:,:)   = getvar3d(cf_ufil, cn_vozocrtx, jpiglo, jpjglo, jpk, ktime=jt )
     !vn(:,:,:)   = getvar3d(cf_vfil, cn_vomecrty, jpiglo, jpjglo, jpk, ktime=jt )
     !wn(:,:,:)   = getvar3d(cf_wfil, cn_vovecrtz, jpiglo, jpjglo, jpk, ktime=jt )
     ub => un
     vb => vn
     !ub(:,:,:)   = getvar3d(cf_ufil, cn_vozocrtx, jpiglo, jpjglo, jpk, ktime=jt-1 )
     !vb(:,:,:)   = getvar3d(cf_vfil, cn_vomecrty, jpiglo, jpjglo, jpk, ktime=jt-1 )


     !-- Advection trends --
     !PRINT *, '-- Compute advection trends --'
     CALL dyn_adv_ubs( jt, jk )

     !-- Construct KE --
     adv_h_ke(:,:) = 0._wp
     adv_z_ke(:,:) = 0._wp
     CALL trd_ken( adv_h_u, adv_h_v, adv_h_ke ) 
     CALL trd_ken( adv_z_u, adv_z_v, adv_z_ke ) 

     !-- output trends --
     !DO jk = 1, jpk
        !- u-mom advection -
        ierr = putvar(ncout_u, id_varout_u(1), adv_h_u(:,:), jk, jpiglo, jpjglo, ktime=jt )
        ierr = putvar(ncout_u, id_varout_u(2), adv_z_u(:,:), jk, jpiglo, jpjglo, ktime=jt )
        !- v-mom advection -
        ierr = putvar(ncout_v, id_varout_v(1), adv_h_v(:,:), jk, jpiglo, jpjglo, ktime=jt )
        ierr = putvar(ncout_v, id_varout_v(2), adv_z_v(:,:), jk, jpiglo, jpjglo, ktime=jt )
        !- KE-mom advection -
        ierr = putvar(ncout_ke, id_varout_ke(1), adv_h_ke(:,:), jk, jpiglo, jpjglo, ktime=jt )
        ierr = putvar(ncout_ke, id_varout_ke(2), adv_z_ke(:,:), jk, jpiglo, jpjglo, ktime=jt )
     !ENDDO

  ENDDO		!jt-loop

  ENDDO		!jk-loop
 
  ierr = closeout(ncout_u)
  ierr = closeout(ncout_v)
  ierr = closeout(ncout_ke)

CONTAINS

   SUBROUTINE dyn_adv_ubs( kt, jk )
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
      ALLOCATE( zfu(jpiglo, jpjglo)       , zfv(jpiglo, jpjglo)       )
      ALLOCATE( zfw(jpiglo, jpjglo, jpkk)                             )
      ALLOCATE( zfu_t(jpiglo, jpjglo)     , zfv_t(jpiglo, jpjglo)     )
      ALLOCATE( zfu_f(jpiglo, jpjglo)     , zfv_f(jpiglo, jpjglo)     )
      ALLOCATE( zfu_uw(jpiglo, jpjglo, jpkk)    , zfv_vw(jpiglo, jpjglo, jpkk)    )
      !
      ALLOCATE( zlu_uu(jpiglo, jpjglo, 2) , zlv_vv(jpiglo, jpjglo, 2) )
      ALLOCATE( zlu_uv(jpiglo, jpjglo, 2) , zlv_vu(jpiglo, jpjglo, 2) )

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
   !   DO jk = 1, jpkm1                       !  Laplacian of the velocity  !
!                                   ! =========================== !
!                                         ! horizontal volume fluxes
         zfu(:,:) = e2u(:,:) * e3u(:,:) * un(:,:,jkk)
         zfv(:,:) = e1v(:,:) * e3v(:,:) * vn(:,:,jkk)
!
         DO jj = 2, jpjm1                          ! laplacian
            DO ji = 2, jpim1   ! vector opt.
!
               zlu_uu(ji,jj,1) = ( ub (ji+1,jj  ,jkk) - 2.*ub (ji,jj,jkk) + ub (ji-1,jj  ,jkk) ) * umask(ji,jj)
               zlv_vv(ji,jj,1) = ( vb (ji  ,jj+1,jkk) - 2.*vb (ji,jj,jkk) + vb (ji  ,jj-1,jkk) ) * vmask(ji,jj)
               zlu_uv(ji,jj,1) = ( ub (ji  ,jj+1,jkk) - ub (ji  ,jj  ,jkk) ) * fmask(ji  ,jj  )   &
                  &               - ( ub (ji  ,jj  ,jkk) - ub (ji  ,jj-1,jkk) ) * fmask(ji  ,jj-1)
               zlv_vu(ji,jj,1) = ( vb (ji+1,jj  ,jkk) - vb (ji  ,jj  ,jkk) ) * fmask(ji  ,jj  )   &
                  &               - ( vb (ji  ,jj  ,jkk) - vb (ji-1,jj  ,jkk) ) * fmask(ji-1,jj  )
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
   !   DO jk = 1, jpkm1                       ! ====================== !
!                                         ! horizontal volume fluxes
         zfu(:,:) = 0.25 * e2u(:,:) * e3u(:,:) * un(:,:,jkk)
         zfv(:,:) = 0.25 * e1v(:,:) * e3v(:,:) * vn(:,:,jkk)
!
         DO jj = 1, jpjm1                          ! horizontal momentum fluxes at T- and F-point
            DO ji = 1, jpim1   ! vector opt.
               zui = ( un(ji,jj,jkk) + un(ji+1,jj  ,jkk) )
               zvj = ( vn(ji,jj,jkk) + vn(ji  ,jj+1,jkk) )
!
               IF (zui > 0) THEN   ;   zl_u = zlu_uu(ji  ,jj,1)
               ELSE                ;   zl_u = zlu_uu(ji+1,jj,1)
               ENDIF
               IF (zvj > 0) THEN   ;   zl_v = zlv_vv(ji,jj  ,1)
               ELSE                ;   zl_v = zlv_vv(ji,jj+1,1)
               ENDIF
!
               zfu_t(ji+1,jj  ) = ( zfu(ji,jj) + zfu(ji+1,jj  )                               &
                  &                    - gamma2 * ( zlu_uu(ji,jj,2) + zlu_uu(ji+1,jj  ,2) )  )   &
                  &                * ( zui - gamma1 * zl_u)
               zfv_t(ji  ,jj+1) = ( zfv(ji,jj) + zfv(ji  ,jj+1)                               &
                  &                    - gamma2 * ( zlv_vv(ji,jj,2) + zlv_vv(ji  ,jj+1,2) )  )   &
                  &                * ( zvj - gamma1 * zl_v)
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
               zfv_f(ji  ,jj  ) = ( zfvi - gamma2 * ( zlv_vu(ji,jj,2) + zlv_vu(ji+1,jj  ,2) )  )   &
                  &                * ( un(ji,jj,jkk) + un(ji  ,jj+1,jkk) - gamma1 * zl_u )
               zfu_f(ji  ,jj  ) = ( zfuj - gamma2 * ( zlu_uv(ji,jj,2) + zlu_uv(ji  ,jj+1,2) )  )   &
                  &                * ( vn(ji,jj,jkk) + vn(ji+1,jj  ,jkk) - gamma1 * zl_v )
            END DO
         END DO
         DO jj = 2, jpjm1                          ! divergence of horizontal momentum fluxes
            DO ji = 2, jpim1   ! vector opt.
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
   !   DO jk = 1, jpkm1                       ! ==================== !
!                                         ! Vertical volume fluxes
         zfw(:,:,jkk) = 0.25 * e1t(:,:) * e2t(:,:) * wn(:,:,jkk)
         zfw(:,:,jkkp1) = 0.25 * e1t(:,:) * e2t(:,:) * wn(:,:,jkkp1)
!
   !      IF( jk == 1 ) THEN                        ! surface/bottom advective fluxes
   !       ... moved after interior fluxes ...
   !      ELSE                                      ! interior fluxes
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  zfu_uw(ji,jj,jkk) = ( zfw(ji,jj,jkk)+ zfw(ji+1,jj  ,jkk) ) * ( un(ji,jj,jkk) + un(ji,jj,jkkm1) )
                  zfv_vw(ji,jj,jkk) = ( zfw(ji,jj,jkk)+ zfw(ji  ,jj+1,jkk) ) * ( vn(ji,jj,jkk) + vn(ji,jj,jkkm1) )
                  zfu_uw(ji,jj,jkkp1) = ( zfw(ji,jj,jkkp1)+ zfw(ji+1,jj  ,jkkp1) ) * ( un(ji,jj,jkkp1) + un(ji,jj,jkk) )
                  zfv_vw(ji,jj,jkkp1) = ( zfw(ji,jj,jkkp1)+ zfw(ji  ,jj+1,jkkp1) ) * ( vn(ji,jj,jkkp1) + vn(ji,jj,jkk) )
               END DO
            END DO
    !     ENDIF
         IF( jk == jpkm1 ) THEN
            zfu_uw(:,:,jkkp1) = 0.e0                      ! Bottom  value : flux set to zero
            zfv_vw(:,:,jkkp1) = 0.e0
         ENDIF
         IF ( jk == 1 ) THEN                        ! Surface value :
!            IF( lk_vvl ) THEN                                ! variable volume : flux set to zero
               zfu_uw(:,:, jkk ) = 0.e0
               zfv_vw(:,:, jkk ) = 0.e0
!            ELSE                                             ! constant volume : advection through the surface
!               ...
!            ENDIF
         ENDIF

   !   END DO
   !   DO jk = 1, jpkm1                             ! divergence of vertical momentum flux divergence
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               adv_z_u(ji,jj) =  - ( zfu_uw(ji,jj,jkk) - zfu_uw(ji,jj,jkkp1) )    &
                  &  / ( e1u(ji,jj) * e2u(ji,jj) * e3u(ji,jj) )
               adv_z_v(ji,jj) =  - ( zfv_vw(ji,jj,jkk) - zfv_vw(ji,jj,jkkp1) )    &
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


   SUBROUTINE trd_ken( putrd, pvtrd, ktrd )
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
      REAL(wp), DIMENSION(:,:), INTENT(in   ) ::   putrd, pvtrd   ! U and V masked trends
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   ktrd           ! KE trend
!
      INTEGER ::   ji, jj, jk       ! dummy loop indices
!
      REAL(wp)                                :: rau0 = 1026._wp    ! volumic mass of reference     [kg/m3] (from phycst.F90)
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   bu, bv   ! volume of u- and v-boxes
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   r1_bt    ! inverse of t-box volume
!!----------------------------------------------------------------------
!
!
      ALLOCATE( bu(jpiglo,jpjglo) , bv(jpiglo,jpjglo) , r1_bt(jpiglo,jpjglo) )
!
!      IF ( lk_vvl .AND. kt /= nkstp ) THEN   ! Variable volume: set box volume at the 1st call of kt time step
!         nkstp = kt
   !      DO jk = 1, jpkm1
            bu   (:,:) =  e1u(:,:) * e2u(:,:) * e3u(:,:)
            bv   (:,:) =  e1v(:,:) * e2v(:,:) * e3v(:,:)
            r1_bt(:,:) = 1._wp / ( e12t(:,:) * e3t(:,:) ) * tmask(:,:)
   !      END DO
!      ENDIF
!
   !   ktrd(:,:,jpk) = 0._wp
      ktrd(1,: ) = 0._wp
      ktrd(:,1 ) = 0._wp
   !   DO jk = 1, jpkm1
         DO jj = 2, jpjglo
            DO ji = 2, jpiglo
               ktrd(ji,jj) = 0.5_wp * rau0 *( un(ji  ,jj,jkk) * putrd(ji  ,jj) * bu(ji  ,jj)  &
                  &                            + un(ji-1,jj,jkk) * putrd(ji-1,jj) * bu(ji-1,jj)  &
                  &                            + vn(ji,jj  ,jkk) * pvtrd(ji,jj  ) * bv(ji,jj  )  &
                  &                            + vn(ji,jj-1,jkk) * pvtrd(ji,jj-1) * bv(ji,jj-1)  ) * r1_bt(ji,jj)
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
    ! define new variables for output
    ipk(:)                        = jpk
    stypvar(1)%ichunk             = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
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
    stypvar(2)%ichunk             = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
    stypvar(2)%cname              = 'advz_u'
    stypvar(2)%cunits             = 'm/s^2'
    stypvar(2)%rmissing_value     = 99999.
    stypvar(2)%valid_min          = -1.
    stypvar(2)%valid_max          = 1.
    stypvar(2)%clong_name         = 'Divergence of vertical u-momentum fluxes'
    stypvar(2)%cshort_name        = 'advz_u'
    stypvar(2)%conline_operation  = 'On u-grid'
    stypvar(2)%caxis              = 'time depthu nav_lon_u nav_lat_u'

    stypvar2(1)%ichunk            = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
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
    stypvar2(2)%ichunk            = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
    stypvar2(2)%cname             = 'advz_v'
    stypvar2(2)%cunits            = 'm/s^2'
    stypvar2(2)%rmissing_value    = 99999.
    stypvar2(2)%valid_min         = -1.
    stypvar2(2)%valid_max         = 1.
    stypvar2(2)%clong_name        = 'Divergence of vertical v-momentum fluxes'
    stypvar2(2)%cshort_name       = 'advz_v'
    stypvar2(2)%conline_operation = 'On v-grid'
    stypvar2(2)%caxis             = 'time depthv nav_lon_v nav_lat_v'

    stypvar3(1)%ichunk            = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
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
    stypvar3(2)%ichunk            = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
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
    ncout_u = create      (cf_out_u, cf_ufil ,  jpiglo, jpjglo, jpk, cdep=cn_vdepthu  , ld_nc4=lnc4 )
    ierr    = createvar   (ncout_u , stypvar ,  pnvarout, ipk , id_varout_u           , ld_nc4=lnc4 )
    ierr    = putheadervar(ncout_u , cf_ufil ,  jpiglo, jpjglo, jpk, nav_lon_u, nav_lat_u, depthu   )

    ncout_v = create      (cf_out_v, cf_vfil ,  jpiglo, jpjglo, jpk, cdep=cn_vdepthv  , ld_nc4=lnc4 )
    ierr    = createvar   (ncout_v , stypvar2,  pnvarout, ipk , id_varout_v           , ld_nc4=lnc4 )
    ierr    = putheadervar(ncout_v , cf_vfil ,  jpiglo, jpjglo, jpk, nav_lon_v, nav_lat_v, depthv   )

    ncout_ke= create      (cf_out_ke, cf_tfil ,  jpiglo, jpjglo, jpk, cdep=cn_vdeptht  , ld_nc4=lnc4 )
    ierr    = createvar   (ncout_ke , stypvar3,  pnvarout, ipk , id_varout_ke          , ld_nc4=lnc4 )
    ierr    = putheadervar(ncout_ke , cf_tfil ,  jpiglo, jpjglo, jpk, nav_lon_t, nav_lat_t, deptht   )


    dtim = getvar1d(cf_ufil , cn_vtimec,   jpt     )
    ierr = putvar1d(ncout_u , dtim,        jpt, 'T')
    ierr = putvar1d(ncout_v , dtim,        jpt, 'T')
    ierr = putvar1d(ncout_ke, dtim,        jpt, 'T')


  END SUBROUTINE CreateOutput


END PROGRAM
