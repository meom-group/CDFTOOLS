PROGRAM cdf_dynspg_ts
  !!======================================================================
  !!                     ***  PROGRAM  cdf_dynspg_ts  ***
  !!=====================================================================
  !!  ** Purpose : Compute the trend due to surface pressure gradient.
  !!               In the variable volume level (vvl) case,
  !!               surface pressure gradients are computed in dynhpg.F90.
  !!               In NEMO, the trends computed by dyn_spg are corrections, i.e.
  !!                  ztrdu = ( ua - ub ) / z2dt - ztrdu
  !!               Here, previous trends (ztrdu in the right hand side)
  !!               are set to zeros (no communication between CDFTOOLS).
  !!               This need to be accounted for in the postprocessing 
  !!               for a full budget analysis.
  !!
  !!  ** Method  : Adapted from dynspg_ts.F90.
  !!
  !! History : 4.0  : 09/2019  : Q. Jamet & J.M. Molines : Original code
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
  INTEGER(KIND=4), PARAMETER                   :: pnvarout = 1             ! number of output variables
  INTEGER(KIND=4)                              :: ji, jj, jk, jt           ! dummy loop indices
  INTEGER(KIND=4)                              :: it                       ! time index for vvl
  INTEGER(KIND=4)                              :: iii                      ! temporaru index
  INTEGER(KIND=4)                              :: ierr                     ! working integer
  INTEGER(KIND=4)                              :: narg, iargc, ijarg       ! 
  INTEGER(KIND=4)                              :: jpiglo, jpjglo, jpk, jpt ! size of the domain
  INTEGER(KIND=4)                              :: jpim1, jpjm1, jpkm1      ! index - 1
  INTEGER(KIND=4)                              :: ncout_u, ncout_v, ncout_ke  ! ncid of output file
  INTEGER(KIND=4)                              :: ncout_ssh                ! ncid of output file
  INTEGER(KIND=4)                              :: icycle                   ! Number of barotropic sub-steps for each internal step nn_baro <= 2.5 nn_baro
  INTEGER(KIND=4), DIMENSION(pnvarout)         :: ipk                      ! level of output variables
  INTEGER(KIND=4), DIMENSION(pnvarout)         :: id_varout_u              ! id of output variables (u-comp)
  INTEGER(KIND=4), DIMENSION(pnvarout)         :: id_varout_v              ! id of output variables (v-comp)
  INTEGER(KIND=4), DIMENSION(pnvarout)         :: id_varout_ke             ! id of output variables (ke-comp)
  INTEGER(KIND=4), DIMENSION(pnvarout)         :: id_varout_ssh            ! id of output variables (ssh, for debugging)
  INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: mbkt, mbku, mbkv         ! vert indices of the shallowest ocean level
  INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: mbathy                   ! number of vertical wet points

  REAL(wp), DIMENSION(:)  ,   ALLOCATABLE      :: dtim                     ! time
  REAL(wp), DIMENSION(:)  ,   ALLOCATABLE      :: deptht, depthu, depthv   ! z-grid (t, u,v)
  REAL(wp), DIMENSION(:)  ,   ALLOCATABLE      :: wgtbtp1, wgtbtp2         ! weights used for time filtering of barotropic variables
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: nav_lon_u, nav_lat_u     ! u-grid hor.
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: nav_lon_v, nav_lat_v     ! v-grid hor.
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: nav_lon_t, nav_lat_t     ! t-grid hor.
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: e1t, e2t                 ! horizontal metrics (t)
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: e12t, r1_e12t            ! horizontal cell surface and inverse at t points
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: e1u, e2u, e1v, e2v       ! horizontal metrics (u,v)
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: r1_e12u, r1_e12v         ! inverse of horizontal cell surface u-, v- points
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: sshn                     ! Sea surface height
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: msl                      ! Atmospheric sea level pressure [Pa]
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: ssh_ib                   ! Inverse barometer sea surface height   [m]
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: emp, rnf                 ! freshwater budget: volume flux ; river runoff  [Kg/m2/s]
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: ht_0, hu_0, hv_0         ! Reference ocean depth at t-, u-, v-points
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: ht, hu, hv               ! Local depth at t-, u-, v-points (meters)
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: hur, hvr                 ! Inverse of local depth at u-, v-points
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: ff                       ! coriolis factor (2.*omega*sin(yphi) ) (s-1)
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: zwz                      ! ff/h at F points
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: ftnw, ftne               ! triad of coriolis parameter
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: ftsw, ftse               ! (only used with een vorticity scheme)
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: bfrua, bfrva             ! bottom friction
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: bfrcoef2d                ! 2D bottom drag coefficient
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: utau, vtau               ! Wind stress
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: tmask, umask, vmask      ! Mask at t-, u- and v-points
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: e3t_0, e3u_0, e3v_0      ! vet. metrics at rest (without vvl)
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: e3u, e3v, e3t            ! vet. metrics, u- v- t- pts
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: un_b, vn_b               ! horizontal barotorpic velocity (n=now)
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: un, vn                   ! horizontal velocity (n=now)
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: ua, va                   ! recomputed trends, ua = (ua-ub) / (2*\Delta t)
  !- before fiel for centered integration -
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: sshb                     ! BEFORE Sea surface height
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: ssh_ibb                  !        Inverse barometer sea surface height   [m]
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: emp_b, rnf_b             !        freshwater budget: volume flux ; river runoff  [Kg/m2/s]
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: utau_b, vtau_b           !        Wind stress
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: ub_b, vb_b               !        Barotropic velocities
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: ehu_b , ehv_b            !        Local depth
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: ehur_b , ehvr_b          !        Inverse of local depth
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: ub, vb                   !        Horizontal velocities 
  !- updated barotropic variables -
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: ssha                     ! updated ssh
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: un_adv, vn_adv           ! Advection vel. at "now" barocl. step
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: ua_b, va_b               ! updated horizontal barotorpic velocity (a=after)
  ! For barotropic loop (befor, now and after variab. need to be allocated)
  REAL(wp)                                     :: rdtbt                    ! Barotropic time step
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: sshbb_e, sshb_e          ! Instantaneous barotropic arrays
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: ubb_e, ub_e              !
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: vbb_e, vb_e              !
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: ua_e, va_e               ! barotropic velocities (after)
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: sshn_e, ssha_e           ! Sea surface height
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: hu_e, hv_e               ! Reference ocean depth at u-, v-points
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: hur_e, hvr_e             ! Inverse of local depth at u-, v-points 
  !-- outputs --
  !REAL(wp), DIMENSION(:,:),   ALLOCATABLE     :: tmpu, tmpv               ! temporary array for debugging
  !REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: tmpu, tmpv               ! temporary array for debugging
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: ssh_bt                   ! temporary array for debugging
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: spgu                     ! Barotropic momentum correction trend (u-mom)
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: spgv                     ! Barotropic momentum correction trend (v-mom)
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: spg_ke                   ! Barotropic momentum correction trend (ke)

  CHARACTER(LEN=256)                           :: cf_tfil                  ! temperature netcdf file name (for mesh only)
  CHARACTER(LEN=256)                           :: cf_ufil                  ! zonal velocity netcdf file name
  CHARACTER(LEN=256)                           :: cf_vfil                  ! meridional velocity netcdf file name
  CHARACTER(LEN=256)                           :: cf_sshfil                ! ssh         netcdf file name (for vvl)
  CHARACTER(LEN=255)                           :: cf_mh                    ! horiz. mesh netcdf file name
  CHARACTER(LEN=255)                           :: cf_mz                    ! mesh       netcdf file name
  CHARACTER(LEN=256)                           :: cf_tauxfil               ! zonal stress netcdf file name
  CHARACTER(LEN=256)                           :: cf_tauyfil               ! merid stress netcdf file name
  CHARACTER(LEN=256)                           :: cf_mslfil                ! atmospheric sea level pressure [Pa] netcdf file name
  CHARACTER(LEN=255)                           :: cf_mask                  ! mask      netcdf file name
  CHARACTER(LEN=255)                           :: cf_bathy                 ! bathymetry netcdf file name
  CHARACTER(LEN=256)                           :: cf_out_u='spg_u.nc'      ! output file name (u-comp)
  CHARACTER(LEN=256)                           :: cf_out_v='spg_v.nc'      ! output file name (v-comp)
  CHARACTER(LEN=256)                           :: cf_out_ke='spg_ke.nc'    ! output file name (ke-comp)
  CHARACTER(LEN=256)                           :: cf_out_ssh='spg_ssh.nc'  ! output file name (ssh, for debugging)
  CHARACTER(LEN=256)                           :: cldum                    ! dummy character variable
  CHARACTER(LEN=256)                           :: cglobal                  ! global attribute

  TYPE (variable), DIMENSION(pnvarout)         :: stypvar                  ! structure for attibutes
  TYPE (variable), DIMENSION(pnvarout)         :: stypvar2                 ! structure for attibutes
  TYPE (variable), DIMENSION(pnvarout)         :: stypvar3                 ! structure for attibutes
  TYPE (variable), DIMENSION(pnvarout)         :: stypvar4                 ! structure for attibutes

  LOGICAL                                      :: l_w   =.FALSE.           ! flag for vertical location of bn2
  LOGICAL                                      :: lchk                     ! check missing files
  LOGICAL                                      :: lfull =.FALSE.           ! full step flag
  LOGICAL                                      :: lnc4  =.FALSE.           ! full step flag

  ! from namelist
  REAL(wp)                                     :: rdt = 80._wp             ! time step [sec] == model output frequency
  REAL(wp)                                     :: rn_pref = 101000._wp     ! reference atmospheric pressure   [N/m2]
  !LOGICAL                                      :: ln_bt_fw = .false.       ! Backward integration of barotropic equations (need u2b and v2b)
  !LOGICAL                                      :: ln_bt_av = .true.        ! Time filtering of barotropic variables
  INTEGER(KIND=4)                              :: nn_baro = 30             ! Number of iterations of barotropic mode
  !INTEGER(KIND=4)                              :: nn_bt_flt = 2            ! Time filter choice (= 2 Boxcar over 2*nn_baro) HARD CODED
  REAL(wp)                                     :: vkarmn = 0.4_wp          ! von Karman constant
  REAL(wp)                                     :: rn_bfri2 = 2.5e-3_wp     ! min. bottom drag coefficient
  REAL(wp)                                     :: rn_bfri2_max = 1.e-1_wp  ! max. bottom drag coefficient
  REAL(wp)                                     :: rn_bfrien = 50._wp       ! local multiplying factor of bfr
  !REAL(wp)                                     :: rn_bfeb2 = 2.5e-3_wp     ! bottom turb. KE background (for tide free run) [m2/s2]
  REAL(wp)                                     :: rn_bfeb2 = 0.0_wp        ! bottom turb. KE background (for tide free run) [m2/s2]
                                                                           !         (=0.0_wp in tides conditions)
  ! from phycst.f90
  REAL(wp)                                     :: rau0 = 1026._wp          ! volumic mass of reference     [kg/m3]
  REAL(wp)                                     :: grav = 9.80665_wp        ! gravity [m/s2]
  REAL(wp)                                     :: r1_grau                  ! = 1.e0 / (grav * rau0)
  REAL(wp)                                     :: rn_bfrz0 = 3.e-3_wp      ! bottom roughness [m] for log. formulation


  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdf_dynspg_ts -t T-file -u U-file -v V-file -ssh SSH-file -tx TAUX-file -ty TAUY-file -msl MSL-file ...'
     PRINT *,'          -mh MESH-file -mz MESZ-file -mask MASK-file -bathy BATHY-file ...'
     PRINT *,'          -o_u OUT-file-U -o_v OUT-file-V -o_ke OUT-file-ke -o_ssh OUT-file-ssh'
     PRINT *,'      '
     PRINT *,'     PURPOSE :Compute the trend due to surface pressure gradient.'
     PRINT *,'              Adapt dynspg_exp.F90.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -t T-file          : netcdf file for temperature (for mesh only)'
     PRINT *,'       -u U-file          : netcdf file for zonal velocity'
     PRINT *,'       -v V-file          : netcdf file for meridional velocity'
     PRINT *,'       -ssh SSH-file      : netcdf file for SSH (for vvl recomputation of vert. grid)'
     PRINT *,'       -tx TAUX-file      : netcdf file for zonal wind stress'
     PRINT *,'       -ty TAUY-file      : netcdf file for meridional wind stress'
     PRINT *,'       -msl MSL-file      : netcdf file for atmospheric sea level pressure [Pa]'
     PRINT *,'       -mh MESH-file      : netcdf file for horizontal mesh'
     PRINT *,'       -mz MESZ-file      : netcdf file for vertical mesh'
     PRINT *,'       -mask MASK-file    : netcdf file for mask'
     PRINT *,'       -bathy BATHY-file  : netcdf file for model bathymetry'
     PRINT *,'       -o_u OUT-file      : netcdf file for surf. press. grad. corr. term for u-momentum'
     PRINT *,'       -o_v OUT-file      : netcdf file for surf. press. grad. corr. term for v-momentum'
     PRINT *,'       -o_ke OUT-file     : netcdf file for surf. press. grad. corr. term for KE'
     PRINT *,'       -o_ssh OUT-file    : netcdf file for SSH associated with surf. press. corr. (for debugging)'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : this is the netcdf file that will be outputed'
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
     CASE ('-t'        ) ; CALL getarg( ijarg, cf_tfil )   ; ijarg=ijarg+1
     CASE ('-u'        ) ; CALL getarg( ijarg, cf_ufil   ) ; ijarg=ijarg+1
     CASE ('-v'        ) ; CALL getarg( ijarg, cf_vfil   ) ; ijarg=ijarg+1
     CASE ('-ssh'      ) ; CALL getarg( ijarg, cf_sshfil ) ; ijarg=ijarg+1
     CASE ('-tx'       ) ; CALL getarg( ijarg, cf_tauxfil) ; ijarg=ijarg+1
     CASE ('-ty'       ) ; CALL getarg( ijarg, cf_tauyfil) ; ijarg=ijarg+1
     CASE ('-msl'      ) ; CALL getarg( ijarg, cf_mslfil)  ; ijarg=ijarg+1
     CASE ('-mh'       ) ; CALL getarg( ijarg, cf_mh     ) ; ijarg=ijarg+1
     CASE ('-mz'       ) ; CALL getarg( ijarg, cf_mz   )   ; ijarg=ijarg+1
     CASE ('-mask'     ) ; CALL getarg( ijarg, cf_mask   ) ; ijarg=ijarg+1
     CASE ('-bathy'    ) ; CALL getarg( ijarg, cf_bathy )  ; ijarg=ijarg+1
        ! options
     CASE ( '-full' ) ; lfull   = .TRUE. ; cglobal = 'full step computation'
     CASE ( '-o_u'    ) ; CALL getarg(ijarg, cf_out_u ) ; ijarg = ijarg + 1
     CASE ( '-o_v'    ) ; CALL getarg(ijarg, cf_out_v ) ; ijarg = ijarg + 1
     CASE ( '-o_ke'   ) ; CALL getarg(ijarg, cf_out_ke) ; ijarg = ijarg + 1
     CASE ( '-o_ssh'  ) ; CALL getarg(ijarg, cf_out_ssh); ijarg = ijarg + 1
     CASE ( '-nc4'  ) ; lnc4    = .TRUE.
     CASE DEFAULT     ; PRINT *,' ERROR : ', TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  END DO

  !-- get dimensions (assuming all files have the same dimension that U-file) --
  jpiglo = getdim (cf_ufil, cn_x)
  jpjglo = getdim (cf_ufil, cn_y)
  jpk    = getdim (cf_ufil, cn_z)
  jpt    = getdim (cf_ufil, cn_t)
  jpim1  = jpiglo-1
  jpjm1  = jpjglo-1
  jpkm1  = jpk-1

  !-- summary --
  PRINT *, 'jpiglo =', jpiglo
  PRINT *, 'jpjglo =', jpjglo
  PRINT *, 'jpt    =', jpt

  !-- Allocate arrays --
  ALLOCATE( deptht(jpk)                 , depthu(jpk)                 , depthv(jpk)                    )
  ALLOCATE( nav_lon_u(jpiglo, jpjglo)   , nav_lat_u(jpiglo, jpjglo)   )
  ALLOCATE( nav_lon_v(jpiglo, jpjglo)   , nav_lat_v(jpiglo, jpjglo)   )
  ALLOCATE( nav_lon_t(jpiglo, jpjglo)   , nav_lat_t(jpiglo, jpjglo)   )
  ALLOCATE( e1t(jpiglo, jpjglo)         , e2t(jpiglo, jpjglo)         )
  ALLOCATE( e12t(jpiglo, jpjglo)        , r1_e12t(jpiglo, jpjglo)     )
  ALLOCATE( e1u(jpiglo, jpjglo)         , e2u(jpiglo, jpjglo)         )
  ALLOCATE( e1v(jpiglo, jpjglo)         , e2v(jpiglo, jpjglo)         )
  ALLOCATE( r1_e12u(jpiglo, jpjglo)     , r1_e12v(jpiglo, jpjglo)     )
  ALLOCATE( sshn(jpiglo, jpjglo)                                      )
  ALLOCATE( msl(jpiglo, jpjglo)         , ssh_ib(jpiglo, jpjglo)      )
  ALLOCATE( emp(jpiglo, jpjglo)         , rnf(jpiglo, jpjglo)         )
  ALLOCATE( ht_0(jpiglo, jpjglo)        , ht(jpiglo, jpjglo)          )
  ALLOCATE( hu_0(jpiglo, jpjglo)        , hv_0(jpiglo, jpjglo)        )
  ALLOCATE( hu(jpiglo, jpjglo)          , hv(jpiglo, jpjglo)          )
  ALLOCATE( hur(jpiglo, jpjglo)         , hvr(jpiglo, jpjglo)         )
  ALLOCATE( ff(jpiglo, jpjglo)                                        )
  ALLOCATE( e3t_0(jpiglo, jpjglo, jpk)  , e3t(jpiglo, jpjglo, jpk)    )
  ALLOCATE( e3u_0(jpiglo, jpjglo, jpk)  , e3v_0(jpiglo, jpjglo, jpk)  )
  ALLOCATE( e3u(jpiglo, jpjglo, jpk)    , e3v(jpiglo, jpjglo, jpk)    )
  ALLOCATE( tmask(jpiglo, jpjglo, jpk)                                )
  ALLOCATE( umask(jpiglo, jpjglo, jpk)  , vmask(jpiglo, jpjglo, jpk)  )
  ALLOCATE( mbathy(jpiglo, jpjglo)                                    )
  ALLOCATE( mbkt(jpiglo, jpjglo)        , mbku(jpiglo, jpjglo)        , mbkv(jpiglo, jpjglo)  )
  ALLOCATE( utau(jpiglo, jpjglo)        , vtau(jpiglo, jpjglo)        )
  !
  ALLOCATE( bfrcoef2d(jpiglo, jpjglo)                                 )
  ALLOCATE( bfrua(jpiglo, jpjglo)       , bfrva(jpiglo, jpjglo)       )
  ALLOCATE( wgtbtp1(3*nn_baro)          , wgtbtp2(3*nn_baro)          )
  ALLOCATE( zwz(jpiglo,jpjglo)                                        )
  ALLOCATE( ftnw(jpiglo,jpjglo)         , ftne(jpiglo,jpjglo)         )
  ALLOCATE( ftsw(jpiglo,jpjglo)         , ftse(jpiglo,jpjglo)         )
  ALLOCATE( un(jpiglo, jpjglo, jpk)     , vn(jpiglo, jpjglo, jpk)     )
  ALLOCATE( ua(jpiglo, jpjglo, jpk)     , va(jpiglo, jpjglo, jpk)     )
  ALLOCATE( un_b(jpiglo, jpjglo)        , vn_b(jpiglo, jpjglo)        )
  !- before fiel for centered integration -
  ALLOCATE( sshb(jpiglo,jpjglo)                                       )
  ALLOCATE( ssh_ibb(jpiglo, jpjglo)                                   )
  ALLOCATE( emp_b(jpiglo, jpjglo)       , rnf_b(jpiglo, jpjglo)       )
  ALLOCATE( utau_b(jpiglo, jpjglo)      , vtau_b(jpiglo, jpjglo)      )
  ALLOCATE( ub_b(jpiglo, jpjglo)        , vb_b(jpiglo, jpjglo)        )
  ALLOCATE( ehu_b(jpiglo, jpjglo)       , ehv_b(jpiglo, jpjglo)       )
  ALLOCATE( ehur_b(jpiglo, jpjglo)      , ehvr_b(jpiglo, jpjglo)      )
  ALLOCATE( ub(jpiglo, jpjglo, jpk)     , vb(jpiglo, jpjglo, jpk)     )
  !- updated barotropic variables -
  ALLOCATE( ssha(jpiglo,jpjglo)                                       )
  ALLOCATE( un_adv(jpiglo,jpjglo)       , vn_adv(jpiglo,jpjglo)       )
  ALLOCATE( ua_b(jpiglo, jpjglo)        , va_b(jpiglo, jpjglo)        )
  ! For barotropic loop (befor, now and after variab. need to be allocated)
  ALLOCATE( sshbb_e(jpiglo,jpjglo)      , sshb_e(jpiglo,jpjglo)       )
  ALLOCATE( ubb_e(jpiglo,jpjglo)        , ub_e(jpiglo,jpjglo)         )
  ALLOCATE( vbb_e(jpiglo,jpjglo)        , vb_e(jpiglo,jpjglo)         )
  ALLOCATE( ua_e(jpiglo,jpjglo)         , va_e(jpiglo,jpjglo)         )
  ALLOCATE( sshn_e(jpiglo,jpjglo)       , ssha_e(jpiglo,jpjglo)       )
  ALLOCATE( hu_e(jpiglo, jpjglo)        , hv_e(jpiglo, jpjglo)        )
  ALLOCATE( hur_e(jpiglo, jpjglo)       , hvr_e(jpiglo, jpjglo)       )
  !-- outputs --
  !ALLOCATE( tmpu(jpiglo,jpjglo)         , tmpv(jpiglo,jpjglo)         )
  !ALLOCATE( tmpu(jpiglo,jpjglo,100)     , tmpv(jpiglo,jpjglo,100)     )
  ALLOCATE( spgu(jpiglo,jpjglo, jpk)    , spgv(jpiglo,jpjglo, jpk)    )
  ALLOCATE( spg_ke(jpiglo,jpjglo, jpk)                                )
  ALLOCATE( ssh_bt(jpiglo,jpjglo, jpk)                                )

  !-- barotropic time step --
  rdtbt = rdt / FLOAT(nn_baro)

  !-- constant definition --
  r1_grau = 1.e0 / (grav * rau0)

  !-- load mesh --
  nav_lon_t    = getvar(cf_tfil, 'nav_lon', 1, jpiglo, jpjglo)
  nav_lat_t    = getvar(cf_tfil, 'nav_lat', 1, jpiglo, jpjglo)
  deptht       = getvar1d(cf_tfil, cn_vdeptht , jpk)
  nav_lon_u    = getvar(cf_ufil, 'nav_lon', 1, jpiglo, jpjglo)
  nav_lat_u    = getvar(cf_ufil, 'nav_lat', 1, jpiglo, jpjglo)
  depthu       = getvar1d(cf_ufil, cn_vdepthu , jpk)
  nav_lon_v    = getvar(cf_vfil, 'nav_lon', 1, jpiglo, jpjglo)
  nav_lat_v    = getvar(cf_vfil, 'nav_lat', 1, jpiglo, jpjglo)
  depthv       = getvar1d(cf_vfil, cn_vdepthv , jpk)
  e1t(:,:)    = getvar(cf_mh  , cn_ve1t  , 1, jpiglo, jpjglo)
  e2t(:,:)    = getvar(cf_mh  , cn_ve2t  , 1, jpiglo, jpjglo)
  e1u(:,:)    = getvar(cf_mh  , cn_ve1u  , 1, jpiglo, jpjglo)
  e2u(:,:)    = getvar(cf_mh  , cn_ve2u  , 1, jpiglo, jpjglo)
  e1v(:,:)    = getvar(cf_mh  , cn_ve1v  , 1, jpiglo, jpjglo)
  e2v(:,:)    = getvar(cf_mh  , cn_ve2v  , 1, jpiglo, jpjglo)
  ff(:,:)     = getvar(cf_mh  , 'ff'     , 1, jpiglo, jpjglo)
  mbathy(:,:) = getvar(cf_mz, cn_mbathy,1, jpiglo, jpjglo )
  ht_0(:,:)   = getvar(cf_bathy, 'gdepw_0', 1, jpiglo, jpjglo )
  hu_0(:,:)   = getvar(cf_bathy, 'hu_0', 1, jpiglo, jpjglo )
  hv_0(:,:)   = getvar(cf_bathy, 'hv_0', 1, jpiglo, jpjglo )
  !- deepest we points (from domzgr.F90) -
  mbkt(:,:) = MAX( mbathy(:,:) , 1 )    ! bottom k-index of T-level (=1 over land)
  !                                     ! bottom k-index of W-level = mbkt+1
  DO jj = 1, jpjm1                      ! bottom k-index of u- (v-) level
     DO ji = 1, jpim1
        mbku(ji,jj) = MIN(  mbkt(ji+1,jj  ) , mbkt(ji,jj)  )
        mbkv(ji,jj) = MIN(  mbkt(ji  ,jj+1) , mbkt(ji,jj)  )
     END DO
  END DO
  !- horizontal cell surface (from domhgr.F90) -
  e12t    (:,:) = e1t(:,:) * e2t(:,:)
  r1_e12t (:,:) = 1._wp    / e12t(:,:)
  r1_e12u (:,:) = 1._wp    / ( e1u(:,:) * e2u(:,:) ) 
  r1_e12v (:,:) = 1._wp    / ( e1v(:,:) * e2v(:,:) )


  !-- load vert. mesh (at rest) and masks (dommsk.f90) --
  e3t_0(:,:,:) = getvar3d(cf_mz  , 'e3t_0' , jpiglo, jpjglo, jpk )
  e3u_0(:,:,:) = e3t_0(:,:,:)
  e3v_0(:,:,:) = e3t_0(:,:,:)
  DO jk = 1,jpk                         ! Computed as the minimum of neighbooring scale factors
     DO jj = 1, jpjm1
        DO ji = 1, jpim1   ! vector opt.
           e3u_0 (ji,jj,jk) = MIN( e3t_0(ji,jj,jk), e3t_0(ji+1,jj,jk) )
           e3v_0 (ji,jj,jk) = MIN( e3t_0(ji,jj,jk), e3t_0(ji,jj+1,jk) )
        END DO
     END DO
  END DO
  tmask(:,:,:) = getvar3d(cf_mask, 'tmask' , jpiglo, jpjglo, jpk )
  umask(:,:,:) = getvar3d(cf_mask, 'umask' , jpiglo, jpjglo, jpk )
  vmask(:,:,:) = getvar3d(cf_mask, 'vmask' , jpiglo, jpjglo, jpk )

  !-- Creat output netcdf files to fill in --
  CALL CreateOutput

  !DO jt = 2,jpt-1
  DO jt = 2,10
     PRINT *, '======= time-step = ', jt


     !-- Before variables and metrics for centered integration need to be recomputed for first time step --
     !-- compute them as now variables and swap right after --
     IF (jt == 2) THEN		! compute before variables and metrics
        !- recompute metric at before time step -
        sshn(:,:)     = getvar(cf_sshfil  , cn_sossheig, 1, jpiglo, jpjglo, ktime=jt-1 )
        !- atmospheric pressure loading, evap, precip, runoff
        msl(:,:)      = getvar(cf_mslfil  , 'msl', 1, jpiglo, jpjglo, ktime=1 )
        ssh_ib(:,:)   = -( msl - rn_pref ) * r1_grau
        emp(:,:)      = 0._wp
        rnf(:,:)      = 0._wp
        !
        e3t(:,:,:)    = 0._wp
        DO jk = 1, jpk
          e3t(:,:,jk)   = e3t_0(:,:,jk) * (1 + sshn/ht_0)
        END DO
        !- at u- and v- pts (domvvl.F90) -
        e3u(:,:,:) = e3u_0(:,:,:)
        e3v(:,:,:) = e3v_0(:,:,:)
        DO jk = 1, jpk
           DO jj = 1, jpjm1
              DO ji = 1, jpim1   ! vector opt.
                 e3u(ji,jj,jk) = e3u_0(ji,jj,jk) + 0.5_wp * umask(ji,jj,jk) * r1_e12u(ji,jj)                 &
                    &                       * (   e12t(ji  ,jj) * ( e3t(ji  ,jj,jk) - e3t_0(ji  ,jj,jk) )    &
                    &                           + e12t(ji+1,jj) * ( e3t(ji+1,jj,jk) - e3t_0(ji+1,jj,jk) ) )
                 e3v(ji,jj,jk) = e3v_0(ji,jj,jk) + 0.5_wp * vmask(ji,jj,jk) * r1_e12v(ji,jj)                 &
                    &                       * (   e12t(ji,jj  ) * ( e3t(ji,jj  ,jk) - e3t_0(ji,jj  ,jk) )    &
                    &                           + e12t(ji,jj+1) * ( e3t(ji,jj+1,jk) - e3t_0(ji,jj+1,jk) ) )
              END DO
           END DO
        END DO
        !- Local depth of the water column at t- points (from domvvl.F90) -
        ht(:,:) = 0._wp
        hu(:,:) = 0._wp
        hv(:,:) = 0._wp
        DO jk = 1, jpkm1
           ht(:,:) = ht(:,:) + e3t(:,:,jk) * tmask(:,:,jk)
           hu(:,:) = hu(:,:) + e3u(:,:,jk) * umask(:,:,jk)
           hv(:,:) = hv(:,:) + e3v(:,:,jk) * vmask(:,:,jk)
        END DO
        !- inverse of local depth -
        hur(:,:) = 1._wp / ( hu(:,:) + 1._wp - umask(:,:,1) ) * umask(:,:,1)      ! to avoid NaNs
        hvr(:,:) = 1._wp / ( hv(:,:) + 1._wp - vmask(:,:,1) ) * vmask(:,:,1)      ! to avoid NaNs
        !- Barotropic velocities -
        un(:,:,:)   = getvar3d(cf_ufil, cn_vozocrtx, jpiglo, jpjglo, jpk, ktime=jt-1   )
        vn(:,:,:)   = getvar3d(cf_vfil, cn_vomecrty, jpiglo, jpjglo, jpk, ktime=jt-1   )
        un_b(:,:) = 0._wp
        vn_b(:,:) = 0._wp
        DO jk = 1, jpk
           DO jj = 1, jpjglo
              DO ji = 1, jpiglo
                 un_b(ji,jj) = un_b(ji,jj) + e3u(ji,jj,jk) * un(ji,jj,jk) * umask(ji,jj,jk)
                 vn_b(ji,jj) = vn_b(ji,jj) + e3v(ji,jj,jk) * vn(ji,jj,jk) * vmask(ji,jj,jk)
              END DO
           END DO
        END DO
        un_b(:,:) = un_b(:,:) * hur(:,:)
        vn_b(:,:) = vn_b(:,:) * hvr(:,:)
        !-- wind stress --
        utau(:,:)     = getvar(cf_tauxfil , 'sozotaux' , 1, jpiglo, jpjglo, ktime=jt-1   )
        vtau(:,:)     = getvar(cf_tauyfil , 'sometauy' , 1, jpiglo, jpjglo, ktime=jt-1   )
     END IF   

     !-- SWAP --
     utau_b(:,:)  = utau(:,:)
     vtau_b(:,:)  = vtau(:,:)
     sshb(:,:)    = sshn(:,:)
     ssh_ibb(:,:) = ssh_ib(:,:)
     emp_b(:,:)   = emp(:,:)
     rnf_b(:,:)   = rnf(:,:)
     ub_b(:,:)    = un_b(:,:)
     vb_b(:,:)    = vn_b(:,:)
     ehu_b(:,:)   = hu(:,:)
     ehv_b(:,:)   = hv(:,:)
     ehur_b(:,:)  = hur(:,:)
     ehvr_b(:,:)  = hvr(:,:)
     ub(:,:,:)    = un(:,:,:)
     vb(:,:,:)    = vn(:,:,:)


     !-- Now variables and metrics --
     !-- recomputed vert. metric due to non-linear free surface (VVL) --
     ! from domvvl.F90 -->> e3t = e3t_0*(1+sshn/ht_0)
     PRINT *, '-- Recompute vert. mesh --'
     sshn(:,:)     = getvar(cf_sshfil  , cn_sossheig, 1, jpiglo, jpjglo, ktime=jt )
     !- atmospheric pressure loading, evap, precip, runoff
     msl(:,:)      = getvar(cf_mslfil  , 'msl', 1, jpiglo, jpjglo, ktime=1 )
     ssh_ib(:,:)   = -( msl - rn_pref ) * r1_grau
     emp(:,:)      = 0._wp
     rnf(:,:)      = 0._wp
     !
     e3t(:,:,:)    = 0._wp
     DO jk = 1, jpk
       e3t(:,:,jk)   = e3t_0(:,:,jk) * (1 + sshn/ht_0)
     END DO
     !- at u- and v- pts (domvvl.F90) -
     DO jk = 1, jpk
        DO jj = 1, jpjm1
           DO ji = 1, jpim1   ! vector opt.
              e3u(ji,jj,jk) = e3u_0(ji,jj,jk) + 0.5_wp * umask(ji,jj,jk) * r1_e12u(ji,jj)                 &
                 &                       * (   e12t(ji  ,jj) * ( e3t(ji  ,jj,jk) - e3t_0(ji  ,jj,jk) )    &
                 &                           + e12t(ji+1,jj) * ( e3t(ji+1,jj,jk) - e3t_0(ji+1,jj,jk) ) )
              e3v(ji,jj,jk) = e3v_0(ji,jj,jk) + 0.5_wp * vmask(ji,jj,jk) * r1_e12v(ji,jj)                 &
                 &                       * (   e12t(ji,jj  ) * ( e3t(ji,jj  ,jk) - e3t_0(ji,jj  ,jk) )    &
                 &                           + e12t(ji,jj+1) * ( e3t(ji,jj+1,jk) - e3t_0(ji,jj+1,jk) ) )
           END DO
        END DO
     END DO
     ! Local depth of the water column at t- points (from domvvl.F90)
     ht(:,:) = 0._wp
     hu(:,:) = 0._wp
     hv(:,:) = 0._wp
     DO jk = 1, jpkm1
        ht(:,:) = ht(:,:) + e3t(:,:,jk) * tmask(:,:,jk)
        hu(:,:) = hu(:,:) + e3u(:,:,jk) * umask(:,:,jk)
        hv(:,:) = hv(:,:) + e3v(:,:,jk) * vmask(:,:,jk)
     END DO
     ! inverse of local depth 
     hur(:,:) = 1._wp / ( hu(:,:) + 1._wp - umask(:,:,1) ) * umask(:,:,1)      ! to avoid NaNs
     hvr(:,:) = 1._wp / ( hv(:,:) + 1._wp - vmask(:,:,1) ) * vmask(:,:,1)      ! to avoid NaNs
     !- Barotropic velocities -
     un(:,:,:)   = getvar3d(cf_ufil, cn_vozocrtx, jpiglo, jpjglo, jpk, ktime=jt   )
     vn(:,:,:)   = getvar3d(cf_vfil, cn_vomecrty, jpiglo, jpjglo, jpk, ktime=jt   )
     un_b(:,:) = 0._wp
     vn_b(:,:) = 0._wp
     DO jk = 1, jpkm1
        DO jj = 1, jpjglo
           DO ji = 1, jpiglo
              un_b(ji,jj) = un_b(ji,jj) + e3u(ji,jj,jk) * un(ji,jj,jk) * umask(ji,jj,jk)
              vn_b(ji,jj) = vn_b(ji,jj) + e3v(ji,jj,jk) * vn(ji,jj,jk) * vmask(ji,jj,jk)
           END DO
        END DO
     END DO
     un_b(:,:) = un_b(:,:) * hur(:,:)
     vn_b(:,:) = vn_b(:,:) * hvr(:,:)
     !-- wind stress --
     utau(:,:)     = getvar(cf_tauxfil , 'sozotaux' , 1, jpiglo, jpjglo, ktime=jt   )
     vtau(:,:)     = getvar(cf_tauyfil , 'sometauy' , 1, jpiglo, jpjglo, ktime=jt   )

     !-- recompute trends; this is an approximation since (ua,va) do not refer to TOTAL trend normally, --
     !-- but only trends du to other dynamical terms (advection, vorticity and pressure gradients) --
     ua(:,:,:)   = getvar3d(cf_ufil, cn_vozocrtx, jpiglo, jpjglo, jpk, ktime=jt+1   )
     va(:,:,:)   = getvar3d(cf_vfil, cn_vomecrty, jpiglo, jpjglo, jpk, ktime=jt+1   )
     ua(:,:,:) = ( ua(:,:,:)-ub(:,:,:) ) / (2*rdt)
     va(:,:,:) = ( va(:,:,:)-vb(:,:,:) ) / (2*rdt)
     !- temporary save for trend computation -
     spgu(:,:,:) = ua(:,:,:)
     spgv(:,:,:) = va(:,:,:)

     !-- apply 2d map (if provided) for enhanced bfr --
     bfrcoef2d(:,:) = 0._wp
     bfrcoef2d(:,:) = rn_bfri2 * ( 1 + rn_bfrien * bfrcoef2d(:,:) )

     !-- compute bfrua/bfrva --
     PRINT *, '-- Compute bottom friction coef. --'
     CALL zdf_bfr( jt )

     !-- compute surface pressure gradient correction --
     PRINT *, '-- Compute surface pressure gradient --'
     ssh_bt(:,:,:) = 0._wp
     CALL dyn_spg_ts( jt )

     !-- compute trends associated with surface pressure correction --
     spgu(:,:,:) = ua(:,:,:) - spgu(:,:,:)
     spgv(:,:,:) = va(:,:,:) - spgv(:,:,:)

     !-- Construct KE --
     spg_ke(:,:,:) = 0._wp
     CALL trd_ken( spgu, spgv, spg_ke)

     !-- save --
     DO jk = 1, jpk
        ierr = putvar(ncout_u, id_varout_u(1), spgu(:,:,jk), jk, jpiglo, jpjglo, ktime=jt )
        ierr = putvar(ncout_v, id_varout_v(1), spgv(:,:,jk), jk, jpiglo, jpjglo, ktime=jt )
        ierr = putvar(ncout_ke, id_varout_ke(1), spg_ke(:,:,jk), jk, jpiglo, jpjglo, ktime=jt )
        ierr = putvar(ncout_ssh, id_varout_ssh(1), ssh_bt(:,:,jk), jk, jpiglo, jpjglo, ktime=jt )
     END DO
     !ierr = putvar(ncout_u, id_varout_u(1), spgu, 1, jpiglo, jpjglo, ktime=jt )
     !ierr = putvar(ncout_v, id_varout_v(1), spgv, 1, jpiglo, jpjglo, ktime=jt )
     !ierr = putvar(ncout_u, id_varout_u(1), tmpu, 1, jpiglo, jpjglo, ktime=jt )
     !ierr = putvar(ncout_v, id_varout_v(1), tmpv, 1, jpiglo, jpjglo, ktime=jt )

  ENDDO		!jt-loop

  ierr = closeout(ncout_u)
  ierr = closeout(ncout_v)
  ierr = closeout(ncout_ke)
  ierr = closeout(ncout_ssh)

CONTAINS

   SUBROUTINE dyn_spg_ts( kt )
!!----------------------------------------------------------------------
!!
!! ** Purpose :
!!      -Compute the now trend due to the explicit time stepping
!!      of the quasi-linear barotropic system.
!!
!! ** Method  :
!!      Barotropic variables are advanced from internal time steps
!!      "n"   to "n+1" if ln_bt_fw=T
!!      or from
!!      "n-1" to "n+1" if ln_bt_fw=F
!!      thanks to a generalized forward-backward time stepping (see ref. below).
!!
!! ** Action :
!!      -Update the filtered free surface at step "n+1"      : ssha
!!      -Update filtered barotropic velocities at step "n+1" : ua_b, va_b
!!      -Compute barotropic advective velocities at step "n" : un_adv, vn_adv
!!      These are used to advect tracers and are compliant with discrete
!!      continuity equation taken at the baroclinic time steps. This
!!      ensures tracers conservation.
!!      -Update 3d trend (ua, va) with barotropic component.
!!
!! References : Shchepetkin, A.F. and J.C. McWilliams, 2005:
!!              The regional oceanic modeling system (ROMS):
!!              a split-explicit, free-surface,
!!              topography-following-coordinate oceanic model.
!!              Ocean Modelling, 9, 347-404.
!!---------------------------------------------------------------------
!
      INTEGER, INTENT(in)  ::   kt   ! ocean time-step index
!
      LOGICAL  ::   ll_fw_start                 ! if true, forward integration
      LOGICAL  ::   ll_init                         ! if true, special startup of 2d equations
      INTEGER  ::   ji, jj, jk, jn              ! dummy loop indices
      INTEGER  ::   ikbu, ikbv, noffset         ! local integers
      !-
      REAL(wp) ::   zraur, z1_2dt_b, z2dt_bf    ! local scalars
      REAL(wp) ::   zx1, zy1, zx2, zy2         !   -      -
      REAL(wp) ::   z1_12, z1_8, z1_4, z1_2       !   -      -
      REAL(wp) ::   zu_spg, zv_spg             !   -      -
      REAL(wp) ::   zhura, zhvra                     !   -      -
      REAL(wp) ::   za0, za1, za2, za3                 !   -      -
!
      REAL(wp), POINTER, DIMENSION(:,:) :: zun_e, zvn_e, zsshp2_e
      REAL(wp), POINTER, DIMENSION(:,:) :: zu_trd, zv_trd, zu_frc, zv_frc, zssh_frc
      REAL(wp), POINTER, DIMENSION(:,:) :: zu_sum, zv_sum, zwx, zwy, zhdiv
      REAL(wp), POINTER, DIMENSION(:,:) :: zhup2_e, zhvp2_e, zhust_e, zhvst_e
      REAL(wp), POINTER, DIMENSION(:,:) :: zsshu_a, zsshv_a
      REAL(wp), POINTER, DIMENSION(:,:) :: zhf
      !QJ: save u,v barotropic trends to force boundaries and avoid instabilities
      REAL(wp), POINTER, DIMENSION(:,:) :: zu_frc2, zv_frc2
!

!!----------------------------------------------------------------------
!
!                                         !* Allocate temporary arrays
      ALLOCATE( zsshp2_e(jpiglo, jpjglo), zhdiv(jpiglo, jpjglo) )
      ALLOCATE( zu_trd(jpiglo, jpjglo), zv_trd(jpiglo, jpjglo), zun_e(jpiglo, jpjglo), zvn_e(jpiglo, jpjglo)  )
      ALLOCATE( zwx(jpiglo, jpjglo), zwy(jpiglo, jpjglo), zu_sum(jpiglo, jpjglo), zv_sum(jpiglo, jpjglo) )
      ALLOCATE( zssh_frc(jpiglo, jpjglo), zu_frc(jpiglo, jpjglo), zv_frc(jpiglo, jpjglo))
      ALLOCATE( zhup2_e(jpiglo, jpjglo), zhvp2_e(jpiglo, jpjglo), zhust_e(jpiglo, jpjglo), zhvst_e(jpiglo, jpjglo))
      ALLOCATE( zsshu_a(jpiglo, jpjglo), zsshv_a(jpiglo, jpjglo)                                   )
      ALLOCATE( zhf(jpiglo, jpjglo) )
      !QJ:
      ALLOCATE( zu_frc2(jpiglo, jpjglo), zv_frc2(jpiglo, jpjglo))
     
!
!                                         !* Local constant initialization
      z1_12 = 1._wp / 12._wp
      z1_8  = 0.125_wp
      z1_4  = 0.25_wp
      z1_2  = 0.5_wp
      zraur = 1._wp / rau0

!      IF( kt == nit000 .AND. neuler == 0 ) THEN         ! reciprocal of baroclinic time step
!        z2dt_bf = rdt
!      ELSE
        z2dt_bf = 2.0_wp * rdt
!      ENDIF
      z1_2dt_b = 1.0_wp / z2dt_bf
!
!      ll_init = ln_bt_av                                ! if no time averaging, then no specific restart
!      ll_fw_start = .FALSE.
!
! time offset in steps for bdy data update
!      IF (.NOT.ln_bt_fw) THEN ; noffset=-nn_baro ; ELSE ;  noffset = 0 ; ENDIF
!
!      IF( kt == nit000 ) THEN                   !* initialisation
!         ....
!
! Set averaging weights and cycle length:
         CALL ts_wgt(icycle, wgtbtp1, wgtbtp2)
!
!
!      ENDIF

!
! Set arrays to remove/compute coriolis trend.
! Do it once at kt=nit000 if volume is fixed, else at each long time step.
! Note that these arrays are also used during barotropic loop. These are however frozen
! although they should be updated in the variable volume case. Not a big approximation.
! To remove this approximation, copy lines below inside barotropic loop
! and update depths at T-F points (ht and zhf resp.) at each barotropic time step
!
!      IF ( kt == nit000 .OR. lk_vvl ) THEN                    ! QJ: lk_vvl=true, thus do it at each time step
!         IF ( ln_dynvor_een_old ) THEN
!            ....
!         ELSE IF ( ln_dynvor_een ) THEN
            DO jj = 1, jpjm1
               DO ji = 1, jpim1
                  zwz(ji,jj) =   ( ht(ji  ,jj+1) + ht(ji+1,jj+1) +                     &
                        &          ht(ji  ,jj  ) + ht(ji+1,jj  )   )                   &
                        &      / ( MAX( 1.0_wp, tmask(ji  ,jj+1, 1) + tmask(ji+1,jj+1, 1) +    &
                        &                       tmask(ji  ,jj  , 1) + tmask(ji+1,jj  , 1) ) )
                  IF( zwz(ji,jj) /= 0._wp )   zwz(ji,jj) = 1._wp / zwz(ji,jj)
               END DO
            END DO
            !CALL lbc_lnk( zwz, 'F', 1._wp )
            zwz(:,:) = ff(:,:) * zwz(:,:)

            ftne(1,:) = 0._wp ; ftnw(1,:) = 0._wp ; ftse(1,:) = 0._wp ; ftsw(1,:) = 0._wp
            DO jj = 2, jpjglo
               DO ji = 2, jpiglo   ! vector opt.
                  ftne(ji,jj) = zwz(ji-1,jj  ) + zwz(ji  ,jj  ) + zwz(ji  ,jj-1)
                  ftnw(ji,jj) = zwz(ji-1,jj-1) + zwz(ji-1,jj  ) + zwz(ji  ,jj  )
                  ftse(ji,jj) = zwz(ji  ,jj  ) + zwz(ji  ,jj-1) + zwz(ji-1,jj-1)
                  ftsw(ji,jj) = zwz(ji  ,jj-1) + zwz(ji-1,jj-1) + zwz(ji-1,jj  )
               END DO
            END DO
!         ELSE
!            ....
!         ENDIF
!      ENDIF
!
! If forward start at previous time step, and centered integration,
! then update averaging weights:
!      IF ((.NOT.ln_bt_fw).AND.((neuler==0).AND.(kt==nit000+1))) THEN
!         ....
!      ENDIF


! -----------------------------------------------------------------------------
!  Phase 1 : Coupling between general trend and barotropic estimates (1st step)
! -----------------------------------------------------------------------------
!
!
!                                   !* e3*d/dt(Ua) (Vertically integrated)
!                                   ! --------------------------------------------------
      zu_frc(:,:) = 0._wp
      zv_frc(:,:) = 0._wp
!
! QJ: additional initialization to avoid boundary instabilities
      zu_trd(:,:)   = 0._wp
      zv_trd(:,:)   = 0._wp
      zhdiv(:,:)    = 0._wp
      zsshu_a(:,:)  = 0._wp
      zsshv_a(:,:)  = 0._wp
      zhust_e(:,:)  = 0._wp
      zhvst_e(:,:)  = 0._wp

      DO jk = 1, jpkm1
         zu_frc(:,:) = zu_frc(:,:) + e3u(:,:,jk) * ua(:,:,jk) * umask(:,:,jk)
         zv_frc(:,:) = zv_frc(:,:) + e3v(:,:,jk) * va(:,:,jk) * vmask(:,:,jk)
      END DO
!
      zu_frc(:,:) = zu_frc(:,:) * hur(:,:)
      zv_frc(:,:) = zv_frc(:,:) * hvr(:,:)
! QJ:
      zu_frc2(:,:) = zu_frc(:,:)
      zv_frc2(:,:) = zv_frc(:,:)


                                   !* baroclinic momentum trend (remove the vertical mean trend)
      DO jk = 1, jpkm1                    ! -----------------------------------------------------------
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               ua(ji,jj,jk) = ua(ji,jj,jk) - zu_frc(ji,jj) * umask(ji,jj,jk)
               va(ji,jj,jk) = va(ji,jj,jk) - zv_frc(ji,jj) * vmask(ji,jj,jk)
            END DO
         END DO
      END DO
                                   !* barotropic Coriolis trends (vorticity scheme dependent)
                                   ! --------------------------------------------------------
      zwx(:,:) = un_b(:,:) * hu(:,:) * e2u(:,:)        ! now fluxes
      zwy(:,:) = vn_b(:,:) * hv(:,:) * e1v(:,:)
!
!      IF( ln_dynvor_ene .OR. ln_dynvor_mix ) THEN      ! energy conserving or mixed scheme
!         ....
!      ELSEIF ( ln_dynvor_ens ) THEN                    ! enstrophy conserving scheme
!         ....
!      ELSEIF ( ln_dynvor_een .or. ln_dynvor_een_old ) THEN  ! enstrophy and energy conserving scheme
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               zu_trd(ji,jj) = + z1_12 / e1u(ji,jj) * (  ftne(ji,jj  ) * zwy(ji  ,jj  ) &
                &                                      + ftnw(ji+1,jj) * zwy(ji+1,jj  ) &
                &                                      + ftse(ji,jj  ) * zwy(ji  ,jj-1) &
                &                                      + ftsw(ji+1,jj) * zwy(ji+1,jj-1) )
               zv_trd(ji,jj) = - z1_12 / e2v(ji,jj) * (  ftsw(ji,jj+1) * zwx(ji-1,jj+1) &
                &                                      + ftse(ji,jj+1) * zwx(ji  ,jj+1) &
                &                                      + ftnw(ji,jj  ) * zwx(ji-1,jj  ) &
                &                                      + ftne(ji,jj  ) * zwx(ji  ,jj  ) )
            END DO
         END DO
!
!      ENDIF
!
!                                   !* Right-Hand-Side of the barotropic momentum equation
!                                   ! ----------------------------------------------------
!      IF( lk_vvl ) THEN                         ! Variable volume : remove surface pressure gradient
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               zu_trd(ji,jj) = zu_trd(ji,jj) - grav * (  sshn(ji+1,jj  ) - sshn(ji  ,jj  )  ) / e1u(ji,jj)
               zv_trd(ji,jj) = zv_trd(ji,jj) - grav * (  sshn(ji  ,jj+1) - sshn(ji  ,jj  )  ) / e2v(ji,jj)
            END DO
         END DO
!      ENDIF

      DO jj = 2, jpjm1                          ! Remove coriolis term (and possibly spg) from barotropic trend
         DO ji = 2, jpim1
             zu_frc(ji,jj) = zu_frc(ji,jj) - zu_trd(ji,jj) * umask(ji,jj,1)
             zv_frc(ji,jj) = zv_frc(ji,jj) - zv_trd(ji,jj) * vmask(ji,jj,1)
          END DO
      END DO
!
!                                               ! Add bottom stress contribution from baroclinic velocities:
!      IF (ln_bt_fw) THEN
!          ....
!      ELSE
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               ikbu = mbku(ji,jj)
               ikbv = mbkv(ji,jj)
               zwx(ji,jj) = ub(ji,jj,ikbu) - ub_b(ji,jj) ! BEFORE bottom baroclinic velocities
               zwy(ji,jj) = vb(ji,jj,ikbv) - vb_b(ji,jj)
            END DO
         END DO
!      ENDIF
!
! Note that the "unclipped" bottom friction parameter is used even with explicit drag
      zu_frc(:,:) = zu_frc(:,:) + hur(:,:) * bfrua(:,:) * zwx(:,:)
      zv_frc(:,:) = zv_frc(:,:) + hvr(:,:) * bfrva(:,:) * zwy(:,:)
!
!      IF (ln_bt_fw) THEN                        ! Add wind forcing
!         .....
!      ELSE
         zu_frc(:,:) =  zu_frc(:,:) + zraur * z1_2 * ( utau_b(:,:) + utau(:,:) ) * hur(:,:)
         zv_frc(:,:) =  zv_frc(:,:) + zraur * z1_2 * ( vtau_b(:,:) + vtau(:,:) ) * hvr(:,:)
!      ENDIF
!
!      IF ( ln_apr_dyn ) THEN                    ! Add atm pressure forcing
!         IF (ln_bt_fw) THEN
!            ....
!         ELSE
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  zu_spg =  grav * z1_2 * (  ssh_ib (ji+1,jj  ) - ssh_ib (ji,jj)    &
                      &                    + ssh_ibb(ji+1,jj  ) - ssh_ibb(ji,jj)  ) /e1u(ji,jj)
                  zv_spg =  grav * z1_2 * (  ssh_ib (ji  ,jj+1) - ssh_ib (ji,jj)    &
                      &                    + ssh_ibb(ji  ,jj+1) - ssh_ibb(ji,jj)  ) /e2v(ji,jj)
                  zu_frc(ji,jj) = zu_frc(ji,jj) + zu_spg
                  zv_frc(ji,jj) = zv_frc(ji,jj) + zv_spg
               END DO
            END DO
!         ENDIF
!      ENDIF
!                                   !* Right-Hand-Side of the barotropic ssh equation
!                                   ! -----------------------------------------------
!                                         ! Surface net water flux and rivers
!      IF (ln_bt_fw) THEN
!         ....
!      ELSE
         zssh_frc(:,:) = zraur * z1_2 * (  emp(:,:) + emp_b(:,:) - rnf(:,:) - rnf_b(:,:)    )
         !       &                        + fwfisf(:,:) + fwfisf_b(:,:)                     )
!      ENDIF

!                                   !* Fill boundary data arrays for AGRIF
!                                   ! ------------------------------------

!
! -----------------------------------------------------------------------
!  Phase 2 : Integration of the barotropic equations
! -----------------------------------------------------------------------
!
!                                             ! ==================== !
!                                             !    Initialisations   !
!                                             ! ==================== !
! Initialize barotropic variables:
!      IF( kt == 0 )THEN
         sshbb_e(:,:) = 0._wp
         ubb_e  (:,:) = 0._wp
         vbb_e  (:,:) = 0._wp
         sshb_e (:,:) = 0._wp
         ub_e   (:,:) = 0._wp
         vb_e   (:,:) = 0._wp
!      ENDIF
!
!      IF (ln_bt_fw) THEN                  ! FORWARD integration: start from NOW fields
!         ...
!      ELSE                                ! CENTRED integration: start from BEFORE fields
         sshn_e(:,:) = sshb (:,:)
         zun_e (:,:) = ub_b (:,:)
         zvn_e (:,:) = vb_b (:,:)
!
         hu_e  (:,:) = ehu_b(:,:)
         hv_e  (:,:) = ehv_b(:,:)
         hur_e (:,:) = ehur_b(:,:)
         hvr_e (:,:) = ehvr_b(:,:)
!      ENDIF
!
!
!
! Initialize sums:
      ua_b  (:,:) = 0._wp       ! After barotropic velocities (or transport if flux form)
      va_b  (:,:) = 0._wp
      ssha  (:,:) = 0._wp       ! Sum for after averaged sea level
      zu_sum(:,:) = 0._wp       ! Sum for now transport issued from ts loop
      zv_sum(:,:) = 0._wp
!                                             ! ==================== !
      DO jn = 1, icycle                             !  sub-time-step loop  !
!                                          ! ==================== !
!                                                !* Update the forcing (BDY and tides)
!                                                !  ------------------
! Update only tidal forcing at open boundaries

!         IF ( lk_bdy .AND. lk_tide ) CALL bdy_dta_tides( kt, kit=jn, time_offset=(noffset+1) )
!         IF ( ln_tide_pot .AND. lk_tide ) CALL upd_tide( kt, kit=jn, time_offset=noffset )

!
! Set extrapolation coefficients for predictor step:
!         IF ((jn<3).AND.ll_init) THEN      ! Forward
!            ...
!         ELSE                              ! AB3-AM4 Coefficients: bet=0.281105
           za1 =  1.781105_wp              ! za1 =   3/2 +   bet
           za2 = -1.06221_wp               ! za2 = -(1/2 + 2*bet)
           za3 =  0.281105_wp              ! za3 = bet
!         ENDIF

! Extrapolate barotropic velocities at step jit+0.5:
         ua_e(:,:) = za1 * zun_e(:,:) + za2 * ub_e(:,:) + za3 * ubb_e(:,:)
         va_e(:,:) = za1 * zvn_e(:,:) + za2 * vb_e(:,:) + za3 * vbb_e(:,:)

!         IF( lk_vvl ) THEN                                !* Update ocean depth (variable volume case only)
!                                             !  ------------------
! Extrapolate Sea Level at step jit+0.5:
            zsshp2_e(:,:) = za1 * sshn_e(:,:)  + za2 * sshb_e(:,:) + za3 * sshbb_e(:,:)
!
            DO jj = 2, jpjm1                                    ! Sea Surface Height at u- & v-points
               DO ji = 2, jpim1   ! Vector opt.
                  zwx(ji,jj) = z1_2 * umask(ji,jj,1)  * r1_e12u(ji,jj)     &
                     &              * ( e12t(ji  ,jj) * zsshp2_e(ji  ,jj)  &
                     &              +   e12t(ji+1,jj) * zsshp2_e(ji+1,jj) )
                  zwy(ji,jj) = z1_2 * vmask(ji,jj,1)  * r1_e12v(ji,jj)     &
                     &              * ( e12t(ji,jj  ) * zsshp2_e(ji,jj  )  &
                     &              +   e12t(ji,jj+1) * zsshp2_e(ji,jj+1) )
               END DO
            END DO
!            CALL lbc_lnk_multi( zwx, 'U', 1._wp, zwy, 'V', 1._wp )
!
            zhup2_e (:,:) = hu_0(:,:) + zwx(:,:)                ! Ocean depth at U- and V-points
            zhvp2_e (:,:) = hv_0(:,:) + zwy(:,:)
!         ELSE
!            ....
!         ENDIF
!                                                !* after ssh
!                                                !  -----------
! One should enforce volume conservation at open boundaries here
! considering fluxes below:
!
         zwx(:,:) = e2u(:,:) * ua_e(:,:) * zhup2_e(:,:)         ! fluxes at jn+0.5
         zwy(:,:) = e1v(:,:) * va_e(:,:) * zhvp2_e(:,:)
!

!
! Sum over sub-time-steps to compute advective velocities
         za2 = wgtbtp2(jn)
         zu_sum  (:,:) = zu_sum  (:,:) + za2 * zwx  (:,:) / e2u  (:,:)
         zv_sum  (:,:) = zv_sum  (:,:) + za2 * zwy  (:,:) / e1v  (:,:)
!
! Set next sea level:
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               zhdiv(ji,jj) = (   zwx(ji,jj) - zwx(ji-1,jj)   &
                  &             + zwy(ji,jj) - zwy(ji,jj-1)   ) * r1_e12t(ji,jj)
            END DO
         END DO
         ssha_e(:,:) = (  sshn_e(:,:) - rdtbt * ( zssh_frc(:,:) + zhdiv(:,:) )  ) * tmask(:,:,1)
!         CALL lbc_lnk( ssha_e, 'T',  1._wp )


! Duplicate sea level across open boundaries (this is only cosmetic if lk_vvl=.false.)
!         IF (lk_bdy) CALL bdy_ssh( ssha_e )


!
! Sea Surface Height at u-,v-points (vvl case only)
!         IF ( lk_vvl ) THEN
            DO jj = 2, jpjm1
               DO ji = 2, jpim1      ! NO Vector Opt.
                  zsshu_a(ji,jj) = z1_2 * umask(ji,jj,1)  * r1_e12u(ji,jj)  &
                     &              * ( e12t(ji  ,jj  ) * ssha_e(ji  ,jj  ) &
                     &              +   e12t(ji+1,jj  ) * ssha_e(ji+1,jj  ) )
                  zsshv_a(ji,jj) = z1_2 * vmask(ji,jj,1)  * r1_e12v(ji,jj)  &
                     &              * ( e12t(ji  ,jj  ) * ssha_e(ji  ,jj  ) &
                     &              +   e12t(ji  ,jj+1) * ssha_e(ji  ,jj+1) )
               END DO
            END DO
!            CALL lbc_lnk_multi( zsshu_a, 'U', 1._wp, zsshv_a, 'V', 1._wp )
!         ENDIF
!
! Half-step back interpolation of SSH for surface pressure computation:
!----------------------------------------------------------------------
!         IF ((jn==1).AND.ll_init) THEN
!            ....
!         ELSEIF ((jn==2).AND.ll_init) THEN  ! AB2-AM3 Coefficients; bet=0 ; gam=-1/6 ; eps=1/12
!            ....
!         ELSE                               ! AB3-AM4 Coefficients; bet=0.281105 ; eps=0.013 ; gam=0.0880
           za0=0.614_wp                     ! za0 = 1/2 +   gam + 2*eps
           za1=0.285_wp                     ! za1 = 1/2 - 2*gam - 3*eps
           za2=0.088_wp                     ! za2 = gam
           za3=0.013_wp                     ! za3 = eps
!         ENDIF

         zsshp2_e(:,:) = za0 *  ssha_e(:,:) + za1 *  sshn_e (:,:) &
          &            + za2 *  sshb_e(:,:) + za3 *  sshbb_e(:,:)

!
! Compute associated depths at U and V points:
!         IF ( lk_vvl.AND.(.NOT.ln_dynadv_vec) ) THEN
!
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  zx1 = z1_2 * umask(ji  ,jj,1) *  r1_e12u(ji  ,jj)    &
                     &      * ( e12t(ji  ,jj  ) * zsshp2_e(ji  ,jj)    &
                     &      +   e12t(ji+1,jj  ) * zsshp2_e(ji+1,jj  ) )
                  zy1 = z1_2 * vmask(ji  ,jj,1) *  r1_e12v(ji  ,jj  )  &
                     &       * ( e12t(ji ,jj  ) * zsshp2_e(ji  ,jj  )  &
                     &       +   e12t(ji ,jj+1) * zsshp2_e(ji  ,jj+1) )
                  zhust_e(ji,jj) = hu_0(ji,jj) + zx1
                  zhvst_e(ji,jj) = hv_0(ji,jj) + zy1
               END DO
            END DO
!         ENDIF
!
! Add Coriolis trend:
! zwz array below or triads normally depend on sea level with 1 and should be updated
! at each time step. We however keep them constant here for optimization.
! Recall that zwx and zwy arrays hold fluxes at this stage:
! zwx(:,:) = e2u(:,:) * ua_e(:,:) * zhup2_e(:,:)   ! fluxes at jn+0.5
! zwy(:,:) = e1v(:,:) * va_e(:,:) * zhvp2_e(:,:)
!
!         IF( ln_dynvor_ene .OR. ln_dynvor_mix ) THEN      !==  energy conserving or mixed scheme  ==!
!            ....
!         ELSEIF ( ln_dynvor_ens ) THEN                    !==  enstrophy conserving scheme  ==!
!            ....
!         ELSEIF ( ln_dynvor_een .or. ln_dynvor_een_old ) THEN !==  energy and enstrophy conserving scheme  ==!
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  zu_trd(ji,jj) = + z1_12 / e1u(ji,jj) * (  ftne(ji,jj  ) * zwy(ji  ,jj  ) &
                     &                                    + ftnw(ji+1,jj) * zwy(ji+1,jj  ) &
                     &                                    + ftse(ji,jj  ) * zwy(ji  ,jj-1) &
                     &                                    + ftsw(ji+1,jj) * zwy(ji+1,jj-1) )
                  zv_trd(ji,jj) = - z1_12 / e2v(ji,jj) * (  ftsw(ji,jj+1) * zwx(ji-1,jj+1) &
                     &                                    + ftse(ji,jj+1) * zwx(ji  ,jj+1) &
                     &                                    + ftnw(ji,jj  ) * zwx(ji-1,jj  ) &
                     &                                    + ftne(ji,jj  ) * zwx(ji  ,jj  ) )
               END DO
            END DO
!
!         ENDIF
!
! Add tidal astronomical forcing if defined
!         IF ( lk_tide.AND.ln_tide_pot ) THEN
!            ....
!         ENDIF
!
! Add bottom stresses:
         zu_trd(:,:) = zu_trd(:,:) + bfrua(:,:) * zun_e(:,:) * hur_e(:,:)
         zv_trd(:,:) = zv_trd(:,:) + bfrva(:,:) * zvn_e(:,:) * hvr_e(:,:)
!
! Surface pressure trend:
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
! Add surface pressure gradient
               zu_spg = - grav * ( zsshp2_e(ji+1,jj) - zsshp2_e(ji,jj) ) / e1u(ji,jj)
               zv_spg = - grav * ( zsshp2_e(ji,jj+1) - zsshp2_e(ji,jj) ) / e2v(ji,jj)
               zwx(ji,jj) = zu_spg
               zwy(ji,jj) = zv_spg
            END DO
         END DO
!
! Set next velocities:
!         IF( ln_dynadv_vec .OR. (.NOT. lk_vvl) ) THEN   ! Vector form
!            ....
!         ELSE                                           ! Flux form
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.

                  zhura = umask(ji,jj,1)/(hu_0(ji,jj) + zsshu_a(ji,jj) + 1._wp - umask(ji,jj,1))
                  zhvra = vmask(ji,jj,1)/(hv_0(ji,jj) + zsshv_a(ji,jj) + 1._wp - vmask(ji,jj,1))

                  ua_e(ji,jj) = (                hu_e(ji,jj)  *  zun_e(ji,jj)   &
                            &     + rdtbt * ( zhust_e(ji,jj)  *    zwx(ji,jj)   &
                            &               + zhup2_e(ji,jj)  * zu_trd(ji,jj)   &
                            &               +      hu(ji,jj)  * zu_frc(ji,jj) ) &
                            &   ) * zhura

                  va_e(ji,jj) = (                hv_e(ji,jj)  *  zvn_e(ji,jj)   &
                            &     + rdtbt * ( zhvst_e(ji,jj)  *    zwy(ji,jj)   &
                            &               + zhvp2_e(ji,jj)  * zv_trd(ji,jj)   &
                            &               +      hv(ji,jj)  * zv_frc(ji,jj) ) &
                            &   ) * zhvra
               END DO
            END DO
!         ENDIF
!
!         IF( lk_vvl ) THEN                             !* Update ocean depth (variable volume case only)
!                                          !  ----------------------------------------------
            hu_e (:,:) = hu_0(:,:) + zsshu_a(:,:)
            hv_e (:,:) = hv_0(:,:) + zsshv_a(:,:)
            hur_e(:,:) = umask(:,:,1) / ( hu_e(:,:) + 1._wp - umask(:,:,1) )
            hvr_e(:,:) = vmask(:,:,1) / ( hv_e(:,:) + 1._wp - vmask(:,:,1) )
!
!         ENDIF
!                                                 !* domain lateral boundary
!                                                 !  -----------------------
!
!         CALL lbc_lnk_multi( ua_e, 'U', -1._wp, va_e , 'V', -1._wp )


! open boundaries
!         IF( lk_bdy ) CALL bdy_dyn2d( jn, ua_e, va_e, zun_e, zvn_e, hur_e, hvr_e, ssha_e )
! QJ: force boundary conditions to first estimates of barotropic trends (zu_frc2, zv_frc2)
!     and now ssh associated with the baroclinic loop
!        ua_e(1     ,:     )   = zu_frc2(1     ,:     )
!        ua_e(jpiglo,:     )   = zu_frc2(jpiglo,:     )
!        ua_e(:     ,1     )   = zu_frc2(:     ,1     )
!        ua_e(:     ,jpjglo)   = zu_frc2(:     ,jpjglo)
!        va_e(1     ,:     )   = zv_frc2(1     ,:     )
!        va_e(jpiglo,:     )   = zv_frc2(jpiglo,:     )
!        va_e(:     ,1     )   = zv_frc2(:     ,1     )
!        va_e(:     ,jpjglo)   = zv_frc2(:     ,jpjglo)
!        ssha_e(1     ,:     ) = sshn(1     ,:     )
!        ssha_e(2     ,:     ) = sshn(2     ,:     )
!        ssha_e(jpiglo,:     ) = sshn(jpiglo,:     )
!        ssha_e(jpiglo-1,:     ) = sshn(jpiglo-1,:     )
!        ssha_e(:     ,1     ) = sshn(:     ,1     )
!        ssha_e(:     ,2     ) = sshn(:     ,2     )
!        ssha_e(:     ,jpjglo) = sshn(:     ,jpjglo)
!        ssha_e(:     ,jpjglo-1) = sshn(:     ,jpjglo-1)


!                                             !* Swap
!                                             !  ----
         ubb_e  (:,:) = ub_e  (:,:)
         ub_e   (:,:) = zun_e (:,:)
         zun_e  (:,:) = ua_e  (:,:)
!
         vbb_e  (:,:) = vb_e  (:,:)
         vb_e   (:,:) = zvn_e (:,:)
         zvn_e  (:,:) = va_e  (:,:)
!
         sshbb_e(:,:) = sshb_e(:,:)
         sshb_e (:,:) = sshn_e(:,:)
         sshn_e (:,:) = ssha_e(:,:)

!                                             !* Sum over whole bt loop
!                                             !  ----------------------
         za1 = wgtbtp1(jn)
!         IF (( ln_dynadv_vec ).OR. (.NOT. lk_vvl)) THEN    ! Sum velocities
!            ....
!        ELSE                                              ! Sum transports
            ua_b  (:,:) = ua_b  (:,:) + za1 * ua_e  (:,:) * hu_e (:,:)
            va_b  (:,:) = va_b  (:,:) + za1 * va_e  (:,:) * hv_e (:,:)
!        ENDIF
!                                                  ! Sum sea level
         ssha(:,:) = ssha(:,:) + za1 * ssha_e(:,:)
!        QJ: debug
         ssh_bt(:,:,jn) = ssha_e(:,:)
         !PRINT *, 'jt, jn, za1 =', jt, jn, za1
!                                                 ! ==================== !
      END DO                                               !        end loop      !
!                                                    ! ==================== !
! -----------------------------------------------------------------------------
! Phase 3. update the general trend with the barotropic trend
! -----------------------------------------------------------------------------
!
! At this stage ssha holds a time averaged value
!                                                ! Sea Surface Height at u-,v- and f-points
!      IF( lk_vvl ) THEN                                ! (required only in 1 case)
         DO jj = 1, jpjm1
            DO ji = 1, jpim1      ! NO Vector Opt.
               zsshu_a(ji,jj) = z1_2 * umask(ji,jj,1)  * r1_e12u(ji,jj) &
                  &              * ( e12t(ji  ,jj) * ssha(ji  ,jj)    &
                  &              +   e12t(ji+1,jj) * ssha(ji+1,jj) )
               zsshv_a(ji,jj) = z1_2 * vmask(ji,jj,1)  * r1_e12v(ji,jj) &
                  &              * ( e12t(ji,jj  ) * ssha(ji,jj  )    &
                  &              +   e12t(ji,jj+1) * ssha(ji,jj+1) )
            END DO
         END DO
!         CALL lbc_lnk_multi( zsshu_a, 'U', 1._wp, zsshv_a, 'V', 1._wp ) ! Boundary conditions
!      ENDIF
!
! Set advection velocity correction:
!      IF (((kt==nit000).AND.(neuler==0)).OR.(.NOT.ln_bt_fw)) THEN
         un_adv(:,:) = zu_sum(:,:)*hur(:,:)
         vn_adv(:,:) = zv_sum(:,:)*hvr(:,:)
!      ELSE
!         ....
!      END IF

!      IF (ln_bt_fw) THEN ! Save integrated transport for next computation
!         ....
!      ENDIF
!
! Update barotropic trend:
!      IF (( ln_dynadv_vec ).OR. (.NOT. lk_vvl)) THEN
!         ....
!      ELSE
!         DO jk=1,jpkm1
!            ua(:,:,jk) = ua(:,:,jk) + hur(:,:) * ( ua_b(:,:) - ub_b(:,:) * ehu_b(:,:) ) * z1_2dt_b
!            va(:,:,jk) = va(:,:,jk) + hvr(:,:) * ( va_b(:,:) - vb_b(:,:) * ehv_b(:,:) ) * z1_2dt_b
            spgu(:,:,jk) = hur(:,:) * ( ua_b(:,:) - ub_b(:,:) * ehu_b(:,:) ) * z1_2dt_b
            spgv(:,:,jk) = hvr(:,:) * ( va_b(:,:) - vb_b(:,:) * ehv_b(:,:) ) * z1_2dt_b
!         END DO
! Save barotropic velocities not transport (this is for advection of tracer, no need):
!         ua_b  (:,:) =  ua_b(:,:) / ( hu_0(:,:) + zsshu_a(:,:) + 1._wp - umask(:,:,1) )
!         va_b  (:,:) =  va_b(:,:) / ( hv_0(:,:) + zsshv_a(:,:) + 1._wp - vmask(:,:,1) )
!      ENDIF
!
      DO jk = 1, jpkm1
! Correct velocities:
!         un(:,:,jk) = ( un(:,:,jk) + un_adv(:,:) - un_b(:,:) )*umask(:,:,jk)
!         vn(:,:,jk) = ( vn(:,:,jk) + vn_adv(:,:) - vn_b(:,:) )*vmask(:,:,jk)
!
      END DO
!
   END SUBROUTINE dyn_spg_ts


   SUBROUTINE ts_wgt( jpit, zwgt1, zwgt2)
!!---------------------------------------------------------------------
!!                   ***  ROUTINE ts_wgt  ***
!!
!! ** Purpose : Set time-splitting weights for temporal averaging (or not)
!!----------------------------------------------------------------------
      !LOGICAL, INTENT(in) ::   ll_av      ! temporal averaging=.true.
      !LOGICAL, INTENT(in) ::   ll_fw      ! forward time splitting =.true.
      INTEGER, INTENT(inout) :: jpit      ! cycle length
      REAL(wp), DIMENSION(3*nn_baro), INTENT(inout) ::   zwgt1, & ! Primary weights
                                                         zwgt2    ! Secondary weights

      INTEGER ::  jic, jn, ji                      ! temporary integers
      REAL(wp) :: za1, za2
!!----------------------------------------------------------------------

      zwgt1(:) = 0._wp
      zwgt2(:) = 0._wp

! Set time index when averaged value is requested
!      IF (ll_fw) THEN
!         jic = nn_baro
!      ELSE
         jic = 2 * nn_baro
!      ENDIF

! Set primary weights:
!      IF (ll_av) THEN
! Define simple boxcar window for primary weights
! (width = nn_baro, centered around jic)
!         SELECT CASE ( nn_bt_flt )
!              CASE( 0 )  ! No averaging
!                 ....
!              CASE( 1 )  ! Boxcar, width = nn_baro
!                 ....
!              CASE( 2 )  ! Boxcar, width = 2 * nn_baro
                 DO jn = 1, 3*nn_baro
                    za1 = ABS(float(jn-jic))/float(nn_baro)
                    IF (za1 < 1._wp) THEN
                      zwgt1(jn) = 1._wp
                      jpit = jn
                    ENDIF
                 ENDDO
!              CASE DEFAULT   ;   CALL ctl_stop( 'unrecognised value for nn_bt_flt' )
!         END SELECT

!      ELSE ! No time averaging
!         ....
!      ENDIF

! Set secondary weights
      DO jn = 1, jpit
        DO ji = jn, jpit
             zwgt2(jn) = zwgt2(jn) + zwgt1(ji)
        END DO
      END DO

! Normalize weigths:
      za1 = 1._wp / SUM(zwgt1(1:jpit))
      za2 = 1._wp / SUM(zwgt2(1:jpit))
      DO jn = 1, jpit
        zwgt1(jn) = zwgt1(jn) * za1
        zwgt2(jn) = zwgt2(jn) * za2
      END DO
!
   END SUBROUTINE ts_wgt




   SUBROUTINE zdf_bfr( kt )
!!----------------------------------------------------------------------
!!                   ***  ROUTINE zdf_bfr  ***
!!
!! ** Purpose :   compute the bottom friction coefficient.
!!
!! ** Method  :   Calculate and store part of the momentum trend due
!!              to bottom friction following the chosen friction type
!!              (free-slip, linear, or quadratic). The component
!!              calculated here is multiplied by the bottom velocity in
!!              dyn_bfr to provide the trend term.
!!                The coefficients are updated at each time step only
!!              in the quadratic case.
!!
!! ** Action  :   bfrua , bfrva   bottom friction coefficients
!!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
!!
      INTEGER  ::   ji, jj                       ! dummy loop indices
      INTEGER  ::   ikbt, ikbu, ikbv             ! local integers
      REAL(wp) ::   zvu, zuv, zecu, zecv, ztmp   ! temporary scalars
      REAL(wp), POINTER, DIMENSION(:,:) ::  zbfrt
!!----------------------------------------------------------------------
!
!      IF( nn_bfr == 2 ) THEN                 ! quadratic bottom friction only
!
         ALLOCATE( zbfrt(jpiglo, jpjglo)    )

!         IF ( ln_loglayer.AND.lk_vvl ) THEN ! "log layer" bottom friction coefficient

            DO jj = 1, jpjglo
               DO ji = 1, jpiglo
                  ikbt = mbkt(ji,jj)
!! JC: possible WAD implementation should modify line below if layers vanish
                  ztmp = tmask(ji,jj,ikbt) * ( vkarmn / LOG( 0.5_wp * e3t(ji,jj,ikbt) / rn_bfrz0 ))**2._wp
                  zbfrt(ji,jj) = MAX(bfrcoef2d(ji,jj), ztmp)
                  zbfrt(ji,jj) = MIN(zbfrt(ji,jj), rn_bfri2_max)
               END DO
            END DO
! (ISF)
!            IF ( ln_isfcav ) THEN
!               ...
!            END IF
!
!         ELSE
!            zbfrt(:,:) = bfrcoef2d(:,:)
!            ztfrt(:,:) = tfrcoef2d(:,:)
!         ENDIF

         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               ikbu = mbku(ji,jj)         ! ocean bottom level at u- and v-points
               ikbv = mbkv(ji,jj)         ! (deepest ocean u- and v-points)
!
               zvu  = 0.25 * (  vn(ji,jj  ,ikbu) + vn(ji+1,jj  ,ikbu)     &
                  &           + vn(ji,jj-1,ikbu) + vn(ji+1,jj-1,ikbu)  )
               zuv  = 0.25 * (  un(ji,jj  ,ikbv) + un(ji-1,jj  ,ikbv)     &
                  &           + un(ji,jj+1,ikbv) + un(ji-1,jj+1,ikbv)  )
!
               zecu = SQRT(  un(ji,jj,ikbu) * un(ji,jj,ikbu) + zvu*zvu + rn_bfeb2  )
               zecv = SQRT(  vn(ji,jj,ikbv) * vn(ji,jj,ikbv) + zuv*zuv + rn_bfeb2  )
!
               bfrua(ji,jj) = - 0.5_wp * ( zbfrt(ji,jj) + zbfrt(ji+1,jj  ) ) * zecu
               bfrva(ji,jj) = - 0.5_wp * ( zbfrt(ji,jj) + zbfrt(ji  ,jj+1) ) * zecv
!
! in case of 2 cell water column, we assume each cell feels the top and bottom friction
!               IF ( ln_isfcav ) THEN
!                  ...
!               END IF
            END DO
         END DO

!         IF ( ln_isfcav ) THEN
!            ...
!         ENDIF
!      ENDIF    ! nn_bfr == 2

         DEALLOCATE( zbfrt )

   END SUBROUTINE zdf_bfr

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
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   putrd, pvtrd   ! U and V masked trends
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   ktrd           ! KE trend
!
      INTEGER ::   ji, jj, jk       ! dummy loop indices
!
      REAL(wp)                                :: rau0 = 1026._wp    ! volumic mass of reference     [kg/m3] (from phycst.F90)
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   bu, bv   ! volume of u- and v-boxes
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   r1_bt    ! inverse of t-box volume
!!----------------------------------------------------------------------
!
!
      ALLOCATE( bu(jpiglo,jpjglo,jpk) , bv(jpiglo,jpjglo,jpk) , r1_bt(jpiglo,jpjglo,jpk) )
!
!      IF ( lk_vvl .AND. kt /= nkstp ) THEN   ! Variable volume: set box volume at the 1st call of kt time step
!         nkstp = kt
         DO jk = 1, jpkm1
            bu   (:,:,jk) =  e1u(:,:) * e2u(:,:) * e3u(:,:,jk)
            bv   (:,:,jk) =  e1v(:,:) * e2v(:,:) * e3v(:,:,jk)
            r1_bt(:,:,jk) = 1._wp / ( e12t(:,:) * e3t(:,:,jk) ) * tmask(:,:,jk)
         END DO
!      ENDIF
!
      ktrd(:,:,jpk) = 0._wp
      ktrd(1,:, : ) = 0._wp
      ktrd(:,1, : ) = 0._wp
      DO jk = 1, jpkm1
         DO jj = 2, jpjglo
            DO ji = 2, jpiglo
               ktrd(ji,jj,jk) = 0.5_wp * rau0 *( un(ji  ,jj,jk) * putrd(ji  ,jj,jk) * bu(ji  ,jj,jk)  &
                  &                            + un(ji-1,jj,jk) * putrd(ji-1,jj,jk) * bu(ji-1,jj,jk)  &
                  &                            + vn(ji,jj  ,jk) * pvtrd(ji,jj  ,jk) * bv(ji,jj  ,jk)  &
                  &                            + vn(ji,jj-1,jk) * pvtrd(ji,jj-1,jk) * bv(ji,jj-1,jk)  ) * r1_bt(ji,jj,jk)
            END DO
         END DO
      END DO

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
    ipk(:)                       = jpk
    stypvar(1)%ichunk            = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
    stypvar(1)%cname             = 'spg_u'
    stypvar(1)%cunits            = 'm/s^2'
    stypvar(1)%rmissing_value    = 99999.
    stypvar(1)%valid_min         = -1.
    stypvar(1)%valid_max         = 1.
    stypvar(1)%clong_name        = 'Surface pressure gradient correction, u-momentum'
    stypvar(1)%cshort_name       = 'spg_u'
    stypvar(1)%conline_operation = 'On u-grid'
    stypvar(1)%caxis             = 'time depthu nav_lon_u nav_lat_u'

    stypvar2(1)%ichunk            = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
    stypvar2(1)%cname             = 'spg_v'
    stypvar2(1)%cunits            = 'm/s^2'
    stypvar2(1)%rmissing_value    = 99999.
    stypvar2(1)%valid_min         = -1.
    stypvar2(1)%valid_max         = 1.
    stypvar2(1)%clong_name        = 'Surface pressure gradient correction, v-momentum'
    stypvar2(1)%cshort_name       = 'spg_v'
    stypvar2(1)%conline_operation = 'On v-grid'
    stypvar2(1)%caxis             = 'time depthv nav_lon_v nav_lat_v'

    stypvar3(1)%ichunk            = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
    stypvar3(1)%cname             = 'spg_ke'
    stypvar3(1)%cunits            = 'm^2/s^3'
    stypvar3(1)%rmissing_value    = 99999.
    stypvar3(1)%valid_min         = -1.
    stypvar3(1)%valid_max         = 1.
    stypvar3(1)%clong_name        = 'Surface pressure gradient correction, KE'
    stypvar3(1)%cshort_name       = 'spg_ke'
    stypvar3(1)%conline_operation = 'On t-grid'
    stypvar3(1)%caxis             = 'time deptht nav_lon_t nav_lat_t'

    stypvar4(1)%ichunk            = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
    stypvar4(1)%cname             = 'spg_ssh_bt'
    stypvar4(1)%cunits            = 'm'
    stypvar4(1)%rmissing_value    = 99999.
    stypvar4(1)%valid_min         = -1.
    stypvar4(1)%valid_max         = 1.
    stypvar4(1)%clong_name        = 'Sea surf heigh associated with barotropic mode of split explicit scheme'
    stypvar4(1)%cshort_name       = 'spg_ssh_bt'
    stypvar4(1)%conline_operation = 'On t-grid'
    stypvar4(1)%caxis             = 'time nn_baro nav_lon_t nav_lat_t'



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

    ncout_ssh= create      (cf_out_ssh, cf_tfil ,  jpiglo, jpjglo, jpk, cdep=cn_vdeptht  , ld_nc4=lnc4 )
    ierr    = createvar   (ncout_ssh , stypvar4,  pnvarout, ipk , id_varout_ssh          , ld_nc4=lnc4 )
    ierr    = putheadervar(ncout_ssh , cf_tfil ,  jpiglo, jpjglo, jpk, nav_lon_t, nav_lat_t, deptht   )

    dtim = getvar1d(cf_ufil , cn_vtimec,   jpt     )
    ierr = putvar1d(ncout_u , dtim,        jpt, 'T')
    ierr = putvar1d(ncout_v , dtim,        jpt, 'T')
    ierr = putvar1d(ncout_ke, dtim,        jpt, 'T')
    ierr = putvar1d(ncout_ssh, dtim,        jpt, 'T')


   END SUBROUTINE CreateOutput

END PROGRAM
