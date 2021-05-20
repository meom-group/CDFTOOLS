PROGRAM cdf_dynzdf_imp
  !!======================================================================
  !!                     ***  PROGRAM  cdf_dynzdf_imp  ***
  !!=====================================================================
  !! ** Purpose :   Compute the trend due to the vert. momentum diffusion
  !!
  !! ** Method  :   Explicit forward time stepping 
  !!      following dyn_zdf_imp NEMO routin (cf bellow). 
  !!      The vertical diffusion of momentum is given by:
  !!         diffu = dz( avmu dz(u) ) = 1/e3u dk+1( avmu/e3uw dk(ub) )
  !!      Surface boundary conditions: wind stress input (averaged over kt-1/2 & kt+1/2)
  !!      NO Bottom boundary conditions : bottom stress is computed separatly in cdf_dynbfr
  !!
  !!      QJ: implement tke_tke and tke_avn (from zdftke.f90) to compute
  !!          vertical eddy viscosity [and diffusivity] coefficients ([avt], avmu, avmu)
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
  INTEGER(KIND=4), PARAMETER                   :: pnvarout = 2             ! number of output variables
  INTEGER(KIND=4), PARAMETER                   :: pnvarout2= 6             ! number of output variables (for restarts)
  INTEGER(KIND=4)                              :: ji, jj, jk, jt           ! dummy loop indices
  INTEGER(KIND=4)                              :: it                       ! time index for vvl
  INTEGER(KIND=4)                              :: ierr                     ! working integer
  INTEGER(KIND=4)                              :: narg, iargc, ijarg       ! 
  INTEGER(KIND=4)                              :: jpiglo, jpjglo, jpk, jpt ! size of the domain
  INTEGER(KIND=4)                              :: jpim1, jpjm1, jpkm1      ! index - 1
  INTEGER(KIND=4)                              :: ncout_u, ncout_v         ! ncid of output file
  INTEGER(KIND=4)                              :: ncout_ke                 ! ncid of output file
  INTEGER(KIND=4)                              :: ncout_rst                ! ncid of output file (restart)
  INTEGER(KIND=4)                              :: jpts = 2                 ! Number of active tracers (=2, i.e. T & S )
  INTEGER(KIND=4)                              :: jp_tem = 1               ! indice for temperature
  INTEGER(KIND=4)                              :: jp_sal = 2               ! indice for salinity
  INTEGER(KIND=4), DIMENSION(pnvarout)         :: ipk                      ! level of output variables
  INTEGER(KIND=4), DIMENSION(pnvarout2)        :: ipk2                     ! level of output variables (for restarts)
  INTEGER(KIND=4), DIMENSION(pnvarout)         :: id_varout_u              ! id of output variables (u-comp)
  INTEGER(KIND=4), DIMENSION(pnvarout)         :: id_varout_v              ! id of output variables (v-comp)
  INTEGER(KIND=4), DIMENSION(pnvarout)         :: id_varout_ke             ! id of output variables (ke-comp)
  INTEGER(KIND=4), DIMENSION(pnvarout2)        :: id_varout_rst            ! id of output variables (restarts)
  INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: mbathy                   ! number of vertical wet points
  INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: mbkt, mbku, mbkv         ! vert indices of the shallowest ocean level

  REAL(wp)                                     :: r_vvl = 1._wp            ! Variable volume indicator
  REAL(wp)                                     :: zcoef                    ! temporary scalar
  REAL(wp)                                     :: rmxl_min                 ! minimum mixing length value 
  REAL(wp), DIMENSION(:)  ,   ALLOCATABLE      :: dtim                     ! time
  REAL(wp), DIMENSION(:)  ,   ALLOCATABLE      :: deptht, depthu, depthv   ! z-grid (t,u,v) -- for output
  REAL(wp), DIMENSION(:)  ,   ALLOCATABLE      :: avmb, avtb               ! background profile of avm and avt
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: gphit                    ! latitude of t- points (for TKE penetration)
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: nav_lon_t, nav_lat_t     ! t-grid hor.
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: nav_lon_u, nav_lat_u     ! u-grid hor.
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: nav_lon_v, nav_lat_v     ! v-grid hor.
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: ht_0                     ! Reference ocean depth at T-points
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: hu, hv, hur, hvr         ! Variable ocean depth at (u,v)-points, and inverse
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: htau                     ! depth of tke penetration (nn_htau=1)
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: sshn                     ! Free surface elevation
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE, TARGET:: utau, vtau, taum       ! wind stress and modulus
  REAL(wp), DIMENSION(:,:),   POINTER          :: utau_b, vtau_b           ! wind stress (b=before)
  !REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: utau, vtau, taum         ! wind stress and modulus
  !REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: utau_b, vtau_b           ! wind stress (b=before)
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: bfrua, bfrva             ! bottom friction
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: bfrcoef2d                ! 2D bottom drag coefficient
  !REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: ua_b, va_b               ! Barotropic velocities at (u,v)-point [m/s]

  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: e1t, e2t                 ! horizontal metric, t-pts
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: e1u, e2u                 ! horizontal metric, u-pts
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: e1v, e2v                 ! horizontal metric, v-pts
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: e12t, r1_e12u, r1_e12v   ! face area at t-pts and inverse at u- v- pts
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: e3t_0, e3w_0             ! vet. metrics, t- w- pts, at rest (without VVL)
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: e3u_0, e3v_0             ! vet. metrics, u- v- pts, at rest (without VVL)
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: e3uw_0, e3vw_0           ! vet. metrics, uw- vw- pts, at rest (without VVL)
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: e3u, e3v, e3t, e3w       ! vet. metrics, u-, v-, t-, w-, pts
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: e3uw_n, e3vw_n        ! vet. metrics, uw-, vw- pts (now)
  REAL(wp), DIMENSION(:,:,:), POINTER          :: e3uw_b, e3vw_b           ! vet. metrics, uw-, vw- pts (before)
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: gdept, gdepw             ! depth of t-, w- pts
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: tmask, wmask             ! Mask at T-, W- points
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: umask, vmask             ! Mask at U-, V- points
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: wumask, wvmask           ! Mask at WU-, WV- points

  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: avt, avm, avmu, avmv     ! eddy diffusion/vicosity coefficients
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: avt_out, avm_out         ! outputed diffusion coef. of t and m
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: avmu_out, avmv_out       ! outputed diffusion coef. of u and v
  !REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: tmpu, tmpv               ! debug
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: avt_k, avm_k             ! Not enhanced Kz 
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: avmu_k, avmv_k           ! Not enhanced Kz 
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: en                       ! now turbulent kinetic energy [m2/s2]
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: rn2                      ! brunt-vaisala frequency**2 [s-2]
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: dissl                    ! mixing lenght of dissipation
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: ztrdu, ztrdv             ! temporary save ua and va
                                                                           ! (needed for the implicit formulation)
                                                                           ! and then trends ouptuts
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: zke                      ! KE       trends due to vert. dissip. (output)


  REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE    :: ts                       ! 4D T-S fields, 1: temperature, 2: salinity
  REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE    :: rab                      ! thermal/haline expansion coef. 
                                                                           !        [Celcius-1,psu-1], 1: temp, 2: sal

  REAL(wp), DIMENSION(:,:,:), POINTER          :: ub, vb                   ! horizontal velocity (b=before, a=after)
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET  :: un, vn               ! horizontal velocity (n=now)
  !REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: ub, vb                   ! horizontal velocity (b=before, a=after)
  !REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: un, vn                   ! horizontal velocity (n=now)
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: ua, va                   ! horizontal velocity (a=now)
                                                                           !    (used as workspace in matrix inv for impl. zdf)

  CHARACTER(LEN=256)                           :: cf_tfil                  ! temperature netcdf file name
  CHARACTER(LEN=256)                           :: cf_sfil                  ! salinity netcdf file name
  CHARACTER(LEN=256)                           :: cf_ufil                  ! zonal velocity netcdf file name
  CHARACTER(LEN=256)                           :: cf_vfil                  ! meridional velocity netcdf file name
  CHARACTER(LEN=256)                           :: cf_sshfil                ! ssh          netcdf file name (for vvl)
  CHARACTER(LEN=256)                           :: cf_tauxfil               ! zonal stress netcdf file name
  CHARACTER(LEN=256)                           :: cf_tauyfil               ! merid stress netcdf file name
  CHARACTER(LEN=255)                           :: cf_mh                    ! horiz. mesh netcdf file name
  CHARACTER(LEN=255)                           :: cf_mz                    ! mesh      netcdf file name
  CHARACTER(LEN=255)                           :: cf_mask                  ! mask      netcdf file name
  CHARACTER(LEN=255)                           :: cf_bathy                 ! bathymetry  netcdf file name
  CHARACTER(LEN=256)                           :: cf_out_u='zdf_u.nc'      ! output file name (u-comp)
  CHARACTER(LEN=256)                           :: cf_out_v='zdf_v.nc'      ! output file name (v-comp)
  CHARACTER(LEN=256)                           :: cf_out_ke='zdf_ke.nc'    ! output file name (ke-comp)
  CHARACTER(LEN=255)                           :: cf_in_rst='no_rst_file.nc'   ! restart file (in)
  CHARACTER(LEN=256)                           :: cf_out_rst='rst_out.nc'  ! restart file (out)
  CHARACTER(LEN=256)                           :: cldum                    ! dummy character variable
  CHARACTER(LEN=256)                           :: cglobal                  ! global attribute

  TYPE (variable), DIMENSION(pnvarout)         :: stypvar                  ! structure for attibutes
  TYPE (variable), DIMENSION(pnvarout)         :: stypvar2                 ! structure for attibutes
  TYPE (variable), DIMENSION(pnvarout)         :: stypvar3                 ! structure for attibutes
  TYPE (variable), DIMENSION(pnvarout2)        :: stypvar4                 ! structure for attibutes

  LOGICAL                                      :: l_w   =.FALSE.           ! flag for vertical location of bn2
  LOGICAL                                      :: lchk                     ! check missing files
  LOGICAL                                      :: lfull =.FALSE.           ! full step flag
  LOGICAL                                      :: lnc4  =.FALSE.           ! full step flag
  LOGICAL                                      :: l_rst                    ! read in restart file 
                                                                           !     (will be re-evaluated depending on the value of cf_in_rst)

  ! from phycst.f90
  REAL(wp)                                     :: rau0 = 1026._wp          ! volumic mass of reference     [kg/m3]
  REAL(wp)                                     :: r1_rau0 = 1. / 1026._wp  ! inverse of volumic mass reference [m3/kg]
  REAL(wp)                                     :: grav  = 9.80665_wp       ! gravity         [m/s2]
  REAL(wp)                                     :: rpi = 3.141592653589793_wp    ! pi
  REAL(wp)                                     :: vkarmn = 0.4_wp          ! von Karman constant
  REAL(wp)                                     :: rsmall = 0.5 * EPSILON( 1.e0 )   ! smallest real computer value
  REAL(wp)                                     :: rn_bfrz0 = 3.e-3_wp      ! bottom roughness [m] for log. formulation

  
  ! from namelist
  LOGICAL                                      :: ln_lc =.TRUE.            ! Langmuir cell parameterisation (Axell 2002)
  LOGICAL                                      :: ln_mxl0 = .true.         ! surface mixing length scale = F(wind stress) (T) or not (F)
  REAL(wp)                                     :: rn_avm0 = 1.e-4_wp       ! vertical (background) eddy viscosity [m2/s]
  REAL(wp)                                     :: rn_avt0 = 1.e-5_wp       ! vertical (background) eddy diffusivity [m2/s]
  REAL(wp)                                     :: rdt = 80._wp             ! time step [sec] == model output frequency
  REAL(wp)                                     :: p2dt                     ! 2* time step [sec] == model output frequency
  REAL(wp)                                     :: rn_ebb = 67.83_wp        ! coef. of the surface input of tke 
                                                                           !        (=67.83 suggested when ln_mxl0=T)
  REAL(wp)                                     :: rn_ediss = 0.7_wp        ! coef. of the Kolmogoroff dissipation
  REAL(wp)                                     :: rn_emin = 1.e-6_wp       ! minimum value of tke [m2/s2]
  REAL(wp)                                     :: rn_emin0 = 1.e-4_wp      ! surface minimum value of tke [m2/s2]
  REAL(wp)                                     :: rn_lc = 0.15_wp          ! coef to compute vert. vel of Langmuir cells
  REAL(wp)                                     :: rn_efr = 0.05_wp         ! fraction of TKE surface value 
                                                                           !          which penetrates in the ocean
  REAL(wp)                                     :: rn_ediff = 0.1_wp        ! coefficient for avt: avt=rn_ediff*mxl*sqrt(e)
  REAL(wp)                                     :: rn_mxl0 = 0.01_wp        ! surface  min value of mixing length 
  REAL(wp)                                     :: rn_bshear = 1.e-20_wp    ! background shear (>0) (numerical threshold)
  REAL(wp)                                     :: rn_bfri2 = 2.5e-3_wp     ! min. bottom drag coefficient
  REAL(wp)                                     :: rn_bfri2_max = 1.e-1_wp  ! max. bottom drag coefficient
  REAL(wp)                                     :: rn_bfrien = 50._wp       ! local multiplying factor of bfr
  REAL(wp)                                     :: rn_bfeb2 = 0.0_wp        ! bottom turb. KE background  [m2/s2]
                                                                           !         (=2.5e-3_wp in no tides conditions)


  !-- for bn2 computation --
  ! TEOS10/EOS80 parameters
   REAL(wp) ::   r1_S0, r1_T0, r1_Z0, rdeltaS
! ALPHA parameters
   REAL(wp) ::   ALP000 , ALP100 , ALP200 , ALP300 , ALP400 , ALP500
   REAL(wp) ::   ALP010 , ALP110 , ALP210 , ALP310 , ALP410
   REAL(wp) ::   ALP020 , ALP120 , ALP220 , ALP320
   REAL(wp) ::   ALP030 , ALP130 , ALP230
   REAL(wp) ::   ALP040 , ALP140
   REAL(wp) ::   ALP050
   REAL(wp) ::   ALP001 , ALP101 , ALP201 , ALP301
   REAL(wp) ::   ALP011 , ALP111 , ALP211
   REAL(wp) ::   ALP021 , ALP121
   REAL(wp) ::   ALP031
   REAL(wp) ::   ALP002 , ALP102
   REAL(wp) ::   ALP012
   REAL(wp) ::   ALP003

! BETA parameters
   REAL(wp) ::   BET000 , BET100 , BET200 , BET300 , BET400 , BET500
   REAL(wp) ::   BET010 , BET110 , BET210 , BET310 , BET410
   REAL(wp) ::   BET020 , BET120 , BET220 , BET320
   REAL(wp) ::   BET030 , BET130 , BET230
   REAL(wp) ::   BET040 , BET140
   REAL(wp) ::   BET050
   REAL(wp) ::   BET001 , BET101 , BET201 , BET301
   REAL(wp) ::   BET011 , BET111 , BET211
   REAL(wp) ::   BET021 , BET121
   REAL(wp) ::   BET031
   REAL(wp) ::   BET002 , BET102
   REAL(wp) ::   BET012
   REAL(wp) ::   BET003


  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdf_dynbfr_imp -t T-file -s S-file  -u U-file -v V-file -ssh SSH-file ...'
     PRINT *,'          -tx TAUX-file -ty TAUY ...'
     PRINT *,'          -mh MESH-file -mz MESZ-file -mask MASK-file -bathy BATHY-file ...'
     PRINT *,'          -o_u OUT-file-U -o_v OUT-file-V -o_ke OUT-file-ke ...'
     PRINT *,'          -o_rst OUT-file-RST -i_rst IN-file-RST'
     PRINT *,'      '
     PRINT *,'     PURPOSE :Compute the trend due to vert. momentum diffusion.'
     PRINT *,'              Adopt dynzdf.F90, including dyn_zdf_exp in dynzdf_exp.F90,'
     PRINT *,'              to compute vert. momentum diffusion trend.'
     PRINT *,'              Explicit formulation with no time-splitting.'
     PRINT *,'              Vertical momentum viscosity coefficients computed following TKE formulation.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -t T-file          : netcdf file for potential temperature'
     PRINT *,'       -s S-file          : netcdf file for salinity'
     PRINT *,'       -u U-file          : netcdf file for zonal velocity'
     PRINT *,'       -v V-file          : netcdf file for meridional velocity'
     PRINT *,'       -ssh SSH-file      : netcdf file for SSH (for vvl recomputation of vert. grid)'
     PRINT *,'       -tx TAUX-file      : netcdf file for zonal wind stress'
     PRINT *,'       -ty TAUY-file      : netcdf file for meridional wind stress'
     PRINT *,'       -mh MESH-file      : netcdf file for horizontal mesh'
     PRINT *,'       -mz MESZ-file      : netcdf file for vertical mesh'
     PRINT *,'       -mask MASK-file    : netcdf file for mask'
     PRINT *,'       -bathy BATHY-file  : netcdf file for model bathymetry'
     PRINT *,'       (OPTIONS)'
     PRINT *,'       -o_u OUT-file-U    : netcdf file for vert. diffusion trend for u-momentum'
     PRINT *,'       -o_v OUT-file-V    : netcdf file for vert. diffusion trend for v-momentum'
     PRINT *,'       -o_ke OUT-file-KE  : netcdf file for vert. diffusion trend for KE'
     PRINT *,'       -i_rst IN-file-RST  : netcdf file for restarts (input)'
     PRINT *,'              (set IN-file-RST to no_rst_file.nc if no restart file available)'
     PRINT *,'       -o_rst OUT-file-RST : netcdf file for restarts (output)'
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
     CASE ('-t'        ) ; CALL getarg( ijarg, cf_tfil   ) ; ijarg=ijarg+1
     CASE ('-s'        ) ; CALL getarg( ijarg, cf_sfil   ) ; ijarg=ijarg+1
     CASE ('-u'        ) ; CALL getarg( ijarg, cf_ufil   ) ; ijarg=ijarg+1
     CASE ('-v'        ) ; CALL getarg( ijarg, cf_vfil   ) ; ijarg=ijarg+1
     CASE ('-ssh'      ) ; CALL getarg( ijarg, cf_sshfil ) ; ijarg=ijarg+1
     CASE ('-tx'       ) ; CALL getarg( ijarg, cf_tauxfil) ; ijarg=ijarg+1
     CASE ('-ty'       ) ; CALL getarg( ijarg, cf_tauyfil) ; ijarg=ijarg+1
     CASE ('-mh'       ) ; CALL getarg( ijarg, cf_mh     ) ; ijarg=ijarg+1
     CASE ('-mz'       ) ; CALL getarg( ijarg, cf_mz     ) ; ijarg=ijarg+1
     CASE ('-mask'     ) ; CALL getarg( ijarg, cf_mask   ) ; ijarg=ijarg+1
     CASE ('-bathy'    ) ; CALL getarg( ijarg, cf_bathy  ) ; ijarg=ijarg+1
        ! options
     CASE ( '-full' ) ; lfull   = .TRUE. ; cglobal = 'full step computation'
     CASE ( '-o_u'    ) ; CALL getarg(ijarg, cf_out_u ) ; ijarg = ijarg + 1
     CASE ( '-o_v'    ) ; CALL getarg(ijarg, cf_out_v ) ; ijarg = ijarg + 1
     CASE ( '-o_ke'   ) ; CALL getarg(ijarg, cf_out_ke) ; ijarg = ijarg + 1
        ! restart files (in/out)
     CASE ( '-i_rst'  ) ; CALL getarg(ijarg, cf_in_rst) ; ijarg = ijarg + 1
     CASE ( '-o_rst'  ) ; CALL getarg(ijarg, cf_out_rst); ijarg = ijarg + 1
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
  PRINT *, 'jpk    =', jpk
  PRINT *, 'jpt    =', jpt


  !-- check for restart file to read in turbulent coefficients --
  IF ( TRIM(cf_in_rst) == 'no_rst_file.nc' ) THEN
     PRINT *,'No input restart file'
     l_rst =.FALSE.
  ELSE
     PRINT *,'Read in restart file:'
     PRINT *, cf_in_rst
     l_rst =.TRUE.
  ENDIF


  !-- Allocate arrays --
  !- mesh  and mask -
  ALLOCATE( deptht(jpk)                                                   )
  ALLOCATE( depthu(jpk)                  , depthv(jpk)                    )
  ALLOCATE( nav_lon_t(jpiglo, jpjglo)    , nav_lat_t(jpiglo, jpjglo)      )
  ALLOCATE( nav_lon_u(jpiglo, jpjglo)    , nav_lat_u(jpiglo, jpjglo)      )
  ALLOCATE( nav_lon_v(jpiglo, jpjglo)    , nav_lat_v(jpiglo, jpjglo)      )
  ALLOCATE( gphit(jpiglo, jpjglo)                                         )
  ALLOCATE( ht_0(jpiglo, jpjglo)                                          )
  ALLOCATE( hu(jpiglo, jpjglo)           , hv(jpiglo, jpjglo)             )
  ALLOCATE( hur(jpiglo, jpjglo)          , hvr(jpiglo, jpjglo)            )
  ALLOCATE( e1t(jpiglo, jpjglo)          , e2t(jpiglo, jpjglo)            )
  ALLOCATE( e1u(jpiglo, jpjglo)          , e2u(jpiglo, jpjglo)            )
  ALLOCATE( e1v(jpiglo, jpjglo)          , e2v(jpiglo, jpjglo)            )
  ALLOCATE( e12t(jpiglo, jpjglo)                                          )
  ALLOCATE( r1_e12u(jpiglo, jpjglo)      , r1_e12v(jpiglo, jpjglo)        )
  ALLOCATE( e3t_0( jpiglo, jpjglo, jpk)  , e3w_0(jpiglo, jpjglo, jpk)     )
  ALLOCATE( e3u_0( jpiglo, jpjglo, jpk)  , e3v_0(jpiglo, jpjglo, jpk)     )
  ALLOCATE( e3uw_0( jpiglo, jpjglo, jpk) , e3vw_0(jpiglo, jpjglo, jpk)    )
  ALLOCATE( e3t( jpiglo, jpjglo, jpk)    , e3w(jpiglo, jpjglo, jpk)       )
  ALLOCATE( e3u( jpiglo, jpjglo, jpk)    , e3v(jpiglo, jpjglo, jpk)       )
  ALLOCATE( e3uw_n( jpiglo, jpjglo, jpk)   , e3vw_n(jpiglo, jpjglo, jpk)      )
  ALLOCATE( gdept( jpiglo, jpjglo, jpk)  , gdepw(jpiglo, jpjglo, jpk)     )
  ALLOCATE( tmask( jpiglo, jpjglo, jpk)  , wmask(jpiglo, jpjglo, jpk)     )
  ALLOCATE( umask( jpiglo, jpjglo, jpk)  , vmask(jpiglo, jpjglo, jpk)     )
  ALLOCATE( wumask( jpiglo, jpjglo, jpk) , wvmask(jpiglo, jpjglo, jpk)    )
  ALLOCATE( mbathy( jpiglo, jpjglo)      , mbkt( jpiglo, jpjglo)          )
  ALLOCATE( mbku(jpiglo, jpjglo)         , mbkv(jpiglo, jpjglo)           )
  !- variables -
  ALLOCATE( avmb(jpk)                    , avtb(jpk)                      )
  ALLOCATE( ub(jpiglo, jpjglo, jpk)      , vb(jpiglo, jpjglo, jpk)        )
  ALLOCATE( un(jpiglo, jpjglo, jpk)      , vn(jpiglo, jpjglo, jpk)        )
  ALLOCATE( ua(jpiglo, jpjglo, jpk)      , va(jpiglo, jpjglo, jpk)        )
  !ALLOCATE( ua_b(jpiglo, jpjglo)         , va_b(jpiglo, jpjglo)           )
  ALLOCATE( avt(jpiglo, jpjglo, jpk)     , avm(jpiglo, jpjglo, jpk)       )
  ALLOCATE( avmu(jpiglo, jpjglo, jpk)    , avmv(jpiglo, jpjglo, jpk)      )
  ALLOCATE( avt_out(jpiglo, jpjglo, jpk) , avm_out(jpiglo, jpjglo, jpk)   )
  ALLOCATE( avmu_out(jpiglo, jpjglo, jpk), avmv_out(jpiglo, jpjglo, jpk)  )
  !ALLOCATE( tmpu(jpiglo, jpjglo, jpk)    , tmpv(jpiglo, jpjglo, jpk)      )
  ALLOCATE( avt_k(jpiglo, jpjglo, jpk)   , avm_k(jpiglo, jpjglo, jpk)     )
  ALLOCATE( avmu_k(jpiglo, jpjglo, jpk)  , avmv_k(jpiglo, jpjglo, jpk)    )
  ALLOCATE( utau(jpiglo, jpjglo)         , vtau(jpiglo, jpjglo)      ,  taum(jpiglo, jpjglo)  )
  ALLOCATE( utau_b(jpiglo, jpjglo)       , vtau_b(jpiglo, jpjglo)         )
  ALLOCATE( bfrua(jpiglo, jpjglo)        , bfrva(jpiglo, jpjglo)          )
  ALLOCATE( bfrcoef2d(jpiglo, jpjglo)                                 )
  ALLOCATE( sshn(jpiglo, jpjglo)                                          )
  ALLOCATE( en(jpiglo, jpjglo, jpk)                                       )
  ALLOCATE( ts(jpiglo, jpjglo, jpk, jpts), rab(jpiglo, jpjglo, jpk, jpts) ) !1: temperature, 2: salinity
  ALLOCATE( rn2(jpiglo, jpjglo, jpk)                                      )
  ALLOCATE( dissl(jpiglo, jpjglo, jpk)                                    )
  ALLOCATE( htau(jpiglo, jpjglo)                                          )
  !
  ALLOCATE( ztrdu(jpiglo, jpjglo, jpk)    , ztrdv(jpiglo, jpjglo, jpk)    )
  ALLOCATE( zke(jpiglo, jpjglo, jpk)                                      )

  !-- define thermal expansion and haline contraction coefficients (polynomial TEOS-10, see eosbn2) --
         rdeltaS = 32._wp
         r1_S0  = 0.875_wp/35.16504_wp
         r1_T0  = 1._wp/40._wp
         r1_Z0  = 1.e-4_wp
!
         ALP000 = -6.5025362670e-01_wp
         ALP100 = 1.6320471316_wp
         ALP200 = -2.0442606277_wp
         ALP300 = 1.4222011580_wp
         ALP400 = -4.4204535284e-01_wp
         ALP500 = 4.7983755487e-02_wp
         ALP010 = 1.8537085209_wp
         ALP110 = -3.0774129064_wp
         ALP210 = 3.0181275751_wp
         ALP310 = -1.4565010626_wp
         ALP410 = 2.7361846370e-01_wp
         ALP020 = -1.6246342147_wp
         ALP120 = 2.5086831352_wp
         ALP220 = -1.4787808849_wp
         ALP320 = 2.3807209899e-01_wp
         ALP030 = 8.3627885467e-01_wp
         ALP130 = -1.1311538584_wp
         ALP230 = 5.3563304045e-01_wp
         ALP040 = -6.7560904739e-02_wp
         ALP140 = -6.0212475204e-02_wp
         ALP050 = 2.8625353333e-02_wp
         ALP001 = 3.3340752782e-01_wp
         ALP101 = 1.1217528644e-01_wp
         ALP201 = -1.2510649515e-01_wp
         ALP301 = 1.6349760916e-02_wp
         ALP011 = -3.3540239802e-01_wp
         ALP111 = -1.7531540640e-01_wp
         ALP211 = 9.3976864981e-02_wp
         ALP021 = 1.8487252150e-01_wp
         ALP121 = 4.1307825959e-02_wp
         ALP031 = -5.5927935970e-02_wp
         ALP002 = -5.1410778748e-02_wp
         ALP102 = 5.3278413794e-03_wp
         ALP012 = 6.2099915132e-02_wp
         ALP003 = -9.4924551138e-03_wp
!
         BET000 = 1.0783203594e+01_wp
         BET100 = -4.4452095908e+01_wp
         BET200 = 7.6048755820e+01_wp
         BET300 = -6.3944280668e+01_wp
         BET400 = 2.6890441098e+01_wp
         BET500 = -4.5221697773_wp
         BET010 = -8.1219372432e-01_wp
         BET110 = 2.0346663041_wp
         BET210 = -2.1232895170_wp
         BET310 = 8.7994140485e-01_wp
         BET410 = -1.1939638360e-01_wp
         BET020 = 7.6574242289e-01_wp
         BET120 = -1.5019813020_wp
         BET220 = 1.0872489522_wp
         BET320 = -2.7233429080e-01_wp
         BET030 = -4.1615152308e-01_wp
         BET130 = 4.9061350869e-01_wp
         BET230 = -1.1847737788e-01_wp
         BET040 = 1.4073062708e-01_wp
         BET140 = -1.3327978879e-01_wp
         BET050 = 5.9929880134e-03_wp
         BET001 = -5.2937873009e-01_wp
         BET101 = 1.2634116779_wp
         BET201 = -1.1547328025_wp
         BET301 = 3.2870876279e-01_wp
         BET011 = -5.5824407214e-02_wp
         BET111 = 1.2451933313e-01_wp
         BET211 = -2.4409539932e-02_wp
         BET021 = 4.3623149752e-02_wp
         BET121 = -4.6767901790e-02_wp
         BET031 = -6.8523260060e-03_wp
         BET002 = -6.1618945251e-02_wp
         BET102 = 6.2255521644e-02_wp
         BET012 = -2.6514181169e-03_wp
         BET003 = -2.3025968587e-04_wp
!

  p2dt = 2. * rdt

  !-- load mesh --
  nav_lon_t    = getvar(cf_tfil, 'nav_lon', 1, jpiglo, jpjglo)
  nav_lat_t    = getvar(cf_tfil, 'nav_lat', 1, jpiglo, jpjglo)
  deptht       = getvar1d(cf_tfil, cn_vdeptht , jpk)
  nav_lon_u    = getvar(cf_ufil, 'nav_lon', 1, jpiglo, jpjglo)
  nav_lat_u    = getvar(cf_ufil, 'nav_lat', 1, jpiglo, jpjglo)
  gphit        = getvar(cf_mh  , 'gphit'  , 1, jpiglo, jpjglo)
  depthu       = getvar1d(cf_ufil, cn_vdepthu , jpk)
  nav_lon_v    = getvar(cf_vfil, 'nav_lon', 1, jpiglo, jpjglo)
  nav_lat_v    = getvar(cf_vfil, 'nav_lat', 1, jpiglo, jpjglo)
  depthv       = getvar1d(cf_vfil, cn_vdepthv , jpk)
  ht_0(:,:)    = getvar(cf_bathy, 'gdepw_0', 1, jpiglo, jpjglo )
  e1t(:,:)     = getvar(cf_mh  , 'e1t'  , 1, jpiglo, jpjglo)
  e2t(:,:)     = getvar(cf_mh  , 'e2t'  , 1, jpiglo, jpjglo)
  e1u(:,:)     = getvar(cf_mh  , 'e1u'  , 1, jpiglo, jpjglo)
  e2u(:,:)     = getvar(cf_mh  , 'e2u'  , 1, jpiglo, jpjglo)
  e1v(:,:)     = getvar(cf_mh  , 'e1v'  , 1, jpiglo, jpjglo)
  e2v(:,:)     = getvar(cf_mh  , 'e2v'  , 1, jpiglo, jpjglo)
  e12t(:,:)    = e1t(:,:) * e2t(:,:)
  r1_e12u(:,:) = 1._wp / (e1u(:,:) * e2u(:,:))
  r1_e12v(:,:) = 1._wp / (e1v(:,:) * e2v(:,:))
  !- deepest we points (from domzgr.F90) -
  mbathy(:,:) = getvar(cf_mz, cn_mbathy,1, jpiglo, jpjglo )
  mbkt(:,:) = MAX( mbathy(:,:) , 1 )    ! bottom k-index of T-level (=1 over land)
  !                                     ! bottom k-index of W-level = mbkt+1
  DO jj = 1, jpjm1                      ! bottom k-index of u- (v-) level
     DO ji = 1, jpim1
        mbku(ji,jj) = MIN(  mbkt(ji+1,jj  ) , mbkt(ji,jj)  )
        mbkv(ji,jj) = MIN(  mbkt(ji  ,jj+1) , mbkt(ji,jj)  )
     END DO
  END DO


  !-- load vert. mesh (at rest) and masks (dommsk.f90) --
  e3t_0(:,:,:) = getvar3d(cf_mz  , 'e3t_0' , jpiglo, jpjglo, jpk )
  e3w_0(:,:,:) = getvar3d(cf_mz  , 'e3w_0' , jpiglo, jpjglo, jpk )
  e3u_0(:,:,:) = e3t_0(:,:,:)
  e3v_0(:,:,:) = e3t_0(:,:,:)
  e3uw_0(:,:,:) = e3w_0(:,:,:)
  e3vw_0(:,:,:) = e3w_0(:,:,:)
  DO jk = 1,jpk                         ! Computed as the minimum of neighbooring scale factors
     DO jj = 1, jpjm1
        DO ji = 1, jpim1   ! vector opt.
           e3u_0 (ji,jj,jk) = MIN( e3t_0(ji,jj,jk), e3t_0(ji+1,jj,jk) )
           e3v_0 (ji,jj,jk) = MIN( e3t_0(ji,jj,jk), e3t_0(ji,jj+1,jk) )
           e3uw_0(ji,jj,jk) = MIN( e3w_0(ji,jj,jk), e3w_0(ji+1,jj,jk) )
           e3vw_0(ji,jj,jk) = MIN( e3w_0(ji,jj,jk), e3w_0(ji,jj+1,jk) )
        END DO
     END DO
  END DO
  tmask(:,:,:) = getvar3d(cf_mask, 'tmask' , jpiglo, jpjglo, jpk )
  umask(:,:,:) = getvar3d(cf_mask, 'umask' , jpiglo, jpjglo, jpk )
  vmask(:,:,:) = getvar3d(cf_mask, 'vmask' , jpiglo, jpjglo, jpk )
  wmask(:,:,1) = tmask(:,:,1)
  wumask(:,:,1)= umask(:,:,1)
  wvmask(:,:,1)= vmask(:,:,1)
  DO jk=2,jpk
    wmask (:,:,jk) = tmask(:,:,jk) * tmask(:,:,jk-1)
    wumask(:,:,jk) = umask(:,:,jk) * umask(:,:,jk-1)
    wvmask(:,:,jk) = vmask(:,:,jk) * vmask(:,:,jk-1)
  END DO

  !-- maximum penetration of TKE --
  htau(:,:) = MAX(  0.5_wp, MIN( 30._wp, 45._wp* ABS( SIN( rpi/180._wp * gphit(:,:) ) ) )   )

  !-- minimum length to recover molecular viscosity --
  rmxl_min = 1.e-6_wp / ( rn_ediff * SQRT( rn_emin ) )  

  !-- set background  dissipation/diffusion coefficients -- 
  avmb(:) = rn_avm0
  avtb(:) = rn_avt0

  !-- Creat output netcdf files to fill in --
  CALL CreateOutput

  DO jt = 1, jpt
     PRINT *, '======= time-step = ', jt

     !-- recomputed vert. metric due to non-linear free surface (VVL) --
     ! from domvvl.F90 -->> e3t = e3t_0*(1+sshn/ht_0)
     PRINT *, '-- Recompute vert. mesh --'
     sshn(:,:)     = getvar(cf_sshfil  , cn_sossheig, 1, jpiglo, jpjglo, ktime=jt )
     DO jk = 1, jpk
       e3t(:,:,jk)   = e3t_0(:,:,jk) * (1 + sshn/ht_0)
       e3w(:,:,jk)   = e3w_0(:,:,jk) * (1 + sshn/ht_0)
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
     !- at uw- and vw- points (domvvl.F90) - 
     e3uw_n(:,:,1) = e3uw_0(:,:,1) + e3u(:,:,1) - e3u_0(:,:,1)
     e3vw_n(:,:,1) = e3vw_0(:,:,1) + e3v(:,:,1) - e3v_0(:,:,1)
     DO jk = 2, jpk
        e3uw_n(:,:,jk) = e3uw_0(:,:,jk) + ( 1.0_wp - 0.5_wp * umask(:,:,jk) ) * ( e3u(:,:,jk-1) - e3u_0(:,:,jk-1) )   &
           &                          +            0.5_wp * umask(:,:,jk)   * ( e3u(:,:,jk  ) - e3u_0(:,:,jk  ) )
        e3vw_n(:,:,jk) = e3vw_0(:,:,jk) + ( 1.0_wp - 0.5_wp * vmask(:,:,jk) ) * ( e3v(:,:,jk-1) - e3v_0(:,:,jk-1) )   &
           &                          +            0.5_wp * vmask(:,:,jk)   * ( e3v(:,:,jk  ) - e3v_0(:,:,jk  ) )
     END DO
     e3uw_b => e3uw_n
     e3vw_b => e3vw_n
     !QJ: to turn e3uw_b and e3vw_b to their actual (before) values, just evaluate them before the now variables.
     ! t- and w- points depth
     ! ----------------------
     ! set the isf depth as it is in the initial step
     gdept(:,:,:) = 0._wp
     gdepw(:,:,:) = 0._wp
     gdept(:,:,1) = 0.5_wp * e3w(:,:,1)
     DO jk = 2, jpk
      DO jj = 1,jpjglo
       DO ji = 1,jpiglo
          zcoef = (tmask(ji,jj,jk) - wmask(ji,jj,jk))   ! 0 everywhere tmask = wmask, ie everywhere expect at jk = mikt
                                                        ! 1 everywhere from mbkt to mikt + 1 or 1 (if no isf)
                                                        ! 0.5 where jk = mikt
          gdepw(ji,jj,jk) = gdepw(ji,jj,jk-1) + e3t(ji,jj,jk-1)
          gdept(ji,jj,jk) =      zcoef  * ( gdepw(ji,jj,jk) + 0.5 * e3w(ji,jj,jk))  &
              &                + (1-zcoef) * ( gdept(ji,jj,jk-1) +       e3w(ji,jj,jk))
       END DO
      END DO
     END DO
     ! Local depth of the water column at t- points (from domvvl.F90)
     hu(:,:) = 0._wp
     hv(:,:) = 0._wp
     DO jk = 1, jpk
        hu(:,:) = hu(:,:) + e3u(:,:,jk) * umask(:,:,jk)
        hv(:,:) = hv(:,:) + e3v(:,:,jk) * vmask(:,:,jk)
     END DO
     ! inverse of local depth
     hur(:,:) = 1._wp / ( hu(:,:) + 1._wp - umask(:,:,1) ) * umask(:,:,1)      ! to avoid NaNs
     hvr(:,:) = 1._wp / ( hv(:,:) + 1._wp - vmask(:,:,1) ) * vmask(:,:,1)      ! to avoid NaNs

     !-- wind stress modulus at T-pts --
     PRINT *, '-- Compute wind stress modulus --'
     utau(:,:)     = getvar(cf_tauxfil , 'sozotaux' , 1, jpiglo, jpjglo, ktime=jt )
     vtau(:,:)     = getvar(cf_tauyfil , 'sometauy' , 1, jpiglo, jpjglo, ktime=jt )
     utau_b => utau
     vtau_b => vtau
     !utau_b(:,:)   = getvar(cf_tauxfil , 'sozotaux' , 1, jpiglo, jpjglo, ktime=jt-1 )
     !vtau_b(:,:)   = getvar(cf_tauyfil , 'sometauy' , 1, jpiglo, jpjglo, ktime=jt-1 )
     DO jj = 2, jpjglo
      DO ji = 2, jpiglo
       taum(ji,jj) = SQRT(  ( 0.5_wp * ( utau(ji, jj) + utau(ji-1, jj  ) ) )**2 &
                          + ( 0.5_wp * ( vtau(ji, jj) + vtau(ji  , jj-1) ) )**2 &
                            ) * tmask(ji, jj, 1)
      END DO
     END DO

     !-- Load variables --
     PRINT *, '-- Load hz. vel --'
     un(:,:,:)   = getvar3d(cf_ufil, cn_vozocrtx, jpiglo, jpjglo, jpk, ktime=jt )
     vn(:,:,:)   = getvar3d(cf_vfil, cn_vomecrty, jpiglo, jpjglo, jpk, ktime=jt )
     ub => un
     vb => vn
     !ub(:,:,:)   = getvar3d(cf_ufil, cn_vozocrtx, jpiglo, jpjglo, jpk, ktime=jt-1 )
     !vb(:,:,:)   = getvar3d(cf_vfil, cn_vomecrty, jpiglo, jpjglo, jpk, ktime=jt-1 )
     ! QJ: The implicit formulation is normally (in NEMO) computed 
     ! based on the update velocities ua=ub+2dt*ua, 
     ! where ua in the rhs refers to the trend due to all other terms (cf step.F90).
     ! Such trends are not available here, thus the computation is based on 
     ! the now velocities directly (using after velocities instead does not improve much).
     !ua(:,:,:)    = getvar3d(cf_ufil, cn_vozocrtx, jpiglo, jpjglo, jpk, ktime=jt+1 )
     !va(:,:,:)    = getvar3d(cf_vfil, cn_vomecrty, jpiglo, jpjglo, jpk, ktime=jt+1 )
     ua(:,:,:)    = getvar3d(cf_ufil, cn_vozocrtx, jpiglo, jpjglo, jpk, ktime=jt )
     va(:,:,:)    = getvar3d(cf_vfil, cn_vomecrty, jpiglo, jpjglo, jpk, ktime=jt )
     ztrdu(:,:,:) = ua(:,:,:)
     ztrdv(:,:,:) = va(:,:,:)
     !- barotropic transports -
     !ua_b(:,:) = 0._wp
     !va_b(:,:) = 0._wp
     !DO jk = 1, jpk
     !   ua_b(:,:) = ua_b(:,:) + ua(:,:,jk) * e3u(:,:,jk) * umask(:,:,jk)
     !   va_b(:,:) = va_b(:,:) + va(:,:,jk) * e3v(:,:,jk) * vmask(:,:,jk)
     !ENDDO
     !ua_b(:,:) = ua_b(:,:) * hur(:,:)
     !va_b(:,:) = va_b(:,:) * hvr(:,:)
     
     !-- compute bfrua/bfrva --
     PRINT *, '-- Compute bottom friction coef. --'
     CALL zdf_bfr( jt )
     ! apply 2d map (if provided) for enhanced bfr 
     bfrcoef2d(:,:) = 0._wp
     bfrcoef2d(:,:) = rn_bfri2 * ( 1 + rn_bfrien * bfrcoef2d(:,:) )

     !-- compute Brunt-Vaisala frequency for Langmuir circulation --
     PRINT *, '-- Compute N2 --'
     ts(:,:,:,jp_tem) = getvar3d(cf_tfil, cn_votemper, jpiglo, jpjglo, jpk, ktime=jt )
     ts(:,:,:,jp_sal) = getvar3d(cf_sfil, cn_vosaline, jpiglo, jpjglo, jpk, ktime=jt )
     CALL rab_3d ( ts, rab )
     CALL bn2    ( ts, rab, rn2 )

     !-- compute vertical momentum dissipation coefficients (using TKE formulation) --
     !-- Set vertical dissipation/diffusion coef. to background values at first iter --
     !-- Should be read in a 'restart file' --
     PRINT *, '-- Compute vert mom dissipation coef --'
     IF ( jt == 1 ) THEN
     !   IF ( l_rst ) THEN
     !      en   (:,:,:) = getvar3d(cf_in_rst, 'en'   , jpiglo, jpjglo, jpk )
     !      dissl(:,:,:) = getvar3d(cf_in_rst, 'dissl', jpiglo, jpjglo, jpk )
     !      avt  (:,:,:) = getvar3d(cf_in_rst, 'avt'  , jpiglo, jpjglo, jpk )
     !      avm  (:,:,:) = getvar3d(cf_in_rst, 'avm'  , jpiglo, jpjglo, jpk )
     !      avmu (:,:,:) = getvar3d(cf_in_rst, 'avmu' , jpiglo, jpjglo, jpk )
     !      avmv (:,:,:) = getvar3d(cf_in_rst, 'avmv' , jpiglo, jpjglo, jpk )
     !   ELSE
           en   (:,:,:) = 0._wp
           dissl(:,:,:) = 0._wp
           avt  (:,:,:) = rn_avt0 * wmask (:,:,:)
           avm  (:,:,:) = rn_avm0 * wmask (:,:,:)
           avmu (:,:,:) = rn_avm0 * wumask(:,:,:)
           avmv (:,:,:) = rn_avm0 * wvmask(:,:,:)
     !   ENDIF
     ENDIF
     CALL zdf_tke( jt )

     !-- compute vertical dissipation of momentum --
     PRINT *, '-- Compute vert. dissipation tendencies --'
     CALL dyn_zdf_imp( jt )	! handle bottom friction

     !-- compute and output trends for momentum --
     DO jk = 1, jpk
        ztrdu(:,:,jk) = ( ua(:,:,jk) - ztrdu(:,:,jk) ) / p2dt 
        ztrdv(:,:,jk) = ( va(:,:,jk) - ztrdv(:,:,jk) ) / p2dt
        ierr = putvar(ncout_u , id_varout_u(1), ztrdu(:,:,jk)    , jk, jpiglo, jpjglo, ktime=jt )
        ierr = putvar(ncout_v , id_varout_v(1), ztrdv(:,:,jk)    , jk, jpiglo, jpjglo, ktime=jt )
        ierr = putvar(ncout_u , id_varout_u(2), avmu_out(:,:,jk) , jk, jpiglo, jpjglo, ktime=jt )
        ierr = putvar(ncout_v , id_varout_v(2), avmv_out(:,:,jk) , jk, jpiglo, jpjglo, ktime=jt )
     ENDDO

     !-- KE trends --
     CALL trd_ken( ztrdu, ztrdv, zke)
     DO jk = 1, jpk
        ierr = putvar(ncout_ke, id_varout_ke(1)  , zke(:,:,jk)      , jk, jpiglo, jpjglo, ktime=jt )
        ierr = putvar(ncout_ke, id_varout_ke(2)  , avm_out(:,:,jk)  , jk, jpiglo, jpjglo, ktime=jt )
     ENDDO

     !-- writte restart file --
     IF ( jt == jpt ) THEN
        PRINT *, 'Write dissipation coefficients for restart'
        !DO jk = 1, jpk  
        !   ierr = putvar(ncout_rst , id_varout_rst(1), en(:,:,jk)    , jk, jpiglo, jpjglo, ktime=1 )
        !   ierr = putvar(ncout_rst , id_varout_rst(2), dissl(:,:,jk) , jk, jpiglo, jpjglo, ktime=1 )
        !   ierr = putvar(ncout_rst , id_varout_rst(3), avt(:,:,jk)   , jk, jpiglo, jpjglo, ktime=1 )
        !   ierr = putvar(ncout_rst , id_varout_rst(4), avm(:,:,jk)   , jk, jpiglo, jpjglo, ktime=1 )
        !   ierr = putvar(ncout_rst , id_varout_rst(5), avmu(:,:,jk)  , jk, jpiglo, jpjglo, ktime=1 )
        !   ierr = putvar(ncout_rst , id_varout_rst(6), avmv(:,:,jk)  , jk, jpiglo, jpjglo, ktime=1 )
        !ENDDO
     END IF


  ENDDO		!jt-loop

  ierr = closeout(ncout_u)
  ierr = closeout(ncout_v)
  ierr = closeout(ncout_ke)
  ierr = closeout(ncout_rst)

CONTAINS


   SUBROUTINE zdf_tke( kt )
!!----------------------------------------------------------------------
!!                   ***  ROUTINE zdf_tke  ***
!!
!! ** Purpose :   Compute the vertical eddy viscosity and diffusivity
!!              coefficients using a turbulent closure scheme (TKE).
!!
!! ** Method  :   The time evolution of the turbulent kinetic energy (tke)
!!              is computed from a prognostic equation :
!!         d(en)/dt = avm (d(u)/dz)**2             ! shear production
!!                  + d( avm d(en)/dz )/dz         ! diffusion of tke
!!                  + avt N^2                      ! stratif. destruc.
!!                  - rn_ediss / emxl en**(2/3)    ! Kolmogoroff dissipation
!!      with the boundary conditions:
!!         surface: en = max( rn_emin0, rn_ebb * taum )
!!         bottom : en = rn_emin
!!      The associated critical Richardson number is: ri_cri = 2/(2+rn_ediss/rn_ediff)
!!
!!        The now Turbulent kinetic energy is computed using the following
!!      time stepping: implicit for vertical diffusion term, linearized semi
!!      implicit for kolmogoroff dissipation term, and explicit forward for
!!      both buoyancy and shear production terms. Therefore a tridiagonal
!!      linear system is solved. Note that buoyancy and shear terms are
!!      discretized in a energy conserving form (Bruchard 2002).
!!
!!        The dissipative and mixing length scale are computed from en and
!!      the stratification (see tke_avn)
!!
!!        The now vertical eddy vicosity and diffusivity coefficients are
!!      given by:
!!              avm = max( avtb, rn_ediff * zmxlm * en^1/2 )
!!              avt = max( avmb, pdl * avm                 )
!!              eav = max( avmb, avm )
!!      where pdl, the inverse of the Prandtl number is 1 if nn_pdl=0 and
!!      given by an empirical funtion of the localRichardson number if nn_pdl=1
!!
!! ** Action  :   compute en (now turbulent kinetic energy)
!!                update avt, avmu, avmv (before vertical eddy coef.)
!!
!! References : Gaspar et al., JGR, 1990,
!!              Blanke and Delecluse, JPO, 1991
!!              Mellor and Blumberg, JPO 2004
!!              Axell, JGR, 2002
!!              Bruchard OM 2002
!!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
!!----------------------------------------------------------------------
!
      PRINT *, 'Compute TKE'
!
      IF( kt /= 1 ) THEN   ! restore before value to compute tke
         avt (:,:,:) = avt_k (:,:,:)
         avm (:,:,:) = avm_k (:,:,:)
         avmu(:,:,:) = avmu_k(:,:,:)
         avmv(:,:,:) = avmv_k(:,:,:)
      ENDIF
!
      CALL tke_tke      ! now tke (en)
!
      PRINT *, 'Compute coefficients'
      CALL tke_avn      ! now avt, avm, avmu, avmv
!
      avt_k (:,:,:) = avt (:,:,:)
      avm_k (:,:,:) = avm (:,:,:)
      avmu_k(:,:,:) = avmu(:,:,:)
      avmv_k(:,:,:) = avmv(:,:,:)
!
   END SUBROUTINE zdf_tke


   SUBROUTINE tke_tke
!!----------------------------------------------------------------------
!!                   ***  ROUTINE tke_tke  ***
!!
!! ** Purpose :   Compute the now Turbulente Kinetic Energy (TKE)
!!
!! ** Method  : - TKE surface boundary condition
!!              - source term due to Langmuir cells (Axell JGR 2002) (ln_lc=T)
!!              - source term due to shear (saved in avmu, avmv arrays)
!!              - Now TKE : resolution of the TKE equation by inverting
!!                a tridiagonal linear system by a "methode de chasse"
!!              - increase TKE due to surface and internal wave breaking
!!
!! ** Action  : - en : now turbulent kinetic energy)
!!              - avmu, avmv : production of TKE by shear at u and v-points
!!                (= Kz dz[Ub] * dz[Un] )
!! ---------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk                      ! dummy loop arguments
!!bfr      INTEGER  ::   ikbu, ikbv, ikbum1, ikbvm1      ! temporary scalar
!!bfr      INTEGER  ::   ikbt, ikbumm1, ikbvmm1          ! temporary scalar
      REAL(wp) ::   zrhoa  = 1.22                   ! Air density kg/m3
      REAL(wp) ::   zcdrag = 1.5e-3                 ! drag coefficient
      REAL(wp) ::   zbbrau, zesh2                   ! temporary scalars
      REAL(wp) ::   zfact1, zfact2, zfact3          !    -         -
      REAL(wp) ::   ztx2  , zty2  , zcof            !    -         -
      REAL(wp) ::   ztau  , zdif                    !    -         -
      REAL(wp) ::   zus   , zwlc  , zind            !    -         -
      REAL(wp) ::   zzd_up, zzd_lw                  !    -         -
!!bfr      REAL(wp) ::   zebot                           !    -         -
      INTEGER , POINTER, DIMENSION(:,:  ) :: imlc
      REAL(wp), POINTER, DIMENSION(:,:  ) :: zhlc
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zpelc, zdiag, zd_up, zd_lw
!!--------------------------------------------------------------------
!
!
      ALLOCATE(imlc(jpiglo,jpjglo)       , zhlc(jpiglo,jpjglo)       )
      ALLOCATE(zpelc(jpiglo,jpjglo,jpk)  , zdiag(jpiglo,jpjglo,jpk)  )
      ALLOCATE(zd_up(jpiglo,jpjglo,jpk)  , zd_lw(jpiglo,jpjglo,jpk)  )

      zbbrau = rn_ebb / rau0       ! Local constant initialisation
      zfact1 = -.5_wp * rdt
      zfact2 = 1.5_wp * rdt * rn_ediss
      zfact3 = 0.5_wp       * rn_ediss
!
!
!                     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!                     !  Surface boundary condition on tke
!                     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!      IF ( ln_isfcav ) THEN
!         ....
!      END IF
      DO jj = 2, jpjm1            ! en(1)   = rn_ebb taum / rau0  (min value rn_emin0)
         DO ji = 2, jpim1   ! vector opt.
            en(ji,jj,1) = MAX( rn_emin0, zbbrau * taum(ji,jj) ) * tmask(ji,jj,1)
         END DO
      END DO

!!bfr   - start commented area
!                     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!                     !  Bottom boundary condition on tke
!                     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Tests to date have found the bottom boundary condition on tke to have very little effect.
! The condition is coded here for completion but commented out until there is proof that the
! computational cost is justified
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!                     en(bot)   = (rn_ebb0/rau0)*0.5*sqrt(u_botfr^2+v_botfr^2) (min value rn_emin)
!CDIR NOVERRCHK
!!    DO jj = 2, jpjm1
!CDIR NOVERRCHK
!!       DO ji = 2, jpim1   ! vector opt.
!!          ztx2 = bfrua(ji-1,jj) * ub(ji-1,jj,mbku(ji-1,jj)) + &
!!                 bfrua(ji  ,jj) * ub(ji  ,jj,mbku(ji  ,jj) )
!!          zty2 = bfrva(ji,jj  ) * vb(ji,jj  ,mbkv(ji,jj  )) + &
!!                 bfrva(ji,jj-1) * vb(ji,jj-1,mbkv(ji,jj-1) )
!!          zebot = 0.001875_wp * SQRT( ztx2 * ztx2 + zty2 * zty2 )   !  where 0.001875 = (rn_ebb0/rau0) * 0.5 = 3.75*0.5/1000.
!!          en (ji,jj,mbkt(ji,jj)+1) = MAX( zebot, rn_emin ) * tmask(ji,jj,1)
!!       END DO
!!    END DO
!!bfr   - end commented area
!
!                     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!      IF( ln_lc ) THEN      !  Langmuir circulation source term added to tke       (Axell JGR 2002)
!                  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!                        !* total energy produce by LC : cumulative sum over jk
         zpelc(:,:,1) =  MAX( rn2(:,:,1), 0._wp ) * gdepw(:,:,1) * e3w(:,:,1)
         DO jk = 2, jpk
            zpelc(:,:,jk)  = zpelc(:,:,jk-1) + MAX( rn2(:,:,jk), 0._wp ) * gdepw(:,:,jk) * e3w(:,:,jk)
         END DO
!                        !* finite Langmuir Circulation depth
         zcof = 0.5 * 0.016 * 0.016 / ( zrhoa * zcdrag )
         imlc(:,:) = mbkt(:,:) + 1       ! Initialization to the number of w ocean point (=2 over land)
         DO jk = jpkm1, 2, -1
            DO jj = 1, jpjglo            ! Last w-level at which zpelc>=0.5*us*us
               DO ji = 1, jpiglo         !      with us=0.016*wind(starting from jpk-1)
                  zus  = zcof * taum(ji,jj)
                  IF( zpelc(ji,jj,jk) > zus )   imlc(ji,jj) = jk
               END DO
            END DO
         END DO
!                               ! finite LC depth
         DO jj = 1, jpjglo
            DO ji = 1, jpiglo
               zhlc(ji,jj) = gdepw(ji,jj,imlc(ji,jj))
            END DO
         END DO
         zcof = 0.016 / SQRT( zrhoa * zcdrag )
!CDIR NOVERRCHK
         DO jk = 2, jpkm1         !* TKE Langmuir circulation source term added to en
!CDIR NOVERRCHK
            DO jj = 2, jpjm1
!CDIR NOVERRCHK
               DO ji = 2, jpim1   ! vector opt.
                  zus  = zcof * SQRT( taum(ji,jj) )           ! Stokes drift
!                                           ! vertical velocity due to LC
                  zind = 0.5 - SIGN( 0.5, gdepw(ji,jj,jk) - zhlc(ji,jj) )
                  zwlc = zind * rn_lc * zus * SIN( rpi * gdepw(ji,jj,jk) / zhlc(ji,jj) )
!                                           ! TKE Langmuir circulation source term
                  ! remove ice fraction fr_i from original formulation
                  en(ji,jj,jk) = en(ji,jj,jk) + rdt * ( zwlc * zwlc * zwlc ) /   &
                     &   zhlc(ji,jj) * wmask(ji,jj,jk) * tmask(ji,jj,1)
               END DO
            END DO
         END DO
!
!      ENDIF
!
!                     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!                     !  Now Turbulent kinetic energy (output in en)
!                     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!                     ! Resolution of a tridiagonal linear system by a "methode de chasse"
!                     ! computation from level 2 to jpkm1  (e(1) already computed and e(jpk)=0 ).
!                     ! zdiag : diagonal zd_up : upper diagonal zd_lw : lower diagonal
!
      DO jk = 2, jpkm1           !* Shear production at uw- and vw-points (energy conserving form)
         DO jj = 1, jpjglo              ! here avmu, avmv used as workspace
            DO ji = 1, jpiglo
               avmu(ji,jj,jk) = avmu(ji,jj,jk) * ( un(ji,jj,jk-1) - un(ji,jj,jk) )   &
                  &                            * ( ub(ji,jj,jk-1) - ub(ji,jj,jk) )   &
                  &                            / (  e3uw_n(ji,jj,jk)               &
                  &                              *  e3uw_b(ji,jj,jk)  )
               avmv(ji,jj,jk) = avmv(ji,jj,jk) * ( vn(ji,jj,jk-1) - vn(ji,jj,jk) )   &
                  &                            * ( vb(ji,jj,jk-1) - vb(ji,jj,jk) )   &
                  &                            / (  e3vw_n(ji,jj,jk)               &
                  &                              *  e3vw_b(ji,jj,jk)  )
            END DO
         END DO
      END DO
!
      DO jk = 2, jpkm1           !* Matrix and right hand side in en
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               zcof   = zfact1 * tmask(ji,jj,jk)

               zzd_up = zcof * ( avm  (ji,jj,jk+1) + avm  (ji,jj,jk  ) )   &  ! upper diagonal
                  &          / ( e3t  (ji,jj,jk  ) * e3w  (ji,jj,jk) )
               zzd_lw = zcof * ( avm  (ji,jj,jk  ) + avm  (ji,jj,jk-1) )   &  ! lower diagonal
                  &          / ( e3t  (ji,jj,jk-1) * e3w  (ji,jj,jk) )

!                                                           ! shear prod. at w-point weightened by mask
               zesh2  =  ( avmu(ji-1,jj,jk) + avmu(ji,jj,jk) ) / MAX( 1._wp , umask(ji-1,jj,jk) + umask(ji,jj,jk) )   &
                  &    + ( avmv(ji,jj-1,jk) + avmv(ji,jj,jk) ) / MAX( 1._wp , vmask(ji,jj-1,jk) + vmask(ji,jj,jk) )
!
               zd_up(ji,jj,jk) = zzd_up            ! Matrix (zdiag, zd_up, zd_lw)
               zd_lw(ji,jj,jk) = zzd_lw
               zdiag(ji,jj,jk) = 1._wp - zzd_lw - zzd_up + zfact2 * dissl(ji,jj,jk) * tmask(ji,jj,jk)
!
!                                   ! right hand side in en
               en(ji,jj,jk) = en(ji,jj,jk) + rdt * (  zesh2  -   avt(ji,jj,jk) * rn2(ji,jj,jk)    &
                  &                                 + zfact3 * dissl(ji,jj,jk) * en (ji,jj,jk)  ) &
                  &                                 * wmask(ji,jj,jk)
            END DO
         END DO
      END DO
!                          !* Matrix inversion from level 2 (tke prescribed at level 1)
      DO jk = 3, jpkm1                             ! First recurrence : Dk = Dk - Lk * Uk-1 / Dk-1
         DO jj = 2, jpjm1
            DO ji = 2, jpim1    ! vector opt.
               zdiag(ji,jj,jk) = zdiag(ji,jj,jk) - zd_lw(ji,jj,jk) * zd_up(ji,jj,jk-1) / zdiag(ji,jj,jk-1)
            END DO
         END DO
      END DO
!
! Second recurrence : Lk = RHSk - Lk / Dk-1 * Lk-1
      DO jj = 2, jpjm1
         DO ji = 2, jpim1   ! vector opt.
            zd_lw(ji,jj,2) = en(ji,jj,2) - zd_lw(ji,jj,2) * en(ji,jj,1)    ! Surface boudary conditions on tke
         END DO
      END DO
      DO jk = 3, jpkm1
         DO jj = 2, jpjm1
            DO ji = 2, jpim1    ! vector opt.
               zd_lw(ji,jj,jk) = en(ji,jj,jk) - zd_lw(ji,jj,jk) / zdiag(ji,jj,jk-1) *zd_lw(ji,jj,jk-1)
            END DO
         END DO
      END DO
!
! thrid recurrence : Ek = ( Lk - Uk * Ek+1 ) / Dk
      DO jj = 2, jpjm1
         DO ji = 2, jpim1   ! vector opt.
            en(ji,jj,jpkm1) = zd_lw(ji,jj,jpkm1) / zdiag(ji,jj,jpkm1)
         END DO
      END DO
      DO jk = jpk-2, 2, -1
         DO jj = 2, jpjm1
            DO ji = 2, jpim1    ! vector opt.
               en(ji,jj,jk) = ( zd_lw(ji,jj,jk) - zd_up(ji,jj,jk) * en(ji,jj,jk+1) ) / zdiag(ji,jj,jk)
            END DO
         END DO
      END DO
      DO jk = 2, jpkm1                             ! set the minimum value of tke
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               en(ji,jj,jk) = MAX( en(ji,jj,jk), rn_emin ) * wmask(ji,jj,jk)
            END DO
         END DO
      END DO

!                            !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!                            !  TKE due to surface and internal wave breaking
!                            !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!      IF( nn_etau == 1 ) THEN           !* penetration below the mixed layer (rn_efr fraction)
         DO jk = 2, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  ! remove ice fraction fr_i from original formulation
                  en(ji,jj,jk) = en(ji,jj,jk) + rn_efr * en(ji,jj,1) * EXP( -gdepw(ji,jj,jk) / htau(ji,jj) )   &
                     &                                 * wmask(ji,jj,jk) * tmask(ji,jj,1)
               END DO
            END DO
         END DO
!      ELSEIF( nn_etau == 2 ) THEN       !* act only at the base of the mixed layer (jk=nmln)  (rn_efr fraction)
!          ....
!      ELSEIF( nn_etau == 3 ) THEN       !* penetration belox the mixed layer (HF variability)
!         ....
!      ENDIF
!
      DEALLOCATE(imlc, zhlc, zpelc, zdiag, zd_up, zd_lw )

   END SUBROUTINE tke_tke



   SUBROUTINE tke_avn
!!----------------------------------------------------------------------
!!                   ***  ROUTINE tke_avn  ***
!!
!! ** Purpose :   Compute the vertical eddy viscosity and diffusivity
!!
!! ** Method  :   At this stage, en, the now TKE, is known (computed in
!!              the tke_tke routine). First, the now mixing lenth is
!!      computed from en and the strafification (N^2), then the mixings
!!      coefficients are computed.
!!              - Mixing length : a first evaluation of the mixing lengh
!!      scales is:
!!                      mxl = sqrt(2*en) / N
!!      where N is the brunt-vaisala frequency, with a minimum value set
!!      to rmxl_min (rn_mxl0) in the interior (surface) ocean.
!!        The mixing and dissipative length scale are bound as follow :
!!         nn_mxl=0 : mxl bounded by the distance to surface and bottom.
!!                        zmxld = zmxlm = mxl
!!         nn_mxl=1 : mxl bounded by the e3w and zmxld = zmxlm = mxl
!!         nn_mxl=2 : mxl bounded such that the vertical derivative of mxl is
!!                    less than 1 (|d/dz(mxl)|<1) and zmxld = zmxlm = mxl
!!         nn_mxl=3 : mxl is bounded from the surface to the bottom usings
!!                    |d/dz(xml)|<1 to obtain lup, and from the bottom to
!!                    the surface to obtain ldown. the resulting length
!!                    scales are:
!!                        zmxld = sqrt( lup * ldown )
!!                        zmxlm = min ( lup , ldown )
!!              - Vertical eddy viscosity and diffusivity:
!!                      avm = max( avtb, rn_ediff * zmxlm * en^1/2 )
!!                      avt = max( avmb, pdlr * avm )
!!      with pdlr=1 if nn_pdl=0, pdlr=1/pdl=F(Ri) otherwise.
!!
!! ** Action  : - avt : now vertical eddy diffusivity (w-point)
!!              - avmu, avmv : now vertical eddy viscosity at uw- and vw-points
!!----------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zrn2, zraug, zcoef, zav     ! local scalars
      REAL(wp) ::   zdku, zpdlr, zri, zsqen     !   -      -
      REAL(wp) ::   zdkv, zemxl, zemlm, zemlp   !   -      -
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zmpdl, zmxlm, zmxld
!!--------------------------------------------------------------------
!
      ALLOCATE( zmpdl(jpiglo,jpjglo,jpk), zmxlm(jpiglo,jpjglo,jpk) )
      ALLOCATE( zmxld(jpiglo,jpjglo,jpk)                           )
!
!                     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!                     !  Mixing length
!                     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!
!                     !* Buoyancy length scale: l=sqrt(2*e/n**2)
!
! initialisation of interior minimum value (avoid a 2d loop with mikt)
      zmxlm(:,:,:)  = rmxl_min
      zmxld(:,:,:)  = rmxl_min
!
!      IF( ln_mxl0 ) THEN            ! surface mixing length = F(stress) : l=vkarmn*2.e5*taum/(rau0*g)
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               zraug = vkarmn * 2.e5_wp / ( rau0 * grav )
               zmxlm(ji,jj,1) = MAX( rn_mxl0, zraug * taum(ji,jj) * tmask(ji,jj,1) )
            END DO
         END DO
!      ELSE
!         zmxlm(:,:,1) = rn_mxl0
!      ENDIF
!
      DO jk = 2, jpkm1              ! interior value : l=sqrt(2*e/n^2)
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               zrn2 = MAX( rn2(ji,jj,jk), rsmall )
               zmxlm(ji,jj,jk) = MAX(  rmxl_min,  SQRT( 2._wp * en(ji,jj,jk) / zrn2 ) )
            END DO
         END DO
      END DO
!
!                     !* Physical limits for the mixing length
!
      zmxld(:,:,1  ) = zmxlm(:,:,1)   ! surface set to the minimum value
      zmxld(:,:,jpk) = rmxl_min       ! last level  set to the minimum value
!
!      SELECT CASE ( nn_mxl )
!      CASE ( 0 )           ! bounded by the distance to surface and bottom
!         ...
!      CASE ( 1 )           ! bounded by the vertical scale factor
!         ...
!      CASE ( 2 )           ! |dk[xml]| bounded by e3t :
!         ...
!      CASE ( 3 )           ! lup and ldown, |dk[xml]| bounded by e3t :
         DO jk = 2, jpkm1         ! from the surface to the bottom : lup
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  zmxld(ji,jj,jk) = MIN( zmxld(ji,jj,jk-1) + e3t(ji,jj,jk-1), zmxlm(ji,jj,jk) )
               END DO
            END DO
         END DO
         DO jk = jpkm1, 2, -1     ! from the bottom to the surface : ldown
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  zmxlm(ji,jj,jk) = MIN( zmxlm(ji,jj,jk+1) + e3t(ji,jj,jk+1), zmxlm(ji,jj,jk) )
               END DO
            END DO
         END DO
         DO jk = 2, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  zemlm = MIN ( zmxld(ji,jj,jk),  zmxlm(ji,jj,jk) )
                  zemlp = SQRT( zmxld(ji,jj,jk) * zmxlm(ji,jj,jk) )
                  zmxlm(ji,jj,jk) = zemlm
                  zmxld(ji,jj,jk) = zemlp
               END DO
            END DO
         END DO


!      END SELECT

!                     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!                     !  Vertical eddy viscosity and diffusivity  (avmu, avmv, avt)
!                     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      DO jk = 1, jpkm1            !* vertical eddy viscosity & diffivity at w-points
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               zsqen = SQRT( en(ji,jj,jk) )
               zav   = rn_ediff * zmxlm(ji,jj,jk) * zsqen
               avm  (ji,jj,jk) = MAX( zav,                  avmb(jk) ) * wmask(ji,jj,jk)
               !QJ: remove the decrease avtb in the equatorial band 15S-15N (avtb_2d)
               avt  (ji,jj,jk) = MAX( zav,                  avtb(jk) ) * wmask(ji,jj,jk)
               dissl(ji,jj,jk) = zsqen / zmxld(ji,jj,jk)
            END DO
         END DO
      END DO
!
      DO jk = 2, jpkm1            !* vertical eddy viscosity at wu- and wv-points
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               avmu(ji,jj,jk) = 0.5 * ( avm(ji,jj,jk) + avm(ji+1,jj  ,jk) ) * wumask(ji,jj,jk)
               avmv(ji,jj,jk) = 0.5 * ( avm(ji,jj,jk) + avm(ji  ,jj+1,jk) ) * wvmask(ji,jj,jk)
            END DO
         END DO
      END DO
!
!      IF( nn_pdl == 1 ) THEN      !* Prandtl number case: update avt
         DO jk = 2, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  zcoef = avm(ji,jj,jk) * 2._wp * e3w(ji,jj,jk) * e3w(ji,jj,jk)
!                                          ! shear
                  zdku = avmu(ji-1,jj,jk) * ( un(ji-1,jj,jk-1) - un(ji-1,jj,jk) ) * ( ub(ji-1,jj,jk-1) - ub(ji-1,jj,jk) )   &
                    &  + avmu(ji  ,jj,jk) * ( un(ji  ,jj,jk-1) - un(ji  ,jj,jk) ) * ( ub(ji  ,jj,jk-1) - ub(ji  ,jj,jk) )
                  zdkv = avmv(ji,jj-1,jk) * ( vn(ji,jj-1,jk-1) - vn(ji,jj-1,jk) ) * ( vb(ji,jj-1,jk-1) - vb(ji,jj-1,jk) )   &
                    &  + avmv(ji,jj  ,jk) * ( vn(ji,jj  ,jk-1) - vn(ji,jj  ,jk) ) * ( vb(ji,jj  ,jk-1) - vb(ji,jj  ,jk) )
!                                          ! local Richardson number
                  zri   = MAX( rn2(ji,jj,jk), 0._wp ) * zcoef / (zdku + zdkv + rn_bshear )
                  zpdlr = MAX(  0.1_wp,  0.2 / MAX( 0.2 , zri )  )
!!gm and even better with the use of the "true" ri_crit=0.22222...  (this change the results!)
!!gm              zpdlr = MAX(  0.1_wp,  ri_crit / MAX( ri_crit , zri )  )
                  !QJ: remove the decrease avtb in the equatorial band 15S-15N (avtb_2d)
                  avt(ji,jj,jk)   = MAX( zpdlr * avt(ji,jj,jk), avtb(jk) ) * wmask(ji,jj,jk)
              END DO
            END DO
         END DO
!      ENDIF
      !QJ: store diffusion coefficients for output
      avm_out(:,:,:) = avm(:,:,:)
      avt_out(:,:,:) = avt(:,:,:)
      avmu_out(:,:,:) = avmu(:,:,:)
      avmv_out(:,:,:) = avmv(:,:,:)
!
      DEALLOCATE( zmpdl, zmxlm, zmxld )


   END SUBROUTINE tke_avn




   SUBROUTINE dyn_zdf_imp( kt )
!!----------------------------------------------------------------------
!!                  ***  ROUTINE dyn_zdf_imp  ***
!!
!! ** Purpose :   Compute the trend due to the vert. momentum diffusion
!!      and the surface forcing, and add it to the general trend of
!!      the momentum equations.
!!
!! ** Method  :   The vertical momentum mixing trend is given by :
!!             dz( avmu dz(u) ) = 1/e3u dk+1( avmu/e3uw dk(ua) )
!!      backward time stepping
!!      Surface boundary conditions: wind stress input (averaged over kt-1/2 & kt+1/2)
!!      Bottom boundary conditions : bottom stress (cf zdfbfr.F)
!!      Add this trend to the general trend ua :
!!         ua = ua + dz( avmu dz(u) )
!!
!! ** Action : - Update (ua,va) arrays with the after vertical diffusive mixing trend.
!!---------------------------------------------------------------------
      INTEGER , INTENT(in) ::  kt     ! ocean time-step index
!!
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      INTEGER  ::   ikbu, ikbv   ! local integers
      REAL(wp) ::   z1_p2dt, zcoef, zzwi, zzws, zrhs   ! local scalars
      REAL(wp) ::   ze3ua, ze3va, zzz
      REAL(wp), POINTER, DIMENSION(:,:)   ::  z2d
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  zwi, zwd, zws
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  z3d
!!----------------------------------------------------------------------

      ALLOCATE( zwi(jpiglo,jpjglo, jpk) )
      ALLOCATE( zwd(jpiglo,jpjglo, jpk) )
      ALLOCATE( zws(jpiglo,jpjglo, jpk) )
!
!      IF( lk_vvl ) THEN   ;    
!      r_vvl = 1._wp       ! Variable volume indicator
!      ELSE                ;    r_vvl = 0._wp
!      ENDIF

! 0. Local constant initialization
! --------------------------------
      z1_p2dt = 1._wp / p2dt      ! inverse of the timestep

! 1. Apply semi-implicit bottom friction
! --------------------------------------
! Only needed for semi-implicit bottom friction setup. The explicit
! bottom friction has been included in "u(v)a" which act as the R.H.S
! column vector of the tri-diagonal matrix equation
!
!      IF( ln_bfrimp ) THEN
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               ikbu = mbku(ji,jj)       ! ocean bottom level at u- and v-points
               ikbv = mbkv(ji,jj)       ! (deepest ocean u- and v-points)
               avmu(ji,jj,ikbu+1) = -bfrua(ji,jj) * e3uw_n(ji,jj,ikbu+1)
               avmv(ji,jj,ikbv+1) = -bfrva(ji,jj) * e3vw_n(ji,jj,ikbv+1)
               !- QJ: add to outputs -
               avmu_out(ji,jj,ikbu+1) = -bfrua(ji,jj) * e3uw_n(ji,jj,ikbu+1)
               avmv_out(ji,jj,ikbv+1) = -bfrva(ji,jj) * e3vw_n(ji,jj,ikbv+1)
            END DO
         END DO
!         IF ( ln_isfcav ) THEN
!            ...
!         ENDIF
!      ENDIF


!      IF( ln_dynadv_vec .OR. .NOT. lk_vvl ) THEN      ! applied on velocity
!         ...
!!      ELSE                                            ! applied on thickness weighted velocity
! QJ: This does not need to be evaluated since (ua, va) as been set to 
!     actual model velocities instead of trends as it is in NEMO at this stage
! 18/02/2021: test this implementation with available momentum trends
!         DO jk = 1, jpkm1
!            ua(:,:,jk) = (          ub(:,:,jk) * e3u(:,:,jk)      &
!               &           + p2dt * ua(:,:,jk) * e3u(:,:,jk)  )   &
!               &                               / e3u(:,:,jk) * umask(:,:,jk)
!            va(:,:,jk) = (          vb(:,:,jk) * e3v(:,:,jk)      &
!               &           + p2dt * va(:,:,jk) * e3v(:,:,jk)  )   &
!               &                               / e3v(:,:,jk) * vmask(:,:,jk)
!         END DO
!!      ENDIF

! QJ: treats barotropic bottom stress separatly -- Not sure how to deal with this
!     Bottom stress due to barotropic component only is included in ua, va
!     since those are model velocities ...
!      IF ( ln_bfrimp .AND.lk_dynspg_ts ) THEN
! remove barotropic velocities:
!         DO jk = 1, jpkm1
!            ua(:,:,jk) = (ua(:,:,jk) - ua_b(:,:)) * umask(:,:,jk)
!            va(:,:,jk) = (va(:,:,jk) - va_b(:,:)) * vmask(:,:,jk)
!         END DO
! Add bottom/top stress due to barotropic component only:
!         DO jj = 2, jpjm1
!            DO ji = 2, jpim1   ! vector opt.
!               ikbu = mbku(ji,jj)         ! ocean bottom level at u- and v-points
!               ikbv = mbkv(ji,jj)         ! (deepest ocean u- and v-points)
!               ze3ua =  ( 1._wp - r_vvl ) * e3u(ji,jj,ikbu) + r_vvl   * e3u(ji,jj,ikbu)
!               ze3va =  ( 1._wp - r_vvl ) * e3v(ji,jj,ikbv) + r_vvl   * e3v(ji,jj,ikbv)
!               ua(ji,jj,ikbu) = ua(ji,jj,ikbu) + p2dt * bfrua(ji,jj) * ua_b(ji,jj) / ze3ua
!               va(ji,jj,ikbv) = va(ji,jj,ikbv) + p2dt * bfrva(ji,jj) * va_b(ji,jj) / ze3va
!            END DO
!         END DO
!         IF ( ln_isfcav ) THEN
!            ...
!         END IF
!      ENDIF



! 2. Vertical diffusion on u
! ---------------------------
! Matrix and second member construction
! bottom boundary condition: both zwi and zws must be masked as avmu can take
! non zero value at the ocean bottom depending on the bottom friction used.
!
      DO jk = 1, jpkm1        ! Matrix
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               ze3ua =  ( 1._wp - r_vvl ) * e3u(ji,jj,jk) + r_vvl   * e3u(ji,jj,jk)   ! after scale factor at T-point
               zcoef = - p2dt / ze3ua
               zzwi          = zcoef * avmu  (ji,jj,jk  ) / e3uw_n(ji,jj,jk)
               zwi(ji,jj,jk) = zzwi  * wumask(ji,jj,jk  )
               zzws          = zcoef * avmu  (ji,jj,jk+1) / e3uw_n(ji,jj,jk+1)
               zws(ji,jj,jk) = zzws  * wumask(ji,jj,jk+1)
               zwd(ji,jj,jk) = 1._wp - zzwi - zzws
            END DO
         END DO
      END DO
      DO jj = 2, jpjm1        ! Surface boundary conditions
         DO ji = 2, jpim1   ! vector opt.
            zwi(ji,jj,1) = 0._wp
            zwd(ji,jj,1) = 1._wp - zws(ji,jj,1)
         END DO
      END DO


! Matrix inversion starting from the first level
!-----------------------------------------------------------------------
!   solve m.x = y  where m is a tri diagonal matrix ( jpk*jpk )
!
!        ( zwd1 zws1   0    0    0  )( zwx1 ) ( zwy1 )
!        ( zwi2 zwd2 zws2   0    0  )( zwx2 ) ( zwy2 )
!        (  0   zwi3 zwd3 zws3   0  )( zwx3 )=( zwy3 )
!        (        ...               )( ...  ) ( ...  )
!        (  0    0    0   zwik zwdk )( zwxk ) ( zwyk )
!
!   m is decomposed in the product of an upper and a lower triangular matrix
!   The 3 diagonal terms are in 2d arrays: zwd, zws, zwi
!   The solution (the after velocity) is in ua
!-----------------------------------------------------------------------
!
!==  First recurrence : Dk = Dk - Lk * Uk-1 / Dk-1   (increasing k)  ==
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               zwd(ji,jj,jk) = zwd(ji,jj,jk) - zwi(ji,jj,jk) * zws(ji,jj,jk-1) / zwd(ji,jj,jk-1)
            END DO
         END DO
      END DO
!
      DO jj = 2, jpjm1        !==  second recurrence:    SOLk = RHSk - Lk / Dk-1  Lk-1  ==
         DO ji = 2, jpim1   ! vector opt.

            ze3ua =  ( 1._wp - r_vvl ) * e3u(ji,jj,1) + r_vvl   * e3u(ji,jj,1)
            ua(ji,jj,1) = ua(ji,jj,1) + p2dt * 0.5_wp * ( utau_b(ji,jj) + utau(ji,jj) )   &
               &                                      / ( ze3ua * rau0 ) * umask(ji,jj,1)
         END DO
      END DO
      !QJ: add surface stress to diffusion coef. for outputs
      DO jj = 1, jpjglo
         DO ji = 1, jpiglo
            avmu_out(ji,jj,1) = 0.5_wp * ( utau_b(ji,jj) + utau(ji,jj) )   &
               &                                      / ( ze3ua * rau0 ) * umask(ji,jj,1)
         ENDDO
      ENDDO
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = 2, jpim1

               zrhs = ua(ji,jj,jk)   ! zrhs=right hand side

               ua(ji,jj,jk) = zrhs - zwi(ji,jj,jk) / zwd(ji,jj,jk-1) * ua(ji,jj,jk-1)
            END DO
         END DO
      END DO
!
      DO jj = 2, jpjm1        !==  thrid recurrence : SOLk = ( Lk - Uk * Ek+1 ) / Dk  ==
         DO ji = 2, jpim1   ! vector opt.
            ua(ji,jj,jpkm1) = ua(ji,jj,jpkm1) / zwd(ji,jj,jpkm1)
         END DO
      END DO
      DO jk = jpk-2, 1, -1
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               ua(ji,jj,jk) = ( ua(ji,jj,jk) - zws(ji,jj,jk) * ua(ji,jj,jk+1) ) / zwd(ji,jj,jk)
            END DO
         END DO
      END DO


! 3. Vertical diffusion on v
! ---------------------------
! Matrix and second member construction
! bottom boundary condition: both zwi and zws must be masked as avmv can take
! non zero value at the ocean bottom depending on the bottom friction used
!
      DO jk = 1, jpkm1        ! Matrix
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               ze3va =  ( 1._wp - r_vvl ) * e3v(ji,jj,jk)  + r_vvl * e3v(ji,jj,jk)   ! after scale factor at T-point
               zcoef = - p2dt / ze3va
               zzwi          = zcoef * avmv (ji,jj,jk  ) / e3vw_n(ji,jj,jk)
               zwi(ji,jj,jk) =  zzwi * wvmask(ji,jj,jk)
               zzws          = zcoef * avmv (ji,jj,jk+1) / e3vw_n(ji,jj,jk+1)
               zws(ji,jj,jk) =  zzws * wvmask(ji,jj,jk+1)
               zwd(ji,jj,jk) = 1._wp - zzwi - zzws
            END DO
         END DO
      END DO
      DO jj = 2, jpjm1        ! Surface boundary conditions
         DO ji = 2, jpim1   ! vector opt.
            zwi(ji,jj,1) = 0._wp
            zwd(ji,jj,1) = 1._wp - zws(ji,jj,1)
         END DO
      END DO

! Matrix inversion
!-----------------------------------------------------------------------
!   solve m.x = y  where m is a tri diagonal matrix ( jpk*jpk )
!
!        ( zwd1 zws1   0    0    0  )( zwx1 ) ( zwy1 )
!        ( zwi2 zwd2 zws2   0    0  )( zwx2 ) ( zwy2 )
!        (  0   zwi3 zwd3 zws3   0  )( zwx3 )=( zwy3 )
!        (        ...               )( ...  ) ( ...  )
!        (  0    0    0   zwik zwdk )( zwxk ) ( zwyk )
!
!   m is decomposed in the product of an upper and lower triangular matrix
!   The 3 diagonal terms are in 2d arrays: zwd, zws, zwi
!   The solution (after velocity) is in 2d array va
!-----------------------------------------------------------------------
!
!==  First recurrence : Dk = Dk - Lk * Uk-1 / Dk-1   (increasing k)  ==
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               zwd(ji,jj,jk) = zwd(ji,jj,jk) - zwi(ji,jj,jk) * zws(ji,jj,jk-1) / zwd(ji,jj,jk-1)
            END DO
         END DO
      END DO
!
      DO jj = 2, jpjm1        !==  second recurrence:    SOLk = RHSk - Lk / Dk-1  Lk-1  ==
         DO ji = 2, jpim1   ! vector opt.

            ze3va =  ( 1._wp - r_vvl ) * e3v(ji,jj,1) + r_vvl   * e3v(ji,jj,1)
            va(ji,jj,1) = va(ji,jj,1) + p2dt * 0.5_wp * ( vtau_b(ji,jj) + vtau(ji,jj) )   &
               &                                      / ( ze3va * rau0 ) * vmask(ji,jj,1)
         END DO
      END DO
      !QJ: add surface stress to diffusion coef. for outputs
      DO jj = 1,jpjglo
         DO ji = 1,jpiglo
            avmv_out(ji,jj,1) = 0.5_wp * ( vtau_b(ji,jj) + vtau(ji,jj) )   &
               &                                      / ( ze3va * rau0 ) * vmask(ji,jj,1)
         ENDDO
      ENDDO
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.

               zrhs = va(ji,jj,jk)   ! zrhs=right hand side

               va(ji,jj,jk) = zrhs - zwi(ji,jj,jk) / zwd(ji,jj,jk-1) * va(ji,jj,jk-1)
            END DO
         END DO
      END DO
      DO jj = 2, jpjm1        !==  third recurrence : SOLk = ( Lk - Uk * SOLk+1 ) / Dk  ==
         DO ji = 2, jpim1   ! vector opt.
            va(ji,jj,jpkm1) = va(ji,jj,jpkm1) / zwd(ji,jj,jpkm1)
         END DO
      END DO
      DO jk = jpk-2, 1, -1
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               va(ji,jj,jk) = ( va(ji,jj,jk) - zws(ji,jj,jk) * va(ji,jj,jk+1) ) / zwd(ji,jj,jk)
            END DO
         END DO
      END DO

!      IF( iom_use( 'dispkevfo' ) ) THEN   ! ocean kinetic energy dissipation per unit area
!         ...
!      ENDIF

!
   END SUBROUTINE dyn_zdf_imp


   SUBROUTINE bn2( pts, pab, pn2 )
!!----------------------------------------------------------------------
!!                  ***  ROUTINE bn2  ***
!!
!! ** Purpose :   Compute the local Brunt-Vaisala frequency at the
!!                time-step of the input arguments
!!
!! ** Method  :   pn2 = grav * (alpha dk[T] + beta dk[S] ) / e3w
!!      where alpha and beta are given in pab, and computed on T-points.
!!      N.B. N^2 is set one for all to zero at jk=1 in istate module.
!!
!! ** Action  :   pn2 : square of the brunt-vaisala frequency at w-point
!!
!!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpiglo,jpjglo,jpk,jpts), INTENT(in   ) :: pts ! pot. temperature and salinity   [Celcius,psu]
      REAL(wp), DIMENSION(jpiglo,jpjglo,jpk,jpts), INTENT(in   ) :: pab ! thermal/haline expansion coef.  [Celcius-1,psu-1]
      REAL(wp), DIMENSION(jpiglo,jpjglo,jpk     ), INTENT(  out) :: pn2 ! Brunt-Vaisala frequency squared [1/s^2]
!
      INTEGER  ::   ji, jj, jk      ! dummy loop indices
      REAL(wp) ::   zaw, zbw, zrw   ! local scalars
!!----------------------------------------------------------------------
!
      DO jk = 2, jpkm1           ! interior points only (2=< jk =< jpkm1 )
         DO jj = 1, jpjglo       ! surface and bottom value set to zero one for all in istate.F90
            DO ji = 1, jpiglo
               zrw =   ( gdepw(ji,jj,jk) - gdept(ji,jj,jk) )   &
                  &  / ( gdept(ji,jj,jk-1) - gdept(ji,jj,jk) )
!
               zaw = pab(ji,jj,jk,jp_tem) * (1. - zrw) + pab(ji,jj,jk-1,jp_tem) * zrw
               zbw = pab(ji,jj,jk,jp_sal) * (1. - zrw) + pab(ji,jj,jk-1,jp_sal) * zrw
!
               pn2(ji,jj,jk) = grav * (  zaw * ( pts(ji,jj,jk-1,jp_tem) - pts(ji,jj,jk,jp_tem) )     &
                  &                    - zbw * ( pts(ji,jj,jk-1,jp_sal) - pts(ji,jj,jk,jp_sal) )  )  &
                  &            / e3w(ji,jj,jk) * tmask(ji,jj,jk)
            END DO
         END DO
      END DO
!
   END SUBROUTINE bn2

   SUBROUTINE rab_3d( pts, pab )
!!----------------------------------------------------------------------
!!                 ***  ROUTINE rab_3d  ***
!!
!! ** Purpose :   Calculates thermal/haline expansion ratio at T-points
!!
!! ** Method  :   calculates alpha / beta at T-points
!!
!! ** Action  : - pab     : thermal/haline expansion ratio at T-points
!!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpiglo,jpjglo,jpk,jpts), INTENT(in   ) ::   pts   ! pot. temperature & salinity
      REAL(wp), DIMENSION(jpiglo,jpjglo,jpk,jpts), INTENT(  out) ::   pab   ! thermal/haline expansion ratio
!
      INTEGER  ::   ji, jj, jk                ! dummy loop indices
      REAL(wp) ::   zt , zh , zs , ztm        ! local scalars
      REAL(wp) ::   zn , zn0, zn1, zn2, zn3   !   -      -
!!----------------------------------------------------------------------
!      SELECT CASE ( nn_eos )
!
!      CASE( -1, 0 )                !==  polynomial TEOS-10 / EOS-80 ==!

         DO jk = 1, jpkm1
            DO jj = 1, jpjglo
               DO ji = 1, jpiglo
!
                  zh  = gdept(ji,jj,jk) * r1_Z0                                ! depth
                  zt  = pts (ji,jj,jk,jp_tem) * r1_T0                           ! temperature
                  zs  = SQRT( ABS( pts(ji,jj,jk,jp_sal) + rdeltaS ) * r1_S0 )   ! square root salinity
                  ztm = tmask(ji,jj,jk)                                         ! tmask
!
! alpha
                  zn3 = ALP003
!
                  zn2 = ALP012*zt + ALP102*zs+ALP002
!
                  zn1 = ((ALP031*zt   &
                     &   + ALP121*zs+ALP021)*zt   &
                     &   + (ALP211*zs+ALP111)*zs+ALP011)*zt   &
                     &   + ((ALP301*zs+ALP201)*zs+ALP101)*zs+ALP001
!
                  zn0 = ((((ALP050*zt   &
                     &   + ALP140*zs+ALP040)*zt   &
                     &   + (ALP230*zs+ALP130)*zs+ALP030)*zt   &
                     &   + ((ALP320*zs+ALP220)*zs+ALP120)*zs+ALP020)*zt   &
                     &   + (((ALP410*zs+ALP310)*zs+ALP210)*zs+ALP110)*zs+ALP010)*zt   &
                     &   + ((((ALP500*zs+ALP400)*zs+ALP300)*zs+ALP200)*zs+ALP100)*zs+ALP000
!
                  zn  = ( ( zn3 * zh + zn2 ) * zh + zn1 ) * zh + zn0
!
                  pab(ji,jj,jk,jp_tem) = zn * r1_rau0 * ztm
!
! beta
                  zn3 = BET003
!
                  zn2 = BET012*zt + BET102*zs+BET002
!
                  zn1 = ((BET031*zt   &
                     &   + BET121*zs+BET021)*zt   &
                     &   + (BET211*zs+BET111)*zs+BET011)*zt   &
                     &   + ((BET301*zs+BET201)*zs+BET101)*zs+BET001
!
                  zn0 = ((((BET050*zt   &
                     &   + BET140*zs+BET040)*zt   &
                     &   + (BET230*zs+BET130)*zs+BET030)*zt   &
                     &   + ((BET320*zs+BET220)*zs+BET120)*zs+BET020)*zt   &
                     &   + (((BET410*zs+BET310)*zs+BET210)*zs+BET110)*zs+BET010)*zt   &
                     &   + ((((BET500*zs+BET400)*zs+BET300)*zs+BET200)*zs+BET100)*zs+BET000
!
                  zn  = ( ( zn3 * zh + zn2 ) * zh + zn1 ) * zh + zn0
!
                  pab(ji,jj,jk,jp_sal) = zn / zs * r1_rau0 * ztm
!
               END DO
            END DO
         END DO




!      CASE( 1 )                  !==  simplified EOS  ==!
!         ....
!      END SELECT



   END SUBROUTINE rab_3d


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
    stypvar(1)%cname             = 'zdf_u'
    stypvar(1)%cunits            = 'm/s^2'
    stypvar(1)%rmissing_value    = 99999.
    stypvar(1)%valid_min         = -1.
    stypvar(1)%valid_max         = 1.
    stypvar(1)%clong_name        = 'Vertical eddy viscosity dissipation from TKE, u-momentum'
    stypvar(1)%cshort_name       = 'zdf_u'
    stypvar(1)%conline_operation = 'On u-grid'
    stypvar(1)%caxis             = 'time depthu nav_lon_u nav_lat_u'
    !
    stypvar(2)%ichunk            = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
    stypvar(2)%cname             = 'avmu'
    stypvar(2)%cunits            = 'm^2/s'
    stypvar(2)%rmissing_value    = 99999.
    stypvar(2)%valid_min         = -1.
    stypvar(2)%valid_max         = 1.
    stypvar(2)%clong_name        = 'Vertical dissipation coefficients, u-momentum'
    stypvar(2)%cshort_name       = 'avmu'
    stypvar(2)%conline_operation = 'On u-grid'
    stypvar(2)%caxis             = 'time depthu nav_lon_u nav_lat_u'

    stypvar2(1)%ichunk            = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
    stypvar2(1)%cname             = 'zdf_v'
    stypvar2(1)%cunits            = 'm/s^2'
    stypvar2(1)%rmissing_value    = 99999.
    stypvar2(1)%valid_min         = -1.
    stypvar2(1)%valid_max         = 1.
    stypvar2(1)%clong_name        = 'Vertical eddy viscosity dissipation from TKE, v-momentum'
    stypvar2(1)%cshort_name       = 'zdf_v'
    stypvar2(1)%conline_operation = 'On v-grid mbkv (see ./domzgr.F90 for definition)'
    stypvar2(1)%caxis             = 'time depthv nav_lon_v nav_lat_v'
    !
    stypvar2(2)%ichunk            = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
    stypvar2(2)%cname             = 'avmv'
    stypvar2(2)%cunits            = 'm^2/s'
    stypvar2(2)%rmissing_value    = 99999.
    stypvar2(2)%valid_min         = -1.
    stypvar2(2)%valid_max         = 1.
    stypvar2(2)%clong_name        = 'Vertical dissipation coefficients, v-momentum'
    stypvar2(2)%cshort_name       = 'avmv'
    stypvar2(2)%conline_operation = 'On v-grid'
    stypvar2(2)%caxis             = 'time depthv nav_lon_v nav_lat_v'

    stypvar3(1)%ichunk            = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
    stypvar3(1)%cname             = 'zdf_ke'
    stypvar3(1)%cunits            = 'm^2/s^3'
    stypvar3(1)%rmissing_value    = 99999.
    stypvar3(1)%valid_min         = -1.
    stypvar3(1)%valid_max         = 1.
    stypvar3(1)%clong_name        = 'KE trend induced by vertical dissipation from TKE'
    stypvar3(1)%cshort_name       = 'zdf_ke'
    stypvar3(1)%conline_operation = 'On t-grid'
    stypvar3(1)%caxis             = 'time deptht nav_lon_t nav_lat_t'
    !
    stypvar3(2)%ichunk            = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
    stypvar3(2)%cname             = 'avm'
    stypvar3(2)%cunits            = 'm^2/s^3'
    stypvar3(2)%rmissing_value    = 99999.
    stypvar3(2)%valid_min         = -1.
    stypvar3(2)%valid_max         = 1.
    stypvar3(2)%clong_name        = 'Vertical dissipation coefficients, t-point'
    stypvar3(2)%cshort_name       = 'avm'
    stypvar3(2)%conline_operation = 'On t-grid'
    stypvar3(2)%caxis             = 'time deptht nav_lon_t nav_lat_t'

    stypvar4(1)%ichunk            = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
    stypvar4(1)%cname             = 'en'
    stypvar4(1)%cunits            = ''
    stypvar4(1)%rmissing_value    = 99999.
    stypvar4(1)%valid_min         = -1.
    stypvar4(1)%valid_max         = 1.
    stypvar4(1)%clong_name        = 'turbulent kinetic energy from TKE -- restart'
    stypvar4(1)%cshort_name       = 'en'
    stypvar4(1)%conline_operation = 'On t-grid'
    stypvar4(1)%caxis             = 'time deptht nav_lon_t nav_lat_t'
    !
    stypvar4(2)%ichunk            = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
    stypvar4(2)%cname             = 'dissl'
    stypvar4(2)%cunits            = ''
    stypvar4(2)%rmissing_value    = 99999.
    stypvar4(2)%valid_min         = -1.
    stypvar4(2)%valid_max         = 1.
    stypvar4(2)%clong_name        = 'turbulent length scale from TKE -- restart'
    stypvar4(2)%cshort_name       = 'dissl'
    stypvar4(2)%conline_operation = 'On t-grid'
    stypvar4(2)%caxis             = 'time deptht nav_lon_t nav_lat_t'
    !
    stypvar4(3)%ichunk            = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
    stypvar4(3)%cname             = 'avt'
    stypvar4(3)%cunits            = ''
    stypvar4(3)%rmissing_value    = 99999.
    stypvar4(3)%valid_min         = -1.
    stypvar4(3)%valid_max         = 1.
    stypvar4(3)%clong_name        = 'Turbulent vertical coefficient for tracer from TKE -- restart'
    stypvar4(3)%cshort_name       = 'avt'
    stypvar4(3)%conline_operation = 'On t-grid'
    stypvar4(3)%caxis             = 'time deptht nav_lon_t nav_lat_t'
    !
    stypvar4(4)%ichunk            = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
    stypvar4(4)%cname             = 'avm'
    stypvar4(4)%cunits            = ''
    stypvar4(4)%rmissing_value    = 99999.
    stypvar4(4)%valid_min         = -1.
    stypvar4(4)%valid_max         = 1.
    stypvar4(4)%clong_name        = 'Turbulent vertical coefficient for momentum from TKE -- restart'
    stypvar4(4)%cshort_name       = 'avm'
    stypvar4(4)%conline_operation = 'On t-grid'
    stypvar4(4)%caxis             = 'time deptht nav_lon_t nav_lat_t'
    !
    stypvar4(5)%ichunk            = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
    stypvar4(5)%cname             = 'avmu'
    stypvar4(5)%cunits            = ''
    stypvar4(5)%rmissing_value    = 99999.
    stypvar4(5)%valid_min         = -1.
    stypvar4(5)%valid_max         = 1.
    stypvar4(5)%clong_name        = 'Turbulent vertical coefficient for momentum (u-comp) from TKE -- restart'
    stypvar4(5)%cshort_name       = 'avmu'
    stypvar4(5)%conline_operation = 'On t-grid'
    stypvar4(5)%caxis             = 'time deptht nav_lon_t nav_lat_t'
    !
    stypvar4(6)%ichunk            = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
    stypvar4(6)%cname             = 'avmv'
    stypvar4(6)%cunits            = ''
    stypvar4(6)%rmissing_value    = 99999.
    stypvar4(6)%valid_min         = -1.
    stypvar4(6)%valid_max         = 1.
    stypvar4(6)%clong_name        = 'Turbulent vertical coefficient for momentum (v-comp) from TKE -- restart'
    stypvar4(6)%cshort_name       = 'avmv'
    stypvar4(6)%conline_operation = 'On t-grid'
    stypvar4(6)%caxis             = 'time deptht nav_lon_t nav_lat_t'

    ! create output fileset
    ncout_u = create      (cf_out_u, cf_ufil ,  jpiglo, jpjglo, jpk, cdep=cn_vdepthu  , ld_nc4=lnc4 )
    ierr    = createvar   (ncout_u , stypvar ,  pnvarout, ipk , id_varout_u           , ld_nc4=lnc4 )
    ierr    = putheadervar(ncout_u , cf_ufil ,  jpiglo, jpjglo, jpk, nav_lon_u, nav_lat_u             )

    ncout_v = create      (cf_out_v, cf_vfil ,  jpiglo, jpjglo, jpk, cdep=cn_vdepthv  , ld_nc4=lnc4 )
    ierr    = createvar   (ncout_v , stypvar2,  pnvarout, ipk,  id_varout_v           , ld_nc4=lnc4 )
    ierr    = putheadervar(ncout_v , cf_vfil ,  jpiglo, jpjglo, jpk, nav_lon_v, nav_lat_v             )

    ncout_ke= create      (cf_out_ke, cf_tfil ,  jpiglo, jpjglo, jpk, cdep=cn_vdeptht  , ld_nc4=lnc4 )
    ierr    = createvar   (ncout_ke , stypvar3,  pnvarout, ipk , id_varout_ke          , ld_nc4=lnc4 )
    ierr    = putheadervar(ncout_ke , cf_tfil ,  jpiglo, jpjglo, jpk, nav_lon_t, nav_lat_t, deptht   )

    ncout_rst= create      (cf_out_rst, cf_tfil ,  jpiglo, jpjglo, jpk, cdep=cn_vdeptht  , ld_nc4=lnc4 )
    ierr     = createvar   (ncout_rst , stypvar4,  pnvarout2, ipk , id_varout_rst        , ld_nc4=lnc4 )
    ierr     = putheadervar(ncout_rst , cf_tfil ,  jpiglo, jpjglo, jpk, nav_lon_t, nav_lat_t, deptht   )

    dtim = getvar1d(cf_ufil, cn_vtimec,   jpt     )
    ierr = putvar1d(ncout_u, dtim,        jpt, 'T')
    ierr = putvar1d(ncout_v, dtim,        jpt, 'T')
    ierr = putvar1d(ncout_ke, dtim,       jpt, 'T')
    ierr = putvar1d(ncout_rst, dtim,      1  , 'T')


   END SUBROUTINE CreateOutput

END PROGRAM
