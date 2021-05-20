PROGRAM cdf_dynhpg 
  !!======================================================================
  !!                     ***  PROGRAM  cdf_dynhpg  ***
  !!=====================================================================
  !!  ** Purpose : Copmpute the horizontal gradient of the hydrostatic pressure
  !!               following the hydrostatic formulation used in eNATL60 
  !!               (should be updated for different hpg_ scheme).
  !!               Follow the s-coordinate standard jacobian formulation of Madec et al. 1996
  !!               Includes the surface pressure gradient computation through the level
  !!               thickness variation in the non-linear freesurface (vvl) formulation.
  !!               Compute hpg trends foru, v and ke.
  !!
  !! History : 4.0  : 09/2019  : Q. Jamet & J.M. Molines : Original code
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames   ! for cdf variable names
  !USE eos           ! to compute in situ density
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class Equation_of_state
  !!----------------------------------------------------------------------
  IMPLICIT NONE


  INTEGER(KIND=4), PARAMETER                   :: wp=8
  INTEGER(KIND=4), PARAMETER                   :: pnvarout = 1             ! number of output variablesi
  INTEGER(KIND=4), PARAMETER                   :: jpts = 2                 ! number of tracers
  INTEGER(KIND=4), PARAMETER                   :: jp_tem = 1               ! indice for temperature
  INTEGER(KIND=4), PARAMETER                   :: jp_sal = 2               ! indice for salinity
  INTEGER(KIND=4), PARAMETER                   :: nn_eos = -1              ! = -1/0 (other not coded) type of eq. of state
  INTEGER(KIND=4)                              :: ji, jj, jk, jt           ! dummy loop indices
  INTEGER(KIND=4)                              :: it                       ! time index for vvl
  INTEGER(KIND=4)                              :: ierr                     ! working integer
  INTEGER(KIND=4)                              :: narg, iargc, ijarg       ! 
  INTEGER(KIND=4)                              :: jpiglo, jpjglo, jpk, jpt ! size of the domain
  INTEGER(KIND=4)                              :: jpim1, jpjm1, jpkm1      ! index - 1
  INTEGER(KIND=4)                              :: jkkm1=1, jkk=2           ! for swapping the levels
  INTEGER(KIND=4)                              :: ncout_u, ncout_v, ncout_ke         ! ncid of output file
  INTEGER(KIND=4), DIMENSION(pnvarout)         :: ipk                      ! level of output variables
  INTEGER(KIND=4), DIMENSION(pnvarout)         :: id_varout_u              ! id of output variables (u-comp)
  INTEGER(KIND=4), DIMENSION(pnvarout)         :: id_varout_v              ! id of output variables (v-comp)
  INTEGER(KIND=4), DIMENSION(pnvarout)         :: id_varout_ke             ! id of output variables (ke-comp)

  REAL(wp)                                     :: zspval                   ! missing value
  REAL(wp)                                     :: zcoef                    ! bottom cell coef for partial step correction
  REAL(wp), DIMENSION(:)      , ALLOCATABLE    :: dtim                     ! time
  REAL(KIND=4), DIMENSION(:)  , ALLOCATABLE    :: deptht, depthu, depthv   ! z-grid (t,u,v) -- for output  
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE    :: nav_lon_t, nav_lat_t     ! t-grid hor.  -- need real*4
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE    :: nav_lon_u, nav_lat_u     ! u-grid hor.  -- need real*4
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE    :: nav_lon_v, nav_lat_v     ! v-grid hor.
  REAL(wp), DIMENSION(:,:)    , ALLOCATABLE    :: ht_0                     ! Reference ocean depth at T-points
  REAL(wp), DIMENSION(:,:)    , ALLOCATABLE    :: sshn                     ! now sea surface height
  REAL(wp), DIMENSION(:,:)    , ALLOCATABLE    :: e1t, e2t                 ! horizontal metric, t-pts
  REAL(wp), DIMENSION(:,:)    , ALLOCATABLE    :: e1u, e2u                 ! horizontal metric, u-pts
  REAL(wp), DIMENSION(:,:)    , ALLOCATABLE    :: e1v, e2v                 ! horizontal metric, v-pts
  REAL(wp), DIMENSION(:,:)    , ALLOCATABLE    :: e12t, r1_e12u, r1_e12v   ! face area at t-pts and inverse at u- v- pts
  REAL(wp), DIMENSION(:,:,:)  , ALLOCATABLE    :: e3t_0, e3u_0, e3v_0, e3w_0      ! vet. metrics at rest (without vvl)
  REAL(wp), DIMENSION(:,:,:)  , ALLOCATABLE    :: e3u, e3v, e3t, e3w       ! vet. metrics, u- v- t- pts
  REAL(wp), DIMENSION(:,:,:)  , ALLOCATABLE    :: tmask, wmask             ! Mask at T- W- points
  REAL(wp), DIMENSION(:,:,:)  , ALLOCATABLE    :: umask, vmask             ! Mask at U- V- points
  REAL(wp), DIMENSION(:,:,:)  , ALLOCATABLE    :: gdep3w                   ! depth of T-points (sum of e3w) (m) 
  REAL(wp), DIMENSION(:,:,:)  , ALLOCATABLE    :: gdept, gdepw             ! depth of the grid cell with NL free surf corr
  REAL(wp), DIMENSION(:,:,:)  , ALLOCATABLE    :: un, vn                   ! 3D hz. velocity (now)
  REAL(wp), DIMENSION(:,:,:)  , ALLOCATABLE    :: rhd                      ! In-situ density
  REAL(wp), DIMENSION(:,:,:)  , ALLOCATABLE    :: zhpi, zhpj, zhpke        ! hydrostatic pressure gradient tendency (outputed)
  REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE    :: ts                       ! temp and salt field (1: temp; 2: salt)

  CHARACTER(LEN=256)                           :: cf_tfil                  ! temperature netcdf file name
  CHARACTER(LEN=256)                           :: cf_sfil                  ! salinity    netcdf file name
  CHARACTER(LEN=256)                           :: cf_ufil                  ! zonal vel   netcdf file name
  CHARACTER(LEN=256)                           :: cf_vfil                  ! merid vel   netcdf file name
  CHARACTER(LEN=256)                           :: cf_sshfil                ! ssh         netcdf file name (for vvl)
  CHARACTER(LEN=255)                           :: cf_mh                    ! horiz. mesh netcdf file name
  CHARACTER(LEN=255)                           :: cf_mz                    ! vert. mesh  netcdf file name
  CHARACTER(LEN=255)                           :: cf_mask                  ! mask        netcdf file name
  CHARACTER(LEN=255)                           :: cf_bathy                 ! bathymetry  netcdf file name
  CHARACTER(LEN=256)                           :: cf_out_u='hpg_u.nc'      ! output file name (u-comp)
  CHARACTER(LEN=256)                           :: cf_out_v='hpg_v.nc'      ! output file name (v-comp)
  CHARACTER(LEN=256)                           :: cf_out_ke='hpg_ke.nc'    ! output file name (ke-comp)
  CHARACTER(LEN=256)                           :: cldum                    ! dummy character variable
  CHARACTER(LEN=256)                           :: cglobal                  ! global attribute

  TYPE (variable), DIMENSION(pnvarout)         :: stypvar                  ! structure for attibutes (u-comp)
  TYPE (variable), DIMENSION(pnvarout)         :: stypvar2                 ! structure for attibutes (v-comp)
  TYPE (variable), DIMENSION(pnvarout)         :: stypvar3                 ! structure for attibutes (ke-comp)

  LOGICAL                                      :: l_w   =.FALSE.           ! flag for vertical location of bn2
  LOGICAL                                      :: lchk                     ! check missing files
  LOGICAL                                      :: lfull =.FALSE.           ! full step flag
  LOGICAL                                      :: lnc4  =.FALSE.           ! full step flag

!-- from phycst.F90 --
  REAL(wp)  :: grav  = 9.80665_wp            ! gravity [m/s2]
  REAL(wp)  :: rau0 = 1026._wp               ! volumic mass of reference     [kg/m3]
  REAL(wp)  :: rcp = 3991.86795711963_wp     ! heat capacity [J/K]
  REAL(wp)  :: r1_rau0                       ! inverse of volumic mass reference [m3/kg]
  REAL(wp)  :: r1_rcp                        ! = 1. / rcp [Kelvin/J]
  REAL(wp)  :: rau0_rcp                      ! = rau0 * rcp
  REAL(wp)  :: r1_rau0_rcp                   ! = 1. / ( rau0 * rcp )

!-- from eosbn2.F90 --
  ! TEOS10/EOS80 parameters
   REAL(wp) ::   r1_S0, r1_T0, r1_Z0, rdeltaS

! EOS parameters
   REAL(wp) ::   EOS000 , EOS100 , EOS200 , EOS300 , EOS400 , EOS500 , EOS600
   REAL(wp) ::   EOS010 , EOS110 , EOS210 , EOS310 , EOS410 , EOS510
   REAL(wp) ::   EOS020 , EOS120 , EOS220 , EOS320 , EOS420
   REAL(wp) ::   EOS030 , EOS130 , EOS230 , EOS330
   REAL(wp) ::   EOS040 , EOS140 , EOS240
   REAL(wp) ::   EOS050 , EOS150
   REAL(wp) ::   EOS060
   REAL(wp) ::   EOS001 , EOS101 , EOS201 , EOS301 , EOS401
   REAL(wp) ::   EOS011 , EOS111 , EOS211 , EOS311
   REAL(wp) ::   EOS021 , EOS121 , EOS221
   REAL(wp) ::   EOS031 , EOS131
   REAL(wp) ::   EOS041
   REAL(wp) ::   EOS002 , EOS102 , EOS202
   REAL(wp) ::   EOS012 , EOS112
   REAL(wp) ::   EOS022
   REAL(wp) ::   EOS003 , EOS103
   REAL(wp) ::   EOS013

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

! PEN parameters
   REAL(wp) ::   PEN000 , PEN100 , PEN200 , PEN300 , PEN400
   REAL(wp) ::   PEN010 , PEN110 , PEN210 , PEN310
   REAL(wp) ::   PEN020 , PEN120 , PEN220
   REAL(wp) ::   PEN030 , PEN130
   REAL(wp) ::   PEN040
   REAL(wp) ::   PEN001 , PEN101 , PEN201
   REAL(wp) ::   PEN011 , PEN111
   REAL(wp) ::   PEN021
   REAL(wp) ::   PEN002 , PEN102
   REAL(wp) ::   PEN012

! ALPHA_PEN parameters
   REAL(wp) ::   APE000 , APE100 , APE200 , APE300
   REAL(wp) ::   APE010 , APE110 , APE210
   REAL(wp) ::   APE020 , APE120
   REAL(wp) ::   APE030
   REAL(wp) ::   APE001 , APE101
   REAL(wp) ::   APE011
   REAL(wp) ::   APE002

! BETA_PEN parameters
   REAL(wp) ::   BPE000 , BPE100 , BPE200 , BPE300
   REAL(wp) ::   BPE010 , BPE110 , BPE210
   REAL(wp) ::   BPE020 , BPE120
   REAL(wp) ::   BPE030
   REAL(wp) ::   BPE001 , BPE101
   REAL(wp) ::   BPE011
   REAL(wp) ::   BPE002

  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdf_dynhpg -t T-file -s S-file -u U-file -v V-file -ssh SSH-file ...'
     PRINT *,'          -mh MESH-file -mz MESZ-file -bathy BATHY-file -mask MASK-file ...'
     PRINT *,'          -o_u OUT-file-u -o_v OUT-file-v -o_ke OUT-file-ke'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'      Compute the now momentum advection trend in flux form'
     PRINT *,'      and the general trend of the momentum equation.'
     PRINT *,'      Compute the associated now KE advection trend.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -t T-file          : netcdf file for potential temperature'
     PRINT *,'       -s S-file          : netcdf file for salinity'
     PRINT *,'       -u U-file          : netcdf file for zonal vel'
     PRINT *,'       -v V-file          : netcdf file for meridional vel'
     PRINT *,'       -ssh SSH-file      : netcdf file for SSH (for vvl recomputation of vert. grid)'
     PRINT *,'       -mh MESH-file      : netcdf file for horizontal mesh'
     PRINT *,'       -mz MESH-file      : netcdf file for vertical mesh'
     PRINT *,'       -mask MASK-file    : netcdf file for mask'
     PRINT *,'       -bathy BATHY-file  : netcdf file for model bathymetry'
     PRINT *,'       -o_u OUT-file      : netcdf file for hydrostatic pressure grad. term for u-momentum'
     PRINT *,'       -o_v OUT-file      : netcdf file for hydrostatic pressure grad. term for v-momentum'
     PRINT *,'       -o_ke OUT-file     : netcdf file for hydrostatic pressure grad. term for KE'
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
     CASE ('-t'        ) ; CALL getarg( ijarg, cf_tfil  ) ; ijarg=ijarg+1
     CASE ('-s'        ) ; CALL getarg( ijarg, cf_sfil  ) ; ijarg=ijarg+1
     CASE ('-u'        ) ; CALL getarg( ijarg, cf_ufil  ) ; ijarg=ijarg+1
     CASE ('-v'        ) ; CALL getarg( ijarg, cf_vfil  ) ; ijarg=ijarg+1
     CASE ('-ssh'      ) ; CALL getarg( ijarg, cf_sshfil) ; ijarg=ijarg+1
     CASE ('-mh'       ) ; CALL getarg( ijarg, cf_mh    ) ; ijarg=ijarg+1
     CASE ('-mz'       ) ; CALL getarg( ijarg, cf_mz    ) ; ijarg=ijarg+1
     CASE ('-mask'     ) ; CALL getarg( ijarg, cf_mask  ) ; ijarg=ijarg+1
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


  !-- get dimensions (assuming all files have the same dimension that U-file) --
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
  PRINT *, 'dyn:hpg_sco : hydrostatic pressure gradient trend'


  !-- Allocate arrays --
  !- mesh -
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
  ALLOCATE( e3t_0(jpiglo, jpjglo, jpk)    , e3t(jpiglo, jpjglo, jpk)       )
  ALLOCATE( e3w_0(jpiglo, jpjglo, jpk)    , e3w(jpiglo, jpjglo, jpk)       )
  ALLOCATE( e3u_0( jpiglo, jpjglo, jpk)   , e3u(jpiglo, jpjglo, jpk)       )
  ALLOCATE( e3v_0(jpiglo, jpjglo, jpk)    , e3v(jpiglo, jpjglo, jpk)       )
  ALLOCATE( tmask(jpiglo, jpjglo, jpk)    , wmask(jpiglo, jpjglo, jpk)    )
  ALLOCATE( umask(jpiglo, jpjglo, jpk)    , vmask(jpiglo, jpjglo, jpk)     )
  ALLOCATE( gdept(jpiglo, jpjglo, jpk)    , gdepw(jpiglo, jpjglo, jpk)     , gdep3w(jpiglo, jpjglo, jpk) )
  !- variables -
  ALLOCATE( un(jpiglo, jpjglo, jpk)       , vn(jpiglo, jpjglo, jpk)         )
  ALLOCATE( rhd(jpiglo, jpjglo, jpk)                                        )
  ALLOCATE( sshn(jpiglo, jpjglo)                                            )
  ALLOCATE( ts(jpiglo, jpjglo, jpk, 2)                                      ) !1: temp, 2: sal
  !- outputs -
  ALLOCATE( zhpi(jpiglo, jpjglo, jpk)     , zhpj(jpiglo, jpjglo, jpk)       )
  ALLOCATE( zhpke(jpiglo, jpjglo, jpk)                                      )


  !-- define thermal expansion and haline contraction coefficients (polynomial TEOS-10, see eosbn2) --
         rdeltaS = 32._wp
         r1_S0  = 0.875_wp/35.16504_wp
         r1_T0  = 1._wp/40._wp
         r1_Z0  = 1.e-4_wp
!
         rau0_rcp    = rau0 * rcp
         r1_rau0     = 1._wp / rau0
         r1_rcp      = 1._wp / rcp
         r1_rau0_rcp = 1._wp / rau0_rcp
!
         EOS000 = 8.0189615746e+02_wp
         EOS100 = 8.6672408165e+02_wp
         EOS200 = -1.7864682637e+03_wp
         EOS300 = 2.0375295546e+03_wp
         EOS400 = -1.2849161071e+03_wp
         EOS500 = 4.3227585684e+02_wp
         EOS600 = -6.0579916612e+01_wp
         EOS010 = 2.6010145068e+01_wp
         EOS110 = -6.5281885265e+01_wp
         EOS210 = 8.1770425108e+01_wp
         EOS310 = -5.6888046321e+01_wp
         EOS410 = 1.7681814114e+01_wp
         EOS510 = -1.9193502195_wp
         EOS020 = -3.7074170417e+01_wp
         EOS120 = 6.1548258127e+01_wp
         EOS220 = -6.0362551501e+01_wp
         EOS320 = 2.9130021253e+01_wp
         EOS420 = -5.4723692739_wp
         EOS030 = 2.1661789529e+01_wp
         EOS130 = -3.3449108469e+01_wp
         EOS230 = 1.9717078466e+01_wp
         EOS330 = -3.1742946532_wp
         EOS040 = -8.3627885467_wp
         EOS140 = 1.1311538584e+01_wp
         EOS240 = -5.3563304045_wp
         EOS050 = 5.4048723791e-01_wp
         EOS150 = 4.8169980163e-01_wp
         EOS060 = -1.9083568888e-01_wp
         EOS001 = 1.9681925209e+01_wp
         EOS101 = -4.2549998214e+01_wp
         EOS201 = 5.0774768218e+01_wp
         EOS301 = -3.0938076334e+01_wp
         EOS401 = 6.6051753097_wp
         EOS011 = -1.3336301113e+01_wp
         EOS111 = -4.4870114575_wp
         EOS211 = 5.0042598061_wp
         EOS311 = -6.5399043664e-01_wp
         EOS021 = 6.7080479603_wp
         EOS121 = 3.5063081279_wp
         EOS221 = -1.8795372996_wp
         EOS031 = -2.4649669534_wp
         EOS131 = -5.5077101279e-01_wp
         EOS041 = 5.5927935970e-01_wp
         EOS002 = 2.0660924175_wp
         EOS102 = -4.9527603989_wp
         EOS202 = 2.5019633244_wp
         EOS012 = 2.0564311499_wp
         EOS112 = -2.1311365518e-01_wp
         EOS022 = -1.2419983026_wp
         EOS003 = -2.3342758797e-02_wp
         EOS103 = -1.8507636718e-02_wp
         EOS013 = 3.7969820455e-01_wp
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
         PEN000 = -9.8409626043_wp
         PEN100 = 2.1274999107e+01_wp
         PEN200 = -2.5387384109e+01_wp
         PEN300 = 1.5469038167e+01_wp
         PEN400 = -3.3025876549_wp
         PEN010 = 6.6681505563_wp
         PEN110 = 2.2435057288_wp
         PEN210 = -2.5021299030_wp
         PEN310 = 3.2699521832e-01_wp
         PEN020 = -3.3540239802_wp
         PEN120 = -1.7531540640_wp
         PEN220 = 9.3976864981e-01_wp
         PEN030 = 1.2324834767_wp
         PEN130 = 2.7538550639e-01_wp
         PEN040 = -2.7963967985e-01_wp
         PEN001 = -1.3773949450_wp
         PEN101 = 3.3018402659_wp
         PEN201 = -1.6679755496_wp
         PEN011 = -1.3709540999_wp
         PEN111 = 1.4207577012e-01_wp
         PEN021 = 8.2799886843e-01_wp
         PEN002 = 1.7507069098e-02_wp
         PEN102 = 1.3880727538e-02_wp
         PEN012 = -2.8477365341e-01_wp
!
         APE000 = -1.6670376391e-01_wp
         APE100 = -5.6087643219e-02_wp
         APE200 = 6.2553247576e-02_wp
         APE300 = -8.1748804580e-03_wp
         APE010 = 1.6770119901e-01_wp
         APE110 = 8.7657703198e-02_wp
         APE210 = -4.6988432490e-02_wp
         APE020 = -9.2436260751e-02_wp
         APE120 = -2.0653912979e-02_wp
         APE030 = 2.7963967985e-02_wp
         APE001 = 3.4273852498e-02_wp
         APE101 = -3.5518942529e-03_wp
         APE011 = -4.1399943421e-02_wp
         APE002 = 7.1193413354e-03_wp
!
         BPE000 = 2.6468936504e-01_wp
         BPE100 = -6.3170583896e-01_wp
         BPE200 = 5.7736640125e-01_wp
         BPE300 = -1.6435438140e-01_wp
         BPE010 = 2.7912203607e-02_wp
         BPE110 = -6.2259666565e-02_wp
         BPE210 = 1.2204769966e-02_wp
         BPE020 = -2.1811574876e-02_wp
         BPE120 = 2.3383950895e-02_wp
         BPE030 = 3.4261630030e-03_wp
         BPE001 = 4.1079296834e-02_wp
         BPE101 = -4.1503681096e-02_wp
         BPE011 = 1.7676120780e-03_wp
         BPE002 = 1.7269476440e-04_wp
!

  !!-- loading --
  PRINT *, '-- LOAD VARIABLES --'
  nav_lon_t    = getvar(cf_tfil, 'nav_lon', 1, jpiglo, jpjglo)
  nav_lat_t    = getvar(cf_tfil, 'nav_lat', 1, jpiglo, jpjglo)
  deptht       = getvar1d(cf_tfil, cn_vdeptht , jpk)
  nav_lon_u = getvar(cf_ufil, 'nav_lon', 1, jpiglo, jpjglo)
  nav_lat_u = getvar(cf_ufil, 'nav_lat', 1, jpiglo, jpjglo) 
  depthu    = getvar1d(cf_ufil, cn_vdepthu , jpk)
  nav_lon_v = getvar(cf_vfil, 'nav_lon', 1, jpiglo, jpjglo)
  nav_lat_v = getvar(cf_vfil, 'nav_lat', 1, jpiglo, jpjglo)
  depthv    = getvar1d(cf_vfil, cn_vdepthv , jpk)
  ht_0(:,:)   = getvar(cf_bathy, 'gdepw_0', 1, jpiglo, jpjglo )
  !- horizontal mesh -
  e1t(:,:)     = getvar(cf_mh  , 'e1t'  , 1, jpiglo, jpjglo)
  e2t(:,:)     = getvar(cf_mh  , 'e2t'  , 1, jpiglo, jpjglo)
  e1u(:,:)     = getvar(cf_mh  , 'e1u'  , 1, jpiglo, jpjglo)
  e2u(:,:)     = getvar(cf_mh  , 'e2u'  , 1, jpiglo, jpjglo)
  e1v(:,:)     = getvar(cf_mh  , 'e1v'  , 1, jpiglo, jpjglo)
  e2v(:,:)     = getvar(cf_mh  , 'e2v'  , 1, jpiglo, jpjglo)
  e12t(:,:)    = e1t(:,:) * e2t(:,:)
  r1_e12u(:,:) = 1._wp / (e1u(:,:) * e2u(:,:))
  r1_e12v(:,:) = 1._wp / (e1v(:,:) * e2v(:,:))
  !-- load vert. mesh (at rest) and masks --
  e3t_0(:,:,:) = getvar3d(cf_mz  , 'e3t_0' , jpiglo, jpjglo, jpk )
  e3w_0(:,:,:) = getvar3d(cf_mz  , 'e3w_0' , jpiglo, jpjglo, jpk )
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
  wmask(:,:,1) = tmask(:,:,1)
  DO jk=2,jpk
    wmask (:,:,jk) = tmask(:,:,jk) * tmask(:,:,jk-1)
  END DO
  
  !-- Creat output netcdf files to fill in --
  CALL CreateOutput


  DO jt = 1,jpt
     PRINT *, '======= time-step = ', jt

     !-- recomputed vert. metric due to non-linear free surface (VVL) --
     ! from domvvl.F90 -->> e3t = e3t_0*(1+sshn/ht_0)
     PRINT *, '-- Recompute vert. mesh --'
     sshn(:,:)     = getvar(cf_sshfil  , cn_sossheig, 1, jpiglo, jpjglo, ktime=jt )
     e3t(:,:,:)   = 0._wp
     e3w(:,:,:)   = 0._wp
     !DO jk = 1, jpk
     !  e3t(:,:,jk)   = e3t_0(:,:,jk) * (1 + sshn/ht_0)
     !  e3w(:,:,jk)   = e3w_0(:,:,jk) * (1 + sshn/ht_0)
     !END DO
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

     ! t- and w- points depth
     ! ----------------------
     ! set the isf depth as it is in the initial step
     gdept(:,:,:) = 0._wp
     gdepw(:,:,:) = 0._wp
     !
     gdept(:,:,1) = 0.5_wp * e3w(:,:,1)
     DO jk = 2, jpk
      DO jj = 1,jpjglo
       DO ji = 1,jpiglo
       !    zcoef = (tmask(ji,jj,jk) - wmask(ji,jj,jk))   ! 0 everywhere tmask = wmask, ie everywhere expect at jk = mikt
                                                          ! 1 everywhere from mbkt to mikt + 1 or 1 (if no isf)
                                                          ! 0.5 where jk = mikt
          zcoef = (tmask(ji,jj,jk) - wmask(ji,jj,jk))
          gdepw(ji,jj,jk) = gdepw(ji,jj,jk-1) + e3t(ji,jj,jk-1)
          gdept(ji,jj,jk) =      zcoef  * ( gdepw(ji,jj,jk) + 0.5 * e3w(ji,jj,jk))  &
              &                + (1-zcoef) * ( gdept(ji,jj,jk-1) +       e3w(ji,jj,jk))
       END DO
      END DO
     END DO
     !- add the SSH up/downlift to free surface here -
     DO jk = 1,jpk
        gdep3w(:,:,jk) = gdept(:,:,jk) - sshn(:,:)
     ENDDO


     !-- load T,S variables and compute in situ density --
     ! in step.F90
     ! CALL eos    ( tsn, rhd, rhop, gdept_n(:,:,:) ) ! now in situ density for hpg computation
     PRINT *, '-- Compute In-Situ Density --'
     ts(:,:,:,jp_tem) = getvar3d(cf_tfil, cn_votemper, jpiglo, jpjglo, jpk, ktime=jt )
     ts(:,:,:,jp_sal) = getvar3d(cf_sfil, cn_vosaline, jpiglo, jpjglo, jpk, ktime=jt )
     CALL eos_insitu    ( ts, rhd, gdept(:,:,:) )	! now in situ density


     !-- compute hydrostatic pressure gradient --
     PRINT *,'-- Compute hydrostatic pressure gradient momentum tendency --'
     CALL hpg_sco    ( jt )      ! s-coordinate (standard jacobian formulation)

     !-- Construct KE --
     PRINT *,'-- Compute KE --'
     un(:,:,:)   = getvar3d(cf_ufil, cn_vozocrtx, jpiglo, jpjglo, jpk, ktime=jt )
     vn(:,:,:)   = getvar3d(cf_vfil, cn_vomecrty, jpiglo, jpjglo, jpk, ktime=jt )
     !
     zhpke(:,:,:) = 0._wp
     CALL trd_ken( zhpi, zhpj, zhpke)

     !- write outputs trend -
     DO jk = 1, jpk
        ierr = putvar(ncout_u , id_varout_u(1) , zhpi(:,:,jk) , jk, jpiglo, jpjglo, ktime=jt )
        ierr = putvar(ncout_v , id_varout_v(1) , zhpj(:,:,jk) , jk, jpiglo, jpjglo, ktime=jt )
        ierr = putvar(ncout_ke, id_varout_ke(1), zhpke(:,:,jk), jk, jpiglo, jpjglo, ktime=jt )
     ENDDO

  ENDDO		!jt-loop

  ierr = closeout(ncout_u)
  ierr = closeout(ncout_v)
  ierr = closeout(ncout_ke)

CONTAINS

   SUBROUTINE hpg_sco( kt )
!!---------------------------------------------------------------------
!!                  ***  ROUTINE hpg_sco  ***
!!
!! ** Method  :   s-coordinate case. Jacobian scheme.
!!      The now hydrostatic pressure gradient at a given level, jk,
!!      is computed by taking the vertical integral of the in-situ
!!      density gradient along the model level from the suface to that
!!      level. s-coordinates (ln_sco): a corrective term is added
!!      to the horizontal pressure gradient :
!!         zhpi = grav .....  + 1/e1u mi(rhd) di[ grav dep3w ]
!!         zhpj = grav .....  + 1/e2v mj(rhd) dj[ grav dep3w ]
!!      add it to the general momentum trend (ua,va).
!!         ua = ua - 1/e1u * zhpi
!!         va = va - 1/e2v * zhpj
!!
!! ** Action : - Update (ua,va) with the now hydrastatic pressure trend
!!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt    ! ocean time-step index
!!
      INTEGER  ::   ji, jj, jk                 ! dummy loop indices
      REAL(wp) ::   zcoef0, zuap, zvap, znad   ! temporary scalars
      REAL(wp), POINTER, DIMENSION(:,:,:)     ::  tmpzhpi, tmpzhpj ! local hpg mom tendency (withou press. grad. corr.)
!!----------------------------------------------------------------------
!
       ALLOCATE( tmpzhpi(jpiglo, jpjglo, jpk)     , tmpzhpj(jpiglo, jpjglo, jpk)       )

! Local constant initialization
      zcoef0 = - grav * 0.5_wp
! To use density and not density anomaly
!      IF ( lk_vvl ) THEN   ;     
         znad = 1._wp          ! Variable volume
!      ELSE                 ;     znad = 0._wp         ! Fixed volume
!      ENDIF

! Surface value
      DO jj = 2, jpjm1
         DO ji = 2, jpim1   ! vector opt.
! hydrostatic pressure gradient along s-surfaces
            tmpzhpi(ji,jj,1) = zcoef0 / e1u(ji,jj) * ( e3w(ji+1,jj,1) * ( znad + rhd(ji+1,jj  ,1) )   &
               &                                  - e3w(ji,jj,1) * ( znad + rhd(ji  ,jj  ,1) ) )
            tmpzhpj(ji,jj,1) = zcoef0 / e2v(ji,jj) * ( e3w(ji,jj+1,1) * ( znad + rhd(ji  ,jj+1,1) )   &
               &                                  - e3w(ji,jj,1) * ( znad + rhd(ji  ,jj  ,1) ) )
! s-coordinate pressure gradient correction
            zuap = -zcoef0 * ( rhd   (ji+1,jj,1) + rhd   (ji,jj,1) + 2._wp * znad )   &
               &           * ( gdep3w(ji+1,jj,1) - gdep3w(ji,jj,1) ) / e1u(ji,jj)
            zvap = -zcoef0 * ( rhd   (ji,jj+1,1) + rhd   (ji,jj,1) + 2._wp * znad )   &
               &           * ( gdep3w(ji,jj+1,1) - gdep3w(ji,jj,1) ) / e2v(ji,jj)
! add to the general momentum trend
            zhpi(ji,jj,1) = ( tmpzhpi(ji,jj,1) + zuap ) * umask(ji,jj,1)
            zhpj(ji,jj,1) = ( tmpzhpj(ji,jj,1) + zvap ) * vmask(ji,jj,1)
         END DO
      END DO

! interior value (2=<jk=<jpkm1)
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
! hydrostatic pressure gradient along s-surfaces
               tmpzhpi(ji,jj,jk) = tmpzhpi(ji,jj,jk-1) + zcoef0 / e1u(ji,jj)   &
                  &           * (  e3w(ji+1,jj,jk) * ( rhd(ji+1,jj,jk) + rhd(ji+1,jj,jk-1) + 2*znad )   &
                  &              - e3w(ji,jj,jk) * ( rhd(ji  ,jj,jk) + rhd(ji  ,jj,jk-1) + 2*znad )  )
               tmpzhpj(ji,jj,jk) = tmpzhpj(ji,jj,jk-1) + zcoef0 / e2v(ji,jj)   &
                  &           * (  e3w(ji,jj+1,jk) * ( rhd(ji,jj+1,jk) + rhd(ji,jj+1,jk-1) + 2*znad )   &
                  &              - e3w(ji,jj,jk) * ( rhd(ji,jj,  jk) + rhd(ji,jj  ,jk-1) + 2*znad )  )
! s-coordinate pressure gradient correction
               zuap = -zcoef0 * ( rhd   (ji+1,jj  ,jk) + rhd   (ji,jj,jk) + 2._wp * znad )   &
                  &           * ( gdep3w(ji+1,jj,jk) - gdep3w(ji,jj,jk) ) / e1u(ji,jj)
               zvap = -zcoef0 * ( rhd   (ji  ,jj+1,jk) + rhd   (ji,jj,jk) + 2._wp * znad )   &
                  &           * ( gdep3w(ji,jj+1,jk) - gdep3w(ji,jj,jk) ) / e2v(ji,jj)
! add to the general momentum trend
               zhpi(ji,jj,jk) = ( tmpzhpi(ji,jj,jk) + zuap ) * umask(ji,jj,jk)
               zhpj(ji,jj,jk) = ( tmpzhpj(ji,jj,jk) + zvap ) * vmask(ji,jj,jk)
            END DO
         END DO
      END DO
!
      DEALLOCATE( tmpzhpi, tmpzhpj )
!
   END SUBROUTINE hpg_sco


   SUBROUTINE eos_insitu( pts, prd, pdep )
!!----------------------------------------------------------------------
!!                   ***  ROUTINE eos_insitu  ***
!!
!! ** Purpose :   Compute the in situ density (ratio rho/rau0) from
!!       potential temperature and salinity using an equation of state
!!       defined through the namelist parameter nn_eos.
!!
!! ** Method  :   prd(t,s,z) = ( rho(t,s,z) - rau0 ) / rau0
!!         with   prd    in situ density anomaly      no units
!!                t      TEOS10: CT or EOS80: PT      Celsius
!!                s      TEOS10: SA or EOS80: SP      TEOS10: g/kg or EOS80: psu
!!                z      depth                        meters
!!                rho    in situ density              kg/m^3
!!                rau0   reference density            kg/m^3
!!
!!     nn_eos = -1 : polynomial TEOS-10 equation of state is used for rho(t,s,z).
!!         Check value: rho = 1028.21993233072 kg/m^3 for z=3000 dbar, ct=3 Celcius, sa=35.5 g/kg
!!
!!     nn_eos =  0 : polynomial EOS-80 equation of state is used for rho(t,s,z).
!!         Check value: rho = 1028.35011066567 kg/m^3 for z=3000 dbar, pt=3 Celcius, sp=35.5 psu
!!
!!     nn_eos =  1 : simplified equation of state
!!              prd(t,s,z) = ( -a0*(1+lambda/2*(T-T0)+mu*z+nu*(S-S0))*(T-T0) + b0*(S-S0) ) / rau0
!!              linear case function of T only: rn_alpha<>0, other coefficients = 0
!!              linear eos function of T and S: rn_alpha and rn_beta<>0, other coefficients=0
!!              Vallis like equation: use default values of coefficients
!!
!! ** Action  :   compute prd , the in situ density (no units)
!!
!! References :   Roquet et al, Ocean Modelling, in preparation (2014)
!!                Vallis, Atmospheric and Oceanic Fluid Dynamics, 2006
!!                TEOS-10 Manual, 2010
!!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpiglo,jpjglo,jpk,jpts), INTENT(in   ) ::   pts   ! 1 : potential temperature  [Celcius]
!                                                               ! 2 : salinity               [psu]
      REAL(wp), DIMENSION(jpiglo,jpjglo,jpk     ), INTENT(  out) ::   prd   ! in situ density            [-]
      REAL(wp), DIMENSION(jpiglo,jpjglo,jpk     ), INTENT(in   ) ::   pdep  ! depth                      [m]
!
      INTEGER  ::   ji, jj, jk                ! dummy loop indices
      REAL(wp) ::   zt , zh , zs , ztm        ! local scalars
      REAL(wp) ::   zn , zn0, zn1, zn2, zn3   !   -      -
!!----------------------------------------------------------------------
!
!      SELECT CASE( nn_eos )
!
!      CASE( -1, 0 )                !==  polynomial TEOS-10 / EOS-80 ==!
!
         DO jk = 1, jpkm1
            DO jj = 1, jpjglo
               DO ji = 1, jpiglo
!
                  zh  = pdep(ji,jj,jk) * r1_Z0                                  ! depth
                  zt  = pts (ji,jj,jk,jp_tem) * r1_T0                           ! temperature
                  zs  = SQRT( ABS( pts(ji,jj,jk,jp_sal) + rdeltaS ) * r1_S0 )   ! square root salinity
                  ztm = tmask(ji,jj,jk)                                         ! tmask
!
                  zn3 = EOS013*zt   &
                     &   + EOS103*zs+EOS003
!
                  zn2 = (EOS022*zt   &
                     &   + EOS112*zs+EOS012)*zt   &
                     &   + (EOS202*zs+EOS102)*zs+EOS002
!
                  zn1 = (((EOS041*zt   &
                     &   + EOS131*zs+EOS031)*zt   &
                     &   + (EOS221*zs+EOS121)*zs+EOS021)*zt   &
                     &   + ((EOS311*zs+EOS211)*zs+EOS111)*zs+EOS011)*zt   &
                     &   + (((EOS401*zs+EOS301)*zs+EOS201)*zs+EOS101)*zs+EOS001
!
                  zn0 = (((((EOS060*zt   &
                     &   + EOS150*zs+EOS050)*zt   &
                     &   + (EOS240*zs+EOS140)*zs+EOS040)*zt   &
                     &   + ((EOS330*zs+EOS230)*zs+EOS130)*zs+EOS030)*zt   &
                     &   + (((EOS420*zs+EOS320)*zs+EOS220)*zs+EOS120)*zs+EOS020)*zt   &
                     &   + ((((EOS510*zs+EOS410)*zs+EOS310)*zs+EOS210)*zs+EOS110)*zs+EOS010)*zt   &
                     &   + (((((EOS600*zs+EOS500)*zs+EOS400)*zs+EOS300)*zs+EOS200)*zs+EOS100)*zs+EOS000
!
                  zn  = ( ( zn3 * zh + zn2 ) * zh + zn1 ) * zh + zn0
                  zn  = ( ( zn3 * zh + zn2 ) * zh + zn1 ) * zh + zn0
!
                  prd(ji,jj,jk) = (  zn * r1_rau0 - 1._wp  ) * ztm  ! density anomaly (masked)
!
               END DO
            END DO
         END DO
!
!      CASE( 1 )                !==  simplified EOS  ==!
!         ....
!      END SELECT
!
   END SUBROUTINE eos_insitu

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
!      REAL(wp)                                :: rau0 = 1026._wp    ! volumic mass of reference     [kg/m3] (from phycst.F90)
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
    stypvar(1)%cname             = 'hpg_u'
    stypvar(1)%cunits            = 'm/s^2'
    stypvar(1)%rmissing_value    = 99999.
    stypvar(1)%valid_min         = -1.
    stypvar(1)%valid_max         = 1.
    stypvar(1)%clong_name        = 'Hydrostatic pressure gradient, u-momentum'
    stypvar(1)%cshort_name       = 'hpg_u'
    stypvar(1)%conline_operation = 'On u-grid'
    stypvar(1)%caxis             = 'time depthu nav_lon_u nav_lat_u'

    stypvar2(1)%ichunk            = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
    stypvar2(1)%cname             = 'hpg_v'
    stypvar2(1)%cunits            = 'm/s^2'
    stypvar2(1)%rmissing_value    = 99999.
    stypvar2(1)%valid_min         = -1.
    stypvar2(1)%valid_max         = 1.
    stypvar2(1)%clong_name        = 'Hydrostatic pressure gradient, v-momentum'
    stypvar2(1)%cshort_name       = 'hpg_v'
    stypvar2(1)%conline_operation = 'On v-grid'
    stypvar2(1)%caxis             = 'time depthv nav_lon_v nav_lat_v'

    stypvar3(1)%ichunk            = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
    stypvar3(1)%cname             = 'hpg_ke'
    stypvar3(1)%cunits            = 'm^2/s^3'
    stypvar3(1)%rmissing_value    = 99999.
    stypvar3(1)%valid_min         = -1.
    stypvar3(1)%valid_max         = 1.
    stypvar3(1)%clong_name        = 'Hydrostatic pressure gradient, KE '
    stypvar3(1)%cshort_name       = 'hpg_ke'
    stypvar3(1)%conline_operation = 'On t-grid'
    stypvar3(1)%caxis             = 'time deptht nav_lon_t nav_lat_t'

    ! create output fileset
    ncout_u = create      (cf_out_u, cf_ufil ,  jpiglo, jpjglo, jpk, cdep=cn_vdepthu  , ld_nc4=lnc4 )
    ierr    = createvar   (ncout_u , stypvar ,  pnvarout, ipk , id_varout_u           , ld_nc4=lnc4 )
    ierr    = putheadervar(ncout_u , cf_ufil ,  jpiglo, jpjglo, jpk, nav_lon_u, nav_lat_u, depthu   )

    ncout_v = create      (cf_out_v, cf_vfil ,  jpiglo, jpjglo, jpk, cdep=cn_vdepthv  , ld_nc4=lnc4 )
    ierr    = createvar   (ncout_v , stypvar2,  pnvarout, ipk,    id_varout_v         , ld_nc4=lnc4 )
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
