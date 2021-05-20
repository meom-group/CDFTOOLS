PROGRAM cdf_dyndt
  !!======================================================================
  !!                     ***  PROGRAM  cdf_dyndt  ***
  !!=====================================================================
  !!  ** Purpose : Compute the time derivate of momentum with a leapfrog 
  !!               time stepping.
  !!
  !!  ** Method  : * Trends are computed using a leap-frog scheme:
  !!                       ztrdu = (ua-ub) / (2*rdt)
  !!                       ztrdv = (va-vb) / (2*rdt)
  !!               Note that with flux form advection and variable volume layer
  !!               (lk_vvl=T), the leap-frog is applied on thickness weighted
  !!               velocity.
  !!               Also note that the trends are computed based on TOTAL model
  !!               output velocities, thus including the Asselin filtering.
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
  INTEGER(KIND=4)                              :: ji, jj, jk, jt           ! dummy loop index
  INTEGER(KIND=4)                              :: jpim1, jpjm1, jpkm1      ! index - 1
  INTEGER(KIND=4)                              :: it                       ! time index for vvl
  INTEGER(KIND=4)                              :: ierr                     ! working integer
  INTEGER(KIND=4)                              :: narg, iargc, ijarg       ! 
  INTEGER(KIND=4)                              :: jpiglo, jpjglo, jpk, jpt ! size of the domain
  INTEGER(KIND=4)                              :: ncout_u, ncout_v, ncout_ke         ! ncid of output file
  INTEGER(KIND=4), DIMENSION(pnvarout)         :: ipk                      ! level of output variables
  INTEGER(KIND=4), DIMENSION(pnvarout)         :: id_varout_u              ! id of output variables (u-comp)
  INTEGER(KIND=4), DIMENSION(pnvarout)         :: id_varout_v              ! id of output variables (v-comp)
  INTEGER(KIND=4), DIMENSION(pnvarout)         :: id_varout_ke             ! id of output variables (ke-comp)

  !REAL(wp)                                     :: rdt=80._wp               ! model 'time step' [sec]
  REAL(wp)                                     :: rdt=3600._wp               ! model 'time step' [sec]
  REAL(wp), DIMENSION(:)    , ALLOCATABLE      :: dtim                     ! time
  REAL(wp), DIMENSION(:)    , ALLOCATABLE      :: deptht, depthu, depthv   ! z-grid (t,u,v)
  REAL(wp), DIMENSION(:,:)  , ALLOCATABLE      :: nav_lon_t, nav_lat_t     ! t-grid hor.
  REAL(wp), DIMENSION(:,:)  , ALLOCATABLE      :: nav_lon_u, nav_lat_u     ! u-grid hor.
  REAL(wp), DIMENSION(:,:)  , ALLOCATABLE      :: nav_lon_v, nav_lat_v     ! v-grid hor.
  !- vertical metrics recomputation -
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: ht_0                     ! Reference ocean depth at T-points
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: e1t, e2t                 ! horizontal metric, t-pts
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: e1u, e2u                 ! horizontal metric, u-pts
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: e1v, e2v                 ! horizontal metric, v-pts
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: e12t, r1_e12u, r1_e12v   ! face area at t-pts and inverse at u- v- pts
  REAL(wp), DIMENSION(:,:),   ALLOCATABLE      :: sshb, sshn, ssha         ! Sea surface height at before, now, after
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: e3t_0, e3u_0, e3v_0      ! vet. metrics at rest (without VVL)
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: e3t_b, e3u_b, e3v_b      ! vet. metrics at BEFORE
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: e3t_n, e3u_n, e3v_n      ! vet. metrics at NOW
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: e3t_a, e3u_a, e3v_a      ! vet. metrics at AFTER
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: tmask, umask, vmask      ! Mask at t-, u-, u- points
  !- velocities -
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: ub, vb                   ! BEFORE velocity
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: un, vn                   ! NOW    velocity
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: ua, va                   ! AFTER  velocity
  !- outputs -
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: dudt, dvdt, dkdt         ! momentum tendency (n=now) 

  CHARACTER(LEN=256)                           :: cf_tfil                  ! temperature netcdf file name (for mesh only)
  CHARACTER(LEN=256)                           :: cf_ufil_b, cf_vfil_b     ! BEFORE vel netcdf file name
  CHARACTER(LEN=256)                           :: cf_ufil_n, cf_vfil_n     ! NOW    vel netcdf file name
  CHARACTER(LEN=256)                           :: cf_ufil_a, cf_vfil_a     ! AFTER  vel netcdf file name
  CHARACTER(LEN=256)                           :: cf_sshfil_b              ! BEFORE ssh netcdf file name (for vvl)
  CHARACTER(LEN=256)                           :: cf_sshfil_n              ! NOW    ssh netcdf file name (for vvl)
  CHARACTER(LEN=256)                           :: cf_sshfil_a              ! AFTER  ssh netcdf file name (for vvl)
  CHARACTER(LEN=255)                           :: cf_mh                    ! horiz. mesh netcdf file name
  CHARACTER(LEN=255)                           :: cf_mz                    ! mesh      netcdf file name
  CHARACTER(LEN=255)                           :: cf_mask                  ! mask      netcdf file name
  CHARACTER(LEN=255)                           :: cf_bathy                 ! bathymetry  netcdf file name
  CHARACTER(LEN=256)                           :: cf_out_u='dt_u.nc'       ! output file name (u-comp)
  CHARACTER(LEN=256)                           :: cf_out_v='dt_v.nc'       ! output file name (v-comp)
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
     PRINT *,' usage : cdf_dyndt -t T-file -ub U-file-b -un U-file-n -ua U-file-a ...'
     PRINT *,'          -vb V-file-b -vn V-file-n -va V-file-a ...'
     PRINT *,'          -sshb SSH-file-b -sshn SSH-file-n -ssha SSH-file-a ...'
     PRINT *,'          -mh MESH-file -mz MESZ-file -mask MASK-file -bathy BATHY-file ...'
     PRINT *,'          -o_u OUT-file-u -o_v OUT-file-v -o_ke OUT-file-ke'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'      Compute the time derivative of momentum as dudt=(ua-ub)/2*dt (leapfrog).'
     PRINT *,'      Then compute KE time rate of change at tracer point'
     PRINT *,'      as rau0/2*(\overline{un*dudt}^i+\overline{vn*dvdt}^j), '
     PRINT *,'      where \overline represent averaging between two u-, v- neigbooring pvelocity points.'
     PRINT *,'      The U/V/SSH-file-b (U/V-file-a) files is used to compute momentum tendency at the '
     PRINT *,'      first (last) time record.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -t T-file          : netcdf file for temperature (for mesh only)'
     PRINT *,'       -ub U-file-b       : netcdf file for before zonal velocity'
     PRINT *,'       -ua U-file-a       : netcdf file for after   zonal velocity'
     PRINT *,'       -vb V-file-b       : netcdf file for before  meridional velocity'
     PRINT *,'       -va V-file-a       : netcdf file for after   meridional velocity'
     PRINT *,'       -sshb SSH-file-b   : netcdf file for before  SSH (for vvl recomputation of vert. grid)'
     PRINT *,'       -sshn SSH-file-n   : netcdf file for now     SSH (for vvl recomputation of vert. grid)'
     PRINT *,'       -ssha SSH-file-a   : netcdf file for after   SSH (for vvl recomputation of vert. grid)'
     PRINT *,'       -mh MESH-file      : netcdf file for horizontal mesh'
     PRINT *,'       -mz MESZ-file      : netcdf file for vertical mesh'
     PRINT *,'       -mask MASK-file    : netcdf file for mask'
     PRINT *,'       -bathy BATHY-file  : netcdf file for model bathymetry'
     PRINT *,'       -o_u OUT-file      : netcdf file for time derivative term for u-momentum'
     PRINT *,'       -o_v OUT-file      : netcdf file for time derivative term for v-momentum'
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
     CASE ('-t'        ) ; CALL getarg( ijarg, cf_tfil )     ; ijarg=ijarg+1
     CASE ('-ub'       ) ; CALL getarg( ijarg, cf_ufil_b   ) ; ijarg=ijarg+1
     CASE ('-un'       ) ; CALL getarg( ijarg, cf_ufil_n   ) ; ijarg=ijarg+1
     CASE ('-ua'       ) ; CALL getarg( ijarg, cf_ufil_a   ) ; ijarg=ijarg+1
     CASE ('-vb'       ) ; CALL getarg( ijarg, cf_vfil_b   ) ; ijarg=ijarg+1
     CASE ('-vn'       ) ; CALL getarg( ijarg, cf_vfil_n   ) ; ijarg=ijarg+1
     CASE ('-va'       ) ; CALL getarg( ijarg, cf_vfil_a   ) ; ijarg=ijarg+1
     CASE ('-sshb'     ) ; CALL getarg( ijarg, cf_sshfil_b ) ; ijarg=ijarg+1
     CASE ('-sshn'     ) ; CALL getarg( ijarg, cf_sshfil_n ) ; ijarg=ijarg+1
     CASE ('-ssha'     ) ; CALL getarg( ijarg, cf_sshfil_a ) ; ijarg=ijarg+1
     CASE ('-mh'       ) ; CALL getarg( ijarg, cf_mh       ) ; ijarg=ijarg+1
     CASE ('-mz'       ) ; CALL getarg( ijarg, cf_mz       ) ; ijarg=ijarg+1
     CASE ('-mask'     ) ; CALL getarg( ijarg, cf_mask     ) ; ijarg=ijarg+1
     CASE ('-bathy'    ) ; CALL getarg( ijarg, cf_bathy    ) ; ijarg=ijarg+1
        ! options
     CASE ( '-full'    ) ; lfull   = .TRUE. ; cglobal = 'full step computation'
     CASE ( '-o_u'     ) ; CALL getarg(ijarg, cf_out_u ) ; ijarg = ijarg + 1
     CASE ( '-o_v'     ) ; CALL getarg(ijarg, cf_out_v ) ; ijarg = ijarg + 1
     CASE ( '-o_ke'    ) ; CALL getarg(ijarg, cf_out_ke) ; ijarg = ijarg + 1
     CASE ( '-nc4'     ) ; lnc4    = .TRUE.
     CASE DEFAULT     ; PRINT *,' ERROR : ', TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  END DO

  !-- get dimensions (all files must have the same dimension that U-file) --
  jpiglo = getdim (cf_ufil_n, cn_x)
  jpjglo = getdim (cf_ufil_n, cn_y)
  jpk    = getdim (cf_ufil_n, cn_z)
  jpt    = getdim (cf_ufil_n, cn_t)
  jpim1  = jpiglo-1
  jpjm1  = jpjglo-1
  jpkm1  = jpk-1
  
  !-- summary --
  PRINT *, 'jpiglo =', jpiglo
  PRINT *, 'jpjglo =', jpjglo
  PRINT *, 'jpk    =', jpk
  PRINT *, 'jpt    =', jpt

  !-- Allocate --
  ALLOCATE( deptht(jpk)                  , depthu(jpk)                  , depthv(jpk)                   )
  ALLOCATE( nav_lon_t(jpiglo, jpjglo)    , nav_lat_t(jpiglo, jpjglo)     )
  ALLOCATE( nav_lon_u(jpiglo, jpjglo)    , nav_lat_u(jpiglo, jpjglo)     )
  ALLOCATE( nav_lon_v(jpiglo, jpjglo)    , nav_lat_v(jpiglo, jpjglo)     )
  !- vertical metrics recomputation -
  ALLOCATE( ht_0(jpiglo, jpjglo)                                         )
  ALLOCATE( e1t(jpiglo, jpjglo)          , e2t(jpiglo, jpjglo)           )
  ALLOCATE( e1u(jpiglo, jpjglo)          , e2u(jpiglo, jpjglo)           )
  ALLOCATE( e1v(jpiglo, jpjglo)          , e2v(jpiglo, jpjglo)           )
  ALLOCATE( e12t(jpiglo, jpjglo)         , r1_e12u(jpiglo, jpjglo)       , r1_e12v(jpiglo, jpjglo)    )
  ALLOCATE( sshb(jpiglo, jpjglo)         , sshn(jpiglo, jpjglo)          , ssha(jpiglo, jpjglo)       )
  ALLOCATE( e3t_0(jpiglo, jpjglo, jpk)   , e3u_0(jpiglo, jpjglo, jpk)    , e3v_0(jpiglo, jpjglo, jpk) )
  ALLOCATE( e3t_b(jpiglo, jpjglo, jpk)   , e3u_b(jpiglo, jpjglo, jpk)    , e3v_b(jpiglo, jpjglo, jpk) )
  ALLOCATE( e3t_n(jpiglo, jpjglo, jpk)   , e3u_n(jpiglo, jpjglo, jpk)    , e3v_n(jpiglo, jpjglo, jpk) )
  ALLOCATE( e3t_a(jpiglo, jpjglo, jpk)   , e3u_a(jpiglo, jpjglo, jpk)    , e3v_a(jpiglo, jpjglo, jpk) )
  ALLOCATE( tmask(jpiglo, jpjglo, jpk)   , umask(jpiglo, jpjglo, jpk)   , vmask(jpiglo, jpjglo, jpk)    )
  !- velocities -
  ALLOCATE( ub(jpiglo, jpjglo, jpk)      , vb(jpiglo, jpjglo, jpk)       )
  ALLOCATE( un(jpiglo, jpjglo, jpk)      , vn(jpiglo, jpjglo, jpk)       )
  ALLOCATE( ua(jpiglo, jpjglo, jpk)      , va(jpiglo, jpjglo, jpk)       )
  !- outputs -
  ALLOCATE( dudt(jpiglo, jpjglo, jpk)    , dvdt(jpiglo, jpjglo, jpk)     , dkdt(jpiglo, jpjglo, jpk) )

  !!-- loading -- 
  PRINT *, '-- LOAD VARIABLES --'
  nav_lon_t    = getvar(cf_tfil, 'nav_lon', 1, jpiglo, jpjglo)
  nav_lat_t    = getvar(cf_tfil, 'nav_lat', 1, jpiglo, jpjglo)
  deptht       = getvar1d(cf_tfil, cn_vdeptht , jpk)
  nav_lon_u = getvar(cf_ufil_n, 'nav_lon', 1, jpiglo, jpjglo)
  nav_lat_u = getvar(cf_ufil_n, 'nav_lat', 1, jpiglo, jpjglo)
  depthu    = getvar1d(cf_ufil_n, cn_vdepthu , jpk)
  nav_lon_v = getvar(cf_vfil_n, 'nav_lon', 1, jpiglo, jpjglo)
  nav_lat_v = getvar(cf_vfil_n, 'nav_lat', 1, jpiglo, jpjglo)
  depthv    = getvar1d(cf_vfil_n, cn_vdepthv , jpk)

  !-- For vertical mesh recomputation --
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

  DO jt = 1, jpt

     PRINT *, '======= time-step = ', jt 

     !-- recomputed vert. metric due to non-linear free surface (VVL) --
     ! from domvvl.F90 -->> e3t = e3t_0*(1+sshn/ht_0)
     PRINT *, '-- Recompute vert. mesh --'
     IF ( jt == 1 ) THEN        ! First time record of a file, look into before file (U/V-file-b)
        sshb(:,:)     = getvar(cf_sshfil_b  , cn_sossheig, 1, jpiglo, jpjglo, ktime=jpt  )
        ssha(:,:)     = getvar(cf_sshfil_n  , cn_sossheig, 1, jpiglo, jpjglo, ktime=jt+1 )
     ELSE IF ( jt == jpt ) THEN ! Last  time record of a file, look into after  file (U/V-file-a)
        sshb(:,:)     = getvar(cf_sshfil_n  , cn_sossheig, 1, jpiglo, jpjglo, ktime=jt-1 )
        ssha(:,:)     = getvar(cf_sshfil_a  , cn_sossheig, 1, jpiglo, jpjglo, ktime=1    )
     ELSE
        sshb(:,:)     = getvar(cf_sshfil_n  , cn_sossheig, 1, jpiglo, jpjglo, ktime=jt-1 )
        ssha(:,:)     = getvar(cf_sshfil_n  , cn_sossheig, 1, jpiglo, jpjglo, ktime=jt+1 )
     END IF
     sshn(:,:)     = getvar(cf_sshfil_n  , cn_sossheig, 1, jpiglo, jpjglo, ktime=jt   )
     DO jk = 1, jpk
       e3t_b(:,:,jk)   = e3t_0(:,:,jk) * (1 + sshb/ht_0)
       e3t_n(:,:,jk)   = e3t_0(:,:,jk) * (1 + sshn/ht_0)
       e3t_a(:,:,jk)   = e3t_0(:,:,jk) * (1 + ssha/ht_0)
     END DO
     !- at u- and v- pts (domvvl.F90) -
     e3u_b(:,:,:) = e3u_0(:,:,:)
     e3u_n(:,:,:) = e3u_0(:,:,:)
     e3u_a(:,:,:) = e3u_0(:,:,:)
     e3v_b(:,:,:) = e3v_0(:,:,:)
     e3v_n(:,:,:) = e3v_0(:,:,:)
     e3v_a(:,:,:) = e3v_0(:,:,:)
     DO jk = 1, jpk
        DO jj = 1, jpjm1
           DO ji = 1, jpim1   ! vector opt.
              e3u_b(ji,jj,jk) = e3u_0(ji,jj,jk) + 0.5_wp * umask(ji,jj,jk) * r1_e12u(ji,jj)                 &
                 &                       * (   e12t(ji  ,jj) * ( e3t_b(ji  ,jj,jk) - e3t_0(ji  ,jj,jk) )    &
                 &                           + e12t(ji+1,jj) * ( e3t_b(ji+1,jj,jk) - e3t_0(ji+1,jj,jk) ) )
              e3u_n(ji,jj,jk) = e3u_0(ji,jj,jk) + 0.5_wp * umask(ji,jj,jk) * r1_e12u(ji,jj)                 &
                 &                       * (   e12t(ji  ,jj) * ( e3t_n(ji  ,jj,jk) - e3t_0(ji  ,jj,jk) )    &
                 &                           + e12t(ji+1,jj) * ( e3t_n(ji+1,jj,jk) - e3t_0(ji+1,jj,jk) ) )
              e3u_a(ji,jj,jk) = e3u_0(ji,jj,jk) + 0.5_wp * umask(ji,jj,jk) * r1_e12u(ji,jj)                 &
                 &                       * (   e12t(ji  ,jj) * ( e3t_a(ji  ,jj,jk) - e3t_0(ji  ,jj,jk) )    &
                 &                           + e12t(ji+1,jj) * ( e3t_a(ji+1,jj,jk) - e3t_0(ji+1,jj,jk) ) )
              !
              e3v_b(ji,jj,jk) = e3v_0(ji,jj,jk) + 0.5_wp * vmask(ji,jj,jk) * r1_e12v(ji,jj)                 &
                 &                       * (   e12t(ji,jj  ) * ( e3t_b(ji,jj  ,jk) - e3t_0(ji,jj  ,jk) )    &
                 &                           + e12t(ji,jj+1) * ( e3t_b(ji,jj+1,jk) - e3t_0(ji,jj+1,jk) ) )
              e3v_n(ji,jj,jk) = e3v_0(ji,jj,jk) + 0.5_wp * vmask(ji,jj,jk) * r1_e12v(ji,jj)                 &
                 &                       * (   e12t(ji,jj  ) * ( e3t_n(ji,jj  ,jk) - e3t_0(ji,jj  ,jk) )    &
                 &                           + e12t(ji,jj+1) * ( e3t_n(ji,jj+1,jk) - e3t_0(ji,jj+1,jk) ) )
              e3v_a(ji,jj,jk) = e3v_0(ji,jj,jk) + 0.5_wp * vmask(ji,jj,jk) * r1_e12v(ji,jj)                 &
                 &                       * (   e12t(ji,jj  ) * ( e3t_a(ji,jj  ,jk) - e3t_0(ji,jj  ,jk) )    &
                 &                           + e12t(ji,jj+1) * ( e3t_a(ji,jj+1,jk) - e3t_0(ji,jj+1,jk) ) )
           END DO
        END DO
     END DO


     !-- Load velocities --
     IF ( jt == 1 ) THEN
        ub(:,:,:)   = getvar3d(cf_ufil_b, cn_vozocrtx, jpiglo, jpjglo, jpk, ktime=jpt  )
        vb(:,:,:)   = getvar3d(cf_vfil_b, cn_vomecrty, jpiglo, jpjglo, jpk, ktime=jpt  )
        ua(:,:,:)   = getvar3d(cf_ufil_n, cn_vozocrtx, jpiglo, jpjglo, jpk, ktime=jt+1 )
        va(:,:,:)   = getvar3d(cf_vfil_n, cn_vomecrty, jpiglo, jpjglo, jpk, ktime=jt+1 )
     ELSEIF( jt == jpt ) THEN    ! Last time record of a file, look into next file (U/V-file-a)
        ub(:,:,:)   = getvar3d(cf_ufil_n, cn_vozocrtx, jpiglo, jpjglo, jpk, ktime=jt-1 )
        vb(:,:,:)   = getvar3d(cf_vfil_n, cn_vomecrty, jpiglo, jpjglo, jpk, ktime=jt-1 )
        ua(:,:,:)   = getvar3d(cf_ufil_a, cn_vozocrtx, jpiglo, jpjglo, jpk, ktime=1    )
        va(:,:,:)   = getvar3d(cf_vfil_a, cn_vomecrty, jpiglo, jpjglo, jpk, ktime=1    ) 
     ELSE                        ! For all other time records of a file
        ub(:,:,:)   = getvar3d(cf_ufil_n, cn_vozocrtx, jpiglo, jpjglo, jpk, ktime=jt-1 )
        vb(:,:,:)   = getvar3d(cf_vfil_n, cn_vomecrty, jpiglo, jpjglo, jpk, ktime=jt-1 )
        ua(:,:,:)   = getvar3d(cf_ufil_n, cn_vozocrtx, jpiglo, jpjglo, jpk, ktime=jt+1 )
        va(:,:,:)   = getvar3d(cf_vfil_n, cn_vomecrty, jpiglo, jpjglo, jpk, ktime=jt+1 )
     ENDIF 
     un(:,:,:)      = getvar3d(cf_ufil_n, cn_vozocrtx, jpiglo, jpjglo, jpk, ktime=jt   )
     vn(:,:,:)      = getvar3d(cf_vfil_n, cn_vomecrty, jpiglo, jpjglo, jpk, ktime=jt   )

     !-- leapfrog time stepping --
     dudt(:,:,:) = ( ua(:,:,:) * e3u_a(:,:,:) - ub(:,:,:) * e3u_b(:,:,:) ) / &
                 &          ( 2 * rdt * e3u_n(:,:,:) )
     dvdt(:,:,:) = ( va(:,:,:) * e3v_a(:,:,:) - vb(:,:,:) * e3v_b(:,:,:) ) / &
                 &          ( 2 * rdt * e3v_n(:,:,:) )

     !-- Construct KE --
     dkdt(:,:,:) = 0._wp
     CALL trd_ken( dudt, dvdt, dkdt)


     !-- output trends --
     DO jk = 1, jpk
        ierr = putvar(ncout_u , id_varout_u(1) , dudt(:,:,jk), jk, jpiglo, jpjglo, ktime=jt )
        ierr = putvar(ncout_v , id_varout_v(1) , dvdt(:,:,jk), jk, jpiglo, jpjglo, ktime=jt )
        ierr = putvar(ncout_ke, id_varout_ke(1), dkdt(:,:,jk), jk, jpiglo, jpjglo, ktime=jt )
     ENDDO

  ENDDO		!jt-loop

  ierr = closeout(ncout_u)
  ierr = closeout(ncout_v)
  ierr = closeout(ncout_ke)

CONTAINS

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
            bu   (:,:,jk) =  e1u(:,:) * e2u(:,:) * e3u_n(:,:,jk)
            bv   (:,:,jk) =  e1v(:,:) * e2v(:,:) * e3v_n(:,:,jk)
            r1_bt(:,:,jk) = 1._wp / ( e12t(:,:) * e3t_n(:,:,jk) ) * tmask(:,:,jk)
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
    ipk(:)                        = jpk
    stypvar(1)%ichunk             = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
    stypvar(1)%cname              = 'dudt'
    stypvar(1)%cunits             = 'm/s^2'
    stypvar(1)%rmissing_value     = 99999.
    stypvar(1)%valid_min          = -1.
    stypvar(1)%valid_max          = 1.
    stypvar(1)%clong_name         = 'Time rate of change of u-momentum'
    stypvar(1)%cshort_name        = 'dudt'
    stypvar(1)%conline_operation  = 'On u-grid'
    stypvar(1)%caxis              = 'time depthu nav_lon_u nav_lat_u'

    stypvar2(1)%ichunk            = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
    stypvar2(1)%cname             = 'dvdt'
    stypvar2(1)%cunits            = 'm/s^2'
    stypvar2(1)%rmissing_value    = 99999.
    stypvar2(1)%valid_min         = -1.
    stypvar2(1)%valid_max         = 1.
    stypvar2(1)%clong_name        = 'Time rate of change of v-momentum'
    stypvar2(1)%cshort_name       = 'dvdt'
    stypvar2(1)%conline_operation = 'On v-grid'
    stypvar2(1)%caxis             = 'time depthv nav_lon_v nav_lat_v'

    stypvar3(1)%ichunk            = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
    stypvar3(1)%cname             = 'dkedt'
    stypvar3(1)%cunits            = 'm^2/s^3'
    stypvar3(1)%rmissing_value    = 99999.
    stypvar3(1)%valid_min         = -1.
    stypvar3(1)%valid_max         = 1.
    stypvar3(1)%clong_name        = 'Time rate of change of KE as rau0/2*(u*dudt + v*dvdt)'
    stypvar3(1)%cshort_name       = 'dkedt'
    stypvar3(1)%conline_operation = 'On t-grid'
    stypvar3(1)%caxis             = 'time deptht nav_lon_t nav_lat_t'


    ! create output fileset
    ncout_u = create      (cf_out_u, cf_ufil_n,  jpiglo, jpjglo, jpk, cdep=cn_vdepthu  , ld_nc4=lnc4 )
    ierr    = createvar   (ncout_u , stypvar  ,  pnvarout, ipk , id_varout_u           , ld_nc4=lnc4 )
    ierr    = putheadervar(ncout_u , cf_ufil_n,  jpiglo, jpjglo, jpk, nav_lon_u, nav_lat_u, depthu   )


    ncout_v = create      (cf_out_v, cf_vfil_n,  jpiglo, jpjglo, jpk, cdep=cn_vdepthv  , ld_nc4=lnc4 )
    ierr    = createvar   (ncout_v , stypvar2 ,  pnvarout, ipk , id_varout_v           , ld_nc4=lnc4 )
    ierr    = putheadervar(ncout_v , cf_vfil_n,  jpiglo, jpjglo, jpk, nav_lon_v, nav_lat_v, depthv   )

    ncout_ke= create      (cf_out_ke, cf_tfil ,  jpiglo, jpjglo, jpk, cdep=cn_vdeptht  , ld_nc4=lnc4 )
    ierr    = createvar   (ncout_ke , stypvar3,  pnvarout, ipk , id_varout_ke          , ld_nc4=lnc4 )
    ierr    = putheadervar(ncout_ke , cf_tfil ,  jpiglo, jpjglo, jpk, nav_lon_t, nav_lat_t, deptht   )

    dtim = getvar1d(cf_ufil_n, cn_vtimec,   jpt     )
    ierr = putvar1d(ncout_u  , dtim,        jpt, 'T')
    ierr = putvar1d(ncout_v  , dtim,        jpt, 'T')
    ierr = putvar1d(ncout_ke , dtim,        jpt, 'T')



  END SUBROUTINE CreateOutput


END PROGRAM
