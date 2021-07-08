PROGRAM cdf_dynadv_ubs 
  !!======================================================================
  !!                     ***  PROGRAM  cdf_dynadv_ubs  ***
  !!=====================================================================
  !!  ** Purpose : Compute momentum and KE advection trends following UBS advection scheme,
  !!               as well as the eddy/mean contributions (optional).
  !!
  !!  ** Method  : Adapt NEMO dynadv_ubs.F90 to CDFTOOLS (cf below for further details)
  !!               following the parameter used in the configuration eNATL60.
  !!
  !! History : 4.0  : 09/2019  : Q. Jamet & J.M. Molines : Original code
  !!                : 05/2021  : Q. Jamet & J.M. Molines : Turn the computation layer per layer
  !!                                                       to avoid memory issues.
  !!                : 06/2021  : Q. Jamet & J.M. Molines : Add eddy/mean decomposition
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
  INTEGER(KIND=4), PARAMETER                   :: jpnvarout1 = 2           ! number of output variables (uv)
  INTEGER(KIND=4)                              :: jpnvarout2               ! number of output variables (ke)
  INTEGER(KIND=4)                              :: ji, jj, jk, jt           ! dummy loop indices
  INTEGER(KIND=4)                              :: it                       ! time index for vvl
  INTEGER(KIND=4)                              :: ierr                     ! working integer
  INTEGER(KIND=4)                              :: narg, iargc, ijarg       ! 
  INTEGER(KIND=4)                              :: npiglo, npjglo, npk, npt ! size of the domain
  INTEGER(KIND=4)                              :: npim1, npjm1, npkm1      ! index for computation of grad/div
  INTEGER(KIND=4)                              :: nkkm1=1, nkk=2, nkkp1=3  ! for swapping the levels
  INTEGER(KIND=4)                              :: npkk=3                   ! number of level to load (nkkm1, nkk, nkkp1)
  INTEGER(KIND=4)                              :: ncout_u, ncout_v         ! ncid of output file
  INTEGER(KIND=4)                              :: ncout_ke                 ! ncid of output file
  INTEGER(KIND=4), DIMENSION(jpnvarout1)       :: ipk1                     ! level of output variables
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE   :: ipk2                     ! level of output variables
  INTEGER(KIND=4), DIMENSION(jpnvarout1)       :: id_varout_u              ! id of output variables (u-comp)
  INTEGER(KIND=4), DIMENSION(jpnvarout1)       :: id_varout_v              ! id of output variables (v-comp)
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE   :: id_varout_ke             ! id of output variables (ke-comp)

  REAL(wp)                                     :: pp_gamma1 = 1._wp/3._wp  ! =0 no dissipation ; =1/4 quick   ; =1/3  3rd order UBS
  REAL(wp), PARAMETER                          :: pp_gamma2 = 1._wp/32._wp ! =0   2nd order  ; =1/32 4th order centred
  REAL(wp), DIMENSION(:)    , ALLOCATABLE      :: deptht, depthu, depthv   ! z-grid (t, u,v)
  REAL(wp), DIMENSION(:,:)  , ALLOCATABLE      :: rlon_t, rlat_t           ! t-grid hor.
  REAL(wp), DIMENSION(:,:)  , ALLOCATABLE      :: rlon_u, rlat_u           ! u-grid hor.
  REAL(wp), DIMENSION(:,:)  , ALLOCATABLE      :: rlon_v, rlat_v           ! v-grid hor.
  REAL(wp), DIMENSION(:,:)  , ALLOCATABLE      :: ht_0                     ! Reference ocean depth at T-points
  REAL(wp), DIMENSION(:,:)  , ALLOCATABLE      :: sshn                     ! now sea surface height
  REAL(wp), DIMENSION(:,:)  , ALLOCATABLE      :: sshnm                    ! now sea surface height - mean
  REAL(wp), DIMENSION(:,:)  , ALLOCATABLE      :: e1t, e2t                 ! horizontal metric, t-pts
  REAL(wp), DIMENSION(:,:)  , ALLOCATABLE      :: e1u, e2u                 ! horizontal metric, u-pts
  REAL(wp), DIMENSION(:,:)  , ALLOCATABLE      :: e1v, e2v                 ! horizontal metric, v-pts
  REAL(wp), DIMENSION(:,:)  , ALLOCATABLE      :: e12t, r1_e12u, r1_e12v   ! face area at t-pts and inverse at u- v- pts
  REAL(wp), DIMENSION(:,:)  , ALLOCATABLE      :: e3t_0, e3u_0, e3v_0      ! vet. metrics at rest (without vvl)
  REAL(wp), DIMENSION(:,:)  , ALLOCATABLE      :: e3u, e3v, e3t            ! vet. metrics, u- v- t- pts
  REAL(wp), DIMENSION(:,:)  , ALLOCATABLE      :: tmask, fmask             ! mesh masks
  REAL(wp), DIMENSION(:,:)  , ALLOCATABLE      :: umask, vmask             ! mesh masks
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: wn                       ! 3D vert. velocity (now)
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: wnm                      ! 3D vert. velocity (now) - mean
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: un, vn                   ! 3D hz. velocity (now)
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: unm, vnm                 ! 3D hz. velocity (now) - mean
  REAL(wp), DIMENSION(:,:)  , ALLOCATABLE      :: adv_h_u, adv_z_u         ! hor. and vert. advection of u-mom. (outputs)
  REAL(wp), DIMENSION(:,:)  , ALLOCATABLE      :: adv_h_v, adv_z_v         ! hor. and vert. advection of v-mom. (outputs)
  REAL(wp), DIMENSION(:,:)  , ALLOCATABLE      :: adv_h_ke, adv_z_ke       ! hor. and vert. advection of KE     (outputs)

  CHARACTER(LEN=256)                           :: cf_tt                    ! temperature netcdf file name (for mesh only)
  CHARACTER(LEN=256)                           :: cf_uu                    ! zonal vel  netcdf file name
  CHARACTER(LEN=256)                           :: cf_um                    ! MEAN zonal vel  netcdf file name
  CHARACTER(LEN=256)                           :: cf_vv                    ! merd vel  netcdf file name
  CHARACTER(LEN=256)                           :: cf_vm                    ! MEAN merd vel  netcdf file name
  CHARACTER(LEN=256)                           :: cf_ww                    ! vert. vel  netcdf file name
  CHARACTER(LEN=256)                           :: cf_wm                    ! MEAN vert. vel  netcdf file name
  CHARACTER(LEN=256)                           :: cf_ssh                   ! Sea surface height
  CHARACTER(LEN=255)                           :: cf_mh                    ! mesh       netcdf file name
  CHARACTER(LEN=255)                           :: cf_mz                    ! mesh       netcdf file name
  CHARACTER(LEN=255)                           :: cf_mask                  ! mask       netcdf file name
  CHARACTER(LEN=255)                           :: cf_bathy                 ! bathymetry netcdf file name
  CHARACTER(LEN=256)                           :: cf_out_u ='adv_u.nc'     ! output file name (u-comp)
  CHARACTER(LEN=256)                           :: cf_out_v ='adv_v.nc'     ! output file name (v-comp)
  CHARACTER(LEN=256)                           :: cf_out_ke='adv_ke.nc'    ! output file name (ke-comp)
  CHARACTER(LEN=256)                           :: cldum                    ! dummy character variable
  CHARACTER(LEN=256)                           :: cglobal                  ! global attribute
  CHARACTER(LEN=256)                           :: eddymean='full'          ! select what to compute

  TYPE (variable), DIMENSION(jpnvarout1)       :: stypvar1                 ! structure for attibutes (u-comp)
  TYPE (variable), DIMENSION(jpnvarout1)       :: stypvar2                 ! structure for attibutes (v-comp)
  TYPE (variable), DIMENSION(:), ALLOCATABLE   :: stypvar3                 ! structure for attibutes (ke-comp)

  LOGICAL                                      :: l_w   =.FALSE.           ! flag for vertical location of bn2
  LOGICAL                                      :: lchk                     ! check missing files
  LOGICAL                                      :: lfull =.FALSE.           ! full step flag
  LOGICAL                                      :: lnc4  =.FALSE.           ! full step flag
  LOGICAL                                      :: nodiss  =.FALSE.         ! to remove dissipation term in UBS scheme
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdf_dynadv_ubs -t T-file -u U-file -v V-file -w W-file SSH-file ...'
     PRINT *,'          [-um Um-file -vm Vm-file -vm Vm-file] ...'
     PRINT *,'          -mh MESH-file -mz MESZ-file -mask MASK-file -bathy BATHY-file ...'
     PRINT *,'          [-o_u OUT-file-u -o_v OUT-file-v -o_ke OUT-file-ke ...'
     PRINT *,'          -em eddy-mean_option -nodiss]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'      Compute the momentum/KE advection trend following UBS advection scheme '
     PRINT *,'      of Shchepetkin & McWilliams (2005).'
     PRINT *,'      Options are provided for computing the eddy/mean decompositions (-em),'
     PRINT *,'      and remove the diffusive term included in this advective scheme (-nodiss).'
     PRINT *,'      The latter is done by turning pp_gamma1 = 0.0.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -t T-file          : netcdf file for      temperature (for mesh only)'
     PRINT *,'       -u U-file          : netcdf file for      zonal velocity'
     PRINT *,'       -v V-file          : netcdf file for      meridional velocity'
     PRINT *,'       -w W-file          : netcdf file for      vertical velocity'
     PRINT *,'       -ssh SSH-file      : netcdf file for      SSH (for vvl recomputation of vert. grid)'
     PRINT *,'       -mh MESH-file      : netcdf file for horizontal mesh'
     PRINT *,'       -mz MESZ-file      : netcdf file for vertical mesh'
     PRINT *,'       -mask MASK-file    : netcdf file for mask'
     PRINT *,'       -bathy BATHY-file  : netcdf file for model bathymetry'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       -em                : Eddy-mean option '
     PRINT *,'                            (full, mean-mean, eddy-mean, mean-eddy, eddy-eddy)'
     PRINT *,'                            (default: '//TRIM(eddymean)//')'
     PRINT *,'             --- For eddy-mean decomposition, no need if the option em=full'
     PRINT *,'       -um Um-file        : netcdf file for MEAN zonal velocity'
     PRINT *,'       -vm Vm-file        : netcdf file for MEAN meridional velocity'
     PRINT *,'       -wm Wm-file        : netcdf file for MEAN vertical velocity'
     PRINT *,'             ---  ' 
     PRINT *,'       -nodiss            : To remove dissipation term in the UBS advection scheme;'
     PRINT *,'                          :   =.TRUE. in case eddy/mean decomposition activated'
     PRINT *,'                          :   to insure balance between full and eddy/mean terms.'
     PRINT *,'       -nc4               : use netcdf4/HDF5 with chunking and deflation on output'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : '
     PRINT *,'       -o_u OUT-file      : netcdf file for advection term for u-momentum '
     PRINT *,'                                          (defaulf: '//TRIM(cf_out_u)//' )'
     PRINT *,'       -o_v OUT-file      : netcdf file for advection term for v-momentum '
     PRINT *,'                                          (defaulf: '//TRIM(cf_out_v)//' )'
     PRINT *,'       -o_ke OUT-file     : netcdf file for advection term for KE         '
     PRINT *,'                           (ubar*putrd + vbar*pvtrd) ; (upr*putrd + vpr*pvtrd)'
     PRINT *,'                                          (defaulf: '//TRIM(cf_out_ke)//')'
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
     CASE ('-t'        ) ; CALL getarg( ijarg, cf_tt )    ; ijarg=ijarg+1
     CASE ('-u'        ) ; CALL getarg( ijarg, cf_uu )    ; ijarg=ijarg+1
     CASE ('-v'        ) ; CALL getarg( ijarg, cf_vv )    ; ijarg=ijarg+1
     CASE ('-w'        ) ; CALL getarg( ijarg, cf_ww )    ; ijarg=ijarg+1
     CASE ('-ssh'      ) ; CALL getarg( ijarg, cf_ssh)    ; ijarg=ijarg+1
     CASE ('-um'       ) ; CALL getarg( ijarg, cf_um )    ; ijarg=ijarg+1
     CASE ('-vm'       ) ; CALL getarg( ijarg, cf_vm )    ; ijarg=ijarg+1
     CASE ('-wm'       ) ; CALL getarg( ijarg, cf_wm )    ; ijarg=ijarg+1
        !
     CASE ('-mh'       ) ; CALL getarg( ijarg, cf_mh   )  ; ijarg=ijarg+1
     CASE ('-mz'       ) ; CALL getarg( ijarg, cf_mz   )  ; ijarg=ijarg+1
     CASE ('-mask'     ) ; CALL getarg( ijarg, cf_mask )  ; ijarg=ijarg+1
     CASE ('-bathy'    ) ; CALL getarg( ijarg, cf_bathy ) ; ijarg=ijarg+1
        ! options
     CASE ( '-full' ) ; lfull   = .TRUE. ; cglobal = 'full step computation'
     CASE ( '-o_u'    ) ; CALL getarg(ijarg, cf_out_u )   ; ijarg = ijarg + 1
     CASE ( '-o_v'    ) ; CALL getarg(ijarg, cf_out_v )   ; ijarg = ijarg + 1
     CASE ( '-o_ke'   ) ; CALL getarg(ijarg, cf_out_ke)   ; ijarg = ijarg + 1
     CASE ( '-em'     ) ; CALL getarg(ijarg, eddymean  )  ; ijarg = ijarg + 1
     CASE ( '-nodiss' ) ; nodiss  = .TRUE.                ; ijarg = ijarg + 1
     CASE ( '-nc4'    ) ; lnc4    = .TRUE.
     CASE DEFAULT     ; PRINT *,' ERROR : ', TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  END DO
  IF ( eddymean .NE. 'full' ) THEN
     jpnvarout2=4
     nodiss  = .TRUE.           ! imposes no dissipation in eddy/mean decomposition
  ELSE
     jpnvarout2=2
  END IF
  !
  IF ( nodiss ) THEN
     pp_gamma1 = 0._wp
  END IF

  !-- get dimensions (all files must have the same dimension that U-file) --
  npiglo = getdim (cf_uu, cn_x)
  npjglo = getdim (cf_uu, cn_y)
  npk    = getdim (cf_uu, cn_z)
  npt    = getdim (cf_uu, cn_t)
  npim1  = npiglo-1
  npjm1  = npjglo-1
  npkm1  = npk-1
  
  !-- summary --
  PRINT *, 'npiglo =', npiglo
  PRINT *, 'npjglo =', npjglo
  PRINT *, 'npk    =', npk
  PRINT *, 'npt    =', npt
  !
  IF ( nodiss ) THEN
          PRINT*, '-- DO NOT compute dissipation --'
          PRINT*, '           -->>  pp_gamma1=  ', pp_gamma1
  ELSE
          PRINT*, '-- Compute dissipation --'
          PRINT*, '           -->>  pp_gamma1=  ', pp_gamma1
  END IF
  !
  SELECT CASE (eddymean)
  CASE ('full' ) ;
          PRINT *, 'Copmute FULL advection'
  CASE ('mean-mean') ;
          PRINT *, 'Compute MEAN-MEAN advection components'
  CASE ('mean-eddy') ;
          PRINT *, 'Compute MEAN-EDDY advection components'
  CASE ('eddy-mean') ;
          PRINT *, 'Compute EDDY-MEAN advection components'
  CASE ('eddy-eddy') ;
          PRINT *, 'Compute EDDY-EDDY advection components'
  CASE DEFAULT ; 
          PRINT *,' ERROR : ', TRIM(eddymean),' : unknown option.' ; STOP 99
  END SELECT

  !-- Allocate --
  ALLOCATE( ipk2(jpnvarout2)                                               )
  ALLOCATE( id_varout_ke(jpnvarout2)                                       )
  ALLOCATE( stypvar3(jpnvarout2)                                           )
  ! mesh
  ALLOCATE( deptht(npk)                   , depthu(npk)                    , depthv(npk)                    )
  ALLOCATE( rlon_t(npiglo, npjglo)        , rlat_t(npiglo, npjglo)         )
  ALLOCATE( rlon_u(npiglo, npjglo)        , rlat_u(npiglo, npjglo)         )
  ALLOCATE( rlon_v(npiglo, npjglo)        , rlat_v(npiglo, npjglo)         )
  ALLOCATE( ht_0(npiglo, npjglo)                                           )
  ALLOCATE( e1t(npiglo, npjglo)           , e2t(npiglo, npjglo)            )
  ALLOCATE( e1u(npiglo, npjglo)           , e2u(npiglo, npjglo)            )
  ALLOCATE( e1v(npiglo, npjglo)           , e2v(npiglo, npjglo)            )
  ALLOCATE( e12t(npiglo, npjglo)                                           )
  ALLOCATE( r1_e12u(npiglo, npjglo)       , r1_e12v(npiglo, npjglo)        )
  ALLOCATE( e3t_0(npiglo, npjglo)         , e3t(npiglo, npjglo)            )
  ALLOCATE( e3u_0( npiglo, npjglo)        , e3v_0(npiglo, npjglo)          )
  ALLOCATE( e3u(npiglo, npjglo)           , e3v(npiglo, npjglo)            )
  ALLOCATE( tmask(npiglo, npjglo)         , fmask(npiglo, npjglo)          )
  ALLOCATE( umask(npiglo, npjglo)         , vmask(npiglo, npjglo)          )
  !! variables
  ALLOCATE( sshn(npiglo, npjglo)                                           )
  ALLOCATE( un(npiglo, npjglo, npkk)      , vn(npiglo, npjglo, npkk)       )
  ALLOCATE( wn(npiglo, npjglo, npkk)                                       )
  IF ( eddymean .NE. 'full' ) THEN
     ALLOCATE( sshnm(npiglo, npjglo)                                       )
     ALLOCATE( unm(npiglo, npjglo, npkk)  , vnm(npiglo, npjglo, npkk)      )
     ALLOCATE( wnm(npiglo, npjglo, npkk)                                   )
  END IF 
  !
  ALLOCATE( adv_h_u(npiglo, npjglo)       , adv_z_u(npiglo, npjglo)        )
  ALLOCATE( adv_h_v(npiglo, npjglo)       , adv_z_v(npiglo, npjglo)        )
  ALLOCATE( adv_h_ke(npiglo, npjglo)      , adv_z_ke(npiglo, npjglo)       )


  !!-- loading -- 
  PRINT *, '-- LOAD VARIABLES --'
  rlon_t       = getvar(cf_tt   , cn_vlon2d, 1, npiglo, npjglo )
  rlat_t       = getvar(cf_tt   , cn_vlat2d, 1, npiglo, npjglo )
  deptht       = getvar1d(cf_tt , cn_vdeptht  , npk            )
  rlon_u       = getvar(cf_uu   , cn_vlon2d, 1, npiglo, npjglo )
  rlat_u       = getvar(cf_uu   , cn_vlat2d, 1, npiglo, npjglo )
  depthu       = getvar1d(cf_uu , cn_vdepthu  , npk            )
  rlon_v       = getvar(cf_vv   , cn_vlon2d, 1, npiglo, npjglo )
  rlat_v       = getvar(cf_vv   , cn_vlat2d, 1, npiglo, npjglo )
  depthv       = getvar1d(cf_vv , cn_vdepthv  , npk            )
  ht_0(:,:)    = getvar(cf_bathy, 'gdepw_0', 1, npiglo, npjglo )
  !- hz. mesh -
  e1t(:,:)     = getvar(cf_mh   , cn_ve1t  , 1, npiglo, npjglo )
  e2t(:,:)     = getvar(cf_mh   , cn_ve2t  , 1, npiglo, npjglo )
  e1u(:,:)     = getvar(cf_mh   , cn_ve1u  , 1, npiglo, npjglo )
  e2u(:,:)     = getvar(cf_mh   , cn_ve2u  , 1, npiglo, npjglo )
  e1v(:,:)     = getvar(cf_mh   , cn_ve1v  , 1, npiglo, npjglo )
  e2v(:,:)     = getvar(cf_mh   , cn_ve2v  , 1, npiglo, npjglo )
  e12t(:,:)    = e1t(:,:) * e2t(:,:)
  r1_e12u(:,:) = 1._wp / (e1u(:,:) * e2u(:,:))
  r1_e12v(:,:) = 1._wp / (e1v(:,:) * e2v(:,:))


  !-- Creat output netcdf files to fill in --
  PRINT *, '-- Creat output --'
  CALL CreateOutput

  DO jk = 1, 1
   PRINT *, '-- klayer: ', jk

   !-- load vert. mesh (at rest) and masks (dommsk.f90) --
   e3t_0(:,:) = getvar(cf_mz  , cn_ve3t0 , jk, npiglo, npjglo )
   e3u_0(:,:) = e3t_0(:,:)
   e3v_0(:,:) = e3t_0(:,:)
   DO jj = 1, npjm1
      DO ji = 1, npim1   ! vector opt.
         e3u_0 (ji,jj) = MIN( e3t_0(ji,jj), e3t_0(ji+1,jj) )
         e3v_0 (ji,jj) = MIN( e3t_0(ji,jj), e3t_0(ji,jj+1) )
      END DO
   END DO
   tmask(:,:) = getvar(cf_mask, cn_tmask , jk, npiglo, npjglo )
   umask(:,:) = getvar(cf_mask, cn_umask , jk, npiglo, npjglo )
   vmask(:,:) = getvar(cf_mask, cn_vmask , jk, npiglo, npjglo )
   !! fmask = 2 on lateral boundaries for no-slip bdy conditions on vorticity !!
   fmask(:,:) = getvar(cf_mask, cn_fmask , jk, npiglo, npjglo )

   DO jt = 1, npt
     sshn(:,:)   = 0._wp
     un(:,:,:)   = 0._wp
     vn(:,:,:)   = 0._wp
     wn(:,:,:)   = 0._wp
     IF ( eddymean .NE. 'full' ) THEN
        unm(:,:,:)   = 0._wp
        vnm(:,:,:)   = 0._wp
        wnm(:,:,:)   = 0._wp
     END IF
     ! 
     sshn(:,:)     = getvar(cf_ssh  , cn_sossheig, 1, npiglo, npjglo, ktime=jt )
     e3t(:,:) = e3t_0(:,:) * (1 + sshn/ht_0)
     !- at u- and v- pts (domvvl.F90) -
     e3u(:,:) = e3u_0(:,:)
     e3v(:,:) = e3v_0(:,:)
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

     !-- Load variables --
     IF ( jk == 1 ) THEN
        !- variable -
        un(:,:,nkk  ) = getvar(cf_uu, cn_vozocrtx, jk  , npiglo, npjglo, ktime=jt )
        vn(:,:,nkk  ) = getvar(cf_vv, cn_vomecrty, jk  , npiglo, npjglo, ktime=jt )
        un(:,:,nkk  ) = getvar(cf_uu, cn_vozocrtx, jk  , npiglo, npjglo, ktime=jt )
        vn(:,:,nkk  ) = getvar(cf_vv, cn_vomecrty, jk  , npiglo, npjglo, ktime=jt )
        wn(:,:,nkk  ) = getvar(cf_ww, cn_vovecrtz, jk  , npiglo, npjglo, ktime=jt )
        un(:,:,nkkp1) = getvar(cf_uu, cn_vozocrtx, jk+1, npiglo, npjglo, ktime=jt )
        vn(:,:,nkkp1) = getvar(cf_vv, cn_vomecrty, jk+1, npiglo, npjglo, ktime=jt )
        wn(:,:,nkkp1) = getvar(cf_ww, cn_vovecrtz, jk+1, npiglo, npjglo, ktime=jt )
        !- ensemble mean -
        IF ( eddymean .NE. 'full' ) THEN
           unm(:,:,nkk  ) = getvar(cf_um, cn_vozocrtx, jk  , npiglo, npjglo, ktime=jt )
           vnm(:,:,nkk  ) = getvar(cf_vm, cn_vomecrty, jk  , npiglo, npjglo, ktime=jt )
           wnm(:,:,nkk  ) = getvar(cf_wm, cn_vovecrtz, jk  , npiglo, npjglo, ktime=jt )
           unm(:,:,nkkp1) = getvar(cf_um, cn_vozocrtx, jk+1, npiglo, npjglo, ktime=jt )
           vnm(:,:,nkkp1) = getvar(cf_vm, cn_vomecrty, jk+1, npiglo, npjglo, ktime=jt )
           wnm(:,:,nkkp1) = getvar(cf_wm, cn_vovecrtz, jk+1, npiglo, npjglo, ktime=jt )
        END IF
     ELSE
        !- variable -
        un(:,:,nkkm1) = getvar(cf_uu, cn_vozocrtx, jk-1, npiglo, npjglo, ktime=jt )
        vn(:,:,nkkm1) = getvar(cf_vv, cn_vomecrty, jk-1, npiglo, npjglo, ktime=jt )
        wn(:,:,nkkm1) = getvar(cf_ww, cn_vovecrtz, jk-1, npiglo, npjglo, ktime=jt )
        un(:,:,nkk  ) = getvar(cf_uu, cn_vozocrtx, jk  , npiglo, npjglo, ktime=jt )
        vn(:,:,nkk  ) = getvar(cf_vv, cn_vomecrty, jk  , npiglo, npjglo, ktime=jt )
        wn(:,:,nkk  ) = getvar(cf_ww, cn_vovecrtz, jk  , npiglo, npjglo, ktime=jt )
        un(:,:,nkkp1) = getvar(cf_uu, cn_vozocrtx, jk+1, npiglo, npjglo, ktime=jt )
        vn(:,:,nkkp1) = getvar(cf_vv, cn_vomecrty, jk+1, npiglo, npjglo, ktime=jt )
        wn(:,:,nkkp1) = getvar(cf_ww, cn_vovecrtz, jk+1, npiglo, npjglo, ktime=jt )
        !- ensemble mean -
        IF ( eddymean .NE. 'full' ) THEN
           unm(:,:,nkkm1) = getvar(cf_um, cn_vozocrtx, jk-1, npiglo, npjglo, ktime=jt )
           vnm(:,:,nkkm1) = getvar(cf_vm, cn_vomecrty, jk-1, npiglo, npjglo, ktime=jt )
           wnm(:,:,nkkm1) = getvar(cf_wm, cn_vovecrtz, jk-1, npiglo, npjglo, ktime=jt )
           unm(:,:,nkk  ) = getvar(cf_um, cn_vozocrtx, jk  , npiglo, npjglo, ktime=jt )
           vnm(:,:,nkk  ) = getvar(cf_vm, cn_vomecrty, jk  , npiglo, npjglo, ktime=jt )
           wnm(:,:,nkk  ) = getvar(cf_wm, cn_vovecrtz, jk  , npiglo, npjglo, ktime=jt )
           unm(:,:,nkkp1) = getvar(cf_um, cn_vozocrtx, jk+1, npiglo, npjglo, ktime=jt )
           vnm(:,:,nkkp1) = getvar(cf_vm, cn_vomecrty, jk+1, npiglo, npjglo, ktime=jt )
           wnm(:,:,nkkp1) = getvar(cf_wm, cn_vovecrtz, jk+1, npiglo, npjglo, ktime=jt )
        END IF
     ENDIF

     !-- Advection and KE trends --
     SELECT CASE (eddymean)
     CASE ('full' ) ;
        CALL dyn_adv_ubs( jt, jk, un, vn, wn, un, vn )
        ierr = putvar(ncout_u , id_varout_u(1) , adv_h_u(:,:) , jk, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout_u , id_varout_u(2) , adv_z_u(:,:) , jk, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout_v , id_varout_v(1) , adv_h_v(:,:) , jk, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout_v , id_varout_v(2) , adv_z_v(:,:) , jk, npiglo, npjglo, ktime=jt )
        ! KE
        CALL trd_ken( adv_h_u, adv_h_v, adv_h_ke, un, vn )
        CALL trd_ken( adv_z_u, adv_z_v, adv_z_ke, un, vn )
        ierr = putvar(ncout_ke, id_varout_ke(1), adv_h_ke(:,:), jk, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout_ke, id_varout_ke(2), adv_z_ke(:,:), jk, npiglo, npjglo, ktime=jt )
     CASE ('mean-mean') ;
        CALL dyn_adv_ubs( jt, jk, unm, vnm, wnm, unm, vnm )
        ierr = putvar(ncout_u , id_varout_u(1) , adv_h_u(:,:) , jk, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout_u , id_varout_u(2) , adv_z_u(:,:) , jk, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout_v , id_varout_v(1) , adv_h_v(:,:) , jk, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout_v , id_varout_v(2) , adv_z_v(:,:) , jk, npiglo, npjglo, ktime=jt )
        ! KE  --  ubar*
        CALL trd_ken( adv_h_u, adv_h_v, adv_h_ke, unm, vnm )
        CALL trd_ken( adv_z_u, adv_z_v, adv_z_ke, unm, vnm )
        ierr = putvar(ncout_ke, id_varout_ke(1), adv_h_ke(:,:), jk, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout_ke, id_varout_ke(2), adv_z_ke(:,:), jk, npiglo, npjglo, ktime=jt )
        ! KE  --  uprime*
        CALL trd_ken( adv_h_u, adv_h_v, adv_h_ke, un(:,:,:)-unm(:,:,:), vn(:,:,:)-vnm(:,:,:) )
        CALL trd_ken( adv_z_u, adv_z_v, adv_z_ke, un(:,:,:)-unm(:,:,:), vn(:,:,:)-vnm(:,:,:) )
        ierr = putvar(ncout_ke, id_varout_ke(3), adv_h_ke(:,:), jk, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout_ke, id_varout_ke(4), adv_z_ke(:,:), jk, npiglo, npjglo, ktime=jt )
     CASE ('mean-eddy') ;
        CALL dyn_adv_ubs( jt, jk, unm, vnm, wnm, un(:,:,:)-unm(:,:,:), vn(:,:,:)-vnm(:,:,:) )
        ierr = putvar(ncout_u , id_varout_u(1) , adv_h_u(:,:) , jk, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout_u , id_varout_u(2) , adv_z_u(:,:) , jk, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout_v , id_varout_v(1) , adv_h_v(:,:) , jk, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout_v , id_varout_v(2) , adv_z_v(:,:) , jk, npiglo, npjglo, ktime=jt )
        ! KE  --  ubar*
        CALL trd_ken( adv_h_u, adv_h_v, adv_h_ke, unm, vnm )
        CALL trd_ken( adv_z_u, adv_z_v, adv_z_ke, unm, vnm )
        ierr = putvar(ncout_ke, id_varout_ke(1), adv_h_ke(:,:), jk, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout_ke, id_varout_ke(2), adv_z_ke(:,:), jk, npiglo, npjglo, ktime=jt )
        ! KE  --  uprime*
        CALL trd_ken( adv_h_u, adv_h_v, adv_h_ke, un(:,:,:)-unm(:,:,:), vn(:,:,:)-vnm(:,:,:) )
        CALL trd_ken( adv_z_u, adv_z_v, adv_z_ke, un(:,:,:)-unm(:,:,:), vn(:,:,:)-vnm(:,:,:) )
        ierr = putvar(ncout_ke, id_varout_ke(3), adv_h_ke(:,:), jk, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout_ke, id_varout_ke(4), adv_z_ke(:,:), jk, npiglo, npjglo, ktime=jt )
     CASE ('eddy-mean') ;
        CALL dyn_adv_ubs( jt, jk, un(:,:,:)-unm(:,:,:), vn(:,:,:)-vnm(:,:,:), &
                   & wn(:,:,:)-wnm(:,:,:), unm, vnm )
        ierr = putvar(ncout_u , id_varout_u(1) , adv_h_u(:,:) , jk, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout_u , id_varout_u(2) , adv_z_u(:,:) , jk, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout_v , id_varout_v(1) , adv_h_v(:,:) , jk, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout_v , id_varout_v(2) , adv_z_v(:,:) , jk, npiglo, npjglo, ktime=jt )
        ! KE  --  ubar*
        CALL trd_ken( adv_h_u, adv_h_v, adv_h_ke, unm, vnm )
        CALL trd_ken( adv_z_u, adv_z_v, adv_z_ke, unm, vnm )
        ierr = putvar(ncout_ke, id_varout_ke(1), adv_h_ke(:,:), jk, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout_ke, id_varout_ke(2), adv_z_ke(:,:), jk, npiglo, npjglo, ktime=jt )
        ! KE  --  uprime*
        CALL trd_ken( adv_h_u, adv_h_v, adv_h_ke, un(:,:,:)-unm(:,:,:), vn(:,:,:)-vnm(:,:,:) )
        CALL trd_ken( adv_z_u, adv_z_v, adv_z_ke, un(:,:,:)-unm(:,:,:), vn(:,:,:)-vnm(:,:,:) )
        ierr = putvar(ncout_ke, id_varout_ke(3), adv_h_ke(:,:), jk, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout_ke, id_varout_ke(4), adv_z_ke(:,:), jk, npiglo, npjglo, ktime=jt )
     CASE ('eddy-eddy') ;
        CALL dyn_adv_ubs( jt, jk, un(:,:,:)-unm(:,:,:), vn(:,:,:)-vnm(:,:,:), wn(:,:,:)-wnm(:,:,:), &
                   & un(:,:,:)-unm(:,:,:), vn(:,:,:)-vnm(:,:,:) )
        ierr = putvar(ncout_u , id_varout_u(1) , adv_h_u(:,:) , jk, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout_u , id_varout_u(2) , adv_z_u(:,:) , jk, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout_v , id_varout_v(1) , adv_h_v(:,:) , jk, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout_v , id_varout_v(2) , adv_z_v(:,:) , jk, npiglo, npjglo, ktime=jt )
        ! KE  --  ubar*
        CALL trd_ken( adv_h_u, adv_h_v, adv_h_ke, unm, vnm )
        CALL trd_ken( adv_z_u, adv_z_v, adv_z_ke, unm, vnm )
        ierr = putvar(ncout_ke, id_varout_ke(1), adv_h_ke(:,:), jk, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout_ke, id_varout_ke(2), adv_z_ke(:,:), jk, npiglo, npjglo, ktime=jt )
        ! KE  --  uprime*
        CALL trd_ken( adv_h_u, adv_h_v, adv_h_ke, un(:,:,:)-unm(:,:,:), vn(:,:,:)-vnm(:,:,:) )
        CALL trd_ken( adv_z_u, adv_z_v, adv_z_ke, un(:,:,:)-unm(:,:,:), vn(:,:,:)-vnm(:,:,:) )
        ierr = putvar(ncout_ke, id_varout_ke(3), adv_h_ke(:,:), jk, npiglo, npjglo, ktime=jt )
        ierr = putvar(ncout_ke, id_varout_ke(4), adv_z_ke(:,:), jk, npiglo, npjglo, ktime=jt )
     END SELECT

   ENDDO       !jt-loop
  ENDDO        !jk-loop
 
  ierr = closeout(ncout_u)
  ierr = closeout(ncout_v)
  ierr = closeout(ncout_ke)

CONTAINS

   SUBROUTINE dyn_adv_ubs( kt, kk, u1, v1, w1, u2, v2 )
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
      INTEGER, INTENT(in) ::   kk     ! ocean vertical level
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   u1, v1, w1   ! U, V and W ocean advecting velocities
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   u2, v2       ! U and V    ocean advected  velocities
!
      INTEGER  ::   ji, jj              ! dummy loop indices
      REAL(wp) ::   zbu, zbv    ! temporary scalars
      REAL(wp) ::   zui, zvj, zfuj, zfvi, zl_u, zl_v   ! temporary scalars
      REAL(wp), POINTER, DIMENSION(:,:    ) ::  zfu, zfv
      REAL(wp), POINTER, DIMENSION(:,:,:  ) ::  zfw
      REAL(wp), POINTER, DIMENSION(:,:    ) ::  zfu_t, zfv_t, zfu_f, zfv_f
      REAL(wp), POINTER, DIMENSION(:,:,:  ) ::  zfu_uw, zfv_vw
      REAL(wp), POINTER, DIMENSION(:,:,:  ) ::  zlu_uu, zlv_vv, zlu_uv, zlv_vu
!!----------------------------------------------------------------------
!
      ALLOCATE( zfu(npiglo, npjglo)          , zfv(npiglo, npjglo)          )
      ALLOCATE( zfw(npiglo, npjglo, npkk)                                   )
      ALLOCATE( zfu_t(npiglo, npjglo)        , zfv_t(npiglo, npjglo)        )
      ALLOCATE( zfu_f(npiglo, npjglo)        , zfv_f(npiglo, npjglo)        )
      ALLOCATE( zfu_uw(npiglo, npjglo, npkk) , zfv_vw(npiglo, npjglo, npkk) )
      !
      ALLOCATE( zlu_uu(npiglo, npjglo, 2)    , zlv_vv(npiglo, npjglo, 2)    )
      ALLOCATE( zlu_uv(npiglo, npjglo, 2)    , zlv_vu(npiglo, npjglo, 2)    )

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
                                       !  Laplacian of the velocity  !
!                                      ! =========================== !
!                                         ! horizontal volume fluxes
         zfu(:,:) = e2u(:,:) * e3u(:,:) * u1(:,:,nkk)
         zfv(:,:) = e1v(:,:) * e3v(:,:) * v1(:,:,nkk)
!
         DO jj = 2, npjm1                          ! laplacian
            DO ji = 2, npim1   ! vector opt.
!
               zlu_uu(ji,jj,1) = ( u2 (ji+1,jj  ,nkk) - 2.*u2 (ji,jj,nkk) + u2 (ji-1,jj  ,nkk) ) * umask(ji,jj)
               zlv_vv(ji,jj,1) = ( v2 (ji  ,jj+1,nkk) - 2.*v2 (ji,jj,nkk) + v2 (ji  ,jj-1,nkk) ) * vmask(ji,jj)
               zlu_uv(ji,jj,1) = ( u2 (ji  ,jj+1,nkk) - u2 (ji  ,jj  ,nkk) ) * fmask(ji  ,jj  )   &
                  &               - ( u2 (ji  ,jj  ,nkk) - u2 (ji  ,jj-1,nkk) ) * fmask(ji  ,jj-1)
               zlv_vu(ji,jj,1) = ( v2 (ji+1,jj  ,nkk) - v2 (ji  ,jj  ,nkk) ) * fmask(ji  ,jj  )   &
                  &               - ( v2 (ji  ,jj  ,nkk) - v2 (ji-1,jj  ,nkk) ) * fmask(ji-1,jj  )
!
               zlu_uu(ji,jj,2) = ( zfu(ji+1,jj  ) - 2.*zfu(ji,jj) + zfu(ji-1,jj  ) ) * umask(ji,jj)
               zlv_vv(ji,jj,2) = ( zfv(ji  ,jj+1) - 2.*zfv(ji,jj) + zfv(ji  ,jj-1) ) * vmask(ji,jj)
               zlu_uv(ji,jj,2) = ( zfu(ji  ,jj+1) - zfu(ji  ,jj  ) ) * fmask(ji  ,jj  )   &
                  &               - ( zfu(ji  ,jj  ) - zfu(ji  ,jj-1) ) * fmask(ji  ,jj-1)
               zlv_vu(ji,jj,2) = ( zfv(ji+1,jj  ) - zfv(ji  ,jj  ) ) * fmask(ji  ,jj  )   &
                  &               - ( zfv(ji  ,jj  ) - zfv(ji-1,jj  ) ) * fmask(ji-1,jj  )
            END DO
         END DO

!                                      ! ====================== !
!                                      !  Horizontal advection  !
                                       ! ====================== !
!                                         ! horizontal volume fluxes
         zfu(:,:) = 0.25 * e2u(:,:) * e3u(:,:) * u1(:,:,nkk)
         zfv(:,:) = 0.25 * e1v(:,:) * e3v(:,:) * v1(:,:,nkk)
!
         DO jj = 1, npjm1                          ! horizontal momentum fluxes at T- and F-point
            DO ji = 1, npim1   ! vector opt.
               zui = ( u2(ji,jj,nkk) + u2(ji+1,jj  ,nkk) )
               zvj = ( v2(ji,jj,nkk) + v2(ji  ,jj+1,nkk) )
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
                  &                * ( u2(ji,jj,nkk) + u2(ji  ,jj+1,nkk) - pp_gamma1 * zl_u )
               zfu_f(ji  ,jj  ) = ( zfuj - pp_gamma2 * ( zlu_uv(ji,jj,2) + zlu_uv(ji  ,jj+1,2) )  )   &
                  &                * ( v2(ji,jj,nkk) + v2(ji+1,jj  ,nkk) - pp_gamma1 * zl_v )
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

!                                      ! ==================== !
!                                      !  Vertical advection  !
                                       ! ==================== !
!                                         ! Vertical volume fluxes
         zfw(:,:,nkk)   = 0.25 * e1t(:,:) * e2t(:,:) * w1(:,:,nkk  )
         zfw(:,:,nkkp1) = 0.25 * e1t(:,:) * e2t(:,:) * w1(:,:,nkkp1)
!
   !      IF( kk == 1 ) THEN                        ! surface/bottom advective fluxes
   !       ... moved after interior fluxes ...
   !      ELSE                                      ! interior fluxes
            DO jj = 2, npjm1
               DO ji = 2, npim1   ! vector opt.
                  zfu_uw(ji,jj,nkk)   = ( zfw(ji,jj,nkk)  + zfw(ji+1,jj  ,nkk) )   * ( u2(ji,jj,nkk)   + u2(ji,jj,nkkm1) )
                  zfv_vw(ji,jj,nkk)   = ( zfw(ji,jj,nkk)  + zfw(ji  ,jj+1,nkk) )   * ( v2(ji,jj,nkk)   + v2(ji,jj,nkkm1) )
                  zfu_uw(ji,jj,nkkp1) = ( zfw(ji,jj,nkkp1)+ zfw(ji+1,jj  ,nkkp1) ) * ( u2(ji,jj,nkkp1) + u2(ji,jj,nkk)   )
                  zfv_vw(ji,jj,nkkp1) = ( zfw(ji,jj,nkkp1)+ zfw(ji  ,jj+1,nkkp1) ) * ( v2(ji,jj,nkkp1) + v2(ji,jj,nkk)   )
               END DO
            END DO
    !     ENDIF
         IF( kk == npkm1 ) THEN
            zfu_uw(:,:,nkkp1) = 0.e0                      ! Bottom  value : flux set to zero
            zfv_vw(:,:,nkkp1) = 0.e0
         ENDIF
         IF ( kk == 1 ) THEN                        ! Surface value :
!            IF( lk_vvl ) THEN                                ! variable volume : flux set to zero
               zfu_uw(:,:, nkk ) = 0.e0
               zfv_vw(:,:, nkk ) = 0.e0
!            ELSE                                             ! constant volume : advection through the surface
!               ...
!            ENDIF
         ENDIF

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
!
      DEALLOCATE( zfu_t , zfv_t  )
      DEALLOCATE( zfu_f , zfv_f  )
      DEALLOCATE( zfu   , zfv    )
      DEALLOCATE( zfw            )
      DEALLOCATE( zlu_uu, zlv_vv )
      DEALLOCATE( zlu_uv, zlv_vu )

   END SUBROUTINE dyn_adv_ubs


   SUBROUTINE trd_ken( putrd, pvtrd, ktrd, u0, v0 )
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
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   u0, v0       ! U and V 
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   ktrd           ! KE trend
!
      INTEGER ::   ji, jj       ! dummy loop indices
!
      REAL(wp)                                :: rau0 = 1026._wp    ! volumic mass of reference     [kg/m3] (from phycst.F90)
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   bu, bv   ! volume of u- and v-boxes
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   r1_bt    ! inverse of t-box volume
!!----------------------------------------------------------------------
!
!
      ALLOCATE( bu(npiglo,npjglo) , bv(npiglo,npjglo) , r1_bt(npiglo,npjglo) )
!
      bu   (:,:) =  e1u(:,:) * e2u(:,:) * e3u(:,:)
      bv   (:,:) =  e1v(:,:) * e2v(:,:) * e3v(:,:)
      r1_bt(:,:) = 1._wp / ( e12t(:,:) * e3t(:,:) ) * tmask(:,:)
!
      ktrd(:,:) = 0._wp
      ktrd(:,:) = 0._wp
      DO jj = 2, npjglo
         DO ji = 2, npiglo
            ktrd(ji,jj) = 0.5_wp * rau0 *( u0(ji  ,jj,nkk) * putrd(ji  ,jj) * bu(ji  ,jj)  &
               &                            + u0(ji-1,jj,nkk) * putrd(ji-1,jj) * bu(ji-1,jj)  &
               &                            + v0(ji,jj  ,nkk) * pvtrd(ji,jj  ) * bv(ji,jj  )  &
               &                            + v0(ji,jj-1,nkk) * pvtrd(ji,jj-1) * bv(ji,jj-1)  ) * r1_bt(ji,jj)
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
    REAL(KIND=8), DIMENSION(npt) :: dltim
    ! define new variables for output
    ipk1(:)                        = npk
    stypvar1(1)%ichunk             = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar1(1)%cname              = 'advh_uu'
    stypvar1(1)%cunits             = 'm/s^2'
    stypvar1(1)%rmissing_value     = 99999.
    stypvar1(1)%valid_min          = -1.
    stypvar1(1)%valid_max          = 1.
    stypvar1(1)%clong_name         = 'Horizontal advection of zonal momentum ('//TRIM(eddymean)//')'
    stypvar1(1)%cshort_name        = 'advh_uu'
    stypvar1(1)%conline_operation  = 'On u-grid'
    stypvar1(1)%caxis              = 'time depthu nav_lon_u nav_lat_u'
    !
    stypvar1(2)%ichunk             = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar1(2)%cname              = 'advz_uu'
    stypvar1(2)%cunits             = 'm/s^2'
    stypvar1(2)%rmissing_value     = 99999.
    stypvar1(2)%valid_min          = -1.
    stypvar1(2)%valid_max          = 1.
    stypvar1(2)%clong_name         = 'Vertical advection of zonal momentum ('//TRIM(eddymean)//')'
    stypvar1(2)%cshort_name        = 'advz_uu'
    stypvar1(2)%conline_operation  = 'On u-grid'
    stypvar1(2)%caxis              = 'time depthu nav_lon_u nav_lat_u'

    stypvar2(1)%ichunk             = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar2(1)%cname              = 'advh_vv'
    stypvar2(1)%cunits             = 'm/s^2'
    stypvar2(1)%rmissing_value     = 99999.
    stypvar2(1)%valid_min          = -1.
    stypvar2(1)%valid_max          = 1.
    stypvar2(1)%clong_name         = 'Horizontal advection of meridional momentum ('//TRIM(eddymean)//')'
    stypvar2(1)%cshort_name        = 'advh_vv'
    stypvar2(1)%conline_operation  = 'On v-grid'
    stypvar2(1)%caxis              = 'time depthv nav_lon_v nav_lat_v'
    !
    stypvar2(2)%ichunk             = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar2(2)%cname              = 'advz_vv'
    stypvar2(2)%cunits             = 'm/s^2'
    stypvar2(2)%rmissing_value     = 99999.
    stypvar2(2)%valid_min          = -1.
    stypvar2(2)%valid_max          = 1.
    stypvar2(2)%clong_name         = 'Vertical advection of meridional momentum ('//TRIM(eddymean)//')'
    stypvar2(2)%cshort_name        = 'advz_vv'
    stypvar2(2)%conline_operation  = 'On v-grid'
    stypvar2(2)%caxis              = 'time depthv nav_lon_v nav_lat_v'

    ipk2(:)                        = npk
    IF ( eddymean .EQ. 'full ') THEN
    stypvar3(1)%ichunk             = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar3(1)%cname              = 'advh_ke'
    stypvar3(1)%cunits             = 'm^2/s^3'
    stypvar3(1)%rmissing_value     = 99999.
    stypvar3(1)%valid_min          = -1.
    stypvar3(1)%valid_max          = 1.
    stypvar3(1)%clong_name         = 'Horizontal advection of Kinetic Energy ('//TRIM(eddymean)//')'
    stypvar3(1)%cshort_name        = 'advh_ke'
    stypvar3(1)%conline_operation  = 'On t-grid'
    stypvar3(1)%caxis              = 'time deptht nav_lon_t nav_lat_t'
    !
    stypvar3(2)%ichunk             = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar3(2)%cname              = 'advz_ke'
    stypvar3(2)%cunits             = 'm^2/s^3'
    stypvar3(2)%rmissing_value     = 99999.
    stypvar3(2)%valid_min          = -1.
    stypvar3(2)%valid_max          = 1.
    stypvar3(2)%clong_name         = 'Vertical advection of Kinetic Energy ('//TRIM(eddymean)//')'
    stypvar3(2)%cshort_name        = 'advz_ke'
    stypvar3(2)%conline_operation  = 'On t-grid'
    stypvar3(2)%caxis              = 'time deptht nav_lon_t nav_lat_t'
    ELSE
    stypvar3(1)%ichunk             = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar3(1)%cname              = 'advh_ke_m'
    stypvar3(1)%cunits             = 'm^2/s^3'
    stypvar3(1)%rmissing_value     = 99999.
    stypvar3(1)%valid_min          = -1.
    stypvar3(1)%valid_max          = 1.
    stypvar3(1)%clong_name         = 'um * advh_uu + vm * advh_vv ('//TRIM(eddymean)//')'
    stypvar3(1)%cshort_name        = 'advh_ke_m'
    stypvar3(1)%conline_operation  = 'On t-grid'
    stypvar3(1)%caxis              = 'time deptht nav_lon_t nav_lat_t'
    !
    stypvar3(2)%ichunk             = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar3(2)%cname              = 'advz_ke_m'
    stypvar3(2)%cunits             = 'm^2/s^3'
    stypvar3(2)%rmissing_value     = 99999.
    stypvar3(2)%valid_min          = -1.
    stypvar3(2)%valid_max          = 1.
    stypvar3(2)%clong_name         = 'um * advz_uu + vm * advz_vv ('//TRIM(eddymean)//')'
    stypvar3(2)%cshort_name        = 'advz_ke_m'
    stypvar3(2)%conline_operation  = 'On t-grid'
    stypvar3(2)%caxis              = 'time deptht nav_lon_t nav_lat_t'
    !
    stypvar3(3)%ichunk             = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar3(3)%cname              = 'advh_ke_pr'
    stypvar3(3)%cunits             = 'm^2/s^3'
    stypvar3(3)%rmissing_value     = 99999.
    stypvar3(3)%valid_min          = -1.
    stypvar3(3)%valid_max          = 1.
    stypvar3(3)%clong_name         = 'upr * advh_uu + vpr * advh_vv ('//TRIM(eddymean)//')'
    stypvar3(3)%cshort_name        = 'advh_ke_pr'
    stypvar3(3)%conline_operation  = 'On t-grid'
    stypvar3(3)%caxis              = 'time deptht nav_lon_t nav_lat_t'
    !
    stypvar3(4)%ichunk             = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar3(4)%cname              = 'advz_ke_pr'
    stypvar3(4)%cunits             = 'm^2/s^3'
    stypvar3(4)%rmissing_value     = 99999.
    stypvar3(4)%valid_min          = -1.
    stypvar3(4)%valid_max          = 1.
    stypvar3(4)%clong_name         = 'upr * advz_uu + vpr * advz_vv ('//TRIM(eddymean)//')'
    stypvar3(4)%cshort_name        = 'advz_ke_pr'
    stypvar3(4)%conline_operation  = 'On t-grid'
    stypvar3(4)%caxis              = 'time deptht nav_lon_t nav_lat_t'
    END IF

    ! create output fileset
    ncout_u = create      (cf_out_u, cf_uu   , npiglo, npjglo, npk  , cdep=cn_vdepthu        , ld_nc4=lnc4 )
    ierr    = createvar   (ncout_u , stypvar1, jpnvarout1    , ipk1 , id_varout_u            , ld_nc4=lnc4 )
    ierr    = putheadervar(ncout_u , cf_uu   , npiglo, npjglo, npk, rlon_u, rlat_u, depthu   )

    ncout_v = create      (cf_out_v, cf_vv   , npiglo, npjglo, npk  , cdep=cn_vdepthv        , ld_nc4=lnc4 )
    ierr    = createvar   (ncout_v , stypvar2, jpnvarout1    , ipk1 , id_varout_v            , ld_nc4=lnc4 )
    ierr    = putheadervar(ncout_v , cf_vv   , npiglo, npjglo, npk, rlon_v, rlat_v, depthv   )

    ncout_ke= create     (cf_out_ke, cf_tt   , npiglo, npjglo, npk  , cdep=cn_vdeptht        , ld_nc4=lnc4 )
    ierr    = createvar   (ncout_ke, stypvar3, jpnvarout2    , ipk2 , id_varout_ke           , ld_nc4=lnc4 )
    ierr    = putheadervar(ncout_ke, cf_tt   , npiglo, npjglo, npk, rlon_t, rlat_t, deptht   )

    dltim = getvar1d(cf_uu  , cn_vtimec,   npt     )
    ierr = putvar1d(ncout_u , dltim,       npt, 'T')
    ierr = putvar1d(ncout_v , dltim,       npt, 'T')
    ierr = putvar1d(ncout_ke, dltim,       npt, 'T')

  END SUBROUTINE CreateOutput


END PROGRAM
