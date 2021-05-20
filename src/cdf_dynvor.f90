PROGRAM cdf_dynvor 
  !!======================================================================
  !!                     ***  PROGRAM  cdf_dynvor  ***
  !!=====================================================================
  !!  ** Purpose : Compute the relative and planetary vorticity terms 
  !!               of the momentum budget
  !!
  !!  ** Method  : Adapt NEMO dynvor.F90 to CDFTOOLS (cf below for further details)
  !!               with the energy and enstrophy conserving scheme (ln_dynvor_een=.true.).
  !!
  !!  ** Comments: Only the part associated with the flux form momentum advection 
  !!               scheme has been coded (12/09/2019)
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
  INTEGER(KIND=4)                              :: ji, jj, jk, jt           ! dummy loop indices
  INTEGER(KIND=4)                              :: it                       ! time index for vvl
  INTEGER(KIND=4)                              :: ierr                     ! working integer
  INTEGER(KIND=4)                              :: narg, iargc, ijarg       ! 
  INTEGER(KIND=4)                              :: jpiglo, jpjglo, jpk, jpt ! size of the domain
  INTEGER(KIND=4)                              :: jpim1, jpjm1, jpkm1      ! indexis -1
  INTEGER(KIND=4)                              :: ncout_fcor = 1           ! coriolis trend output number
  INTEGER(KIND=4)                              :: ncout_metr = 2           ! metric trend output number
  INTEGER(KIND=4)                              :: ncout_u, ncout_v, ncout_ke         ! ncid of output file
  INTEGER(KIND=4), DIMENSION(pnvarout)         :: ipk                      ! level of output variables
  INTEGER(KIND=4), DIMENSION(pnvarout)         :: id_varout_u              ! id of output variables (u-comp)
  INTEGER(KIND=4), DIMENSION(pnvarout)         :: id_varout_v              ! id of output variables (v-comp)
  INTEGER(KIND=4), DIMENSION(pnvarout)         :: id_varout_ke             ! id of output variables (ke-comp)

  REAL(wp), DIMENSION(:)    , ALLOCATABLE      :: dtim                     ! time
  REAL(wp), DIMENSION(:)    , ALLOCATABLE      :: deptht, depthu, depthv   ! z-grid (t,u,v)
  REAL(wp), DIMENSION(:,:)  , ALLOCATABLE      :: ff                       ! coriolis factor (2.*omega*sin(yphi) ) (s-1)
  REAL(wp), DIMENSION(:,:)  , ALLOCATABLE      :: nav_lon_t, nav_lat_t     ! t-grid hor.
  REAL(wp), DIMENSION(:,:)  , ALLOCATABLE      :: nav_lon_u, nav_lat_u     ! u-grid hor.
  REAL(wp), DIMENSION(:,:)  , ALLOCATABLE      :: nav_lon_v, nav_lat_v     ! v-grid hor.
  REAL(wp), DIMENSION(:,:)  , ALLOCATABLE      :: ht_0                     ! Reference ocean depth at T-points
  REAL(wp), DIMENSION(:,:)  , ALLOCATABLE      :: sshn                     ! now sea surface height
  REAL(wp), DIMENSION(:,:)  , ALLOCATABLE      :: e1t, e2t                 ! horizontal metric, t-pts
  REAL(wp), DIMENSION(:,:)  , ALLOCATABLE      :: e1u, e2u                 ! horizontal metric, u-pts
  REAL(wp), DIMENSION(:,:)  , ALLOCATABLE      :: e1v, e2v                 ! horizontal metric, v-pts
  REAL(wp), DIMENSION(:,:)  , ALLOCATABLE      :: e1f, e2f                 ! horizontal metrics
  REAL(wp), DIMENSION(:,:)  , ALLOCATABLE      :: e12t, r1_e12u, r1_e12v   ! face area at t-pts and inverse at u- v- pts
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: e3t_0, e3w_0             ! vet. metrics, t- w- pts, at rest (without VVL)
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: e3u_0, e3v_0             ! vet. metrics, u- v- pts, at rest (without VVL)
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: e3u, e3v, e3t            ! vet. metrics, u-, v-, t- pts
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: umask, vmask, tmask      ! Mask at U- V- points
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: un, vn                   ! velocity field
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: zua, zva                 ! u- v-mom trends due to coriolis/metric (output)
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE      :: zke                      ! KE       trends due to coriolis/metric (output)

  CHARACTER(LEN=256)                           :: cf_tfil                  ! temperature   netcdf file name (for mesh only)
  CHARACTER(LEN=256)                           :: cf_ufil                  ! zonal vel   netcdf file name
  CHARACTER(LEN=256)                           :: cf_vfil                  ! merid vel   netcdf file name
  CHARACTER(LEN=256)                           :: cf_sshfil                ! ssh         netcdf file name (for vvl)
  CHARACTER(LEN=255)                           :: cf_mh                    ! horiz. mesh netcdf file nam
  CHARACTER(LEN=255)                           :: cf_mz                    ! vert. mesh  netcdf file nam
  CHARACTER(LEN=255)                           :: cf_mask                  ! mask        netcdf file name
  CHARACTER(LEN=255)                           :: cf_bathy                 ! bathymetry  netcdf file name
  CHARACTER(LEN=256)                           :: cf_out_u='vor_u.nc'      ! output file name (u-comp)
  CHARACTER(LEN=256)                           :: cf_out_v='vor_v.nc'      ! output file name (v-comp)
  CHARACTER(LEN=256)                           :: cf_out_ke='vor_ke.nc'      ! output file name (ke-comp)
  CHARACTER(LEN=256)                           :: cldum                    ! dummy character variable
  CHARACTER(LEN=256)                           :: cglobal                  ! global attribute

  TYPE (variable), DIMENSION(pnvarout)         :: stypvar                  ! structure for attibutes (u-comp)
  TYPE (variable), DIMENSION(pnvarout)         :: stypvar2                 ! structure for attibutes (v-comp)
  TYPE (variable), DIMENSION(pnvarout)         :: stypvar3                 ! structure for attibutes (ke-comp)

  LOGICAL                                      :: lchk                     ! check missing files
  LOGICAL                                      :: lfull =.FALSE.           ! full step flag
  LOGICAL                                      :: lnc4  =.FALSE.           ! full step flag
  
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdf_dynvor -t T-file -u U-file -v V-file -ssh SSH-file ...'
     PRINT *,'          -mh MESH-file -mz MESZ-file -mask MASK-file -bathy BATHY-file ...'
     PRINT *,'          -o_u OUT-file-u -o_v OUT-file-v -o_ke OUT-file-ke'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'      Compute the relative and planetary vorticity terms'
     PRINT *,'      of the momentum budget.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -t T-file          : netcdf file for temperature (for mesh only)'
     PRINT *,'       -u U-file          : netcdf file for zonal vel'
     PRINT *,'       -v V-file          : netcdf file for meridional vel'
     PRINT *,'       -ssh SSH-file      : netcdf file for SSH (for vvl recomputation of vert. grid)'
     PRINT *,'       -mh MESH-file      : netcdf file for horizontal mesh'
     PRINT *,'       -mz MESH-file      : netcdf file for vertical mesh'
     PRINT *,'       -mask MASK-file    : netcdf file for mask'
     PRINT *,'       -bathy BATHY-file  : netcdf file for model bathymetry'
     PRINT *,'       -o_u OUT-file      : netcdf file for vortivity/metric term for u-momentum'
     PRINT *,'       -o_v OUT-file      : netcdf file for vortivity/metric term for v-momentum'
     PRINT *,'       -o_ke OUT-file     : netcdf file for vortivity/metric term for KE'
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
     CASE ('-t'        ) ; CALL getarg( ijarg, cf_tfil ) ; ijarg=ijarg+1
     CASE ('-u'        ) ; CALL getarg( ijarg, cf_ufil ) ; ijarg=ijarg+1
     CASE ('-v'        ) ; CALL getarg( ijarg, cf_vfil ) ; ijarg=ijarg+1
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

  !-- get dimensions (assuming all files have the same dimension that U-file) --
  jpiglo = getdim (cf_ufil, cn_x)
  jpjglo = getdim (cf_ufil, cn_y)
  jpk    = getdim (cf_ufil, cn_z)
  jpt    = getdim (cf_ufil, cn_t)
  jpim1 = jpiglo-1
  jpjm1 = jpjglo-1
  jpkm1 = jpk-1

  PRINT *, 'jpiglo =', jpiglo
  PRINT *, 'jpjglo =', jpjglo
  PRINT *, 'jpk    =', jpk
  PRINT *, 'jpt    =', jpt

  !-- Allocate arrays --
  ALLOCATE( ff(jpiglo, jpjglo)                                           )
  ALLOCATE( deptht(jpk)                   , depthu(jpk)                  , depthv(jpk)                )
  ALLOCATE( nav_lon_t(jpiglo, jpjglo)     , nav_lat_t(jpiglo, jpjglo)      )
  ALLOCATE( nav_lon_u(jpiglo, jpjglo)     , nav_lat_u(jpiglo, jpjglo)    )
  ALLOCATE( nav_lon_v(jpiglo, jpjglo)     , nav_lat_v(jpiglo, jpjglo)    )
  ALLOCATE( ht_0(jpiglo, jpjglo)                                         )
  ALLOCATE( sshn(jpiglo, jpjglo)                                         )
  ALLOCATE( e1t(jpiglo, jpjglo)           , e2t(jpiglo, jpjglo)          )
  ALLOCATE( e1u(jpiglo, jpjglo)           , e2u(jpiglo, jpjglo)          )
  ALLOCATE( e1v(jpiglo, jpjglo)           , e2v(jpiglo, jpjglo)          )
  ALLOCATE( e1f(jpiglo, jpjglo)           , e2f(jpiglo, jpjglo)          )
  ALLOCATE( e12t(jpiglo, jpjglo)                                         )
  ALLOCATE( r1_e12u(jpiglo, jpjglo)       , r1_e12v(jpiglo, jpjglo)      )
  ALLOCATE( e3u_0( jpiglo, jpjglo, jpk)   , e3v_0(jpiglo, jpjglo, jpk)   )
  ALLOCATE( e3u( jpiglo, jpjglo, jpk)     , e3v(jpiglo, jpjglo, jpk)     )
  ALLOCATE( e3t_0(jpiglo, jpjglo, jpk)    , e3t(jpiglo, jpjglo, jpk)     )
  ALLOCATE( tmask(jpiglo, jpjglo, jpk)    , umask(jpiglo, jpjglo, jpk)   , vmask(jpiglo, jpjglo, jpk) )
  ALLOCATE( un(jpiglo, jpjglo, jpk)       , vn(jpiglo, jpjglo, jpk)      ) 
  ALLOCATE( zua(jpiglo, jpjglo, jpk)      , zva(jpiglo, jpjglo, jpk)     ) 
  ALLOCATE( zke(jpiglo, jpjglo, jpk)                                     ) 

  !!-- loading -- 
  PRINT *, '-- LOAD VARIABLES --'
  !
  nav_lon_t    = getvar(cf_tfil, 'nav_lon', 1, jpiglo, jpjglo)
  nav_lat_t    = getvar(cf_tfil, 'nav_lat', 1, jpiglo, jpjglo)
  deptht       = getvar1d(cf_tfil, cn_vdeptht , jpk)
  nav_lon_u    = getvar(cf_ufil, 'nav_lon', 1, jpiglo, jpjglo)
  nav_lat_u    = getvar(cf_ufil, 'nav_lat', 1, jpiglo, jpjglo)
  depthu       = getvar1d(cf_ufil, cn_vdepthu , jpk)
  nav_lon_v    = getvar(cf_vfil, 'nav_lon', 1, jpiglo, jpjglo)
  nav_lat_v    = getvar(cf_vfil, 'nav_lat', 1, jpiglo, jpjglo)
  depthv       = getvar1d(cf_vfil, cn_vdepthv , jpk)
  !
  ht_0(:,:)    = getvar(cf_bathy, 'gdepw_0', 1, jpiglo, jpjglo )
  !ht_0(:,:)    = getvar(cf_mz    , 'hdepw', 1, jpiglo, jpjglo )
  e1t(:,:)     = getvar(cf_mh  , 'e1t'  , 1, jpiglo, jpjglo)
  e2t(:,:)     = getvar(cf_mh  , 'e2t'  , 1, jpiglo, jpjglo)
  e1u(:,:)     = getvar(cf_mh  , 'e1u'  , 1, jpiglo, jpjglo)
  e2u(:,:)     = getvar(cf_mh  , 'e2u'  , 1, jpiglo, jpjglo)
  e1v(:,:)     = getvar(cf_mh  , 'e1v'  , 1, jpiglo, jpjglo)
  e2v(:,:)     = getvar(cf_mh  , 'e2v'  , 1, jpiglo, jpjglo)
  e1f(:,:)     = getvar(cf_mh, cn_ve1f,  1, jpiglo, jpjglo )
  e2f(:,:)     = getvar(cf_mh, cn_ve2f,  1, jpiglo, jpjglo )
  e12t(:,:)    = e1t(:,:) * e2t(:,:)
  r1_e12u(:,:) = 1._wp / (e1u(:,:) * e2u(:,:))
  r1_e12v(:,:) = 1._wp / (e1v(:,:) * e2v(:,:))
  ff(:,:)      = getvar(cf_mh, 'ff'   ,  1, jpiglo, jpjglo ) 
  !
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

  DO jt = 1,jpt
     PRINT *, '======= time-step = ', jt

     !-- recomputed vert. metric due to non-linear free surface (VVL) --
     ! from domvvl.F90 -->> e3t = e3t_0*(1+sshn/ht_0)
     PRINT *, '-- Recompute vert. mesh --'
     sshn(:,:)     = getvar(cf_sshfil  , cn_sossheig, 1, jpiglo, jpjglo, ktime=jt )
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
   
     PRINT*, '-- Load velocities --'
     un(:,:,:)   = getvar3d(cf_ufil, cn_vozocrtx, jpiglo, jpjglo, jpk, ktime=jt )
     vn(:,:,:)   = getvar3d(cf_vfil, cn_vomecrty, jpiglo, jpjglo, jpk, ktime=jt )

     !-- compute momentum trend due to coriolis term --
     zua(:,:,:) = 0._wp
     zva(:,:,:) = 0._wp
     zke(:,:,:) = 0._wp
     CALL vor_een( jt, 1 ) !planetary vorticity (Coriolis)
     CALL trd_ken( zua, zva, zke)
     DO jk = 1,jpk
        ierr = putvar(ncout_u , id_varout_u(ncout_fcor) , zua(:,:,jk), jk, jpiglo, jpjglo, ktime=jt )
        ierr = putvar(ncout_v , id_varout_v(ncout_fcor) , zva(:,:,jk), jk, jpiglo, jpjglo, ktime=jt )
        ierr = putvar(ncout_ke, id_varout_ke(ncout_fcor), zke(:,:,jk), jk, jpiglo, jpjglo, ktime=jt )
     ENDDO

     !-- compute momentum trend due to metric term --
     zua(:,:,:) = 0._wp
     zva(:,:,:) = 0._wp
     zke(:,:,:) = 0._wp
     CALL vor_een( jt, 3 ) !metric term
     CALL trd_ken( zua, zva, zke)
     DO jk = 1,jpk
        ierr = putvar(ncout_u , id_varout_u(ncout_metr) , zua(:,:,jk), jk, jpiglo, jpjglo, ktime=jt )
        ierr = putvar(ncout_v , id_varout_v(ncout_metr) , zva(:,:,jk), jk, jpiglo, jpjglo, ktime=jt )
        ierr = putvar(ncout_ke, id_varout_ke(ncout_metr), zke(:,:,jk), jk, jpiglo, jpjglo, ktime=jt )
     ENDDO

  ENDDO		!jt-loop

  ierr = closeout(ncout_u)
  ierr = closeout(ncout_v)
  ierr = closeout(ncout_ke)


CONTAINS

   SUBROUTINE vor_een( kt, kvor )
!!----------------------------------------------------------------------
!!                ***  ROUTINE vor_een  ***
!!
!! ** Purpose :   Compute the now total vorticity trend and add it to
!!      the general trend of the momentum equation.
!!
!! ** Method  :   Trend evaluated using now fields (centered in time)
!!      and the Arakawa and Lamb (1980) flux form formulation : conserves
!!      both the horizontal kinetic energy and the potential enstrophy
!!      when horizontal divergence is zero (see the NEMO documentation)
!!      Add this trend to the general momentum trend (ua,va).
!!
!! ** Action : - Update (ua,va) with the now vorticity term trend
!!
!! References : Arakawa and Lamb 1980, Mon. Wea. Rev., 109, 18-36
!!----------------------------------------------------------------------
!
      INTEGER , INTENT(in   )                         ::   kt     ! ocean time-step index
      INTEGER , INTENT(in   )                         ::   kvor   ! =ncor (planetary) ; =ntot (total) ;
!                                                           ! =nrvm (relative vorticity or metric)
      INTEGER  ::   ji, jj, jk                                    ! dummy loop indices
      INTEGER  ::   ierr                                          ! local integer
      REAL(wp) ::   zfac12                                        ! local scalars
      REAL(wp) ::   zmsk, ze3                                     ! local scalars
!                                                           !  3D workspace
      REAL(wp), POINTER    , DIMENSION(:,:  )         :: zwx, zwy, zwz
      REAL(wp), POINTER    , DIMENSION(:,:  )         :: ztnw, ztne, ztsw, ztse

      REAL(wp), POINTER    , DIMENSION(:,:,:)         :: ze3f     !  3D workspace (lk_vvl=T)

!!----------------------------------------------------------------------
!
!      IF( kt == nit000 .OR. lk_vvl ) THEN      ! reciprocal of e3 at F-point (masked averaging of e3t over ocean points)
      ALLOCATE( zwx(jpiglo, jpjglo) , zwy(jpiglo, jpjglo), zwz(jpiglo, jpjglo)   )
      ALLOCATE( ztnw(jpiglo, jpjglo), ztne(jpiglo, jpjglo)                       )
      ALLOCATE( ztsw(jpiglo, jpjglo), ztse(jpiglo, jpjglo)                       )
      ALLOCATE( ze3f(jpiglo, jpjglo, jpk)                                        )

!         IF( ln_dynvor_een_old ) THEN ! original formulation
!              ....
!         ELSE ! new formulation from NEMO 3.6
            DO jk = 1, jpk
               DO jj = 1, jpjm1
                  DO ji = 1, jpim1
                     ze3  = ( e3t(ji,jj+1,jk)*tmask(ji,jj+1,jk) + e3t(ji+1,jj+1,jk)*tmask(ji+1,jj+1,jk)   &
                        &   + e3t(ji,jj  ,jk)*tmask(ji,jj  ,jk) + e3t(ji+1,jj  ,jk)*tmask(ji+1,jj  ,jk) )
                     zmsk = (                   tmask(ji,jj+1,jk) +                     tmask(ji+1,jj+1,jk)   &
                        &                     + tmask(ji,jj  ,jk) +                     tmask(ji+1,jj  ,jk) )
                     IF   ( ze3 /= 0._wp ) THEN ;   ze3f(ji,jj,jk) = zmsk / ze3
                     ELSE                       ;   ze3f(ji,jj,jk) = 0.0_wp
                     ENDIF
                  END DO
               END DO
            END DO
!         ENDIF
!      ENDIF
      zfac12 = 1._wp / 12._wp    ! Local constant initialization
!
!                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
!                                             ! ===============

! Potential vorticity and horizontal fluxes
! -----------------------------------------
         SELECT CASE( kvor )      ! vorticity considered
         CASE ( 1 )                                                ! planetary vorticity (Coriolis)
            zwz(:,:) = ff(:,:)      * ze3f(:,:,jk)
!         CASE ( 2 )                                                ! relative  vorticity
!            ....
         CASE ( 3 )                                                ! metric term
            DO jj = 1, jpjm1
               DO ji = 1, jpim1   ! vector opt.
                  zwz(ji,jj) = (   ( vn(ji+1,jj  ,jk) + vn (ji,jj,jk) ) * ( e2v(ji+1,jj  ) - e2v(ji,jj) )       &
                       &         - ( un(ji  ,jj+1,jk) + un (ji,jj,jk) ) * ( e1u(ji  ,jj+1) - e1u(ji,jj) )   )   &
                       &     * 0.5 / ( e1f(ji,jj) * e2f(ji,jj) ) * ze3f(ji,jj,jk)
               END DO
            END DO
!            CALL lbc_lnk( zwz, 'F', 1. )
!        CASE ( 4 )                                                ! total (relative + planetary vorticity)
!            zwz(:,:) = ( rotn(:,:,jk) + ff(:,:) ) * ze3f(:,:,jk)
!         CASE ( 5 )                                                ! total (coriolis + metric)
!             ....
         END SELECT

         zwx(:,:) = e2u(:,:) * e3u(:,:,jk) * un(:,:,jk)
         zwy(:,:) = e1v(:,:) * e3v(:,:,jk) * vn(:,:,jk)

! Compute and add the vorticity term trend
! ----------------------------------------
         jj = 2
         ztne(1,:) = 0   ;   ztnw(1,:) = 0   ;   ztse(1,:) = 0   ;   ztsw(1,:) = 0
         DO ji = 2, jpiglo
               ztne(ji,jj) = zwz(ji-1,jj  ) + zwz(ji  ,jj  ) + zwz(ji  ,jj-1)
               ztnw(ji,jj) = zwz(ji-1,jj-1) + zwz(ji-1,jj  ) + zwz(ji  ,jj  )
               ztse(ji,jj) = zwz(ji  ,jj  ) + zwz(ji  ,jj-1) + zwz(ji-1,jj-1)
               ztsw(ji,jj) = zwz(ji  ,jj-1) + zwz(ji-1,jj-1) + zwz(ji-1,jj  )
         END DO
         DO jj = 3, jpjglo
            DO ji = 2, jpiglo   ! vector opt. ok because we start at jj = 3
               ztne(ji,jj) = zwz(ji-1,jj  ) + zwz(ji  ,jj  ) + zwz(ji  ,jj-1)
               ztnw(ji,jj) = zwz(ji-1,jj-1) + zwz(ji-1,jj  ) + zwz(ji  ,jj  )
               ztse(ji,jj) = zwz(ji  ,jj  ) + zwz(ji  ,jj-1) + zwz(ji-1,jj-1)
               ztsw(ji,jj) = zwz(ji  ,jj-1) + zwz(ji-1,jj-1) + zwz(ji-1,jj  )
            END DO
         END DO
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               zua(ji,jj,jk) = + zfac12 / e1u(ji,jj) * (  ztne(ji,jj  ) * zwy(ji  ,jj  ) + ztnw(ji+1,jj) * zwy(ji+1,jj  )   &
                  &                           + ztse(ji,jj  ) * zwy(ji  ,jj-1) + ztsw(ji+1,jj) * zwy(ji+1,jj-1) )
               zva(ji,jj,jk) = - zfac12 / e2v(ji,jj) * (  ztsw(ji,jj+1) * zwx(ji-1,jj+1) + ztse(ji,jj+1) * zwx(ji  ,jj+1)   &
                  &                           + ztnw(ji,jj  ) * zwx(ji-1,jj  ) + ztne(ji,jj  ) * zwx(ji  ,jj  ) )
               zua(ji,jj,jk) = zua(ji,jj,jk) * umask(ji,jj,jk)
               zva(ji,jj,jk) = zva(ji,jj,jk) * vmask(ji,jj,jk)
            END DO
         END DO
!                                             ! ===============
      END DO                                           !   End of slab
!                                                ! ===============
      DEALLOCATE( zwx , zwy, zwz, ztnw, ztne, ztsw, ztse, ze3f )

   END SUBROUTINE vor_een

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
    ipk(:)                        = jpk
    stypvar(1)%ichunk             = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
    stypvar(1)%cname              = 'fvor_u'
    stypvar(1)%cunits             = 'm/s^2'
    stypvar(1)%rmissing_value     = 99999.
    stypvar(1)%valid_min          = -1.
    stypvar(1)%valid_max          = 1.
    stypvar(1)%clong_name         = 'U-momentum trend induced by the planetary vorticity (Coriolis)'
    stypvar(1)%cshort_name        = 'fvor_u'
    stypvar(1)%conline_operation  = 'On u-grid'
    stypvar(1)%caxis              = 'time depthu nav_lon_u nav_lat_u'
    !
    stypvar(2)%ichunk             = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
    stypvar(2)%cname              = 'metric_u'
    stypvar(2)%cunits             = 'm/s^2'
    stypvar(2)%rmissing_value     = 99999.
    stypvar(2)%valid_min          = -1.
    stypvar(2)%valid_max          = 1.
    stypvar(2)%clong_name         = 'U-momentum trend induced by metric term (flux formulation)'
    stypvar(2)%cshort_name        = 'metric_u'
    stypvar(2)%conline_operation  = 'On u-grid'
    stypvar(2)%caxis              = 'time depthu nav_lon_u nav_lat_u'

    stypvar2(1)%ichunk            = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
    stypvar2(1)%cname             = 'fvor_v'
    stypvar2(1)%cunits            = 'm/s^2'
    stypvar2(1)%rmissing_value    = 99999.
    stypvar2(1)%valid_min         = -1.
    stypvar2(1)%valid_max         = 1.
    stypvar2(1)%clong_name        = 'V-momentum trend induced by the planetary vorticity (Coriolis)'
    stypvar2(1)%cshort_name       = 'fvor_v'
    stypvar2(1)%conline_operation = 'On v-grid'
    stypvar2(1)%caxis             = 'time depthv nav_lon_v nav_lat_v'
    !
    stypvar2(2)%ichunk            = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
    stypvar2(2)%cname             = 'metric_v'
    stypvar2(2)%cunits            = 'm/s^2'
    stypvar2(2)%rmissing_value    = 99999.
    stypvar2(2)%valid_min         = -1.
    stypvar2(2)%valid_max         = 1.
    stypvar2(2)%clong_name        = 'V-momentum trend induced by metric term (flux formulation)'
    stypvar2(2)%cshort_name       = 'metric_v'
    stypvar2(2)%conline_operation = 'On v-grid'
    stypvar2(2)%caxis             = 'time depthv nav_lon_v nav_lat_v'

    stypvar3(1)%ichunk            = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
    stypvar3(1)%cname             = 'fvor_ke'
    stypvar3(1)%cunits            = 'm^2/s^3'
    stypvar3(1)%rmissing_value    = 99999.
    stypvar3(1)%valid_min         = -1.
    stypvar3(1)%valid_max         = 1.
    stypvar3(1)%clong_name        = 'KE trend induced by the planetary vorticity (Coriolis)'
    stypvar3(1)%cshort_name       = 'fvor_ke'
    stypvar3(1)%conline_operation = 'On t-grid'
    stypvar3(1)%caxis             = 'time deptht nav_lon_t nav_lat_t'
    !
    stypvar3(2)%ichunk            = (/jpiglo,MAX(1,jpjglo/30),1,1 /)
    stypvar3(2)%cname             = 'metric_ke'
    stypvar3(2)%cunits            = 'm^2/s^3'
    stypvar3(2)%rmissing_value    = 99999.
    stypvar3(2)%valid_min         = -1.
    stypvar3(2)%valid_max         = 1.
    stypvar3(2)%clong_name        = 'KE trend induced by metric term (flux formulation)'
    stypvar3(2)%cshort_name       = 'metric_ke'
    stypvar3(2)%conline_operation = 'On t-grid'
    stypvar3(2)%caxis             = 'time deptht nav_lon_t nav_lat_t'


    ! create output fileset
    ncout_u = create      (cf_out_u, cf_ufil ,  jpiglo, jpjglo, jpk, cdep=cn_vdepthu  , ld_nc4=lnc4 )
    ierr    = createvar   (ncout_u , stypvar ,  pnvarout, ipk , id_varout_u           , ld_nc4=lnc4 )
    ierr    = putheadervar(ncout_u , cf_ufil ,  jpiglo, jpjglo, jpk, nav_lon_u, nav_lat_u, depthu   )

    ncout_v = create      (cf_out_v, cf_vfil ,  jpiglo, jpjglo, jpk, cdep=cn_vdepthv  , ld_nc4=lnc4 )
    ierr    = createvar   (ncout_v , stypvar2,  pnvarout, ipk ,  id_varout_v          , ld_nc4=lnc4 )
    ierr    = putheadervar(ncout_v , cf_vfil ,  jpiglo, jpjglo, jpk, nav_lon_v, nav_lat_v, depthv   )

    ncout_ke= create      (cf_out_ke, cf_tfil ,  jpiglo, jpjglo, jpk, cdep=cn_vdeptht  , ld_nc4=lnc4 )
    ierr    = createvar   (ncout_ke , stypvar3,  pnvarout, ipk , id_varout_ke          , ld_nc4=lnc4 )
    ierr    = putheadervar(ncout_ke , cf_tfil ,  jpiglo, jpjglo, jpk, nav_lon_t, nav_lat_t, deptht   )

    dtim = getvar1d(cf_ufil , cn_vtimec,   jpt     )
    ierr = putvar1d(ncout_u , dtim,        jpt, 'T')
    ierr = putvar1d(ncout_v , dtim,        jpt, 'T')
    ierr = putvar1d(ncout_ke, dtim,       jpt, 'T')


  END SUBROUTINE CreateOutput

END PROGRAM
