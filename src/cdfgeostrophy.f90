PROGRAM cdfgeostrophy
  !!======================================================================
  !!                     ***  PROGRAM  cdfgeostrophy  ***
  !!=====================================================================
  !!  ** Purpose : Compute the ug and vg component of the geostrophic 
  !!               velocity from ssh and density field
  !!
  !!  ** Method  : * Integrate pressure from surface to current level
  !!                 P(n) = rho_insitu(1) * g * ssh 
  !!                      + sum( rho_insitu(k) * g * h(k) ) k=1,n-1
  !!                      + rho_insitu(n) * g * h(n) / 2
  !!
  !!                 h = level thickness
  !!
  !!               * Interpolation of Pressure on F points
  !!                 values on F-point are given
  !!                 by the demi-sum of X points (on the diagonal)
  !!
  !!                 p1 = 0.5 * ( A + B )
  !!                 p1 = 0.5 * ( B + C )
  !!                 p1 = 0.5 * ( C + D )
  !!
  !!                 F--------F--------F--------F
  !!                 |        |        |        |
  !!                 |        |   B    |        |
  !!                 |        |        |        |
  !!                 F--------p1--V----p2-------F
  !!                 |        |        |        |
  !!                 |   A    |  Pi,j  U   C    |
  !!                 |        |        |        |
  !!                 F--------F--------p3-------F
  !!                 |        |        |        |
  !!                 |        |   D    |        |
  !!                 |        |        |        |
  !!                 F--------F--------F--------F
  !!
  !!               * Compute local coriolis parameter at U and V point
  !!
  !!                 F--------F1--V----F2-------F
  !!                 |        |        |        |
  !!                 |        |   Pij  U        |
  !!                 |        |        |        |
  !!                 F--------F--------F3-------F
  !!
  !!                 Vg computation : Fij_v = 0.5 * ( F1 + F2 )
  !!                 Ug computation : Fij_u = 0.5 * ( F2 + F3 )
  !!
  !!               * Compute geostrophic balance
  !!
  !!                 Vg(i,j) = +1 * ( 1 / rho0 * Fij_v ) * ( p2 - p1 ) / e1v(i,j)
  !!                 Ug(i,j) = -1 * ( 1 / rho0 * Fij_u ) * ( p2 - p3 ) / e2u(i,j)
  !!               
  !!               * Masking :
  !!
  !!                 - if A,B or C are land points -> Vg = 0
  !!                 - if B,C or D are land points -> Ug = 0
  !!                 - multiplied by umask and vmask
  !!                 - if f < 1e-5, we mask
  !!
  !!  **  Note : Ug is located on a U grid point
  !!             Vg                 V grid point
  !!
  !!
  !! History : 3.0  : 01/2011  : R.Dussin : original code
  !!----------------------------------------------------------------------
  USE cdfio
  USE eos
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2011, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                           :: ji, jj, jt     ! dummy loop index
  INTEGER(KIND=4)                           :: npiglo, npjglo ! size of the domain
  INTEGER(KIND=4)                           :: npk, npt       ! size of the domain
  INTEGER(KIND=4)                           :: jk             ! vertical index
  INTEGER(KIND=4)                           :: narg, iargc    ! browse line
  INTEGER(KIND=4)                           :: ncoutu         ! ncid for ugeo file
  INTEGER(KIND=4)                           :: ncoutv         ! ncid for vgeo file
  INTEGER(KIND=4)                           :: ierr           ! error status
  INTEGER(KIND=4), DIMENSION(1)             :: ipk            ! levels of output vars
  INTEGER(KIND=4), DIMENSION(1)             :: id_varoutu     ! varid for ugeo
  INTEGER(KIND=4), DIMENSION(1)             :: id_varoutv     ! varid for vgeo

  REAL(KIND=4)                              :: grav           ! gravity
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: tim            ! time counter
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: deptht, depthw
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e1u, e2v, ff   ! horiz metrics, coriolis (f-point)
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e1v, e2u       ! horiz metrics
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e3             ! vertic metrics
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: glamu, gphiu   ! longitude latitude u-point
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: glamv, gphiv   ! longitude latitude v-point
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: umask, vmask   ! mask at u and v points
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: tmask          ! mask at t points
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zsigsurf       ! density at first level (used for zpsurf)
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zsiglevel      ! density at current level (used for zplevel/zphalflevel)
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zt, zsal       ! temporary arrays for temperature and salinity

  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: un, vn         ! velocity components
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: zsshn          ! ssh
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: zpupper        ! total pressure above current level
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: zphalflevel    ! pressure at T-point of current level
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: zplevel        ! pressure at bottom W-point of current level
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: zpsurf         ! pressure due to SSH
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: zptot          ! total pressure at current level

  REAL(KIND=4)                              :: zhlevel        ! thickness of current level
  REAL(KIND=4)                              :: zhhalflevel    ! thickness of half the current level
  REAL(KIND=4)                              :: zrho0          ! reference density in geos balance
  REAL(KIND=4)                              :: zohr0          ! reference density in geos balance
  REAL(KIND=4)                              :: zffu, zffv     ! local coriolis parameter
  REAL(KIND=4)                              :: zp1, zp2       ! dummy for pressure interp
  REAL(KIND=4)                              :: zp3, zp4       ! dummy for pressure interp
  REAL(KIND=4)                              :: zumask, zvmask ! dummy for mask 

  CHARACTER(LEN=256)                        :: cf_tfil        ! input file name
  CHARACTER(LEN=256)                        :: cf_uout='ugeo.nc' 
  CHARACTER(LEN=256)                        :: cf_vout='vgeo.nc'

  TYPE(variable), DIMENSION(1)              :: stypvaru       ! attributes for ugeo
  TYPE(variable), DIMENSION(1)              :: stypvarv       ! attributes for vgeo

  LOGICAL                                   :: lchk           ! file existence flag
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  grav = 9.81  ! gravity
  zrho0 = 1025 ! reference density
  zohr0 = 1. / zrho0

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfgeostrophy T-file'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the geostrophic velocity component from the pressure gradient '
     PRINT *,'       computed from SSH and in-situ density (T,S of input file) '
     PRINT *,'      '
     PRINT *,'     WARNING : USE AT YOUR OWN RISKS'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       T-file : netcdf file with SSH, T and S.' 
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        ',TRIM(cn_fmsk),' ',TRIM(cn_fhgr),' and ',TRIM(cn_fzgr)
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       - netcdf file : ', TRIM(cf_uout) 
     PRINT *,'           variables : ', TRIM(cn_vozocrtx)
     PRINT *,'       - netcdf file : ', TRIM(cf_vout) 
     PRINT *,'           variables : ', TRIM(cn_vomecrty)
     STOP
  ENDIF

  CALL getarg(1, cf_tfil)

  lchk = chkfile(cn_fhgr)
  lchk = chkfile(cn_fzgr) .OR. lchk
  lchk = chkfile(cn_fmsk) .OR. lchk
  lchk = chkfile(cf_tfil) .OR. lchk
  IF ( lchk ) STOP 99 ! missing file

  npiglo = getdim(cf_tfil, cn_x)
  npjglo = getdim(cf_tfil, cn_y)
  npk    = getdim(cf_tfil, cn_z) 
  npt    = getdim(cf_tfil, cn_t) 

  PRINT *, 'npiglo =', npiglo
  PRINT *, 'npjglo =', npjglo
  PRINT *, 'npk =', npk

  ipk(1)                        = npk
  stypvaru(1)%cname             = TRIM(cn_vozocrtx)
  stypvaru(1)%cunits            = 'm/s'
  stypvaru(1)%rmissing_value    = 0.
  stypvaru(1)%valid_min         = -20.
  stypvaru(1)%valid_max         = 20.
  stypvaru(1)%clong_name        = 'Zonal_Geostrophic_Velocity'
  stypvaru(1)%cshort_name       = TRIM(cn_vozocrtx)
  stypvaru(1)%conline_operation = 'N/A'
  stypvaru(1)%caxis             = 'TZYX'

  stypvarv(1)%cname             = TRIM(cn_vomecrty)
  stypvarv(1)%cunits            = 'm/s'
  stypvarv(1)%rmissing_value    = 0.
  stypvarv(1)%valid_min         = -20.
  stypvarv(1)%valid_max         = 20.
  stypvarv(1)%clong_name        = 'Meridional_Geostrophic_Velocity'
  stypvarv(1)%cshort_name       = TRIM(cn_vomecrty)
  stypvarv(1)%conline_operation = 'N/A'
  stypvarv(1)%caxis             = 'TZYX'

  ! Allocate the memory
  ALLOCATE ( e1v(npiglo,npjglo), e2u(npiglo,npjglo) )
  ALLOCATE ( ff(npiglo,npjglo) )
  ALLOCATE ( glamu(npiglo,npjglo), gphiu(npiglo,npjglo)  )
  ALLOCATE ( glamv(npiglo,npjglo), gphiv(npiglo,npjglo)  )
  ALLOCATE ( un(npiglo,npjglo), vn(npiglo,npjglo)  )
  ALLOCATE ( zsshn(npiglo,npjglo) )
  ALLOCATE ( umask(npiglo,npjglo), vmask(npiglo,npjglo) )
  ALLOCATE ( tmask(npiglo,npjglo) )
  ALLOCATE ( zpupper(npiglo,npjglo), zpsurf(npiglo,npjglo) )
  ALLOCATE ( zphalflevel(npiglo,npjglo), zplevel(npiglo,npjglo) )
  ALLOCATE ( zptot(npiglo,npjglo) )
  ALLOCATE ( zt(npiglo,npjglo), zsal(npiglo,npjglo) )
  ALLOCATE ( deptht(npk), depthw(npk) )
  ALLOCATE ( zsigsurf(npiglo,npjglo) , zsiglevel(npiglo,npjglo) )
  ALLOCATE ( e3(npiglo,npjglo) )
  ALLOCATE ( tim(npt) )

  ! Read the metrics from the mesh_hgr file
  e2u   = getvar(cn_fhgr, cn_ve2u,  1, npiglo, npjglo)
  e1v   = getvar(cn_fhgr, cn_ve1v,  1, npiglo, npjglo)
  ff    = getvar(cn_fhgr, cn_vff,   1, npiglo, npjglo) 

  glamu = getvar(cn_fhgr, cn_glamu, 1, npiglo, npjglo)
  gphiu = getvar(cn_fhgr, cn_gphiu, 1, npiglo, npjglo)
  glamv = getvar(cn_fhgr, cn_glamv, 1, npiglo, npjglo)
  gphiv = getvar(cn_fhgr, cn_gphiv, 1, npiglo, npjglo)

  deptht(:) = getvar1d(cf_tfil, cn_vdeptht, npk    )

  ! create output filesets
  ! U geo  
  ncoutu = create      (cf_uout, cf_tfil,  npiglo, npjglo, npk                            )
  ierr   = createvar   (ncoutu,  stypvaru, 1,      ipk,    id_varoutu                     )
  ierr   = putheadervar(ncoutu,  cf_tfil,  npiglo, npjglo, npk, pnavlon=glamu, pnavlat=gphiu)
  tim  = getvar1d(cf_tfil, cn_vtimec, npt     )
  ierr = putvar1d(ncoutu,  tim,       npt, 'T')

  ! V geo 
  ncoutv = create      (cf_vout, cf_tfil,  npiglo, npjglo, npk                            )
  ierr   = createvar   (ncoutv,  stypvarv, 1,      ipk,    id_varoutv                     )
  ierr   = putheadervar(ncoutv,  cf_tfil,  npiglo, npjglo, npk, pnavlon=glamv, pnavlat=gphiv)

  tim  = getvar1d(cf_tfil, cn_vtimec, npt     )
  ierr = putvar1d(ncoutv,  tim,       npt, 'T')

  ! time loop
  DO jt=1,npt

     ! Read ssh
     zsshn = getvar(cf_tfil, cn_sossheig, 1, npiglo, npjglo, ktime=jt)
     ! Read temperature and salinity
     zt   = getvar(cf_tfil, cn_votemper, 1, npiglo, npjglo, ktime=jt)
     zsal = getvar(cf_tfil, cn_vosaline, 1, npiglo, npjglo, ktime=jt)
     ! Compute density at first level
     zsigsurf(:,:) = 1000. + sigmai ( zt,zsal,deptht(1),npiglo,npjglo )
     ! Compute psurf (pressure due to SSH)
     zpsurf(:,:) = zsigsurf * grav * zsshn

     zpupper(:,:) = 0.d0

     DO jk=1,npk

         tmask = getvar(cn_fmsk, 'tmask',  jk, npiglo, npjglo)
         umask = getvar(cn_fmsk, 'umask',  jk, npiglo, npjglo)
         vmask = getvar(cn_fmsk, 'vmask',  jk, npiglo, npjglo)

         PRINT *,'Working on level ', jk
         !! 1. First we compute integrated pressure from the surface to current level

         ! Thickness
         e3  = getvar(cn_fzgr, cn_ve3t1d, jk, npiglo, npjglo)
         ! MAXVAL is used to avoid partial steps
         zhlevel     = MAXVAL(e3)
         zhhalflevel = 0.5 * MAXVAL(e3)
         ! 
         !PRINT *,' At level ', jk, ' thickness is ', zhlevel
         ! Read temperature and salinity at current level
         zt   = getvar(cf_tfil, cn_votemper, jk, npiglo, npjglo, ktime=jt)
         zsal = getvar(cf_tfil, cn_vosaline, jk, npiglo, npjglo, ktime=jt)
         ! Compute density of this level
         zsiglevel(:,:) = 1000. + sigmai ( zt,zsal,deptht(jk),npiglo,npjglo )
         ! Compute the pressure at T-point 
         zphalflevel(:,:) = zsiglevel * grav * zhhalflevel
         ! Compute the pressure at bottom W-point
         zplevel(:,:) =  zsiglevel * grav * zhlevel
         ! Compute the total pression -> This one is used in the geostrophic balance !
         zptot(:,:) = zpsurf(:,:) + zpupper(:,:) + zphalflevel(:,:)
         ! update zpupper for next level
         zpupper(:,:) = zpupper(:,:) + zplevel(:,:)

         !! 2. We compute the velocities from geostrophic balance

         un(:,:) = 0.d0
         vn(:,:) = 0.d0

         DO jj=2,npjglo-1
            DO ji=2,npiglo-1

               ! local coriolis parameter
               zffu = 0.5 * ( ff(ji,jj) + ff(ji,jj-1) )
               zffv = 0.5 * ( ff(ji,jj) + ff(ji-1,jj) )

               ! interp on F points 
               zp1 = 0.5 * ( zptot(ji-1,jj) + zptot(ji,jj+1) ) 
               zp2 = 0.5 * ( zptot(ji+1,jj) + zptot(ji,jj+1) )
               zp3 = 0.5 * ( zptot(ji,jj-1) + zptot(ji+1,jj) )

               zumask = tmask(ji,jj-1) * tmask(ji+1,jj) * tmask(ji,jj+1)
               zvmask = tmask(ji-1,jj) * tmask(ji,jj+1) * tmask(ji+1,jj)

               ! geostrophic balance
               vn(ji,jj) = +1 * ( zohr0 / zffv ) * ( zp2 - zp1 ) / e1v(ji,jj)
               un(ji,jj) = -1 * ( zohr0 / zffu ) * ( zp2 - zp3 ) / e2u(ji,jj)

               vn(ji,jj) = vn(ji,jj) * zvmask
               un(ji,jj) = un(ji,jj) * zumask

            ENDDO
         ENDDO

     WHERE ( ABS(ff) < 1.e-5 ) un(:,:) = 0.d0
     WHERE ( ABS(ff) < 1.e-5 ) vn(:,:) = 0.d0
 
     un(:,:) = un(:,:) * umask(:,:)
     vn(:,:) = vn(:,:) * vmask(:,:)

     ! write un and vn  ...
     ierr = putvar(ncoutu, id_varoutu(1), un(:,:), jk, npiglo, npjglo, ktime=jt)
     ierr = putvar(ncoutv, id_varoutv(1), vn(:,:), jk, npiglo, npjglo, ktime=jt)


    ENDDO ! vertical loop

  END DO  ! time loop

  ierr = closeout(ncoutu)
  ierr = closeout(ncoutv)

END PROGRAM cdfgeostrophy
