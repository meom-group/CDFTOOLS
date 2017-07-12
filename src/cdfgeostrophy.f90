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
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  USE cdfio
  USE eos
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class derived_fields
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                           :: ji, jj, jt     ! dummy loop index
  INTEGER(KIND=4)                           :: it             ! time index for vvl
  INTEGER(KIND=4)                           :: npiglo, npjglo ! size of the domain
  INTEGER(KIND=4)                           :: npk, npt       ! size of the domain
  INTEGER(KIND=4)                           :: jk             ! vertical index
  INTEGER(KIND=4)                           :: narg, iargc    ! browse line
  INTEGER(KIND=4)                           :: ijarg          ! browse line
  INTEGER(KIND=4)                           :: ncoutu         ! ncid for ugeo file
  INTEGER(KIND=4)                           :: ncoutv         ! ncid for vgeo file
  INTEGER(KIND=4)                           :: ierr           ! error status
  INTEGER(KIND=4), DIMENSION(1)             :: ipk            ! levels of output vars
  INTEGER(KIND=4), DIMENSION(1)             :: id_varoutu     ! varid for ugeo
  INTEGER(KIND=4), DIMENSION(1)             :: id_varoutv     ! varid for vgeo

  REAL(KIND=4)                              :: grav           ! gravity
  REAL(KIND=4)                              :: zhlevel        ! thickness of current level
  REAL(KIND=4)                              :: zhhalflevel    ! thickness of half the current level
  REAL(KIND=4)                              :: zrho0          ! reference density in geos balance
  REAL(KIND=4)                              :: zohr0          ! reference density in geos balance
  REAL(KIND=4)                              :: zffu, zffv     ! local coriolis parameter
  REAL(KIND=4)                              :: zumask, zvmask ! dummy for mask 
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: deptht, depthw
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e1u, e2v, ff   ! horiz metrics, coriolis (f-point)
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e1v, e2u       ! horiz metrics
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e3             ! vertic metrics
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: glamu, gphiu   ! longitude latitude u-point
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: glamv, gphiv   ! longitude latitude v-point
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: umask, vmask   ! mask at u and v points
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: tmask          ! mask at t points
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zt, zsal       ! temporary arrays for temperature and salinity

  REAL(KIND=8)                              :: dp1, dp2, dp3  ! dummy for pressure interp
  REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: dtim           ! time counter
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: dsigsurf       ! density at first level (used for dpsurf)
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: dsiglevel      ! density at current level (used for dplevel/dphalflevel)
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dun, dvn       ! velocity components
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dsshn          ! ssh
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dpupper        ! total pressure above current level
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dphalflevel    ! pressure at T-point of current level
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dplevel        ! pressure at bottom W-point of current level
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dpsurf         ! pressure due to SSH
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dptot          ! total pressure at current level

  CHARACTER(LEN=256)                        :: cf_tfil        ! input file name
  CHARACTER(LEN=256)                        :: cf_uout='ugeo.nc' 
  CHARACTER(LEN=256)                        :: cf_vout='vgeo.nc'
  CHARACTER(LEN=256)                        :: cldum          ! working char variable

  TYPE(variable), DIMENSION(1)              :: stypvaru       ! attributes for ugeo
  TYPE(variable), DIMENSION(1)              :: stypvarv       ! attributes for vgeo

  LOGICAL                                   :: lchk           ! file existence flag
  LOGICAL                                   :: lnc4 = .FALSE. ! Use nc4 with chunking and deflation
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  grav  = 9.81  ! m/s2   gravity
  zrho0 = 1025  ! kg/m3  reference density
  zohr0 = 1. / zrho0

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfgeostrophy -f T-file [-o OUT-ufile OUT-vfile] [-nc4] [-vvl]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the geostrophic velocity components from the pressure gradient'
     PRINT *,'       induced by SSH and in-situ density (T,S of input file), using the '
     PRINT *,'       thermal wind equation.'
     PRINT *,'      '
     PRINT *,'     WARNING : USE AT YOUR OWN RISKS. '
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f T-file : netcdf file with SSH, T and S.' 
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-o OUT-ufile OUT-vfile] : Specify output files name''s  instead of '
     PRINT *,'            ',TRIM(cf_uout),' and ', TRIM(cf_vout)
     PRINT *,'       [-nc4 ]:  Use netcdf4 output with chunking and deflation level 1.'
     PRINT *,'                 This option is effective only if cdftools are compiled with'
     PRINT *,'                 a netcdf library supporting chunking and deflation.'
     PRINT *,'       [-vvl] : use time varying vertical metrics.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        ',TRIM(cn_fmsk),' ',TRIM(cn_fhgr),' and ',TRIM(cn_fzgr)
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       - netcdf file : ', TRIM(cf_uout) 
     PRINT *,'           variables : ', TRIM(cn_vozocrtx)
     PRINT *,'       - netcdf file : ', TRIM(cf_vout) 
     PRINT *,'           variables : ', TRIM(cn_vomecrty)
     PRINT *,'      '
     STOP 
  ENDIF

  ijarg=1
  DO WHILE (ijarg <= narg ) 
     CALL getarg(ijarg, cldum) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-f '  ) ; CALL getarg(ijarg, cf_tfil) ; ijarg=ijarg+1
        ! option
     CASE ( '-o '  ) ; CALL getarg(ijarg, cf_uout) ; ijarg=ijarg+1
        ;              CALL getarg(ijarg, cf_vout) ; ijarg=ijarg+1
     CASE ( '-nc4' ) ; lnc4   = .TRUE.
     CASE ( '-vvl' ) ; lg_vvl = .TRUE.
     CASE DEFAULT    ; PRINT *,' ERROR : ', TRIM(cldum), ' : unknown option.' ; STOP 99
     END SELECT
  ENDDO

  lchk = chkfile(cn_fhgr)
  lchk = chkfile(cn_fzgr) .OR. lchk
  lchk = chkfile(cn_fmsk) .OR. lchk
  lchk = chkfile(cf_tfil) .OR. lchk
  IF ( lchk ) STOP 99 ! missing file

  IF ( lg_vvl ) THEN
     cn_fe3t = cf_tfil
     cn_ve3t = cn_ve3tvvl
  ENDIF

  npiglo = getdim(cf_tfil, cn_x)
  npjglo = getdim(cf_tfil, cn_y)
  npk    = getdim(cf_tfil, cn_z) 
  npt    = getdim(cf_tfil, cn_t) 

  PRINT *, 'npiglo =', npiglo
  PRINT *, 'npjglo =', npjglo
  PRINT *, 'npk =', npk

  ! Allocate the memory
  ALLOCATE ( e1v(npiglo,npjglo), e2u(npiglo,npjglo) )
  ALLOCATE ( ff(npiglo,npjglo) )
  ALLOCATE ( glamu(npiglo,npjglo), gphiu(npiglo,npjglo)  )
  ALLOCATE ( glamv(npiglo,npjglo), gphiv(npiglo,npjglo)  )
  ALLOCATE ( dun(npiglo,npjglo),     dvn(npiglo,npjglo)  )
  ALLOCATE ( dsshn(npiglo,npjglo) )
  ALLOCATE ( umask(npiglo,npjglo),  vmask(npiglo,npjglo) )
  ALLOCATE ( tmask(npiglo,npjglo) )
  ALLOCATE ( dpupper(npiglo,npjglo),      dpsurf(npiglo,npjglo) )
  ALLOCATE ( dphalflevel(npiglo,npjglo), dplevel(npiglo,npjglo) )
  ALLOCATE ( dptot(npiglo,npjglo) )
  ALLOCATE ( zt(npiglo,npjglo), zsal(npiglo,npjglo) )
  ALLOCATE ( deptht(npk), depthw(npk) )
  ALLOCATE ( dsigsurf(npiglo,npjglo) , dsiglevel(npiglo,npjglo) )
  ALLOCATE ( e3(npiglo,npjglo) )
  ALLOCATE ( dtim(npt) )

  ! Read the metrics from the mesh_hgr file
  e2u   = getvar(cn_fhgr, cn_ve2u,  1, npiglo, npjglo)
  e1v   = getvar(cn_fhgr, cn_ve1v,  1, npiglo, npjglo)
  ff    = getvar(cn_fhgr, cn_vff,   1, npiglo, npjglo) 

  glamu = getvar(cn_fhgr, cn_glamu, 1, npiglo, npjglo)
  gphiu = getvar(cn_fhgr, cn_gphiu, 1, npiglo, npjglo)
  glamv = getvar(cn_fhgr, cn_glamv, 1, npiglo, npjglo)
  gphiv = getvar(cn_fhgr, cn_gphiv, 1, npiglo, npjglo)

  deptht(:) = getvar1d(cf_tfil, cn_vdeptht, npk    )

  CALL CreateOutputUV

  ! time loop
  DO jt=1,npt
     IF ( lg_vvl ) THEN  ; it=jt
     ELSE                ; it=1
     ENDIF

     ! Read ssh
     dsshn = getvar(cf_tfil, cn_sossheig, 1, npiglo, npjglo, ktime=jt)*1.d0
     ! Read temperature and salinity
     zt   = getvar(cf_tfil, cn_votemper, 1, npiglo, npjglo, ktime=jt)
     zsal = getvar(cf_tfil, cn_vosaline, 1, npiglo, npjglo, ktime=jt)
     ! Compute density at first level
     dsigsurf(:,:) = 1000.d0 + sigmai ( zt,zsal,deptht(1),npiglo,npjglo )
     ! Compute psurf (pressure due to SSH)
     dpsurf(:,:)  = dsigsurf * grav * dsshn
     dpupper(:,:) = 0.d0

     DO jk=1,npk
        tmask = getvar(cn_fmsk, cn_tmask,  jk, npiglo, npjglo)
        umask = getvar(cn_fmsk, cn_umask,  jk, npiglo, npjglo)
        vmask = getvar(cn_fmsk, cn_vmask,  jk, npiglo, npjglo)

        PRINT *,'Working on level ', jk
        !! 1. First we compute integrated pressure from the surface to current level

        ! Thickness
        e3(:,:)  = getvar(cn_fe3t, cn_ve3t, jk, npiglo, npjglo, ktime=it, ldiom=.NOT.lg_vvl )
        ! MAXVAL is used to avoid partial steps 
        ! JMM : this is probably not valid for vvl ...
        zhlevel     = MAXVAL(e3)
        zhhalflevel = 0.5 * MAXVAL(e3)
        ! 
        ! Read temperature and salinity at current level
        zt   = getvar(cf_tfil, cn_votemper, jk, npiglo, npjglo, ktime=jt)
        zsal = getvar(cf_tfil, cn_vosaline, jk, npiglo, npjglo, ktime=jt)
        ! Compute density of this level
        dsiglevel(:,:) = 1000.d0 + sigmai ( zt,zsal,deptht(jk),npiglo,npjglo )
        ! Compute the pressure at T-point 
        dphalflevel(:,:) = dsiglevel * grav * zhhalflevel
        ! Compute the pressure at bottom W-point
        dplevel(:,:) =  dsiglevel * grav * zhlevel
        ! Compute the total pression -> This one is used in the geostrophic balance !
        dptot(:,:) = dpsurf(:,:) + dpupper(:,:) + dphalflevel(:,:)
        ! update dpupper for next level
        dpupper(:,:) = dpupper(:,:) + dplevel(:,:)

        !! 2. We compute the velocities from geostrophic balance
        dun(:,:) = 0.d0
        dvn(:,:) = 0.d0

        DO jj=2,npjglo-1
           DO ji=2,npiglo-1

              ! local coriolis parameter
              zffu = 0.5 * ( ff(ji,jj) + ff(ji,jj-1) )
              zffv = 0.5 * ( ff(ji,jj) + ff(ji-1,jj) )

              ! interp on F points 
              dp1 = 0.5 * ( dptot(ji-1,jj) + dptot(ji,jj+1) ) 
              dp2 = 0.5 * ( dptot(ji+1,jj) + dptot(ji,jj+1) )
              dp3 = 0.5 * ( dptot(ji,jj-1) + dptot(ji+1,jj) )

              zumask = tmask(ji,jj-1) * tmask(ji+1,jj) * tmask(ji,jj+1)
              zvmask = tmask(ji-1,jj) * tmask(ji,jj+1) * tmask(ji+1,jj)

              ! geostrophic balance
              dvn(ji,jj) = +1 * ( zohr0 / zffv ) * ( dp2 - dp1 ) / e1v(ji,jj)
              dun(ji,jj) = -1 * ( zohr0 / zffu ) * ( dp2 - dp3 ) / e2u(ji,jj)

              dvn(ji,jj) = dvn(ji,jj) * zvmask
              dun(ji,jj) = dun(ji,jj) * zumask

           ENDDO
        ENDDO

        WHERE ( ABS(ff) < 1.e-5 ) dun(:,:) = 0.d0
        WHERE ( ABS(ff) < 1.e-5 ) dvn(:,:) = 0.d0

        dun(:,:) = dun(:,:) * umask(:,:)
        dvn(:,:) = dvn(:,:) * vmask(:,:)

        ! write dun and dvn  ...
        ierr = putvar(ncoutu, id_varoutu(1), dun(:,:), jk, npiglo, npjglo, ktime=jt)
        ierr = putvar(ncoutv, id_varoutv(1), dvn(:,:), jk, npiglo, npjglo, ktime=jt)
     ENDDO ! vertical loop
  END DO  ! time loop

  ierr = closeout(ncoutu)
  ierr = closeout(ncoutv)

CONTAINS

  SUBROUTINE CreateOutputUV
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutputUV  ***
    !!
    !! ** Purpose :  Create netcdf output file(s)  ugeo and vgeo
    !!
    !! ** Method  :  Use stypvar global description of variables
    !!
    !!----------------------------------------------------------------------
    ipk(1)                        = npk
    stypvaru(1)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvaru(1)%cname             = TRIM(cn_vozocrtx)
    stypvaru(1)%cunits            = 'm/s'
    stypvaru(1)%rmissing_value    = 0.
    stypvaru(1)%valid_min         = -20.
    stypvaru(1)%valid_max         = 20.
    stypvaru(1)%clong_name        = 'Zonal_Geostrophic_Velocity'
    stypvaru(1)%cshort_name       = TRIM(cn_vozocrtx)
    stypvaru(1)%conline_operation = 'N/A'
    stypvaru(1)%caxis             = 'TZYX'

    stypvarv(1)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvarv(1)%cname             = TRIM(cn_vomecrty)
    stypvarv(1)%cunits            = 'm/s'
    stypvarv(1)%rmissing_value    = 0.
    stypvarv(1)%valid_min         = -20.
    stypvarv(1)%valid_max         = 20.
    stypvarv(1)%clong_name        = 'Meridional_Geostrophic_Velocity'
    stypvarv(1)%cshort_name       = TRIM(cn_vomecrty)
    stypvarv(1)%conline_operation = 'N/A'
    stypvarv(1)%caxis             = 'TZYX'

    ! create output filesets
    ! U geo  
    ncoutu = create      (cf_uout, cf_tfil,  npiglo, npjglo, npk        , ld_nc4=lnc4         )
    ierr   = createvar   (ncoutu,  stypvaru, 1,      ipk,    id_varoutu , ld_nc4=lnc4         )
    ierr   = putheadervar(ncoutu,  cf_tfil,  npiglo, npjglo, npk, pnavlon=glamu, pnavlat=gphiu)

    dtim = getvar1d(cf_tfil, cn_vtimec, npt     )
    ierr = putvar1d(ncoutu,  dtim,      npt, 'T')

    ! V geo 
    ncoutv = create      (cf_vout, cf_tfil,  npiglo, npjglo, npk        , ld_nc4=lnc4         )
    ierr   = createvar   (ncoutv,  stypvarv, 1,      ipk,    id_varoutv , ld_nc4=lnc4         )
    ierr   = putheadervar(ncoutv,  cf_tfil,  npiglo, npjglo, npk, pnavlon=glamv, pnavlat=gphiv)

    dtim = getvar1d(cf_tfil, cn_vtimec, npt     )
    ierr = putvar1d(ncoutv, dtim,       npt, 'T')

  END SUBROUTINE CreateOutputUV

END PROGRAM cdfgeostrophy
