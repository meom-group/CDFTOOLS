PROGRAM cdfmeshmask
  !!======================================================================
  !!                     ***  PROGRAM  cdfmeshmask  ***
  !!=====================================================================
  !!  ** Purpose : build mesh mask file from bathymetry
  !!
  !!  ** Method  : use nemo3.6 simplified zgr_xxx routines for zps
  !!               In this tools, for similarities with NEMO, NF90
  !!               functions are used in this program.
  !!
  !! History : 3.0  : 10/2014  : J.M. Molines 
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!----------------------------------------------------------------------
  USE netcdf
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class mask
  !!----------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER, PARAMETER :: wp=8

  INTEGER :: narg, iargc, ijarg
  INTEGER :: npiglo, npjglo, nkmax
  INTEGER :: nbloc_sz=193  , nbloc=-1
  INTEGER, DIMENSION(:), ALLOCATABLE :: nbloc_szt, njbloct
  INTEGER, DIMENSION(:), ALLOCATABLE :: nblock_szt, nkbloct
  INTEGER :: jbloc, ij, i1, i2, ji     !
  INTEGER :: nperio = 0    ! closed boundary only for the time being
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: mbathy

  INTEGER :: ierr
  INTEGER :: nczgr, nchgr, ncmsk, nccoo  ! ncdif of output files
  INTEGER :: id_dept1d, id_depw1d, id_e3t1d, id_e3w1d, id_dept, id_depw, id_e3t, id_e3w, id_mbat  ! zgr
  INTEGER :: id_e3u, id_e3v, id_hdept, id_hdepw, id_navlat, id_navlon, id_navlev, id_time
  INTEGER :: idx, idy, idz, idt, idvar
  INTEGER :: id_tmsk, id_umsk, id_vmsk, id_fmsk
  INTEGER :: id_tmsku, id_umsku, id_vmsku, id_fmsku

  INTEGER :: id_lamt, id_lamu, id_lamv, id_lamf    ! hgr  copy of coordinates  !!!
  INTEGER :: id_phit, id_phiu, id_phiv, id_phif 

  INTEGER :: id_tmask, id_umask, id_vmask, id_fmask  ! mask file

  REAL(wp) :: rhmin
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: bathy, rlon, rlat

  REAL(wp), DIMENSION(:),   ALLOCATABLE :: gdepw_1d, gdept_1d, e3w_1d, e3t_1d
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: gdept_0  ! npiglo, jbloc, jpk
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: gdepw_0  ! npiglo, jbloc, jpk
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: e3t_0    ! npiglo, jbloc, jpk  then npiglo, npjglo, kbloc
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: e3w_0    ! npiglo, jbloc, jpk  then npiglo, npjglo, kbloc
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: e3u_0    ! npiglo, jbloc, jpk  then npiglo, npjglo, kbloc
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: e3v_0    ! npiglo, jbloc, jpk  then npiglo, npjglo, kbloc
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: tmask    ! npiglo, npjglo, kbloc
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: umask    ! npiglo, npjglo, kbloc
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: vmask    ! npiglo, npjglo, kbloc
  REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: fmask    ! npiglo, npjglo, kbloc

  LOGICAL :: lchk=.FALSE.

  CHARACTER(LEN=80) :: cldum
  CHARACTER(LEN=80) :: cf_nam
  CHARACTER(LEN=80) :: cf_bat
  CHARACTER(LEN=80) :: cf_coo
  CHARACTER(LEN=80) :: cf_zgr = 'cdf_mesh_zgr.nc'
  CHARACTER(LEN=80) :: cf_hgr = 'cdf_mesh_hgr.nc'
  CHARACTER(LEN=80) :: cf_msk = 'cdf_mask.nc'

  INTEGER ::    jpk     , jpkm1
  INTEGER  ::   nn_bathy        !: = 0/1 ,compute/read the bathymetry file
  REAL(wp) ::   rn_bathy        !: depth of flat bottom (active if nn_bathy=0; if =0 depth=jpkm1)
  REAL(wp) ::   rn_hmin         !: minimum ocean depth (>0) or minimum number of ocean levels (<0)
  REAL(wp) ::   rn_e3zps_min    !: miminum thickness for partial steps (meters)
  REAL(wp) ::   rn_e3zps_rat    !: minimum thickness ration for partial steps
  INTEGER  ::   nn_msh          !: = 1 create a mesh-mask file
  INTEGER  ::   nn_acc          !: = 0/1 use of the acceleration of convergence technique
  REAL(wp) ::   rn_atfp         !: asselin time filter parameter
  REAL(wp) ::   rn_rdt          !: time step for the dynamics (and tracer if nacc=0)
  REAL(wp) ::   rn_rdtmin       !: minimum time step on tracers
  REAL(wp) ::   rn_rdtmax       !: maximum time step on tracers
  REAL(wp) ::   rn_rdth         !: depth variation of tracer step
  INTEGER  ::   nn_closea       !: =0 suppress closed sea/lake from the ORCA domain or not (=1)
  INTEGER  ::   nn_euler        !: =0 start with forward time step or not (=1)
  LOGICAL  ::   ln_crs          !: Apply grid coarsening to dynamical model output or online passive tracers

  INTEGER       ::   jphgr_msh        !: type of horizontal mesh
  !                                       !  = 0 curvilinear coordinate on the sphere read in coordinate.nc
  !                                       !  = 1 geographical mesh on the sphere with regular grid-spacing
  !                                       !  = 2 f-plane with regular grid-spacing
  !                                       !  = 3 beta-plane with regular grid-spacing
  !                                       !  = 4 Mercator grid with T/U point at the equator

  REAL(wp)      ::   ppglam0              !: longitude of first raw and column T-point (jphgr_msh = 1)
  REAL(wp)      ::   ppgphi0              !: latitude  of first raw and column T-point (jphgr_msh = 1)
  !                                                        !  used for Coriolis & Beta parameters (jphgr_msh = 2 or 3)
  REAL(wp)      ::   ppe1_deg             !: zonal      grid-spacing (degrees)
  REAL(wp)      ::   ppe2_deg             !: meridional grid-spacing (degrees)
  REAL(wp)      ::   ppe1_m               !: zonal      grid-spacing (degrees)
  REAL(wp)      ::   ppe2_m               !: meridional grid-spacing (degrees)

  REAL(wp)      ::   ppsur                !: ORCA r4, r2 and r05 coefficients
  REAL(wp)      ::   ppa0                 !: (default coefficients)
  REAL(wp)      ::   ppa1                 !:
  REAL(wp)      ::   ppkth                !:
  REAL(wp)      ::   ppacr                !:
  !
  !  If both ppa0 ppa1 and ppsur are specified to 0, then
  !  they are computed from ppdzmin, pphmax , ppkth, ppacr in dom_zgr
  REAL(wp)      ::   ppdzmin              !: Minimum vertical spacing
  REAL(wp)      ::   pphmax               !: Maximum depth
  !
  LOGICAL       ::   ldbletanh            !: Use/do not use double tanf function for vertical coordinates
  REAL(wp)      ::   ppa2                 !: Double tanh function parameters
  REAL(wp)      ::   ppkth2               !:
  REAL(wp)      ::   ppacr2               !:
  !-----------------------------------------------------------------------------
  narg=iargc()

  IF ( narg == 0 ) THEN
     PRINT *,' usage :  cdfmeshmask -n NAMELIST-file -b BATHY-file  -c COORD-file ...'
     PRINT *,'         ... [-njbloc nbloc] '
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Create mesh_mask from bathymetry and namdom information (namelist)' 
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -n NAMELIST-file : name of the namelist file (with NEMO namdom block)' 
     PRINT *,'       -b BATHY-file : name of bathymetry (meters)'
     PRINT *,'       -c COORD-file : name of coordinates file'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-njbloc nbloc]: number of j-bloc of rows to treat together. Increasing'
     PRINT *,'             nbloc decreases memory usage but increases writing time.'
     PRINT *,'            default : nbloc = npjglo (worst condition) '
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        namelist, bathymetry and coordinated passed on the command line' 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf files : ', TRIM(cf_zgr),' ', TRIM(cf_hgr),' and ', TRIM(cf_msk)
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'       '
     PRINT *,'      '
     STOP 
  ENDIF

  ijarg=1
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ('-n '    ) ; CALL getarg(ijarg,cf_nam) ; ijarg=ijarg+1
     CASE ('-b '    ) ; CALL getarg(ijarg,cf_bat) ; ijarg=ijarg+1
     CASE ('-c '    ) ; CALL getarg(ijarg,cf_coo) ; ijarg=ijarg+1
     CASE ('-njbloc') ; CALL getarg(ijarg,cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) nbloc
     CASE DEFAULT     ; PRINT *,' ERROR : ',TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  ENDDO
  lchk = lchk .OR. chkfile(cf_nam) 
  lchk = lchk .OR. chkfile(cf_bat) 
  lchk = lchk .OR. chkfile(cf_coo) 
  IF ( lchk ) STOP 99 ! missing files

  npiglo = getdim(cf_coo, cn_x)
  npjglo = getdim(cf_coo, cn_y)
  IF ( nbloc == -1 ) nbloc = npjglo
  ! build nbloc_szt and njbloct
  ALLOCATE(  nbloc_szt(nbloc), njbloct(nbloc) )
  i1 = npjglo/nbloc
  i2 = MOD(npjglo, nbloc )
  nbloc_szt(1         :nbloc-i2) = i1
  nbloc_szt(nbloc-i2+1:nbloc   ) = i1+1
  njbloct(1) = 1
  DO ji=2, nbloc
     njbloct(ji) = njbloct(ji-1) + nbloc_szt(ji-1)
  ENDDO

  ALLOCATE (bathy(npiglo, npjglo ), rlon(npiglo, npjglo )  , rlat(npiglo, npjglo )  ) 
  ALLOCATE (mbathy(npiglo, npjglo )) 

  CALL zgr_z
  bathy(:,:) = getvar(cf_bat, cn_bathymet, 1, npiglo, npjglo)
  ! set minimum value here :
  IF( rn_hmin < 0._wp ) THEN ; nkmax = - INT( rn_hmin )                                       ! from a nb of level
  ELSE                       ; nkmax = MINLOC( gdepw_1d, mask = gdepw_1d > rn_hmin, dim = 1 ) ! from a depth
  ENDIF

  rhmin = gdepw_1d(nkmax+1)                                        ! minimum depth = ik+1 w-levels
  WHERE( bathy(:,:) <= 0._wp ) ; bathy(:,:) = 0._wp                         ! min=0     over the lands
ELSE WHERE                   ; bathy(:,:) = MAX(  rhmin , bathy(:,:)  )   ! min=rhmin over the oceans
END WHERE

PRINT *,  'Minimum ocean depth: ', rhmin, ' minimum number of ocean levels : ', nkmax

CALL zgr_zps

CONTAINS
SUBROUTINE zgr_z
  !!---------------------------------------------------------------------
  !!                  ***  ROUTINE zgr_z  ***
  !!
  !! ** Purpose :    set the depth of model levels and the resulting 
  !!         vertical scale factors.
  !!
  !! ** Method  :   from NEMO 3.6
  !!
  !! References :  
  !!----------------------------------------------------------------------
 INTEGER  ::   jk                     ! dummy loop indices
 INTEGER  ::   inum = 27              !
 REAL(wp) ::   zt, zw                 ! temporary scalars
 REAL(wp) ::   zsur, za0, za1, zkth   ! Values set from parameters in
 REAL(wp) ::   zacr, zdzmin, zhmax    ! par_CONFIG_Rxx.h90
 REAL(wp) ::   zrefdep                ! depth of the reference level (~10m)
 REAL(wp) ::   za2, zkth2, zacr2      ! Values for optional double tanh function set from parameters 

 NAMELIST/namdom/ jpk, nn_bathy , rn_bathy, rn_e3zps_min, rn_e3zps_rat, nn_msh    , rn_hmin, &
      &             nn_acc   , rn_atfp     , rn_rdt      , rn_rdtmin,           &
      &             rn_rdtmax, rn_rdth     , nn_closea , ln_crs,    &
      &             jphgr_msh, &
      &             ppglam0, ppgphi0, ppe1_deg, ppe2_deg, ppe1_m, ppe2_m, &
      &             ppsur, ppa0, ppa1, ppkth, ppacr, ppdzmin, pphmax, ldbletanh, &
      &             ppa2, ppkth2, ppacr2
 !! -----------------------------------------------------------------------------------------
 OPEN(inum, FILE=cf_nam)
 READ(inum, namdom)
 CLOSE(inum)
 jpkm1=jpk-1

 ! build nbloc_szt and njbloct
 ALLOCATE(  nblock_szt(nbloc), nkbloct(nbloc) )
 i1 = jpk/nbloc
 i2 = MOD(jpk, nbloc )
 nblock_szt(1         :nbloc-i2) = i1
 nblock_szt(nbloc-i2+1:nbloc   ) = i1+1
 nkbloct(1) = 1
 DO ji=2, nbloc
    nkbloct(ji) = nkbloct(ji-1) + nblock_szt(ji-1)
 ENDDO

 ALLOCATE (gdept_1d(jpk), gdepw_1d(jpk), e3t_1d(jpk), e3w_1d(jpk) )

 ! Set variables from parameters
 ! ------------------------------
 zkth = ppkth       ;   zacr = ppacr
 zdzmin = ppdzmin   ;   zhmax = pphmax
 zkth2 = ppkth2     ;   zacr2 = ppacr2   ! optional (ldbletanh=T) double tanh parameters

 IF(   ppa1  == 999    .AND.  &
      &  ppa0  == 999  .AND.  &
      &  ppsur == 999         ) THEN
    !
    za1  = (  ppdzmin - pphmax / FLOAT(jpkm1)  )                                                      &
         & / ( TANH((1-ppkth)/ppacr) - ppacr/FLOAT(jpk-1) * (  LOG( COSH( (jpk - ppkth) / ppacr) )      &
         &                                                   - LOG( COSH( ( 1  - ppkth) / ppacr) )  )  )
    za0  = ppdzmin - za1 *              TANH( (1-ppkth) / ppacr )
    zsur =   - za0 - za1 * ppacr * LOG( COSH( (1-ppkth) / ppacr )  )
 ELSE
    za1 = ppa1 ;       za0 = ppa0 ;          zsur = ppsur
    za2 = ppa2                            ! optional (ldbletanh=T) double tanh parameter
 ENDIF

 IF( .NOT. ldbletanh ) THEN
    DO jk = 1, jpk
       zw = REAL( jk , wp )
       zt = REAL( jk , wp ) + 0.5_wp
       gdepw_1d(jk) = ( zsur + za0 * zw + za1 * zacr * LOG ( COSH( (zw-zkth) / zacr ) )  )
       gdept_1d(jk) = ( zsur + za0 * zt + za1 * zacr * LOG ( COSH( (zt-zkth) / zacr ) )  )
       e3w_1d  (jk) =          za0      + za1        * TANH(       (zw-zkth) / zacr   )
       e3t_1d  (jk) =          za0      + za1        * TANH(       (zt-zkth) / zacr   )
    END DO
 ELSE
    DO jk = 1, jpk
       zw = FLOAT( jk )
       zt = FLOAT( jk ) + 0.5_wp
       ! Double tanh function
       gdepw_1d(jk) = ( zsur + za0 * zw + za1 * zacr * LOG ( COSH( (zw-zkth ) / zacr  ) )    &
            &                             + za2 * zacr2* LOG ( COSH( (zw-zkth2) / zacr2 ) )  )
       gdept_1d(jk) = ( zsur + za0 * zt + za1 * zacr * LOG ( COSH( (zt-zkth ) / zacr  ) )    &
            &                             + za2 * zacr2* LOG ( COSH( (zt-zkth2) / zacr2 ) )  )
       e3w_1d  (jk) =          za0      + za1        * TANH(       (zw-zkth ) / zacr  )      &
            &                             + za2        * TANH(       (zw-zkth2) / zacr2 )
       e3t_1d  (jk) =          za0      + za1        * TANH(       (zt-zkth ) / zacr  )      &
            &                             + za2        * TANH(       (zt-zkth2) / zacr2 )
    END DO
 ENDIF
 gdepw_1d(1) = 0._wp                    ! force first w-level to be exactly at zero
 CALL CreateMeshZgrFile
 CALL CreateMaskFile

END SUBROUTINE zgr_z

SUBROUTINE zgr_zps
  !!----------------------------------------------------------------------
  !!                  ***  ROUTINE zgr_zps  ***
  !!                     
  !! ** Purpose :   the depth and vertical scale factor in partial step
  !!      z-coordinate case
  !!
  !! ** Method  :   Partial steps : computes the 3D vertical scale factors
  !!      of T-, U-, V-, W-, UW-, VW and F-points that are associated with
  !!      a partial step representation of bottom topography.
  !!----------------------------------------------------------------------
 INTEGER  ::   ji, jj, jk       ! dummy loop indices
 INTEGER  ::   ik, it           ! temporary integers
 LOGICAL  ::   ll_print         ! Allow  control print for debugging
 REAL(wp) ::   ze3tp , ze3wp    ! Last ocean level thickness at T- and W-points
 REAL(wp) ::   zdepwp, zdepth   ! Ajusted ocean depth to avoid too small e3t
 REAL(wp) ::   zmax             ! Maximum depth
 REAL(wp) ::   zdiff            ! temporary scalar
 REAL(wp) ::   zrefdep          ! temporary scalar
 REAL(wp), DIMENSION(:,:), ALLOCATABLE :: zdep
 !!----------------------------------------------------------------------


 zmax = gdepw_1d(jpk) + e3t_1d(jpk)        ! maximum depth (i.e. the last ocean level thickness <= 2*e3t_1d(jpkm1) )
 bathy(:,:) = MIN( zmax ,  bathy(:,:) )    ! bounded value of bathy (min already set at the end of zgr_bat)

 WHERE( bathy(:,:) == 0._wp )   
    mbathy(:,:) = 0       ! land  : set mbathy to 0
 ELSE WHERE                     
    mbathy(:,:) = jpkm1   ! ocean : initialize mbathy to the max ocean level
 END WHERE

 ! Compute mbathy for ocean points (i.e. the number of ocean levels)
 ! find the number of ocean levels such that the last level thickness
 ! is larger than the minimum of e3zps_min and e3zps_rat * e3t_1d (where
 ! e3t_1d is the reference level thickness
 DO jk = jpkm1, 1, -1
    zdepth = gdepw_1d(jk) + MIN( rn_e3zps_min, e3t_1d(jk)*rn_e3zps_rat )
    WHERE( 0._wp < bathy(:,:) .AND. bathy(:,:) <= zdepth )   mbathy(:,:) = jk-1
 END DO

 CALL zgr_bat_ctl

 ! compute mask from mbathy : for writing performance and memory management, each mask is computed
 !  and written on file sequentially. Thus repeating  tmask determination (cheap)
 DO jbloc=1, nbloc
    nbloc_sz = nblock_szt(jbloc)
    ALLOCATE ( tmask(npiglo, npjglo,nbloc_sz) )
    tmask(:,:,:) = 0._wp
    DO jk = 1, nbloc_sz
       ik = nkbloct(jbloc) +jk -1
       PRINT *, 'T Masks for jk = ', ik
       DO jj = 1, npjglo
          DO ji = 1, npiglo
             IF( REAL( mbathy(ji,jj) - ik, wp ) + 0.1_wp >= 0._wp )   tmask(ji,jj,jk) = 1._wp
          END DO
       END DO
       IF ( ik == 1 ) THEN
          ierr = NF90_PUT_VAR( ncmsk, id_tmsku, tmask(:,:,1),    start=(/1,1,1/), count=(/npiglo,npjglo,1/) )
          IF ( ierr /= NF90_NOERR ) THEN  ; PRINT *, NF90_STRERROR(ierr) ; STOP 99 ;
          ENDIF
       ENDIF
    ENDDO
    ierr = NF90_PUT_VAR( ncmsk, id_tmsk, tmask,    start=(/1,1,nkbloct(jbloc),1/), count=(/npiglo,npjglo,nbloc_sz,1/) )
    IF ( ierr /= NF90_NOERR ) THEN  ; PRINT *, NF90_STRERROR(ierr) ; STOP 99 ; 
    ENDIF
    DEALLOCATE ( tmask )
 ENDDO ! bloc

 DO jbloc=1, nbloc
    nbloc_sz = nblock_szt(jbloc)
    ALLOCATE ( tmask(npiglo, npjglo,nbloc_sz) )
    ALLOCATE ( umask(npiglo, npjglo,nbloc_sz) )
    tmask(:,:,:) = 0._wp
    umask(:,:,:) = 0._wp
    DO jk = 1, nbloc_sz
       ik = nkbloct(jbloc) +jk -1
       PRINT *, 'U Masks for jk = ', ik
       DO jj = 1, npjglo
          DO ji = 1, npiglo
             IF( REAL( mbathy(ji,jj) - ik, wp ) + 0.1_wp >= 0._wp )   tmask(ji,jj,jk) = 1._wp
          END DO
       END DO
       DO jj = 1, npjglo -1
          DO ji = 1, npiglo -1
             umask(ji,jj,jk) = tmask(ji,jj,jk) * tmask(ji+1,jj,jk  )
          ENDDO
       ENDDO
       IF ( ik == 1 ) THEN
          ierr = NF90_PUT_VAR( ncmsk, id_umsku, umask(:,:,1),    start=(/1,1,1/), count=(/npiglo,npjglo,1/) )
          IF ( ierr /= NF90_NOERR ) THEN  ; PRINT *, NF90_STRERROR(ierr) ; STOP 99 ;
          ENDIF
       ENDIF
    ENDDO
    ierr = NF90_PUT_VAR( ncmsk, id_umsk, umask,    start=(/1,1,nkbloct(jbloc),1/), count=(/npiglo,npjglo,nbloc_sz,1/) )
    IF ( ierr /= NF90_NOERR ) THEN  ; PRINT *, NF90_STRERROR(ierr) ; STOP 99 ;
    ENDIF
    DEALLOCATE ( tmask , umask)
 ENDDO
 DO jbloc=1, nbloc
    nbloc_sz = nblock_szt(jbloc)
    ALLOCATE ( tmask(npiglo, npjglo,nbloc_sz) )
    ALLOCATE ( vmask(npiglo, npjglo,nbloc_sz) )
    tmask(:,:,:) = 0._wp
    vmask(:,:,:) = 0._wp

    DO jk = 1, nbloc_sz
       ik = nkbloct(jbloc) +jk -1
       PRINT *, 'V Masks for jk = ', ik
       DO jj = 1, npjglo
          DO ji = 1, npiglo
             IF( REAL( mbathy(ji,jj) - ik, wp ) + 0.1_wp >= 0._wp )   tmask(ji,jj,jk) = 1._wp
          END DO
       END DO
       DO jj = 1, npjglo -1
          DO ji = 1, npiglo -1
             vmask(ji,jj,jk) = tmask(ji,jj,jk) * tmask(ji  ,jj+1,jk)
          ENDDO
       ENDDO
       IF ( ik == 1 ) THEN
          ierr = NF90_PUT_VAR( ncmsk, id_vmsku, vmask(:,:,1),    start=(/1,1,1/), count=(/npiglo,npjglo,1/) )
          IF ( ierr /= NF90_NOERR ) THEN  ; PRINT *, NF90_STRERROR(ierr) ; STOP 99 ;
          ENDIF
       ENDIF
    ENDDO
    ierr = NF90_PUT_VAR( ncmsk, id_vmsk, vmask,    start=(/1,1,nkbloct(jbloc),1/), count=(/npiglo,npjglo,nbloc_sz,1/) )
    IF ( ierr /= NF90_NOERR ) THEN  ; PRINT *, NF90_STRERROR(ierr) ; STOP 99 ;
    ENDIF
    DEALLOCATE ( tmask , vmask)
 ENDDO

 DO jbloc=1, nbloc
    nbloc_sz = nblock_szt(jbloc)
    ALLOCATE ( tmask(npiglo, npjglo,nbloc_sz) )
    ALLOCATE ( fmask(npiglo, npjglo,nbloc_sz) )
    tmask(:,:,:) = 0._wp
    fmask(:,:,:) = 0._wp

    DO jk = 1, nbloc_sz
       ik = nkbloct(jbloc) +jk -1
       PRINT *, 'F Masks for jk = ', ik
       DO jj = 1, npjglo
          DO ji = 1, npiglo
             IF( REAL( mbathy(ji,jj) - ik, wp ) + 0.1_wp >= 0._wp )   tmask(ji,jj,jk) = 1._wp
          END DO
       END DO
       DO jj = 1, npjglo -1
          DO ji = 1, npiglo -1
             fmask(ji,jj,jk) = tmask(ji,jj,jk  ) * tmask(ji+1,jj,jk  )   &
                  &  * tmask(ji,jj+1,jk) * tmask(ji+1,jj+1,jk)
          ENDDO
       ENDDO
       IF ( ik == 1 ) THEN
          ierr = NF90_PUT_VAR( ncmsk, id_fmsku, fmask(:,:,1),    start=(/1,1,1/), count=(/npiglo,npjglo,1/) )
          IF ( ierr /= NF90_NOERR ) THEN  ; PRINT *, NF90_STRERROR(ierr) ; STOP 99 ;
          ENDIF
       ENDIF
    ENDDO
    ierr = NF90_PUT_VAR( ncmsk, id_fmsk, fmask,    start=(/1,1,nkbloct(jbloc),1/), count=(/npiglo,npjglo,nbloc_sz,1/) )
    IF ( ierr /= NF90_NOERR ) THEN  ; PRINT *, NF90_STRERROR(ierr) ; STOP 99 ;
    ENDIF
    DEALLOCATE ( tmask , fmask)
 ENDDO

 ierr = NF90_CLOSE(ncmsk)

 ! Now work for vertical metrics
 DO jbloc=1,nbloc
    PRINT *, ' WORKING for jbloc-slab ', jbloc
    nbloc_sz = nbloc_szt(jbloc)
    ALLOCATE ( gdept_0(npiglo,nbloc_sz, jpk) )
    ALLOCATE ( gdepw_0(npiglo,nbloc_sz, jpk) )
    ALLOCATE ( zdep(npiglo,nbloc_sz) )
    ALLOCATE ( e3t_0(npiglo,nbloc_sz, jpk) )
    ALLOCATE ( e3w_0(npiglo,nbloc_sz, jpk) )
    DO jk=1, jpk
       gdept_0(:,:,jk) = gdept_1d(jk)
       gdepw_0(:,:,jk) = gdepw_1d(jk)
       e3t_0  (:,:,jk) = e3t_1d  (jk)
       e3w_0  (:,:,jk) = e3w_1d  (jk)
    ENDDO
    DO jj=1, nbloc_sz
       ij=njbloct(jbloc) + jj - 1
       DO ji=1, npiglo
          ik = mbathy(ji,ij)
          IF( ik > 0 ) THEN               ! ocean point only
             ! max ocean level case
             IF( ik == jpkm1 ) THEN
                zdepwp = bathy(ji,ij)
                ze3tp  = bathy(ji,ij) - gdepw_1d(ik)
                ze3wp = 0.5_wp * e3w_1d(ik) * ( 1._wp + ( ze3tp/e3t_1d(ik) ) )
                e3t_0(ji,jj,ik  ) = ze3tp
                e3t_0(ji,jj,ik+1) = ze3tp
                e3w_0(ji,jj,ik  ) = ze3wp
                e3w_0(ji,jj,ik+1) = ze3tp
                gdepw_0(ji,jj,ik+1) = zdepwp
                gdept_0(ji,jj,ik  ) = gdept_1d(ik-1) + ze3wp
                gdept_0(ji,jj,ik+1) = gdept_0(ji,jj,ik) + ze3tp
                !
             ELSE                         ! standard case
                IF( bathy(ji,ij) <= gdepw_1d(ik+1) ) THEN  ;   gdepw_0(ji,jj,ik+1) = bathy(ji,ij)
                ELSE                                       ;   gdepw_0(ji,jj,ik+1) = gdepw_1d(ik+1)
                ENDIF
                !       ... on ik
                gdept_0(ji,jj,ik) = gdepw_1d(ik) + ( gdepw_0   (ji,jj,ik+1) - gdepw_1d(ik) )   &
                     &                          * ((gdept_1d(     ik  ) - gdepw_1d(ik) )   &
                     &                           /( gdepw_1d(     ik+1) - gdepw_1d(ik) ))
                e3t_0(ji,jj,ik) = e3t_1d (ik)    * ( gdepw_0   (ji,jj,ik+1) - gdepw_1d(ik) )   &
                     &                          / ( gdepw_1d(     ik+1) - gdepw_1d(ik) )
                e3w_0(ji,jj,ik) = 0.5_wp * ( gdepw_0(ji,jj,ik+1) + gdepw_1d(ik+1) - 2._wp * gdepw_1d(ik) )   &
                     &                  * ( e3w_1d(ik) / ( gdepw_1d(ik+1) - gdepw_1d(ik) ) )
                !       ... on ik+1
                e3w_0  (ji,jj,ik+1) = e3t_0  (ji,jj,ik)
                e3t_0  (ji,jj,ik+1) = e3t_0  (ji,jj,ik)
                gdept_0(ji,jj,ik+1) = gdept_0(ji,jj,ik) + e3t_0(ji,jj,ik)
             ENDIF
          ENDIF
       END DO
    ENDDO

    ! write mesh_zgr here with gdept_0, gdepw_0, e3t_0, e3w_0 only
    ! do it quick and dirty at this time ...
    ! write vertical slab (maybe slow but save lot of memory
    DO jj=1, nbloc_sz
       !         ij=(jbloc -1 )*nbloc_sz +jj
       ij=njbloct(jbloc) + jj - 1
       DO ji = 1, npiglo
          ik=MAX( mbathy(ji,ij),1)
          zdep(ji,jj)=gdept_0(ji,jj,ik) 
       ENDDO
    ENDDO
    ierr = NF90_PUT_VAR( nczgr, id_hdept, zdep,    start=(/1,(jbloc-1)*nbloc_sz+1  ,1/), count=(/npiglo,nbloc_sz    ,1/) )

    DO jj=1, nbloc_sz
       !         ij=(jbloc -1 )*nbloc_sz +jj
       ij=njbloct(jbloc) + jj - 1
       DO ji = 1, npiglo
          ik=MAX( mbathy(ji,ij),1)
          zdep(ji,jj)=gdepw_0(ji,jj,ik+1) 
       ENDDO
    ENDDO
    ierr = NF90_PUT_VAR( nczgr, id_hdepw, zdep,    start=(/1,(jbloc-1)*nbloc_sz+1  ,1/), count=(/npiglo,nbloc_sz    ,1/) )

    ierr = NF90_PUT_VAR( nczgr, id_e3t, e3t_0,    start=(/1,(jbloc-1)*nbloc_sz+1,1,1/), count=(/npiglo,nbloc_sz,jpk,1/) )
    ierr = NF90_PUT_VAR( nczgr, id_e3w, e3w_0,    start=(/1,(jbloc-1)*nbloc_sz+1,1,1/), count=(/npiglo,nbloc_sz,jpk,1/) )

    DEALLOCATE ( gdept_0,  gdepw_0, zdep, e3t_0,  e3w_0 )
 END DO
 ierr = NF90_CLOSE(nczgr )

 ! re open mesh_zgr in order to read e3t, e3w horizontally and compute e3u, e3v level by level
 !  use rlon, rlat as temporary arrays for e3t, e3w, e3u, e3v ...
 ierr = NF90_OPEN(cf_zgr,NF90_WRITE,nczgr )
 ierr = NF90_INQ_VARID(nczgr,'e3t',id_e3t )
 ierr = NF90_INQ_VARID(nczgr,'e3w',id_e3w ) 
 ierr = NF90_INQ_VARID(nczgr,'e3u',id_e3u )
 ierr = NF90_INQ_VARID(nczgr,'e3v',id_e3v )

 DO jbloc=1, nbloc
    nbloc_sz = nblock_szt(jbloc)
    ALLOCATE ( e3t_0(npiglo, npjglo,nbloc_sz) )
    ALLOCATE ( e3u_0 (npiglo, npjglo,nbloc_sz) )
    ierr= NF90_GET_VAR( nczgr, id_e3t, e3t_0, start=(/1,1,nkbloct(jbloc),1/), count=(/npiglo,npjglo,nbloc_sz,1/) )

    DO jk = 1, nbloc_sz
       ik = nkbloct(jbloc) +jk -1
       PRINT *,' e3u at level ', ik
       e3u_0(:,:,jk) =  e3t_1d(ik)
       DO jj=1,npjglo - 1
          DO ji=1,npiglo -1
             e3u_0(ji,jj,jk) = MIN( e3t_0(ji,jj,jk), e3t_0(ji+1,jj,jk) )
          ENDDO
       ENDDO
       !         ierr = NF90_PUT_VAR(nczgr, id_e3u, e3u_0, start=(/1,1,ik,1/), count=(/npiglo,npjglo,1,1/) )
    ENDDO
    ierr = NF90_PUT_VAR(nczgr, id_e3u, e3u_0, start=(/1,1,nkbloct(jbloc),1/), count=(/npiglo,npjglo,nbloc_sz,1/) )
    DEALLOCATE (e3t_0, e3u_0 )
 ENDDO

 DO jbloc=1, nbloc
    nbloc_sz = nblock_szt(jbloc)
    ALLOCATE ( e3t_0(npiglo, npjglo,nbloc_sz) )
    ALLOCATE ( e3v_0 (npiglo, npjglo,nbloc_sz) )
    ierr= NF90_GET_VAR( nczgr, id_e3t, e3t_0, start=(/1,1,nkbloct(jbloc),1/), count=(/npiglo,npjglo,nbloc_sz,1/) )

    DO jk = 1, nbloc_sz
       ik = nkbloct(jbloc) +jk -1
       PRINT *,' e3v at level ', ik
       e3v_0(:,:,jk) =  e3t_1d(ik)
       DO jj=1,npjglo - 1
          DO ji=1,npiglo -1
             e3v_0(ji,jj,jk) = MIN( e3t_0(ji,jj,jk), e3t_0(ji,jj+1,jk) )
          ENDDO
       ENDDO
       PRINT * , 'e3v_0 (81,2)  e3t_0(81,2), (81,3)', e3v_0(81,2,jk), e3t_0(81,2,jk) ,e3t_0(81,3,jk) 
       !          ierr = NF90_PUT_VAR(nczgr, id_e3v, e3v_0, start=(/1,1,ik,1/), count=(/npiglo,npjglo,1,1/) )
    ENDDO
    ierr = NF90_PUT_VAR(nczgr, id_e3v, e3v_0, start=(/1,1,nkbloct(jbloc),1/), count=(/npiglo,npjglo,nbloc_sz,1/) )
    DEALLOCATE (e3t_0, e3v_0 )
 ENDDO

 ierr = NF90_CLOSE(nczgr )

END SUBROUTINE zgr_zps

SUBROUTINE zgr_bat_ctl
  !!----------------------------------------------------------------------
  !!                    ***  ROUTINE zgr_bat_ctl  ***
  !!
  !! ** Purpose :   check the bathymetry in levels
  !!
  !! ** Method  :   The array mbathy is checked to verified its consistency
  !!      with the model options. in particular:
  !!            mbathy must have at least 1 land grid-points (mbathy<=0)
  !!                  along closed boundary.
  !!            mbathy must be cyclic IF jperio=1.
  !!            mbathy must be lower or equal to jpk-1.
  !!            isolated ocean grid points are suppressed from mbathy
  !!                  since they are only connected to remaining
  !!                  ocean through vertical diffusion.
  !!      C A U T I O N : mbathy will be modified during the initializa-
  !!      tion phase to become the number of non-zero w-levels of a water
  !!      column, with a minimum value of 1.
  !!
  !! ** Action  : - update mbathy: level bathymetry (in level index)
  !!              - update bathy : meter bathymetry (in meters)
  !!----------------------------------------------------------------------
  !!
 INTEGER ::   ji, jj, jl                    ! dummy loop indices
 INTEGER ::   icompt, ibtest, ikmax         ! temporary integers
 !!----------------------------------------------------------------------
 icompt = 0
 DO jl = 1, 2
    IF( nperio == 1 .OR. nperio  ==  4 .OR. nperio  ==  6 ) THEN
       mbathy( 1 ,:) = mbathy(npiglo-1,:)           ! local domain is cyclic east-west
       mbathy(npiglo,:) = mbathy(  2  ,:)
    ENDIF
    DO jj = 2, npjglo-1
       DO ji = 2, npiglo-1
          ibtest = MAX(  mbathy(ji-1,jj), mbathy(ji+1,jj),   &
               &           mbathy(ji,jj-1), mbathy(ji,jj+1)  )
          IF( ibtest < mbathy(ji,jj) ) THEN
             PRINT *, ' the number of ocean level at ',   &
                  &   'grid-point (i,j) =  ',ji,jj,' is changed from ', mbathy(ji,jj),' to ', ibtest
             mbathy(ji,jj) = ibtest
             icompt = icompt + 1
          ENDIF
       END DO
    END DO
 END DO
 PRINT *, icompt,' ocean grid points suppressed'

 !                                          ! East-west cyclic boundary conditions
 IF( nperio == 0 ) THEN
    PRINT *, ' mbathy set to 0 along east and west boundary: nperio = ', nperio
    mbathy( 1 ,:) = 0
    mbathy(npiglo,:) = 0
 ELSEIF( nperio == 1 .OR. nperio == 4 .OR. nperio ==  6 ) THEN
    PRINT *, ' east-west cyclic boundary conditions on mbathy: nperio = ', nperio
    mbathy( 1 ,:) = mbathy(npiglo-1,:)
    mbathy(npiglo,:) = mbathy(  2  ,:)
 ELSEIF( nperio == 2 ) THEN
    PRINT *, '   equatorial boundary conditions on mbathy: nperio = ', nperio
 ELSE
    PRINT *, '    e r r o r'
    PRINT *, '    parameter , nperio = ', nperio
    STOP 99 
 ENDIF

 ! write mbathy to file mesh_zgr
 ierr = NF90_PUT_VAR( nczgr, id_mbat, mbathy, start=(/1,1,1/), count=(/npiglo,npjglo,1/) )

END SUBROUTINE zgr_bat_ctl

SUBROUTINE CreateMeshZgrFile
  !netcdf ORCA12.L75-MAL83_mesh_zgr {
  !dimensions:
  !	x = 4322 ;
  !	y = 3059 ;
  !	z = 75 ;
  !	t = UNLIMITED ; // (1 currently)
  !variables:
  !	float nav_lon(y, x) ;
  !	float nav_lat(y, x) ;
  !	float nav_lev(z) ;
  !	float time_counter(t) ;
  !	float gdept_0(t, z) ;
  !	float gdepw_0(t, z) ;
  !	float e3t_0(t, z) ;
  !	float e3w_0(t, z) ;
  !	float mbathy(t, y, x) ;
  !	float hdept(t, y, x) ;
  !	float hdepw(t, y, x) ;
  !	float e3t(t, z, y, x) ;
  !	float e3u(t, z, y, x) ;
  !	float e3v(t, z, y, x) ;
  !	float e3w(t, z, y, x) ;
  !  ierr= NF90_CREATE(cf_zgr, or(NF90_CLOBBER,NF90_64BIT_OFFSET), nczgr) 

 ierr= NF90_CREATE(cf_zgr, or(NF90_CLOBBER,NF90_NETCDF4), nczgr) 
 ierr= NF90_DEF_DIM(nczgr, 'x', npiglo, idx)
 ierr= NF90_DEF_DIM(nczgr, 'y', npjglo, idy)
 ierr= NF90_DEF_DIM(nczgr, 'z', jpk,    idz)
 ierr= NF90_DEF_DIM(nczgr, 't', NF90_UNLIMITED,  idt)

 ierr=NF90_DEF_VAR(nczgr, 'nav_lon',       NF90_FLOAT, (/idx,idy/), id_navlon )
 ierr=NF90_DEF_VAR(nczgr, 'nav_lat',       NF90_FLOAT, (/idx,idy/), id_navlat )
 ierr=NF90_DEF_VAR(nczgr, 'nav_lev',       NF90_FLOAT, (/idz/)    , id_navlev )
 ierr=NF90_DEF_VAR(nczgr, 'time_counter',  NF90_FLOAT, (/idt/)    , id_time   )

 ierr=NF90_DEF_VAR(nczgr, 'gdept_0',       NF90_FLOAT, (/idz,idt/), id_dept1d )
 ierr=NF90_DEF_VAR(nczgr, 'gdepw_0',       NF90_FLOAT, (/idz,idt/), id_depw1d )
 ierr=NF90_DEF_VAR(nczgr, 'e3t_0',         NF90_FLOAT, (/idz,idt/), id_e3t1d  )
 ierr=NF90_DEF_VAR(nczgr, 'e3w_0',         NF90_FLOAT, (/idz,idt/), id_e3w1d  )

 ierr=NF90_DEF_VAR(nczgr, 'mbathy',        NF90_FLOAT, (/idx,idy,idt/), id_mbat )
 ierr=NF90_DEF_VAR(nczgr, 'hdept',         NF90_FLOAT, (/idx,idy,idt/), id_hdept )
 ierr=NF90_DEF_VAR(nczgr, 'hdepw',         NF90_FLOAT, (/idx,idy,idt/), id_hdepw )

 ierr=NF90_DEF_VAR(nczgr, 'e3t',           NF90_FLOAT, (/idx,idy,idz,idt/), id_e3t )
 ierr=NF90_DEF_VAR(nczgr, 'e3w',           NF90_FLOAT, (/idx,idy,idz,idt/), id_e3w )
 ierr=NF90_DEF_VAR(nczgr, 'e3u',           NF90_FLOAT, (/idx,idy,idz,idt/), id_e3u )
 ierr=NF90_DEF_VAR(nczgr, 'e3v',           NF90_FLOAT, (/idx,idy,idz,idt/), id_e3v )

 ierr = NF90_ENDDEF(nczgr)
 ! put dimension related variables
 ierr = NF90_OPEN(cf_coo, NF90_NOCLOBBER, nccoo )
 ierr = NF90_INQ_VARID(nccoo, cn_glamt, idvar ) ; ierr=NF90_GET_VAR(nccoo, idvar, rlon )
 ierr = NF90_INQ_VARID(nccoo, cn_gphit, idvar ) ; ierr=NF90_GET_VAR(nccoo, idvar, rlat )
 ierr = NF90_CLOSE(nccoo )

 ierr = NF90_PUT_VAR(nczgr, id_navlon, rlon)
 ierr = NF90_PUT_VAR(nczgr, id_navlat, rlat)
 ierr = NF90_PUT_VAR(nczgr, id_navlev, gdept_1d )
 ierr = NF90_PUT_VAR(nczgr, id_time  , (/0./) )

 ! put 1D vertical variables (reference depth)
 ierr = NF90_PUT_VAR(nczgr, id_dept1d, gdept_1d, start=(/1,1/), count=(/jpk,1/) )
 ierr = NF90_PUT_VAR(nczgr, id_depw1d, gdepw_1d, start=(/1,1/), count=(/jpk,1/) )
 ierr = NF90_PUT_VAR(nczgr, id_e3t1d,  e3t_1d,   start=(/1,1/), count=(/jpk,1/) )
 ierr = NF90_PUT_VAR(nczgr, id_e3w1d,  e3w_1d,   start=(/1,1/), count=(/jpk,1/) )

END SUBROUTINE CreateMeshZgrFile

SUBROUTINE CreateMaskFile
  !netcdf ORCA025.L75-MJM101.1_byte_mask {
  !dimensions:
  !	x = 1442 ;
  !	y = 1021 ;
  !	z = 75 ;
  !	t = UNLIMITED ; // (1 currently)
  !variables:
  !	float nav_lon(y, x) ;
  !	float nav_lat(y, x) ;
  !	float nav_lev(z) ;
  !	float time_counter(t) ;
  !	byte tmaskutil(t, y, x) ;
  !	byte umaskutil(t, y, x) ;
  !	byte vmaskutil(t, y, x) ;
  !	byte fmaskutil(t, y, x) ;
  !	byte tmask(t, z, y, x) ;
  !	byte umask(t, z, y, x) ;
  !	byte vmask(t, z, y, x) ;
  !	byte fmask(t, z, y, x) ;

 ierr= NF90_CREATE(cf_msk, or(NF90_CLOBBER,NF90_NETCDF4), ncmsk) 
 IF ( ierr /= NF90_NOERR ) THEN  ; PRINT *, NF90_STRERROR(ierr) ; STOP 99 ;
 ENDIF
 ierr= NF90_DEF_DIM(ncmsk, 'x', npiglo, idx)
 IF ( ierr /= NF90_NOERR ) THEN  ; PRINT *, NF90_STRERROR(ierr) ; STOP 99 ;
 ENDIF
 ierr= NF90_DEF_DIM(ncmsk, 'y', npjglo, idy)
 IF ( ierr /= NF90_NOERR ) THEN  ; PRINT *, NF90_STRERROR(ierr) ; STOP 99 ;
 ENDIF
 ierr= NF90_DEF_DIM(ncmsk, 'z', jpk,    idz)
 IF ( ierr /= NF90_NOERR ) THEN  ; PRINT *, NF90_STRERROR(ierr) ; STOP 99 ;
 ENDIF
 ierr= NF90_DEF_DIM(ncmsk, 't', NF90_UNLIMITED,  idt)
 IF ( ierr /= NF90_NOERR ) THEN  ; PRINT *, NF90_STRERROR(ierr) ; STOP 99 ;
 ENDIF

 ierr=NF90_DEF_VAR(ncmsk, 'nav_lon',       NF90_FLOAT, (/idx,idy/), id_navlon )
 IF ( ierr /= NF90_NOERR ) THEN  ; PRINT *, NF90_STRERROR(ierr) ; STOP 99 ;
 ENDIF
 ierr=NF90_DEF_VAR(ncmsk, 'nav_lat',       NF90_FLOAT, (/idx,idy/), id_navlat )
 IF ( ierr /= NF90_NOERR ) THEN  ; PRINT *, NF90_STRERROR(ierr) ; STOP 99 ;
 ENDIF
 ierr=NF90_DEF_VAR(ncmsk, 'nav_lev',       NF90_FLOAT, (/idz/)    , id_navlev )
 IF ( ierr /= NF90_NOERR ) THEN  ; PRINT *, NF90_STRERROR(ierr) ; STOP 99 ;
 ENDIF
 ierr=NF90_DEF_VAR(ncmsk, 'time_counter',  NF90_FLOAT, (/idt/)    , id_time   )
 IF ( ierr /= NF90_NOERR ) THEN  ; PRINT *, NF90_STRERROR(ierr) ; STOP 99 ;
 ENDIF

 ierr=NF90_DEF_VAR(ncmsk, 'tmaskutil',    NF90_BYTE, (/idx,idy,idt/), id_tmsku )
 IF ( ierr /= NF90_NOERR ) THEN  ; PRINT *, NF90_STRERROR(ierr) ; STOP 99 ;
 ENDIF
 ierr=NF90_DEF_VAR(ncmsk, 'umaskutil',    NF90_BYTE, (/idx,idy,idt/), id_umsku )
 IF ( ierr /= NF90_NOERR ) THEN  ; PRINT *, NF90_STRERROR(ierr) ; STOP 99 ;
 ENDIF
 ierr=NF90_DEF_VAR(ncmsk, 'vmaskutil',    NF90_BYTE, (/idx,idy,idt/), id_vmsku )
 IF ( ierr /= NF90_NOERR ) THEN  ; PRINT *, NF90_STRERROR(ierr) ; STOP 99 ;
 ENDIF
 ierr=NF90_DEF_VAR(ncmsk, 'fmaskutil',    NF90_BYTE, (/idx,idy,idt/), id_fmsku )
 IF ( ierr /= NF90_NOERR ) THEN  ; PRINT *, NF90_STRERROR(ierr) ; STOP 99 ;
 ENDIF

 ierr=NF90_DEF_VAR(ncmsk, 'tmask',         NF90_BYTE, (/idx,idy,idz,idt/), id_tmsk )
 IF ( ierr /= NF90_NOERR ) THEN  ; PRINT *, NF90_STRERROR(ierr) ; STOP 99 ; 
 ENDIF
 ierr=NF90_DEF_VAR(ncmsk, 'umask',         NF90_BYTE, (/idx,idy,idz,idt/), id_umsk )
 IF ( ierr /= NF90_NOERR ) THEN  ; PRINT *, NF90_STRERROR(ierr) ; STOP 99 ;
 ENDIF
 ierr=NF90_DEF_VAR(ncmsk, 'vmask',         NF90_BYTE, (/idx,idy,idz,idt/), id_vmsk )
 IF ( ierr /= NF90_NOERR ) THEN  ; PRINT *, NF90_STRERROR(ierr) ; STOP 99 ;
 ENDIF
 ierr=NF90_DEF_VAR(ncmsk, 'fmask',         NF90_BYTE, (/idx,idy,idz,idt/), id_fmsk )

 ierr = NF90_ENDDEF(ncmsk)
 IF ( ierr /= NF90_NOERR ) THEN  ; PRINT *, NF90_STRERROR(ierr) ; STOP 99 ; 
 ENDIF
 ! put dimension related variables

 ierr = NF90_PUT_VAR(ncmsk, id_navlon, rlon)
 IF ( ierr /= NF90_NOERR ) THEN  ; PRINT *, NF90_STRERROR(ierr) ; STOP 99 ;
 ENDIF
 ierr = NF90_PUT_VAR(ncmsk, id_navlat, rlat)
 IF ( ierr /= NF90_NOERR ) THEN  ; PRINT *, NF90_STRERROR(ierr) ; STOP 99 ;
 ENDIF
 ierr = NF90_PUT_VAR(ncmsk, id_navlev, gdept_1d )
 IF ( ierr /= NF90_NOERR ) THEN  ; PRINT *, NF90_STRERROR(ierr) ; STOP 99 ;
 ENDIF
 ierr = NF90_PUT_VAR(ncmsk, id_time  , (/0./) )
 IF ( ierr /= NF90_NOERR ) THEN  ; PRINT *, NF90_STRERROR(ierr) ; STOP 99 ;
 ENDIF

END SUBROUTINE CreateMaskFile


END PROGRAM
