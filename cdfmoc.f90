PROGRAM cdfmoc
  !!======================================================================
  !!                     ***  PROGRAM  cdfmoc  ***
  !!=====================================================================
  !!  ** Purpose : Compute the Meridional Overturning Cell (MOC)
  !!
  !!  ** Method  : The MOC is computed from the V velocity field, integrated
  !!               from the bottom to the surface, then zonally averaged with
  !!               eventual masking for oceanic basins.
  !!               The program looks for the file "new_maskglo.nc". If it 
  !!               does not exist, only the calculation over all the domain
  !!               is performed (this is adequate for a basin configuration).
  !!               In new_maskglo.nc the masking corresponds to the global
  !!               configuration. MOC for Global, Atlantic, Indo-Pacific, 
  !!               Indian, Pacific ocean.
  !!               Results are saved on moc.nc file with variables name 
  !!               respectively zomsfglo, zomsfatl, zomsfinp, zomsfind, zomsfpac
  !!
  !! History : 2.1  : 07/2005  : J.M. Molines  : Original code
  !!                : 04/2006  : A.M. Treguier : Adaptation to NATL4 case
  !!                : 09/2007  : G. Smith      : MOC decomposition
  !!                : 01/2008  : A. Lecointre  : MOC decomposition adaptation 
  !!           3.0  : 03/2011  : J.M. Molines  : Merge all MOC prog, Doctor norm + Lic.
  !!
  !! References :  For MOC decomposition : Lee & Marotzke (1998), 
  !!               Baehr, Hirschi, Beismann &  Marotzke (2004),
  !!               Cabanes, Lee, & Fu (2007),  Koehl & Stammer (2007).
  !!               See also the powerpoint presentation by Tony Lee at the third 
  !!               CLIVAR-GSOP intercomparison  available at : 
  !!    http://www.clivar.org/organization/gsop/synthesis/mit/talks/lee_MOC_comparison.ppt
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  USE eos
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2011, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=2), DIMENSION(:,:,:), ALLOCATABLE :: ibmask       !  nbasins x npiglo x npjglo
  INTEGER(KIND=2), DIMENSION(:,:), ALLOCATABLE :: ivmask         ! ivmask (used to mask e3v)

  INTEGER(KIND=4)                             :: npglo, npatl, npinp
  INTEGER(KIND=4)                             :: npind, nppac
  INTEGER(KIND=4)                             :: jbasin, jj, jk  ! dummy loop index
  INTEGER(KIND=4)                             :: ji, jt          ! dummy loop index
  INTEGER(KIND=4)                             :: nbasins, ibasin ! number of sub basins
  INTEGER(KIND=4)                             :: ierr            ! working integer
  INTEGER(KIND=4)                             :: narg, iargc     ! command line browser
  INTEGER(KIND=4)                             :: ijarg, ii       !  "             "
  INTEGER(KIND=4)                             :: npiglo,npjglo   ! size of the domain
  INTEGER(KIND=4)                             :: npk, npt        ! size of the domain
  INTEGER(KIND=4)                             :: ncout           ! out put file id
  INTEGER(KIND=4)                             :: nvarout         ! number of output variables
  INTEGER(KIND=4)                             :: ijvar           ! index for output variable
  INTEGER(KIND=4), DIMENSION(:),  ALLOCATABLE :: ipk, id_varout  ! output variables info
  INTEGER(KIND=4), DIMENSION(2)               :: iloc            ! working integer array

  REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE :: e3v             ! Vertical e3v masked by vmask
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: e1v,      gphiv ! metrics, velocity
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: zv              ! meridional velocity
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: rdumlon         ! dummy longitude = 0.
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: rdumlat         ! latitude for i = north pole
  REAL(KIND=4), DIMENSION(:),     ALLOCATABLE :: gdepw           ! depthw
  REAL(KIND=4), DIMENSION(:),     ALLOCATABLE :: gdept           ! deptht
  REAL(KIND=4), DIMENSION(:),     ALLOCATABLE :: e31d            ! e3 1D : used if full step
  REAL(KIND=4), DIMENSION(:),     ALLOCATABLE :: tim             ! time counter array

  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: dmoc            ! nbasins x npjglo x npk

  CHARACTER(LEN=256)                          :: cf_vfil         ! meridional velocity file
  CHARACTER(LEN=256)                          :: cf_moc = 'moc.nc'  ! output file name
  CHARACTER(LEN=256)                          :: cglobal         ! Global attribute for output file
  CHARACTER(LEN=256)                          :: cldum           ! dummy char variable

  TYPE(variable) ,DIMENSION(:),   ALLOCATABLE :: stypvar         ! structure for attribute

  LOGICAL                                     :: lbas  = .FALSE. ! new_maskglo.nc file  flag
  LOGICAL                                     :: lfull = .FALSE. ! full step flag
  LOGICAL                                     :: lchk  = .FALSE. ! check for missing files
  LOGICAL                                     :: ldec  = .FALSE. ! check for missing files

  ! Variables used only when MOC decomposition is requested
  INTEGER(KIND=2), DIMENSION(:,:), ALLOCATABLE :: iumask         ! iumask (used if decomposition)
  INTEGER(KIND=2), DIMENSION(:,:), ALLOCATABLE :: itmask         ! itmask (used if decomposition)

  INTEGER(KIND=4)                             :: itmp, iup, ido  ! up and down index for work

  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: e1u             ! used if ldec
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: hdep            ! total depth at v point
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: zcoef           ! coefficient for geostrophic calc
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: ztemp           ! temperature
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: zsal            ! salinity
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: zsig0           ! density
  REAL(KIND=4)                                :: zmsv
  REAL(KIND=4)                                :: rpi             ! pi 
  REAL(KIND=4)                                :: grav = 9.81     ! gravity
  REAL(KIND=4)                                :: rau0 = 1025.    ! mean density

  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: dmoc_sh         ! nbasins x npjglo x npk
  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: dmoc_bt         ! nbasins x npjglo x npk
  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: dmoc_btw        ! nbasins x npjglo x npk
  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: dmoc_ag         ! nbasins x npjglo x npk
  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: dvgeo           !  npiglo x npjglo x 2
  REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: dvbt            ! Barotropic velocity
  REAL(KIND=8)                                :: dgeo            ! Barotropic velocity

  CHARACTER(LEN=256)                          :: cf_tfil         ! Grid T file (case of decomposition)
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfmoc  V_file [-full] [-decomp ] [T_file] '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Computes the MOC for oceanic sub basins as described '
     PRINT *,'       in ',TRIM(cn_fbasins)
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       V_file : file with meridional velocity component.'
     PRINT *,'       T_file : file with temperature and salinity'
     PRINT *,'               (required only for -decomp option).'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-full ] : use full step instead of default partial step' 
     PRINT *,'       [-decomp ] : decompose MOC in 3 components: Geostrophic,'
     PRINT *,'                 Barotropic,  Ageostrophic). For this option a '
     PRINT *,'                 gridT file is required.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       Files ',TRIM(cn_fhgr),' ', TRIM(cn_fhgr),' and ', TRIM(cn_fmsk)
     PRINT *,'       File ',TRIM(cn_fbasins),'. If this latter file is not available '
     PRINT *,'             only the MOC for the global domain is computed'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_moc)
     PRINT *,'       variables ',TRIM( cn_zomsfglo),' : Global ocean '
     PRINT *,'       variables ',TRIM( cn_zomsfatl),' : Atlantic Ocean '
     PRINT *,'       variables ',TRIM( cn_zomsfinp),' : Indo Pacific '
     PRINT *,'       variables ',TRIM( cn_zomsfind),' : Indian Ocean alone'
     PRINT *,'       variables ',TRIM( cn_zomsfpac),' : Pacific Ocean alone'
     PRINT *,'      '
     PRINT *,'       If decomposition is required , ( option -decomp ) add 3 additional'
     PRINT *,'       variables per basin with suffixes _sh, _bt, _ag.'
     STOP
  ENDIF

  cglobal = 'Partial step computation'
  ijarg = 1 ; ii = 0
  DO WHILE ( ijarg <= narg )
     CALL getarg (ijarg, cldum) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ('-full') 
        lfull   = .TRUE.
        cglobal = 'Full step computation'
     CASE ('-decomp') 
        ldec    = .TRUE.
     CASE DEFAULT
        ii=ii+1
        SELECT CASE (ii)
        CASE ( 1 ) ; cf_vfil = cldum
        CASE ( 2 ) ; cf_tfil = cldum
        CASE DEFAULT
           PRINT*, 'ERROR : Too many arguments ...'
           STOP
        END SELECT
     END SELECT
  END DO

  lchk = lchk .OR. chkfile ( cn_fhgr )
  lchk = lchk .OR. chkfile ( cn_fzgr )
  lchk = lchk .OR. chkfile ( cn_fmsk )
  lchk = lchk .OR. chkfile ( cf_vfil )
  IF ( ldec ) lchk = lchk .OR. chkfile ( TRIM(cf_tfil) ) 
  IF ( lchk ) STOP  ! missing file(s)

  npiglo = getdim (cf_vfil,cn_x)
  npjglo = getdim (cf_vfil,cn_y)
  npk    = getdim (cf_vfil,cn_z)
  npt    = getdim (cf_vfil,cn_t)

  PRINT *, 'Working with cdfmoc ...'
  PRINT *, '  npiglo =', npiglo
  PRINT *, '  npjglo =', npjglo
  PRINT *, '  npk    =', npk
  PRINT *, '  npt    =', npt

  !  Detects newmaskglo file 
  lbas = .NOT. chkfile (cn_fbasins )

  IF (lbas) THEN
     nbasins = 5
  ELSE
     nbasins = 1
  ENDIF

  IF ( ldec ) THEN
     nvarout=nbasins * 4   ! total, _sh, _bt, _ag
  ELSE
     nvarout=nbasins       ! total
  ENDIF

  ALLOCATE ( stypvar(nvarout), ipk(nvarout), id_varout(nvarout) )

  ! define new variables for output 
  !    all variables
  stypvar%cunits            = 'Sverdrup'
  stypvar%rmissing_value    = 99999.
  stypvar%valid_min         = -1000.
  stypvar%valid_max         =  1000.
  stypvar%scale_factor      = 1.
  stypvar%add_offset        = 0.
  stypvar%savelog10         = 0.
  stypvar%conline_operation = 'N/A'
  stypvar%caxis             = 'TZY'
  ipk(:) = npk  ! All variables are vertical slices 1 x npjglo x npk

  ii=1 ; ibasin=1
  PRINT *, 'Variable ',ii,' is zomsfglo'
  npglo=ibasin  ; ibasin = ibasin + 1 
  stypvar(ii)%cname          = TRIM(cn_zomsfglo)
  stypvar(ii)%clong_name     = 'Meridional_Overt.Cell_Global'
  stypvar(ii)%cshort_name    = TRIM(cn_zomsfglo)
  ii=ii+1

  IF ( ldec ) THEN
     PRINT *, 'Variable ',ii,' is zomsfglo_sh'
     stypvar(ii)%cname       = TRIM(cn_zomsfglo)//'_sh'
     stypvar(ii)%clong_name  = 'GeoShear_Merid_StreamFunction'
     stypvar(ii)%cshort_name = TRIM(cn_zomsfglo)//'_sh'
     ii= ii+1
     PRINT *, 'Variable ',ii,' is zomsfglo_bt'
     stypvar(ii)%cname       = TRIM(cn_zomsfglo)//'_bt'
     stypvar(ii)%clong_name  = 'Barotropic_Merid_StreamFunction'
     stypvar(ii)%cshort_name = TRIM(cn_zomsfglo)//'_bt'
     ii= ii+1
     PRINT *, 'Variable ',ii,' is zomsfglo_ag'
     stypvar(ii)%cname       = TRIM(cn_zomsfglo)//'_ag'
     stypvar(ii)%clong_name  = 'Ageostoph_Merid_StreamFunction'
     stypvar(ii)%cshort_name = TRIM(cn_zomsfglo)//'_ag'
     ii= ii+1
  ENDIF

  IF (lbas) THEN
     npatl=ibasin  ; ibasin = ibasin + 1
     PRINT *, 'Variable ',ii,' is zomsfatl'
     stypvar(ii)%cname       = TRIM(cn_zomsfatl)
     stypvar(ii)%clong_name  = 'Meridional_Overt.Cell_Atlantic'
     stypvar(ii)%cshort_name = TRIM(cn_zomsfatl)
     ii= ii+1
     IF ( ldec ) THEN 
        PRINT *, 'Variable ',ii,' is zomsfatl_sh'
        stypvar(ii)%cname       = TRIM(cn_zomsfatl)//'_sh'
        stypvar(ii)%clong_name  = 'GeoShear_Merid_StreamFunction_Atlantic'
        stypvar(ii)%cshort_name = TRIM(cn_zomsfatl)//'_sh'
        ii= ii+1
        PRINT *, 'Variable ',ii,' is zomsfatl_bt'
        stypvar(ii)%cname       = TRIM(cn_zomsfatl)//'_bt'
        stypvar(ii)%clong_name  = 'Barotropic_Merid_StreamFunction_Atlantic'
        stypvar(ii)%cshort_name = TRIM(cn_zomsfatl)//'_bt'
        ii= ii+1
        PRINT *, 'Variable ',ii,' is zomsfatl_ag'
        stypvar(ii)%cname       = TRIM(cn_zomsfatl)//'_ag'
        stypvar(ii)%clong_name  = 'Ageostroph_Merid_StreamFunction_Atlantic'
        stypvar(ii)%cshort_name = TRIM(cn_zomsfatl)//'_ag'
        ii= ii+1
     ENDIF

     npinp=ibasin  ; ibasin = ibasin + 1
     PRINT *, 'Variable ',ii,' is zomsfinp'
     stypvar(ii)%cname       = TRIM(cn_zomsfinp)
     stypvar(ii)%clong_name  = 'Meridional_Overt.Cell_IndoPacif'
     stypvar(ii)%cshort_name = TRIM(cn_zomsfinp)
     ii= ii+1

     IF ( ldec ) THEN
        PRINT *, 'Variable ',ii,' is zomsfinp_sh'
        stypvar(ii)%cname       = TRIM(cn_zomsfinp)//'_sh'
        stypvar(ii)%clong_name  = 'GeoShear_Merid_StreamFunction_IndoPacif'
        stypvar(ii)%cshort_name = TRIM(cn_zomsfinp)//'_sh'
        ii= ii+1
        PRINT *, 'Variable ',ii,' is zomsfinp_bt'
        stypvar(ii)%cname       = TRIM(cn_zomsfinp)//'_bt'
        stypvar(ii)%clong_name  = 'Barotropic_Merid_StreamFunction_IndoPacif'
        stypvar(ii)%cshort_name = TRIM(cn_zomsfinp)//'_bt'
        ii= ii+1
        PRINT *, 'Variable ',ii,' is zomsfinp_ag'
        stypvar(ii)%cname       = TRIM(cn_zomsfinp)//'_ag'
        stypvar(ii)%clong_name  = 'Ageostroph_Merid_StreamFunction_IndoPacif'
        stypvar(ii)%cshort_name = TRIM(cn_zomsfinp)//'_ag'
        ii= ii+1
     ENDIF

     npind=ibasin  ; ibasin = ibasin + 1
     PRINT *, 'Variable ',ii,' is zomsfind'
     stypvar(ii)%cname       = TRIM(cn_zomsfind)
     stypvar(ii)%clong_name  = 'Meridional_Overt.Cell_Indian'
     stypvar(ii)%cshort_name = TRIM(cn_zomsfind)
     ii= ii+1

     IF ( ldec ) THEN
        PRINT *, 'Variable ',ii,' is zomsfind_sh'
        stypvar(ii)%cname       = TRIM(cn_zomsfind)//'_sh'
        stypvar(ii)%clong_name  = 'GeoShear_Merid_StreamFunction_Indian'
        stypvar(ii)%cshort_name = TRIM(cn_zomsfind)//'_sh'
        ii= ii+1
        PRINT *, 'Variable ',ii,' is zomsfind_bt'
        stypvar(ii)%cname       = TRIM(cn_zomsfind)//'_bt'
        stypvar(ii)%clong_name  = 'Barotropic_Merid_StreamFunction_Indian'
        stypvar(ii)%cshort_name = TRIM(cn_zomsfind)//'_bt'
        ii= ii+1
        PRINT *, 'Variable ',ii,' is zomsfind_ag'
        stypvar(ii)%cname       = TRIM(cn_zomsfind)//'_ag'
        stypvar(ii)%clong_name  = 'Ageostroph_Merid_StreamFunction_Indian'
        stypvar(ii)%cshort_name = TRIM(cn_zomsfind)//'_ag'
        ii= ii+1
     ENDIF

     nppac=ibasin  ; ibasin = ibasin + 1
     PRINT *, 'Variable ',ii,' is zomsfpac'
     stypvar(ii)%cname       = TRIM(cn_zomsfpac)
     stypvar(ii)%clong_name  = 'Meridional_Overt.Cell_pacif'
     stypvar(ii)%cshort_name = TRIM(cn_zomsfpac)
     ii= ii+1

     IF ( ldec ) THEN
        PRINT *, 'Variable ',ii,' is zomsfpac_sh'
        stypvar(ii)%cname       = TRIM(cn_zomsfpac)//'_sh'
        stypvar(ii)%clong_name  = 'GeoShear_Merid_StreamFunction_Pacif'
        stypvar(ii)%cshort_name = TRIM(cn_zomsfpac)//'_sh'
        ii= ii+1
        PRINT *, 'Variable ',ii,' is zomsfpac_bt'
        stypvar(ii)%cname       = TRIM(cn_zomsfpac)//'_bt'
        stypvar(ii)%clong_name  = 'Barotropic_Merid_StreamFunction_Pacif'
        stypvar(ii)%cshort_name = TRIM(cn_zomsfpac)//'_bt'
        ii= ii+1
        PRINT *, 'Variable ',ii,' is zomsfpac_ag'
        stypvar(ii)%cname       = TRIM(cn_zomsfpac)//'_ag'
        stypvar(ii)%clong_name  = 'Ageostroph_Merid_StreamFunction_Pacif'
        stypvar(ii)%cshort_name = TRIM(cn_zomsfpac)//'_ag'
     ENDIF
  ENDIF

  ! Allocate arrays
  ALLOCATE ( ibmask(nbasins, npiglo, npjglo) )
  ALLOCATE ( zv(npiglo, npjglo), e1v(npiglo,npjglo), e3v(npiglo,npjglo,npk) )
  ALLOCATE ( gphiv(npiglo,npjglo) )
  ALLOCATE ( rdumlon(1,npjglo), rdumlat(1,npjglo))
  ALLOCATE ( gdepw(npk), gdept(npk), e31d(npk) )
  ALLOCATE ( tim(npt) )
  ALLOCATE ( dmoc( nbasins, npjglo, npk   ) )
  ALLOCATE ( ivmask(npiglo, npjglo) )
  IF ( ldec ) THEN 
     ALLOCATE ( iumask(npiglo, npjglo) )
     ALLOCATE ( itmask(npiglo, npjglo) )
     ALLOCATE ( ztemp(npiglo, npjglo) )
     ALLOCATE ( zsal(npiglo, npjglo) )
     ALLOCATE ( zsig0(npiglo, npjglo) )
     ALLOCATE ( e1u(npiglo, npjglo) )
     ALLOCATE ( zcoef(npiglo, npjglo) )
     ALLOCATE ( dvbt(npiglo, npjglo), hdep(npiglo,npjglo) )
     ALLOCATE ( dmoc_sh(nbasins, npjglo, npk) )
     ALLOCATE ( dmoc_bt(nbasins, npjglo, npk) )
     ALLOCATE ( dmoc_btw(nbasins, npjglo, npk) )
     ALLOCATE ( dmoc_ag(nbasins, npjglo, npk)  )
     ALLOCATE ( dvgeo(npiglo, npjglo, 2 ) )
  ENDIF

  e1v(:,:)   = getvar  (cn_fhgr, cn_ve1v,  1, npiglo,npjglo) 
  gphiv(:,:) = getvar  (cn_fhgr, cn_gphiv, 1, npiglo,npjglo)
  gdepw(:)   = getvare3(cn_fzgr, cn_gdepw, npk             )
  gdepw(:)   = -1.* gdepw(:)
  DO jk= 1, npk
     ! save e3v masked with vmask as 3d array
     e3v(:,:,jk) = get_e3v(jk)
  END DO

  IF ( ldec  ) gdept(:) = getvare3(cn_fzgr, cn_gdept, npk             )
  IF ( ldec  ) e1u(:,:) = getvar  (cn_fhgr, cn_ve1u,  1, npiglo,npjglo)
  IF ( lfull ) e31d(:)  = getvare3(cn_fzgr, cn_ve3t, npk)

  iloc=MAXLOC(gphiv)
  rdumlat(1,:) = gphiv(iloc(1),:)
  rdumlon(:,:) = 0.   ! set the dummy longitude to 0

  ! create output fileset
  ncout = create      ( cf_moc,  'none',    1, npjglo, npk, cdep=cn_vdepthw )
  ierr  = createvar   ( ncout,   stypvar,   nvarout,   ipk, id_varout, cdglobal=TRIM(cglobal)           )
  ierr  = putheadervar( ncout,   cf_vfil,   1, npjglo, npk, pnavlon=rdumlon, pnavlat=rdumlat, pdep=gdepw)
  tim   = getvar1d    ( cf_vfil, cn_vtimec, npt                    )
  ierr  = putvar1d    ( ncout,   tim,       npt, 'T')
 
  ! 1 : global ; 2 : Atlantic ; 3 : Indo-Pacif ; 4 : Indian ; 5 : Pacif
     ibmask(npglo,:,:) = getvar(cn_fmsk,   'vmask', 1, npiglo, npjglo)
  IF ( lbas ) THEN
     ibmask(npatl,:,:) = getvar(cn_fbasins, 'tmaskatl', 1, npiglo, npjglo)
     ibmask(npind,:,:) = getvar(cn_fbasins, 'tmaskind', 1, npiglo, npjglo)
     ibmask(nppac,:,:) = getvar(cn_fbasins, 'tmaskpac', 1, npiglo, npjglo)
     ibmask(npinp,:,:) = ibmask(nppac,:,:) + ibmask(npind,:,:)  ! indo pacific mask
     ! ensure that there are no overlapping on the masks
     WHERE(ibmask(npinp,:,:) > 0 ) ibmask(npinp,:,:) = 1
     ! change global mask for GLOBAL periodic condition
     ibmask(1,1,     :) = 0.
     ibmask(1,npiglo,:) = 0.
  ENDIF

  DO jt = 1, npt
     ! --------------------------
     ! 1) Compute total MOC: dmoc
     ! --------------------------
     dmoc(:,:,:) = 0.d0        ! initialize moc to 0
     IF ( ldec) THEN ; dvbt=0.d0 ; hdep=0.0 ; dmoc_bt=0.d0 ; ENDIF
     DO jk = 1, npk-1
        ! Get velocities v at jk, time = jt
        zv(:,:)= getvar(cf_vfil, cn_vomecrty,  jk, npiglo, npjglo, ktime=jt)

        IF ( ldec ) THEN
        ! compute barotropic component when requested
        ! this contribution is computed here in order to use zv(jk)
           dvbt(:,:)    = dvbt(:,:) + e3v(:,:,jk)*zv(:,:)*1.d0
           hdep(:,:)    = hdep(:,:) + e3v(:,:,jk)
        ENDIF

        ! integrates 'zonally' (along i-coordinate)
        DO ji=1,npiglo
           ! For all basins 
           DO jbasin = 1, nbasins
              DO jj=1,npjglo
                 dmoc(jbasin,jj,jk)=dmoc(jbasin,jj,jk) -  &
                      &             e1v(ji,jj)*e3v(ji,jj,jk)* ibmask(jbasin,ji,jj)*zv(ji,jj)*1.d0
              ENDDO
           END DO
        END DO
     END DO

     ! integrates vertically from bottom to surface
     DO jk = npk-1, 1, -1
        dmoc(:,:,jk)    = dmoc(:,:,jk+1)    + dmoc(:,:,jk)/1.d6
     END DO  

     IF ( ldec ) THEN
     !--------------------------------------------------
     ! 2) compute extra term if decomposition requested
     !--------------------------------------------------
     !  2.1 : Barotropic MOC : dmoc_bt
     !  """"""""""""""""""""
        ! compute vertical mean of the meridional velocity
        WHERE ( hdep /= 0 )
          dvbt(:,:) = dvbt(:,:) / hdep(:,:)
        ELSEWHERE
          dvbt(:,:) = 0.d0
        ENDWHERE

        DO jk=1, npk-1

           ! integrates 'zonally' (along i-coordinate)
           DO ji=1,npiglo
              ! For all basins
              DO jbasin = 1, nbasins
                 DO jj=1,npjglo
                    dmoc_bt(jbasin,jj,jk)=dmoc_bt(jbasin,jj,jk) -  &
                         &    e1v(ji,jj)*e3v(ji,jj,jk)* ibmask(jbasin,ji,jj)*dvbt(ji,jj)
                 ENDDO
              END DO
           END DO
        END DO
        ! integrates vertically   from bottom to surface
        DO jk = npk-1, 1, -1
           dmoc_bt(:,:,jk) = dmoc_bt(:,:,jk+1) + dmoc_bt(:,:,jk)/1.d6
        END DO  

     !  2.2 : Geostrophic Shear MOC : dmoc_sh
     !  """""""""""""""""""""""""""""
        ! using equation 2.7 of Lecointre (2008 
        ! f. Dv/Dz = -g/rau0. Drho/Dx 
        rau0 = 1025.0
        grav = 9.81
        rpi  = ACOS( -1.)
        zcoef(:,:) =  2*2*rpi/( 24.0 * 3600. )* SIN ( rpi * gphiv(:,:) /180.0) ! f at v point
        WHERE ( zcoef /= 0 ) 
           zcoef(:,:) = -grav/ rau0 / zcoef(:,:)
        ELSEWHERE
           zcoef(:,:) = 0.
        END WHERE

        dvgeo(:,:,:) = 0.0
        dvbt(:,:)    = 0.d0
        iup = 1 ; ido = 2
        DO jk=npk-1, 1, -1
           iumask(:,:) = getvar(cn_fmsk, 'umask', jk, npiglo, npjglo)
           itmask(:,:) = getvar(cn_fmsk, 'tmask', jk, npiglo, npjglo)
           ztemp(:,:)  = getvar(cf_tfil, cn_votemper, jk, npiglo, npjglo, ktime=jt )
           zsal(:,:)   = getvar(cf_tfil, cn_vosaline, jk, npiglo, npjglo, ktime=jt )
           zsig0(:,:)  = sigmai (ztemp, zsal, gdept(jk), npiglo, npjglo )* itmask(:,:)

           ! dgeo is Drho/dx at V point ( average on the 4 neighbours U points)
           ! thus, dgeo is -f.rau0/g. Dv/Dz 
           DO jj = 2, npjglo -1
              DO ji = 2, npiglo -1
                 zmsv =  1. / MAX (1, iumask(ji-1,jj+1)+iumask(ji,jj+1)+iumask(ji-1,jj)+iumask(ji,jj) )
                 dgeo = ( ( zsig0(ji,  jj+1) - zsig0(ji-1,jj+1) ) * iumask(ji-1, jj+1) / e1u(ji-1, jj+1) &
                      &  +( zsig0(ji+1,jj+1) - zsig0(ji  ,jj+1) ) * iumask(ji,   jj+1) / e1u(ji,   jj+1) &
                      &  +( zsig0(ji,  jj  ) - zsig0(ji-1,jj  ) ) * iumask(ji-1, jj  ) / e1u(ji-1, jj  ) &
                      &  +( zsig0(ji+1,jj  ) - zsig0(ji,  jj  ) ) * iumask(ji,   jj  ) / e1u(ji,   jj  )  )*1.d0
                 ! 
                 ! dvgeo is the geostrophic velocity at w point(jk) obtained by vertical integration of Dv/Dz
                 ! between bottom and jk
                 dvgeo(ji,jj,iup) = dvgeo(ji,jj,ido) + zcoef(ji,jj) * dgeo * zmsv * ibmask(npglo,ji,jj) *e3v(ji,jj,jk)
                 ! zv is the geostrophic velocity located at v-level (jk)
                 zv(ji,jj) = 0.5 *( dvgeo(ji,jj,iup) + dvgeo(ji,jj,ido) )
              ENDDO
           ENDDO
           ! compute the vertical mean of geostrophic velocity
           ! for memory management purpose we re-use dvbt which is not used any longer.
           dvbt(:,:) = dvbt(:,:) + e3v(:,:,jk)*zv(:,:)*1.d0

           ! integrates 'zonally' (along i-coordinate)
           DO ji=1,npiglo
              ! For all basins
              DO jbasin = 1, nbasins
                 DO jj=1,npjglo
                    dmoc_sh(jbasin,jj,jk)=dmoc_sh(jbasin,jj,jk) -  &
                         &             e1v(ji,jj)*e3v(ji,jj,jk)* ibmask(jbasin,ji,jj)*zv(ji,jj)*1.d0
                 ENDDO
              END DO
           END DO
           ! swap up and down for next level computation
           itmp=iup ;  iup = ido ; ido = itmp
        ENDDO   ! end of level loop

        WHERE ( hdep /=0 )
          dvbt(:,:) = dvbt(:,:) / hdep(:,:)
        ELSEWHERE
          dvbt(:,:) = 0.d0
        END WHERE

     !  2.2.1 : Barotropic Geostrophic Shear MOC : dmoc_btw
     !  """"""""""""""""""""""""""""""""""""""""""
        ! compute corresponding MOC for this unwanted pseudo barotropic contribution
        dmoc_btw(:,:,:) = 0.d0
        DO jk=1, npk-1

           ! integrates 'zonally' (along i-coordinate)
           DO ji=1,npiglo
              ! For all basins
              DO jbasin = 1, nbasins
                 DO jj=1,npjglo
                    dmoc_btw(jbasin,jj,jk)=dmoc_btw(jbasin,jj,jk) -  &
                         &         e1v(ji,jj)*e3v(ji,jj,jk)* ibmask(jbasin,ji,jj)*dvbt(ji,jj)
                 ENDDO
              END DO
           END DO
        END DO

        ! apply correction to dmoc_sh
        dmoc_sh(:,:,:) = dmoc_sh(:,:,:) - dmoc_btw(:,:,:)

        ! integrates vertically   from bottom to surface
        DO jk = npk-1, 1, -1
           dmoc_sh(:,:,jk) = dmoc_sh(:,:,jk+1) + dmoc_sh(:,:,jk)/1.e6
        END DO  ! 

     !  2.3 : Barotropic Geostrophic Shear MOC : dmoc_ag
     ! ----------------------------------------
        ! compute ageostrophic component 
        !  AGEO        =   MOC total    Geo-Shear        Barotropic
        dmoc_ag(:,:,:) = dmoc(:,:,:) - dmoc_sh(:,:,:) - dmoc_bt(:,:,:)
     ENDIF

     ! netcdf output
     ijvar=1
     DO jbasin = 1, nbasins
        DO jk = 1, npk
           ierr = putvar (ncout, id_varout(ijvar), REAL(dmoc(jbasin,:,jk)), jk, 1, npjglo, ktime=jt)
        END DO
        ijvar = ijvar + 1
        IF ( ldec ) THEN
!           print *, dmoc_sh(jbasin,60,10)
            DO jk = 1, npk 
              ierr = putvar (ncout, id_varout(ijvar), REAL(dmoc_sh(jbasin,:,jk)), jk, 1, npjglo, ktime=jt)
            END DO
!           print *, dmoc_bt(jbasin,60,10)
            ijvar = ijvar + 1 
            DO jk = 1, npk 
              ierr = putvar (ncout, id_varout(ijvar), REAL(dmoc_bt(jbasin,:,jk)), jk, 1, npjglo, ktime=jt)
            END DO
!           print *, dmoc_ag(jbasin,60,10)
            ijvar = ijvar + 1 
            DO jk = 1, npk 
              ierr = putvar (ncout, id_varout(ijvar), REAL(dmoc_ag(jbasin,:,jk)), jk, 1, npjglo, ktime=jt)
            END DO
            ijvar = ijvar + 1 
         ENDIF
     END DO
  ENDDO  ! time loop

  ierr = closeout(ncout)
CONTAINS
   FUNCTION get_e3v(kk)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION get_e3  ***
    !!
    !! ** Purpose : Send back e3v at level kk selecting
    !!              full step or partial step case
    !!
    !! ** Method  :  check for global flag lfull  
    !!
    !!----------------------------------------------------------------------
    INTEGER(KIND=4), INTENT(in)  :: kk  ! level to work with
    REAL(KIND=4), DIMENSION(npiglo,npjglo) :: get_e3v

    ivmask(:,:) = getvar(cn_fmsk, 'vmask', jk, npiglo, npjglo)
    IF ( lfull ) THEN
        get_e3v(:,:) = e31d(jk)
    ELSE
        get_e3v(:,:) = getvar(cn_fzgr, 'e3v_ps', jk, npiglo, npjglo, ldiom=.TRUE.)
    ENDIF
        get_e3v(:,:) = get_e3v(:,:) * ivmask(:,:)
   
   END FUNCTION get_e3v
 


END PROGRAM cdfmoc
