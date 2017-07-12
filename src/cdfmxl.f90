PROGRAM cdfmxl
  !!======================================================================
  !!                     ***  PROGRAM  cdfmxl  ***
  !!=====================================================================
  !!  ** Purpose : Compute mixed layer depth
  !!
  !!  ** Method  : - compute surface properties
  !!               - initialize depths and model levels number
  !!               - from bottom to top compute rho and
  !!                 check if rho > rho_surf +rho_c
  !!                 where rho_c is a density criteria given as argument
  !!
  !! History : 2.1  : 10/2005  : J.M. Molines : Original code
  !!           3.0  : 01/2011  : J.M. Molines : Doctor norm + Lic.
  !!           3.0  : 07/2012  : F. Hernandez: Optional S-FILE input
  !!           3.0  : 07/2012  : F. Hernandez: Add new MLD computation for GSOP/GODAE
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
  !! @class mixed_layer
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4), PARAMETER                   :: jp_varout = 7  ! number of output variables 
  INTEGER(KIND=4)                              :: ji, jj, jk, jt ! dummy loop index
  INTEGER(KIND=4)                              :: ik1, ik2, ikt  ! k vertical index of mixed layers 
  INTEGER(KIND=4), DIMENSION(1)                :: nkref10        ! vertical index for 10m depth T layer  
  INTEGER(KIND=4)                              :: narg, iargc    ! browse line
  INTEGER(KIND=4)                              :: ijarg, ixtra   ! browse line
  INTEGER(KIND=4)                              :: npiglo, npjglo ! domain size
  INTEGER(KIND=4)                              :: npk, npt       ! domain size
  INTEGER(KIND=4)                              :: ncout, ierr    ! ncid of output file, error status
  INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: mbathy         ! number of w levels in water <= npk
  INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: nmln1          ! last level where rho > rho + rho_c1
  INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: nmln2          ! last level where rho > rho + rho_c2
  INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: nmln3          ! last level where rho > rho10 + rho_c2 
  INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: nmln4          ! last level where rho > rho10 + rho_c3 
  INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: nmlnt          ! last level where T - SST > temp_c
  INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: nmlnt2         ! last level where T-T10 > temp_c  
  INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: nmlnt3         ! last level where T-T10 > temp_c2 
  INTEGER(KIND=4), DIMENSION(jp_varout)        :: ipk, id_varout ! levels and varid's of output vars

  REAL(KIND=4)                                 :: rmisval=32767. ! Missing value of Mercator fields 
  REAL(KIND=4)                                 :: rr1,rr2        ! Coef for T(z=10m) interp. 
  REAL(KIND=4)                                 :: rho_c1=0.01    ! 1rst density criterion
  REAL(KIND=4)                                 :: rho_c2=0.03    ! 2nd density criterion
  REAL(KIND=4)                                 :: rho_c3=0.125   ! 3rd density criterion 
  REAL(KIND=4)                                 :: temp_c=-0.2    ! temperature criterion
  REAL(KIND=4)                                 :: temp_c2=-0.5   ! 2nd temperature criterion 
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: rtem           ! temperature
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: rtem10         ! 10m depth temperature 
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: rsal           ! salinity
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: rsal10         ! 10m depth salinity 
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: rho            ! density
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: rho10          ! 10m depth density 
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: rho_surf       ! surface density
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: tem_surf       ! surface temperature
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: tmask_surf     ! surface tmask
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: tmask_10       ! 10m-depth tmask 
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: tmask          ! level tmask
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: hmlp1          ! mxl depth based on density criterion 1
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: hmlp2          ! mxl depth based on density criterion 2
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: hmlp3          ! mxl depth based on density criterion 2 and 10m
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: hmlp4          ! mxl depth based on density criterion 3 and 10m
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: hmlt           ! mxl depth based on temperature criterion
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: hmlt2          ! mxl depth based on temperature criterion and 10m
  REAL(KIND=4), DIMENSION(:,:),    ALLOCATABLE :: hmlt3          ! mxl depth based on temperature criterion 2 and 10m
  REAL(KIND=4), DIMENSION(:),      ALLOCATABLE :: gdepw          ! depth of w levels
  REAL(KIND=4), DIMENSION(:),      ALLOCATABLE :: gdept          ! depth of T levels 
  REAL(KIND=4), DIMENSION(1)                   :: rdep           ! dummy depth for output

  REAL(KIND=8), DIMENSION(:),      ALLOCATABLE :: dtim           ! time counter

  CHARACTER(LEN=256)                           :: cf_tfil        ! input T file
  CHARACTER(LEN=256)                           :: cf_sfil        ! input S file (F.Hernandez)
  CHARACTER(LEN=256)                           :: cf_out='mxl.nc'! output file name
  CHARACTER(LEN=256)                           :: cldum          ! dummy character variable

  TYPE(variable), DIMENSION(jp_varout)         :: stypvar        ! structure for attributes 

  LOGICAL                                      :: lexist         ! flag for existence of bathy_level file
  LOGICAL                                      :: lnc4=.FALSE.   ! flag for netcdf4 output with chunking and deflation
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfmxl -t T-file [-s S-file] [-o OUT-file] [-nc4]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute 7 estimates of the mixed layer depth from temperature and '
     PRINT *,'       salinity given in the input file, based on 7 different criteria:'
     PRINT *,'       1- Density criterion (0.01 kg/m3 difference between surface and MLD)' 
     PRINT *,'       2- Density criterion (0.03 kg/m3 difference between surface and MLD)' 
     PRINT *,'       3- Temperature criterion (0.2 C absolute difference between surface '
     PRINT *,'          and MLD)'
     PRINT *,'       4- Temperature criterion (0.2 C absolute difference between T at 10m '
     PRINT *,'          and MLD)'
     PRINT *,'       5- Temperature criterion (0.5 C absolute difference between T at 10m '
     PRINT *,'          and MLD)'
     PRINT *,'       6- Density criterion (0.03 kg/m3 difference between rho at 10m and MLD) '
     PRINT *,'       7- Density criterion (0.125 kg/m3 difference between rho at 10m and MLD) '
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -t T-file   : input netcdf file (gridT)' 
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-s S-file] : input netcdf file (gridS), if vosaline not in T-file' 
     PRINT *,'       [-o OUT-file] : specify the name of output file instead of ', TRIM(cf_out)
     PRINT *,'       [-nc4] : use netcdf4 chunking and deflation on output '
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        ',TRIM(cn_fzgr)
     PRINT *,'         In case of FULL STEP configuration, ',TRIM(cn_fbathylev),' is also required.'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : somxl010    = mld on density criterion 0.01 ref. surf.'
     PRINT *,'                     somxl030    = mld on density criterion 0.03 ref. surf.'
     PRINT *,'                     somxlt02    = mld on temperature criterion -0.2 ref. surf.'
     PRINT *,'                     somxlt02z10 = mld on temperature criterion -0.2 ref. 10m'
     PRINT *,'                     somxlt05z10 = mld on temperature criterion -0.5 ref. 10m'
     PRINT *,'                     somxl030z10 = mld on density criterion 0.03 ref. 10m'
     PRINT *,'                     somxl125z10 = mld on density criterion 0.125 ref. 10m'
     PRINT *,'      '
     STOP 
  ENDIF

  ijarg = 1  ; cf_sfil='none'

  DO WHILE ( ijarg <= narg )
     CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1
     SELECT CASE (cldum )
     CASE ( '-t'   ) ; CALL getarg (ijarg,cf_tfil ) ; ijarg = ijarg + 1
     CASE ( '-s'   ) ; CALL getarg (ijarg,cf_sfil ) ; ijarg = ijarg + 1
        ! options
     CASE ( '-nc4' ) ; lnc4 = .TRUE.
     CASE ( '-o'   ) ; CALL getarg (ijarg,cf_out  ) ; ijarg = ijarg + 1
     CASE DEFAULT    ; PRINT *,' ERROR : ',TRIM(cldum) ,' :  unknown option.' ; STOP 99
     END SELECT
  ENDDO

  IF ( cf_sfil == 'none' ) cf_sfil=cf_tfil

  IF ( chkfile(cf_tfil) .OR. chkfile(cn_fzgr) .OR. chkfile(cf_sfil)  ) STOP 99 ! missing file

  ! read dimensions 
  npiglo = getdim (cf_tfil,cn_x)
  npjglo = getdim (cf_tfil,cn_y)
  npk    = getdim (cf_tfil,cn_z)
  npt    = getdim (cf_tfil,cn_t)

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  ALLOCATE (rtem  (npiglo,npjglo), rsal  (npiglo,npjglo), rho  (npiglo,npjglo)    )
  ALLOCATE (rtem10(npiglo,npjglo), rsal10(npiglo,npjglo), rho10(npiglo,npjglo)    )

  ALLOCATE (hmlp1(npiglo,npjglo), hmlp2(npiglo,npjglo), hmlt(npiglo,npjglo)       )
  ALLOCATE (hmlp3(npiglo,npjglo), hmlp4(npiglo,npjglo)                            ) 
  ALLOCATE (hmlt2(npiglo,npjglo), hmlt3(npiglo,npjglo)                            )

  ALLOCATE (nmln1 (npiglo,npjglo), nmln2 (npiglo,npjglo), nmlnt(npiglo,npjglo)    )
  ALLOCATE (nmln3 (npiglo,npjglo), nmln4 (npiglo,npjglo)                          )
  ALLOCATE (nmlnt2(npiglo,npjglo), nmlnt3(npiglo,npjglo)                          )  

  ALLOCATE (tmask(npiglo,npjglo), tmask_surf(npiglo,npjglo), tmask_10(npiglo,npjglo))
  ALLOCATE (rho_surf(npiglo,npjglo), tem_surf(npiglo,npjglo)                )
  ALLOCATE (mbathy(npiglo,npjglo)                                           )
  ALLOCATE (gdepw(0:npk), gdept(npk), dtim(npt)                              )

  ! read mbathy and gdepw use real rtem(:,:) as template (getvar is used for real only)
  IF ( chkfile( cn_fbathylev,ld_verbose=.FALSE.)  ) THEN
     PRINT *, 'Read mbathy in ', TRIM(cn_fzgr),' ...'
     rtem(:,:) = getvar(cn_fzgr,      'mbathy',    1, npiglo, npjglo)
  ELSE
     PRINT *, 'Read mbathy in ', TRIM(cn_bathylev),' ...'
     rtem(:,:) = getvar(cn_fbathylev, cn_bathylev, 1, npiglo, npjglo)
  ENDIF

  mbathy(:,:)  = rtem(:,:)
  gdepw(0)     = 99999. ! dummy value, always masked -but eventually accessed on land-
  gdepw(1:npk) = getvare3(cn_fzgr, cn_gdepw, npk)
  gdept(:)     = getvare3(cn_fzgr, cn_gdept, npk)

  ! find the T-reference level for 10m (F.Hernandez)
  nkref10 = MINLOC(gdept,gdept>=10.) - 1 ;  IF ( nkref10(1) < 1 ) nkref10(1)=1

  ! coef for linear interpolation of T at 10m between nkref10 and nkref10+1
  rr1 = (10. - gdept(nkref10(1)+1) ) / (gdept(nkref10(1))-gdept(nkref10(1)+1))
  rr2 = (gdept(nkref10(1)) - 10.   ) / (gdept(nkref10(1))-gdept(nkref10(1)+1))

  ! find W levels for later computation
  nkref10 = MINLOC(gdepw(1:npk),gdepw(1:npk)>=10)-1 ;  IF ( nkref10(1) < 1 ) nkref10(1)=1

  CALL CreateOutput

  DO jt=1,npt

     ! read T/S levels around 10m and interpolate 
     rtem  (:,:) = getvar(cf_tfil, cn_votemper, nkref10(1),   npiglo, npjglo, ktime=jt )
     rtem10(:,:) = getvar(cf_tfil, cn_votemper, nkref10(1)+1, npiglo, npjglo, ktime=jt )
     WHERE ( rtem == rmisval ) rtem10 = rmisval
     WHERE ( .NOT. (rtem10 == rmisval) ) rtem10 = rtem*rr1 + rtem10*rr2

     rsal  (:,:) = getvar(cf_sfil, cn_vosaline, nkref10(1),   npiglo, npjglo, ktime=jt )
     rsal10(:,:) = getvar(cf_sfil, cn_vosaline, nkref10(1)+1, npiglo, npjglo, ktime=jt )
     WHERE ( rsal == rmisval ) rsal10 = rmisval
     WHERE ( .NOT. (rsal10 == rmisval) ) rsal10 = rsal*rr1 + rsal10*rr2

     ! read surface T/S
     rtem(:,:) = getvar(cf_tfil, cn_votemper, 1, npiglo, npjglo, ktime=jt )
     rsal(:,:) = getvar(cf_sfil, cn_vosaline, 1, npiglo, npjglo, ktime=jt )

     ! .. and deduce land-mask from salinity 
     ! ... modified to take into account fill_value = 32767 F.Hernandez 
     IF (jt == 1 ) THEN
        ! For surface criteria
        tmask(:,:) = 1.
        WHERE ( rsal == 0. .OR. rsal == rmisval .OR. rtem == rmisval ) tmask = 0.
        tmask_surf(:,:) = tmask(:,:)

        ! For 10m depth criteria (F. Hernandez)
        tmask(:,:) = 1.
        WHERE ( rsal10 == 0. .OR. rsal10 == rmisval .OR. rtem10 == rmisval) tmask = 0.
        tmask_10(:,:) = tmask(:,:)
     ENDIF

     ! compute rho_surf
     rho_surf(:,:) = sigma0 (rtem, rsal, npiglo, npjglo )* tmask_surf(:,:)
     tem_surf(:,:) = rtem(:,:)

     ! compute rho at 10m-depth
     rho10(:,:) = sigma0 (rtem10, rsal10, npiglo, npjglo )* tmask_10(:,:)

     ! Initialization to the number of w ocean point mbathy
     nmln1(:,:) = mbathy(:,:)
     nmln2(:,:) = mbathy(:,:)
     nmln3(:,:) = mbathy(:,:)  
     nmln4(:,:) = mbathy(:,:)  
     nmlnt(:,:) = mbathy(:,:)
     nmlnt2(:,:) = mbathy(:,:) 
     nmlnt3(:,:) = mbathy(:,:) 

     ! compute mixed layer depth
     ! Last w-level at which rhop>=rho surf+rho_c (starting from jpk-1)
     ! (rhop defined at t-point, thus jk-1 for w-level just above)
     DO jk = npk-1, 2, -1
        rtem (:,:) = getvar(cf_tfil, cn_votemper, jk ,npiglo, npjglo, ktime=jt)
        rsal (:,:) = getvar(cf_sfil, cn_vosaline, jk ,npiglo, npjglo, ktime=jt)
        tmask(:,:) = 1. ! take into account missing values 32767 
        WHERE ( rsal == 0. .OR. rsal >= rmisval .OR. rtem == rmisval ) tmask = 0.
        rho  (:,:) = sigma0 (rtem, rsal, npiglo, npjglo )* tmask(:,:)

        DO jj = 1, npjglo
           DO ji = 1, npiglo
              IF( rho(ji,jj)  > rho_surf(ji,jj) + rho_c1 )   nmln1(ji,jj) = jk
              IF( rho(ji,jj)  > rho_surf(ji,jj) + rho_c2 )   nmln2(ji,jj) = jk
              IF( ABS(rtem(ji,jj) - tem_surf(ji,jj)) > ABS( temp_c)  )   nmlnt(ji,jj) = jk
           END DO
        END DO

        ! Compute with the 10m depth reference: stop if level < nkref10+1 (F.Hernandez)
        IF ( jk > nkref10(1) ) THEN
           DO jj = 1, npjglo
              DO ji = 1, npiglo
                 IF( rho(ji,jj)  > rho10(ji,jj) + rho_c2 )   nmln3(ji,jj) = jk
                 IF( rho(ji,jj)  > rho10(ji,jj) + rho_c3 )   nmln4(ji,jj) = jk
                 IF( ABS(rtem(ji,jj) - rtem10(ji,jj)) > ABS( temp_c)  )   nmlnt2(ji,jj) = jk
                 IF( ABS(rtem(ji,jj) - rtem10(ji,jj)) > ABS( temp_c2)  )   nmlnt3(ji,jj) = jk
              END DO
           END DO
        ENDIF

     END DO

     ! Mixed layer depth
     DO jj = 1, npjglo
        DO ji = 1, npiglo
           ik1 = nmln1(ji,jj) ; ik2 = nmln2(ji,jj) ; ikt = nmlnt(ji,jj)
           hmlp1 (ji,jj) = gdepw(ik1) * tmask_surf(ji,jj)
           hmlp2 (ji,jj) = gdepw(ik2) * tmask_surf(ji,jj)
           hmlp3 (ji,jj) = gdepw(nmln3(ji,jj)) * tmask_10(ji,jj)  
           hmlp4 (ji,jj) = gdepw(nmln4(ji,jj)) * tmask_10(ji,jj)            
           hmlt (ji,jj)  = gdepw(ikt) * tmask_surf(ji,jj)
           hmlt2 (ji,jj) = gdepw(nmlnt2(ji,jj)) * tmask_10(ji,jj) 
           hmlt3 (ji,jj) = gdepw(nmlnt3(ji,jj)) * tmask_10(ji,jj) 
        END DO
     END DO

     ! Correct for missing values = 32767 
     WHERE ( tmask_surf == 0. )
        hmlp1 = rmisval ; hmlp2 = rmisval ; hmlt = rmisval
     END WHERE
     WHERE ( tmask_10 == 0. )
        hmlp3 = rmisval ; hmlp4 = rmisval ; hmlt2 = rmisval ; hmlt3 = rmisval
     END WHERE

     ierr = putvar(ncout, id_varout(1), hmlp1, 1, npiglo, npjglo, ktime=jt)
     ierr = putvar(ncout, id_varout(2), hmlp2, 1, npiglo, npjglo, ktime=jt)
     ierr = putvar(ncout, id_varout(3), hmlt , 1, npiglo, npjglo, ktime=jt)
     ierr = putvar(ncout, id_varout(4), hmlt2, 1, npiglo, npjglo, ktime=jt) 
     ierr = putvar(ncout, id_varout(5), hmlt3, 1, npiglo, npjglo, ktime=jt) 
     ierr = putvar(ncout, id_varout(6), hmlp3, 1, npiglo, npjglo, ktime=jt) 
     ierr = putvar(ncout, id_varout(7), hmlp4, 1, npiglo, npjglo, ktime=jt) 

  END DO ! time loop

  ierr = closeout(ncout)

CONTAINS

  SUBROUTINE CreateOutput
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutput  ***
    !!
    !! ** Purpose :  Create netcdf output file(s) 
    !!
    !! ** Method  :  Use stypvar global description of variables
    !!
    !!----------------------------------------------------------------------
    rdep(1) = 0.

    DO ji = 1, jp_varout
       stypvar(ji)%ichunk = (/npiglo, MAX(1,npjglo/30), 1, 1 /)
    ENDDO

    ipk(:)                    = 1
    stypvar(1)%cname          = 'somxl010'
    stypvar(2)%cname          = 'somxl030'
    stypvar(3)%cname          = 'somxlt02'
    stypvar(4)%cname          = 'somxlt02z10' 
    stypvar(5)%cname          = 'somxlt05z10'
    stypvar(6)%cname          = 'somxl030z10'
    stypvar(7)%cname          = 'somxl125z10'
    stypvar%cunits            = 'm'
    stypvar%rmissing_value    = rmisval  ! to be compliant with Mercator standards
    stypvar%valid_min         = 0.
    stypvar%valid_max         = 7000.
    stypvar(1)%clong_name     = 'Mixed_Layer_Depth_on_0.01_rho_crit'
    stypvar(2)%clong_name     = 'Mixed_Layer_Depth_on_0.03_rho_crit'
    stypvar(3)%clong_name     = 'Mixed_Layer_Depth_on_-0.2_temp_crit'
    stypvar(4)%clong_name     = 'Mixed_Layer_Depth_on_-0.2_temp_crit ref. 10m'
    stypvar(5)%clong_name     = 'Mixed_Layer_Depth_on_-0.5_temp_crit ref. 10m'
    stypvar(6)%clong_name     = 'Mixed_Layer_Depth_on_0.03_rho_crit ref. 10m'
    stypvar(7)%clong_name     = 'Mixed_Layer_Depth_on_0.125_rho_crit ref. 10m'
    stypvar(1)%cshort_name    = 'somxl010'
    stypvar(2)%cshort_name    = 'somxl030'
    stypvar(3)%cshort_name    = 'somxlt02'
    stypvar(4)%cshort_name    = 'ILD02z10'
    stypvar(5)%cshort_name    = 'ILD05z10'
    stypvar(6)%cshort_name    = 'MLD030z10'
    stypvar(7)%cshort_name    = 'MLD125z10'
    stypvar%conline_operation = 'N/A'
    stypvar%caxis             = 'TYX'

    ncout = create      (cf_out, cf_tfil, npiglo, npjglo, 1             , ld_nc4=lnc4 )
    ierr  = createvar   (ncout,  stypvar, jp_varout,      ipk, id_varout, ld_nc4=lnc4 )
    ierr  = putheadervar(ncout,  cf_tfil, npiglo, npjglo, 1, pdep=rdep)

    dtim = getvar1d(cf_tfil, cn_vtimec, npt     )
    ierr = putvar1d(ncout,  dtim,       npt, 'T')

  END SUBROUTINE CreateOutput

END PROGRAM cdfmxl
