PROGRAM cdfocape
  !!======================================================================
  !!                     ***  PROGRAM cdfocape   ***
  !!=====================================================================
  !!  ** Purpose : Compute ocean convective potential energy
  !!
  !!  ** Method  : Paul Myers
  !!
  !! History :  4.0  : 06/2019  : P. Verezemskaya : J.M. Molines :  
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!----------------------------------------------------------------------
  USE cdfio
  USE eos
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017
  !! $Id$
  !! Copyright (c) 2012, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE
  
  INTEGER(KIND=4)                           :: jk, jti, it, jt    ! dummy loop index
  INTEGER(KIND=4)                           :: iimin=0, iimax=0   ! domain limitation for computation
  INTEGER(KIND=4)                           :: ijmin=0, ijmax=0   ! domain limitation for computation
  INTEGER(KIND=4)                           :: ierr               ! error status
  INTEGER(KIND=4)                           :: narg, iargc, ijarg ! browse command line
  INTEGER(KIND=4)                           :: npiglo, npjglo     ! size of the domain
  INTEGER(KIND=4)                           :: npiglo_fi, npjglo_fi ! size of the zoomed domain
  INTEGER(KIND=4)                           :: npk, npt           ! size of the domain
  INTEGER(KIND=4)                           :: ncout, ncout1d     ! ncid of output file
  INTEGER(KIND=4), DIMENSION(2)             :: ipk, ipk1d         ! level 
  INTEGER(KIND=4), DIMENSION(2)             :: id_varout          ! varid's
  INTEGER(KIND=4), DIMENSION(2)             :: id_varout1d        ! varid's

  REAL(KIND=4)                              :: grvt=9.8           ! gravity constant (m/s2)
  REAL(KIND=4), DIMENSION(1)                :: zarea              ! area of a domain points
  REAL(KIND=4), DIMENSION(1)                :: zcape_int          ! area of a domain points
  REAL(KIND=4)                              :: zspval             ! missing value
  REAL(KIND=4)                              :: ref_dep            ! reference depth
  REAL(KIND=4), DIMENSION(:), ALLOCATABLE   :: gdep               ! depth array
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e1t                ! e1t - dz for integration
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e2t                ! e2t - dz for integration
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e3t                ! e3t - dz for integration
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: hdept              ! real depth
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zmask              ! mask 
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ztemp              ! temperature
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zsal               ! salinity

  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: zsigh              ! sigma at h level
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: zsigz              ! sigma at z level
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: zsigvi             ! vertically integrated sigma
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: zcape              ! convective energy
  REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: dtim               ! time counter

  CHARACTER(LEN=256)                        :: cf_tfil            ! input filename
  CHARACTER(LEN=256)                        :: cf_sfil            ! salinity file (option)
  CHARACTER(LEN=256)                        :: cf_root            ! readable name
  CHARACTER(LEN=256)                        :: cf_out='ocape.nc'  ! output file name
  CHARACTER(LEN=256)                        :: cf_out_int='intocape.nc' ! integral output file name
  CHARACTER(LEN=256)                        :: cv_out='vocape'         ! 2d ocape name
  CHARACTER(LEN=256)                        :: cv_out_name='vocape'         ! 2d ocape name
  CHARACTER(LEN=256)                        :: cv_out_int='iocape'     ! integral ocape name
  CHARACTER(LEN=256)                        :: cv_out_name_int='iocape'     ! integral ocape name
  CHARACTER(LEN=256)                        :: cldum              ! dummy string
  CHARACTER(LEN=256)                        :: cv_msk='tmask'     ! mask

  TYPE (variable), DIMENSION(2)             :: stypvar            ! structure for attributes
  TYPE (variable), DIMENSION(2)             :: stypvar1d          ! structure for attributes int

  LOGICAL                                   :: lchk  = .FALSE.    ! flag for missing files
  LOGICAL                                   :: lnc4  = .FALSE.    ! flag for missing files
  !!--------------------------------------------------------------------------------------
  CALL ReadCdfNames()
  !!  Read command line and output usage message if not compliant.
  
  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfocape -dep REF-dep -t T-file [-s S-file] [-w imin imax jmin jmax] [-nc4] [-o OUT-file] '  
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'        Compute OCAPE from TS and vertical integration '
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'        -dep REF-dep : reference depth to which you want to ventilate the ocean '
     PRINT *,'        -t T-file : name of the file holding temperature and salinity '
     PRINT *,'      '
     PRINT *,'     OPTIONS :'  
     PRINT *,'        -s S-file : Specify salinity file if not T-file. ' 
     PRINT *,'        -w imin imax jmin jmax : spatial window where mean value'
     PRINT *,'           is computed:'
     PRINT *,'                  if imin = 0 then ALL i are taken'
     PRINT *,'                  if jmin = 0 then ALL j are taken'
     PRINT *,'        -o : Specify output file name instead of ocape.nc ' 
     PRINT *,'        -nc4 : Use netcdf4 output with chunking and deflation level 1 ' 
     PRINT *,'         This option is effective only if cdftools are compiled with' 
     PRINT *,'         a netcdf library supporting chunking and deflation.' 
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'     mesh_hgr.nc, mesh_zgr.nc, mask.nc '
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ocape.nc '           , TRIM(cf_out) 
     PRINT *,'         variables : 2d ocape '           , TRIM(cv_out),' ( Joules/m3 )'
     PRINT *,'       netcdf file : intocape.nc '        , TRIM(cf_out_int) 
     PRINT *,'         variables : weighted mean ocape ', TRIM(cv_out_int),' ( Joules/m3 )'
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'       cdfsigi, cdfmxl '
     PRINT *,'      '
     STOP
  ENDIF

  cf_tfil='none'
  cf_sfil='none'

  ijarg = 1
  cf_sfil='none' 
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-t'   ) ; CALL getarg(ijarg, cf_tfil ) ; ijarg=ijarg+1
     CASE ( '-dep' ) ; CALL getarg(ijarg, cldum )   ; ijarg=ijarg+1 ; READ(cldum,*) ref_dep
        ! option
     CASE ( '-s'   ) ; CALL getarg (ijarg, cf_sfil) ; ijarg=ijarg+1
     CASE ( '-w'   ) ; CALL getarg(ijarg, cldum   ) ; ijarg=ijarg+1 ; READ(cldum,*) iimin
                     ; CALL getarg(ijarg, cldum   ) ; ijarg=ijarg+1 ; READ(cldum,*) iimax
                     ; CALL getarg(ijarg, cldum   ) ; ijarg=ijarg+1 ; READ(cldum,*) ijmin
                     ; CALL getarg(ijarg, cldum   ) ; ijarg=ijarg+1 ; READ(cldum,*) ijmax
     CASE ( '-o'   ) ; CALL getarg(ijarg, cf_root ) ; ijarg=ijarg+1
     CASE ( '-nc4' ) ; lnc4 = .TRUE.
     CASE DEFAULT    ; PRINT *, ' ERROR : ', TRIM(cldum),' : unknown option.'; STOP 1
     END SELECT
  ENDDO

 WRITE(cldum, '(I4.4)') INT(ref_dep)

 ! Define file names according to the one specified without .nc
 cf_out=TRIM(cf_root )//'ocape'//TRIM(cldum)//'.nc'
 cf_out_int=TRIM(cf_root)//'iocape'//TRIM(cldum)//'.nc'

 PRINT *, ref_dep

  IF ( cf_sfil == 'none' ) cf_sfil=cf_tfil
  ! check for missing files - obligatory
  ! Check file existence mask, hgr, zgr:
  lchk = lchk .OR. chkfile( cn_fhgr )
  lchk = lchk .OR. chkfile( cn_fzgr )
  lchk = lchk .OR. chkfile( cn_fmsk )
  lchk = lchk .OR. chkfile( cf_tfil )   ! chkfile is returning F when input file name is 'none'
  lchk = lchk .OR. chkfile( cf_sfil )   ! chkfile is returning F when input file name is 'none'
  IF ( lchk ) STOP 99 ! missing files

  ! Read the domain dimensions - obgigatory
  npiglo = getdim (cf_tfil, cn_x)
  npjglo = getdim (cf_tfil, cn_y)
  npk    = getdim (cf_tfil, cn_z)
  npt    = getdim (cf_tfil, cn_t)

  ! save original npiglo, npiglo - for zoom
  npiglo_fi = npiglo
  npjglo_fi = npjglo
  ! Set new npiglo, npjglo according to zoom
  IF (iimin /= 0 ) THEN ; npiglo = iimax -iimin + 1;  ELSE ; iimin=1 ;
  ENDIF
  IF (ijmin /= 0 ) THEN ; npjglo = ijmax -ijmin + 1;  ELSE ; ijmin=1 ;
  ENDIF

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt
  
  ! Allocate aarays for read, compute and write + time
  ALLOCATE (ztemp(npiglo,npjglo), zsal (npiglo,npjglo))
  ALLOCATE (zsigh(npiglo,npjglo), zsigz (npiglo,npjglo))
  ALLOCATE (zsigvi(npiglo,npjglo), e3t(npiglo,npjglo))
  ALLOCATE ( e1t(npiglo,npjglo), e2t(npiglo,npjglo)) 
  ALLOCATE (zcape(npiglo,npjglo), zmask (npiglo,npjglo))
  ALLOCATE (hdept(npiglo,npjglo))
  ALLOCATE (dtim(npt) , gdep(npk))

  ! Prepare output netcdf files 

  ! The strange thing is that the main part starts here!
  CALL CreateOutput
  
  ! first read missing vaue for salinity
  zspval= getatt(cf_sfil, cn_vosaline, cn_missing_value)
  gdep(:) = getvare3(cn_fzgr, cn_gdept, npk)
  hdept(:,:) = getvar(cn_fzgr, cn_hdept, 1, npiglo, npjglo)
  e1t(:,:) =  getvar(cn_fhgr, cn_ve1t, 1, npiglo, npjglo)
  e2t(:,:) =  getvar(cn_fhgr, cn_ve2t, 1, npiglo, npjglo)
  ! loop over time
  it = 1.
  DO jt=1, npt
     ! loop over depth
     jk = 1
     DO WHILE (gdep(jk) <= ref_dep)
        zmask(:,:) = getvar(cn_fmsk, cv_msk, jk, npiglo, npjglo, kimin=iimin, kjmin=ijmin, ktime=jt)

        ztemp(:,:) = getvar( cf_tfil, cn_votemper,  jk, npiglo, npjglo, kimin=iimin, kjmin=ijmin, ktime=jt )
        zsal( :,:) = getvar( cf_sfil, cn_vosaline,  jk, npiglo, npjglo, kimin=iimin, kjmin=ijmin, ktime=jt )

        e3t(:,:) = getvar(cn_fzgr, cn_ve3t, jk, npiglo, npjglo, kimin=iimin, kjmin=ijmin, ktime=it, ldiom=.NOT.lg_vvl )

        ! we compute sigma using sigmai from eos and multiply it by 0-1 mask
        zsigz(:,:) = sigmai(ztemp, zsal, gdep(jk) , npiglo, npjglo )* zmask(:,:)
        
        ! vertical integration
        zsigvi(:,:) = zsigvi(:,:) + zsigz(:,:)*e3t(:,:)
        jk = jk+1
     END DO
     ! And after this we can read zsigh
     zsigh(:,:) = sigmai(ztemp, zsal, ref_dep, npiglo, npjglo )* zmask(:,:)
     ! and compute 2d convective energy
     zcape(:,:) = (ref_dep*zsigh(:,:) - zsigvi(:,:)) * zmask(:,:)
     ! integrate by area:
     zarea = SUM (e1t(:,:)*e2t(:,:))
     zcape_int = SUM( zcape(:,:)*e1t(:,:)*e2t(:,:) ) * grvt/zarea
     
     zcape = zcape * grvt
 
     ierr = putvar(ncout, id_varout(1), zcape, 1, npiglo, npjglo, ktime=jt)
     ierr = putvar(ncout1d, id_varout1d(1), zcape_int, 1, 1, 1, ktime=jt)
  END DO

  ierr = closeout(ncout)
  ierr = closeout(ncout1d)
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
    ipk(:)= npk  ! all variables (input and output are 3D)
    stypvar(1)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(1)%cname             = cv_out_name ! Where do I get this? 
    stypvar(1)%cunits            = 'J/m3'
    stypvar(1)%rmissing_value    = 0.
    stypvar(1)%valid_min         = 0.001
    stypvar(1)%valid_max         = 5000.
    stypvar(1)%clong_name        = 'Convective energy:refered to '//TRIM(cldum)//' m'
    stypvar(1)%cshort_name       = cv_out 
    stypvar(1)%conline_operation = 'N/A'
    stypvar(1)%caxis             = 'TYX'

    ! create output fileset
    ncout = create      (cf_out, cf_tfil, npiglo, npjglo, 1  , ld_nc4=lnc4 )
    ierr  = createvar   (ncout,  stypvar, 1,   ipk,    id_varout , ld_nc4=lnc4 )
    ierr  = putheadervar(ncout,  cf_tfil, npiglo, npjglo, 1)

    dtim  = getvar1d(cf_tfil, cn_vtimec, npt     )
    ierr  = putvar1d(ncout,   dtim,      npt, 'T')
    
   ! Create 1d output file additionally
   ipk1d(:)=1
   stypvar1d(1)%cname             = cv_out_name_int
   stypvar1d(1)%cunits            = 'J/m3'
   stypvar1d(1)%rmissing_value    = 0.
   stypvar1d(1)%valid_min         = 0.001
   stypvar1d(1)%valid_max         = 5000.
   stypvar1d(1)%clong_name        = 'Weighed mean OCAPE :refered to '//TRIM(cldum)//' m'
   stypvar1d(1)%cshort_name       = cv_out_int
   stypvar1d(1)%conline_operation = 'N/A'
   stypvar1d(1)%caxis             = 'T'

   ! Creating output fileset for int
   ncout1d = create      (cf_out_int, cf_tfil, 1, 1, 1  , ld_nc4=lnc4 )
   ierr    = createvar   (ncout1d,  stypvar1d, 1,   ipk1d,  id_varout1d , ld_nc4=lnc4 )
   ierr    = putheadervar(ncout1d,  cf_tfil, 1, 1, 1)

   dtim    = getvar1d(cf_tfil, cn_vtimec, npt     )
   ierr    = putvar1d(ncout1d,   dtim,      npt, 'T')
  END SUBROUTINE CreateOutput
  
END PROGRAM cdfocape
