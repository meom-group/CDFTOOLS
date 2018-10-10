PROGRAM cdfmxlsaltc
  !!======================================================================
  !!                     ***  PROGRAM  cdfmxlsaltc  ***
  !!=====================================================================
  !!  ** Purpose : Compute the salt content in the mixed layer. Work for
  !!               partial steps (default) or full step (-full option)
  !!
  !!  ** Method  : compute the sum ( rho S  * e1 *e2 * e3 *mask )
  !!               for the mixed layer stored into gridT file
  !!
  !! History : 2.1  : 04/2006  : J.M. Molines : Original code
  !!           3.0  : 01/2011  : J.M. Molines : Doctor norm + Lic.
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  USE modutils
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class mixed_layer
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                               :: jk, jt              ! dummy loop index
  INTEGER(KIND=4)                               :: it                  ! time index for vvl
  INTEGER(KIND=4)                               :: narg, iargc, ijarg  ! command line 
  INTEGER(KIND=4)                               :: npiglo, npjglo      ! size of the domain
  INTEGER(KIND=4)                               :: npk, npt            ! size of the domain,
  INTEGER(KIND=4)                               :: ncout, ierr         ! ncid and error status
  INTEGER(KIND=4), DIMENSION(1)                 :: ipk, id_varout      ! levels and varid's of output vars

  REAL(KIND=4), PARAMETER                       :: rprho0=1020.        ! rho reference density
  REAL(KIND=4), PARAMETER                       :: rpcp=4000.          ! calorific capacity
  REAL(KIND=4), DIMENSION(1)                    :: rdep                ! dummy depth output
  REAL(KIND=4), DIMENSION(:),       ALLOCATABLE :: gdepw               ! vertical levels
  REAL(KIND=4), DIMENSION(:),       ALLOCATABLE :: e31d                ! vertical metric full
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: e3                  ! metrics
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: zs                  ! temperature in the MXL
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: zmxl                ! depth of the MXL
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: zmask               ! mask

  REAL(KIND=8)                                  :: dvol                ! total volume
  REAL(KIND=8), DIMENSION(:),       ALLOCATABLE :: dtim                ! time counter
  REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: dmxlsaltc           ! heat content

  CHARACTER(LEN=256)                            :: cf_tfil             ! input file name
  CHARACTER(LEN=256)                            :: cf_sfil             ! input file name with salinity (option)
  CHARACTER(LEN=256)                            :: cf_mfil             ! input file name with mld (option)
  CHARACTER(LEN=256)                            :: cf_out='mxlsaltc.nc'! output file
  CHARACTER(LEN=256)                            :: cv_out='somxlsaltc' ! input file name
  CHARACTER(LEN=256)                            :: cglobal             ! global attribute
  CHARACTER(LEN=256)                            :: cldum               ! dummy string

  TYPE(variable), DIMENSION(1)                  :: stypvar             ! stucture for attributes (output)

  LOGICAL                                       :: lfull=.FALSE.       ! full step flag
  LOGICAL                                       :: lnc4 =.FALSE.       ! netcdf4  flag
  LOGICAL                                       :: lchk                ! file existence flag (true if missing)
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfmxlsaltc -t T-file [-s S-file] [-mxlf MXL-file] [-full] [-vvl] ...'
     PRINT *,'                 ...   [-o OUT-file] [-nc4] [-vvl]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the salt content in the mixed layer.' 
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -t T-file : netcdf input file with salinity and mld (gridT).' 
     PRINT *,'         If salinity not in T-file use -s option.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-s S-file     ]: specify salinity file if not T-file.'
     PRINT *,'       [-mxlf MXL-file]: netcdf input file with mld if not in T-file.' 
     PRINT *,'       [-full ] : for full step configurations, default is partial step.' 
     PRINT *,'       [-o OUT-file ] : specify output file instead of ',TRIM(cf_out) 
     PRINT *,'       [-nc4 ] : Use netcdf4 output with chunking and deflation level 1.'
     PRINT *,'                 This option is effective only if cdftools are compiled with'
     PRINT *,'                 a netcdf library supporting chunking and deflation.'
     PRINT *,'       [-vvl ] : use time-varying vertical metrics.'
     PRINT *,'      '
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       ', TRIM(cn_fzgr),' and ', TRIM(cn_fmsk) 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out),' unless option -o is used.'
     PRINT *,'         variables : ', TRIM(cv_out),' (Kg/m2)'
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'       cdfmxl, cdfmxlhcsc, cdfmxlheatc ' 
     PRINT *,'      '
     STOP 
  ENDIF

  cf_mfil='none'
  cf_sfil='none'
  ijarg = 1
  DO WHILE ( ijarg <= narg )
    CALL getarg (ijarg, cldum   ) ; ijarg = ijarg + 1 
    SELECT CASE ( cldum )
    CASE ( '-t','-f'  ) ; CALL getarg (ijarg, cf_tfil  ) ; ijarg = ijarg + 1
    ! options
    CASE ( '-s'       ) ; CALL getarg (ijarg, cf_sfil  ) ; ijarg = ijarg + 1
    CASE ( '-mxlf'    ) ; CALL getarg (ijarg, cf_mfil  ) ; ijarg = ijarg + 1
    CASE ( '-full'    ) ; lfull  = .TRUE.
    CASE ( '-o'       ) ; CALL getarg (ijarg, cf_out ) ; ijarg = ijarg + 1
    CASE ( '-nc4'     ) ; lnc4   = .TRUE.
    CASE ( '-vvl'     ) ; lg_vvl = .TRUE.
    CASE DEFAULT  ; PRINT *, 'ERROR: ', TRIM(cldum),' : unknown option.' ; STOP 99
    END SELECT
  END DO

  IF ( cf_mfil == 'none' ) cf_mfil = cf_tfil 
  IF ( cf_sfil == 'none' ) cf_sfil = cf_tfil 

  lchk = chkfile (cn_fzgr  )
  lchk = chkfile (cn_fmsk  ) .OR. lchk
  lchk = chkfile (cf_tfil  ) .OR. lchk
  lchk = chkfile (cf_sfil  ) .OR. lchk
  lchk = chkfile (cf_mfil  ) .OR. lchk
  IF ( lchk ) STOP 99 ! missing files
  IF ( lg_vvl ) THEN 
     cn_fe3t = cf_tfil
     cn_ve3t = cn_ve3tvvl
  ENDIF

  CALL SetGlobalAtt( cglobal )

  npiglo = getdim (cf_tfil,cn_x)
  npjglo = getdim (cf_tfil,cn_y)
  npk    = getdim (cf_tfil,cn_z)
  npt    = getdim (cf_tfil,cn_t)


  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  ! Allocate arrays
  ALLOCATE ( zmask(npiglo,npjglo), dmxlsaltc(npiglo, npjglo) )
  ALLOCATE ( zs(npiglo,npjglo), zmxl(npiglo,npjglo)          )
  ALLOCATE ( e3(npiglo,npjglo)                               )
  ALLOCATE ( gdepw(npk), dtim(npt)                           )
  
  CALL CreateOutput

  IF ( lfull ) ALLOCATE ( e31d(npk)                          )
               gdepw(:) = getvare3(cn_fzgr, cn_gdepw, npk)
  IF ( lfull ) e31d( :) = getvare3(cn_fzgr, cn_ve3t1d,  npk)

  DO jt=1,npt
     IF ( lg_vvl ) THEN ; it=jt
     ELSE               ; it=1
     ENDIF
     dmxlsaltc(:,:) = 0.d0
     dvol           = 0.d0
     zmxl( :,:) = getvar(cf_mfil, cn_somxl010, 1,  npiglo, npjglo, ktime=jt)

     DO jk = 1, npk
        zs(   :,:) = getvar(cf_sfil, cn_vosaline, jk, npiglo, npjglo, ktime=jt)
        zmask(:,:) = getvar(cn_fmsk, cn_tmask,    jk, npiglo, npjglo          )

        ! get e3 at level jk ( ps...)
        IF ( lfull ) THEN
           e3(:,:) = e31d(jk)
        ELSE
           e3(:,:) = getvar(cn_fe3t, cn_ve3t, jk, npiglo, npjglo, ktime=it, ldiom=.NOT.lg_vvl )
        ENDIF

        !  e3 is used as a flag for the mixed layer; It is 0 outside the mixed layer
        e3(:,:)=MAX ( 0., MIN(e3, zmxl-gdepw(jk) ) )
        WHERE ( e3 == 0 ) zmask = 0.

        dvol      = SUM( DBLE(e3 * zmask) )
        dmxlsaltc = zs * e3 * zmask * 1.d0 + dmxlsaltc

        IF (dvol /= 0 )THEN
           !   go on !
        ELSE
           !   no more layer below !
           EXIT   ! get out of the jk loop
        ENDIF
     END DO

     ! Output to netcdf file : Kg/m2
     dmxlsaltc = rprho0*dmxlsaltc
     ierr = putvar(ncout, id_varout(1), REAL(dmxlsaltc), 1, npiglo, npjglo, ktime=jt)
  END DO

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
  rdep(1)                      = 0.

  ipk(:)                       = 1
  stypvar(1)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
  stypvar(1)%cname             = cv_out
  stypvar(1)%cunits            = 'kg/m2'
  stypvar(1)%rmissing_value    = 0.
  stypvar(1)%valid_min         = 0
  stypvar(1)%valid_max         =  1.e9
  stypvar(1)%clong_name        = 'Mixed_Layer_Salt_Content'
  stypvar(1)%cshort_name       = cv_out
  stypvar(1)%conline_operation = 'N/A'
  stypvar(1)%caxis             = 'TYX'

  ! Initialize output file
  ncout = create      (cf_out, cf_tfil, npiglo, npjglo, 1                          ,ld_nc4=lnc4 )
  ierr  = createvar   (ncout,  stypvar, 1,      ipk,    id_varout, cdglobal=cglobal,ld_nc4=lnc4 )
  ierr  = putheadervar(ncout,  cf_tfil, npiglo, npjglo, 1, pdep=rdep                )

  dtim = getvar1d(cf_tfil, cn_vtimec, npt     )
  ierr = putvar1d(ncout,  dtim,       npt, 'T')

  END SUBROUTINE CreateOutput

END PROGRAM cdfmxlsaltc
