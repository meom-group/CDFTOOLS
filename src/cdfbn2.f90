PROGRAM cdfbn2
  !!======================================================================
  !!                     ***  PROGRAM  cdfbn2  ***
  !!=====================================================================
  !!  ** Purpose : Compute the Brunt Vaissala frequency
  !!               using same algoritm than NEMO
  !!
  !!  ** Method  : Try to avoid 3 d arrays : work with 2 levels a a time
  !!              The brunt-vaisala frequency is computed using the
  !!              polynomial expression of McDougall (1987):
  !!              N^2 = grav * beta * ( alpha/beta*dk[ t ] - dk[ s ] )/e3w
  !!              N2 is then insterpolated at T levels
  !!
  !! History : 2.0  : 11/2004  : J.M. Molines : Original code
  !!           2.1  : 04/2005  : J.M. Molines : use cdfio
  !!           3.0  : 12/2010  : J.M. Molines : Doctor norm + Lic.
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames   ! for cdf variable names
  USE eos
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class Equation_of_state
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                              :: jk, jt                   ! dummy loop index
  INTEGER(KIND=4)                              :: it                       ! time index for vvl
  INTEGER(KIND=4)                              :: ierr                     ! working integer
  INTEGER(KIND=4)                              :: narg, iargc, ijarg       ! 
  INTEGER(KIND=4)                              :: npiglo, npjglo, npk, npt ! size of the domain
  INTEGER(KIND=4)                              :: iup = 1, idown = 2, itmp ! for swapping the levels
  INTEGER(KIND=4)                              :: ncout                    ! ncid of output file
  INTEGER(KIND=4), DIMENSION(1)                :: ipk, id_varout           ! level and id of output variables

  REAL(KIND=4)                                 :: zsps                     ! Missing value for salinity
  REAL(KIND=4), DIMENSION (:,:,:), ALLOCATABLE :: ztemp, zsal, zwk         ! Array to read 2 layer of data
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE :: zn2                      ! Brunt Vaissala Frequency (N2)
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE :: zmask, e3w               ! mask and metric
  REAL(KIND=4), DIMENSION (:),     ALLOCATABLE :: gdep,  e3w1d             ! depth and time

  REAL(KIND=8), DIMENSION (:),     ALLOCATABLE :: dtim                     ! depth and time

  CHARACTER(LEN=256)                           :: cf_tfil, cldum, cv_dep   ! input file name, ...
  CHARACTER(LEN=256)                           :: cf_sfil                  !
  CHARACTER(LEN=256)                           :: cf_out = 'bn2.nc'        ! output file name
  CHARACTER(LEN=256)                           :: cglobal                  ! global attribute
  CHARACTER(LEN=80)                            :: cf_e3w                   ! file with e3w in case of vvl
  CHARACTER(LEN=80)                            :: cv_bn2  = 'vobn2'        ! cdf variable name for N2

  TYPE(variable), DIMENSION(1)                 :: stypvar                  ! variable attribute

  LOGICAL                                      :: l_w   =.FALSE.           ! flag for vertical location of bn2
  LOGICAL                                      :: lchk                     ! check missing files
  LOGICAL                                      :: lfull =.FALSE.           ! full step flag
  LOGICAL                                      :: lnc4  =.FALSE.           ! full step flag
  LOGICAL                                      :: ll_teos10  = .FALSE.     ! teos10 flag

  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  cf_sfil='none'
  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfbn2  -t T-file [-s S-file] [-W] [-full] [-o OUT-file]...'
     PRINT *,'               ... [-nc4] [-vvl W-file] [-teos10] '
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the Brunt-Vaissala frequency (N2) according to temperature and' 
     PRINT *,'       salinity given in the input file.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -t T-file : netcdf input gridT file for temperature and salinity.' 
     PRINT *,'             If salinity is not in T-file, use -s option...' 
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-s S-file ] : specify salinity file name if not T-file' 
     PRINT *,'       [-W ] : keep N2 at W points. Default is to interpolate N2 at T point on' 
     PRINT *,'             the vertical.'
     PRINT *,'       [-full ] : indicate a full step configuration instead of the default'
     PRINT *,'             partial steps.'
     PRINT *,'       [-o OUT-file ] : specify output file name instead of ',TRIM(cf_out),'.'
     PRINT *,'       [-nc4 ]  : Use netcdf4 output with chunking and deflation level 1.'
     PRINT *,'                 This option is effective only if cdftools are compiled with'
     PRINT *,'                 a netcdf library supporting chunking and deflation.'
     PRINT *,'       [-vvl W-file ] : use time-varying vertical metrics, W-file is a file'
     PRINT *,'                 holding e3w(t) for vvl.'
     PRINT *,'       [-teos10] : use TEOS10 equation of state instead of default EOS80'
     PRINT *,'                 Temperature should be conservative temperature (CT) in deg C.'
     PRINT *,'                 Salinity should be absolute salinity (SA) in g/kg.'
     PRINT *,'      '
     PRINT *,'     OPENMP SUPPORT : yes'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       ',TRIM(cn_fzgr),' is needed for this program.' 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) ,' unless -o option specified'
     PRINT *,'         variables : ', TRIM(cv_bn2)
     PRINT *,'      '
     PRINT *,'    SEE ALSO :'
     PRINT *,'       cdfsig0, cdfsigi, cdfsiginsitu, cdfsigntr '
     STOP 
  ENDIF

  cglobal = 'Partial step computation'

  ijarg = 1
  DO WHILE ( ijarg <= narg ) 
     CALL getarg(ijarg, cldum) ; ijarg = ijarg + 1
     SELECT CASE (cldum)
     CASE ( '-t'      ) ; CALL getarg(ijarg, cf_tfil) ; ijarg = ijarg + 1
        ! options
     CASE ( '-s'      ) ; CALL getarg(ijarg, cf_sfil) ; ijarg = ijarg + 1
     CASE ( '-W'      ) ; l_w       = .TRUE.
     CASE ( '-full'   ) ; lfull     = .TRUE. ; cglobal = 'full step computation'
     CASE ( '-o'      ) ; CALL getarg(ijarg, cf_out ) ; ijarg = ijarg + 1
     CASE ( '-nc4'    ) ; lnc4      = .TRUE.
     CASE ( '-vvl'    ) ; lg_vvl    = .TRUE. 
        ; CALL getarg(ijarg, cf_e3w ) ; ijarg = ijarg + 1
     CASE ( '-teos10' ) ; ll_teos10 = .TRUE. 
     CASE DEFAULT     ; PRINT *,' ERROR : ', TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  END DO

  CALL eos_init(ll_teos10)

  IF ( cf_sfil == 'none' ) cf_sfil=cf_tfil

  lchk = chkfile (cn_fzgr )
  lchk = lchk .OR. chkfile (cf_tfil  )
  lchk = lchk .OR. chkfile (cf_sfil  )
  IF (lg_vvl ) lchk = lchk .OR. chkfile (cf_e3w) 
  IF ( lchk  ) STOP 99  ! missing files 
  ! Look for missing value for salinity
  zsps = getspval(cf_sfil, cn_vosaline)

  IF ( lg_vvl ) THEN
    cn_fe3w = cf_e3w
    cn_ve3w = cn_ve3wvvl
  ENDIF

  npiglo = getdim (cf_tfil, cn_x)
  npjglo = getdim (cf_tfil, cn_y)
  npk    = getdim (cf_tfil, cn_z)
  npt    = getdim (cf_tfil, cn_t)


  PRINT *, 'npiglo =', npiglo
  PRINT *, 'npjglo =', npjglo
  PRINT *, 'npk    =', npk
  PRINT *, 'npt    =', npt

  ALLOCATE (ztemp(npiglo,npjglo,2), zsal(npiglo,npjglo,2) )
  ALLOCATE (zwk(npiglo,npjglo,2), zmask(npiglo,npjglo)    )
  ALLOCATE (zn2(npiglo,npjglo), e3w(npiglo,npjglo)        )
  ALLOCATE (gdep(npk), dtim(npt)                          )
  IF ( lfull ) ALLOCATE (e3w1d(npk) )

  cv_dep=cn_gdept
  IF (l_w) cv_dep=cn_gdepw

  gdep(:) = getvare3(cn_fzgr, cv_dep, npk)

  CALL CreateOutput

  IF ( lfull )  e3w1d(:) = getvare3(cn_fzgr, cn_ve3w1d, npk)

  gdep(:) = getvare3(cn_fzgr, cn_gdepw, npk)
  DO jt=1,npt
     IF ( lg_vvl ) THEN ; it=jt
     ELSE               ; it=1
     ENDIF
     !  2 levels of T and S are required : iup,idown (with respect to W level)
     !  Compute from bottom to top (for vertical integration)
     ztemp(:,:,idown) = getvar(cf_tfil, cn_votemper,  npk-1, npiglo, npjglo, ktime=jt)
     zsal( :,:,idown) = getvar(cf_sfil, cn_vosaline,  npk-1, npiglo, npjglo, ktime=jt)
     zwk(:,:,idown)   = 0.0

     DO jk = npk-1, 2, -1 
        PRINT *,'level ',jk
        zmask(:,:)=1.
        ztemp(:,:,iup)= getvar(cf_tfil, cn_votemper, jk-1, npiglo, npjglo, ktime=jt)
        zsal(:,:,iup) = getvar(cf_sfil, cn_vosaline, jk-1, npiglo, npjglo, ktime=jt)
        WHERE(zsal(:,:,idown) == zsps ) zmask = 0

        IF ( lfull ) THEN ; e3w(:,:) = e3w1d(jk)
        ELSE              ; e3w(:,:) = getvar(cn_fe3w, cn_ve3w , jk, npiglo, npjglo, ktime=it, ldiom=.NOT.lg_vvl )
        ENDIF

        IF ( ll_teos10 ) THEN
           zwk(:,:,iup) = eosbn2(ztemp, zsal, gdep(jk), e3w, npiglo, npjglo ,iup, idown )* zmask(:,:)
        ELSE
           zwk(:,:,iup) = eosbn2(ztemp, zsal, gdep(jk), e3w, npiglo, npjglo ,iup, idown, .TRUE. )* zmask(:,:)
        ENDIF

        IF ( .NOT. l_w ) THEN
           ! now put zn2 at T level (k )
           WHERE ( zwk(:,:,idown) == 0 ) ; zn2(:,:) =  zwk(:,:,iup)
           ELSEWHERE                     ; zn2(:,:) = 0.5 * ( zwk(:,:,iup) + zwk(:,:,idown) ) * zmask(:,:)
           END WHERE
        ELSE
          zn2(:,:) = zwk(:,:,iup)
        ENDIF

        WHERE ( zn2 == 0 ) zn2 = -1000.
        ierr = putvar(ncout, id_varout(1), zn2, jk, npiglo, npjglo, ktime=jt )
        itmp = idown ; idown = iup ; iup = itmp

     END DO  ! loop to next level
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
    ipk(1)                       = npk  !  3D
    stypvar(1)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(1)%cname             = cv_bn2
    stypvar(1)%cunits            = 's-1'
    stypvar(1)%rmissing_value    = -1000.
    stypvar(1)%valid_min         = 0.
    stypvar(1)%valid_max         = 50000.
    stypvar(1)%clong_name        = 'Brunt_Vaissala_Frequency'
    stypvar(1)%cshort_name       = cv_bn2
    stypvar(1)%conline_operation = 'N/A'
    stypvar(1)%caxis             = 'TZYX'

    ! create output fileset
    ncout = create      (cf_out,   cf_tfil,  npiglo, npjglo, npk,                               ld_nc4=lnc4 )
    ierr  = createvar   (ncout ,   stypvar,  1,      ipk,    id_varout, cdglobal=TRIM(cglobal), ld_nc4=lnc4 )
    ierr  = putheadervar(ncout,    cf_tfil,  npiglo, npjglo, npk, pdep=gdep)

    dtim = getvar1d(cf_tfil, cn_vtimec, npt    )
    ierr = putvar1d(ncout,  dtim,       npt,'T')

  END SUBROUTINE CreateOutput

END PROGRAM
