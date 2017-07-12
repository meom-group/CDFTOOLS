PROGRAM cdfimprovechk
  !!======================================================================
  !!                     ***  PROGRAM  cdfimprovechk  ***
  !!=====================================================================
  !!  ** Purpose : Estimate the improvement/deterioration of a test run, 
  !!               compared with a reference run relative to some observations
  !!
  !!  ** Method  : Given zobs (observed field), zref (reference run field)
  !!               and ztst (test run field), compute zchk as the ratio:
  !!                 zchk=(zref - ztst) / (zref - zobs )
  !!
  !!               Where  0 < zchk <=1 correction act in the right direction
  !!               Where  1 < zchk     correction is too strong, in the right way
  !!               Where  zchk < 0     correction is in the wrong way (deterioration)
  !!
  !! History : 2.1  : 11/2005  : J.M. Molines : Original code
  !!           3.0  : 01/2011  : J.M. Molines : Doctor norm + Lic.
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class miscellaneous
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                           :: jk, jt            ! dummy loop index
  INTEGER(KIND=4)                           :: ierr              ! working integer
  INTEGER(KIND=4)                           :: narg, iargc, ijarg! browse line
  INTEGER(KIND=4)                           :: npiglo, npjglo    ! size of the domain
  INTEGER(KIND=4)                           :: npk, npt          ! size of the domain
  INTEGER(KIND=4)                           :: nvpk              ! dim of the working variable
  INTEGER(KIND=4)                           :: ncout             ! ncid of output file
  INTEGER(KIND=4), DIMENSION(1)             :: ipk, id_varout    ! levels and varid of output vars

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zobs              ! observation array
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zref              ! reference array
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ztst              ! test array
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zmask             ! 2D mask at surface
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zchk              ! check index output

  REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: dtim              ! time counter

  CHARACTER(LEN=256)                        :: cf_obs            ! observation-file name
  CHARACTER(LEN=256)                        :: cf_ref            ! reference-file name
  CHARACTER(LEN=256)                        :: cf_tst            ! test-file name
  CHARACTER(LEN=256)                        :: cv_in             ! cdf variable name
  CHARACTER(LEN=256)                        :: cf_out='chk.nc'   ! output filename
  CHARACTER(LEN=256)                        :: cldum             ! working char variable

  TYPE (variable), DIMENSION(1)             :: stypvar           ! structure for attributes

  LOGICAL                                   :: lnc4 = .FALSE.    ! Use nc4 with chunking and deflation
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfimprovechk -v IN-var -obs OBS-file -r REF-file -t TST-file ...'
     PRINT *,'                ... [-o OUT-file] [-nc4]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'        Estimates the improvement/deterioration of a test run, compared with a'
     PRINT *,'        reference run relative to some observations.'
     PRINT *,'        This program computes the field zchk= ( REF - TST )/(REF - OBS).'
     PRINT *,'        Where 0 < zchk <= 1, the TST is better than the reference'
     PRINT *,'        Where 1 < zchk, the TST  was corrected in the right sense but too much'
     PRINT *,'        Where  zchk < 0, the TST  was corrected was corrected in the wrong way.'
     PRINT *,'        Although not very much used, this program is maintained as one of the'
     PRINT *,'        first CDFTOOLS.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'        -v IN-var    : netcdf input variable'
     PRINT *,'        -obs OBS-file: netcdf observation file'
     PRINT *,'        -r REF-file  : netcdf reference file'
     PRINT *,'        -t TST-file  : netcdf test file'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'        [-o OUT-file] : specifiy the output file name instead of ',TRIM(cf_out)
     PRINT *,'        [-nc4]    : Use netcdf4 output with chunking and deflation level 1.'
     PRINT *,'               This option is effective only if cdftools are compiled with'
     PRINT *,'               a netcdf library supporting chunking and deflation.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       none' 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : same as input variable.'
     PRINT *,'      '
     STOP 
  ENDIF

  ijarg=1
  DO WHILE (ijarg <= narg )
     CALL getarg(ijarg,cldum) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-v'  ) ; CALL getarg(ijarg,cv_in ) ; ijarg=ijarg+1
     CASE ( '-obs') ; CALL getarg(ijarg,cf_obs) ; ijarg=ijarg+1
     CASE ( '-r'  ) ; CALL getarg(ijarg,cf_ref) ; ijarg=ijarg+1
     CASE ( '-t'  ) ; CALL getarg(ijarg,cf_tst) ; ijarg=ijarg+1
        ! options
     CASE ( '-o'  ) ; CALL getarg(ijarg,cf_out) ; ijarg=ijarg+1
     CASE ( '-nc4') ; lnc4 = .TRUE.
     CASE DEFAULT   ; PRINT *,' ERROR : ',TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  ENDDO

  IF ( chkfile(cf_obs) .OR. chkfile(cf_ref) .OR. chkfile(cf_tst) ) STOP 99 ! missing files

  npiglo = getdim(cf_ref, cn_x)
  npjglo = getdim(cf_ref, cn_y)
  npk    = getdim(cf_ref, cn_z)
  npt    = getdim(cf_ref, cn_t)

  nvpk   = getvdim(cf_ref, cv_in)
  IF (nvpk == 2 ) nvpk = 1
  IF (nvpk == 3 ) nvpk = npk

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  ALLOCATE (zobs(npiglo,npjglo), zref(npiglo,npjglo), ztst(npiglo,npjglo), zmask(npiglo,npjglo))
  ALLOCATE (zchk(npiglo,npjglo), dtim(npt) )

  CALL CreateOutput
  zref  = 0. ;   zobs  = 0. ;  zmask = 1.

  DO jt = 1,npt
     DO jk = 1,npk
        PRINT *,'level ',jk
        zchk = 0.
        zobs(:,:) = getvar(cf_obs, cv_in, jk ,npiglo, npjglo, ktime=jt)
        zref(:,:) = getvar(cf_ref, cv_in, jk ,npiglo, npjglo, ktime=jt)
        ztst(:,:) = getvar(cf_tst, cv_in, jk ,npiglo, npjglo, ktime=jt)

        IF (jk == 1  )  THEN
           WHERE( zref == 0. ) zmask = 0.
        END IF
        WHERE  ( (zref - zobs ) /= 0 ) 
           zchk = (zref - ztst ) / ( zref - zobs) * zmask
        END WHERE
        ierr = putvar(ncout, id_varout(1), zchk, jk, npiglo, npjglo, ktime=jt)
     END DO
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
    ipk(:)                       = nvpk  ! all variables 
    stypvar(1)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(1)%cname             = TRIM(cv_in)
    stypvar(1)%cunits            = '%'
    stypvar(1)%rmissing_value    = 0.
    stypvar(1)%valid_min         = 0.
    stypvar(1)%valid_max         = 100.
    stypvar(1)%clong_name        = 'Checking ratio for'//TRIM(cv_in)
    stypvar(1)%cshort_name       = cv_in
    stypvar(1)%conline_operation = 'N/A'

    IF (nvpk == npk ) stypvar(1)%caxis='TZYX'
    IF (nvpk == 1   ) stypvar(1)%caxis='TYX'

    ! create output fileset

    ncout = create      (cf_out, cf_ref,  npiglo, npjglo, npk       , ld_nc4=lnc4 )
    ierr  = createvar   (ncout,  stypvar, 1,      ipk,    id_varout , ld_nc4=lnc4 )
    ierr  = putheadervar(ncout,  cf_ref,  npiglo, npjglo, npk       )

    dtim = getvar1d(cf_ref, cn_vtimec, npt     )
    ierr = putvar1d(ncout,  dtim,      npt, 'T')

  END SUBROUTINE CreateOutput

END PROGRAM cdfimprovechk
