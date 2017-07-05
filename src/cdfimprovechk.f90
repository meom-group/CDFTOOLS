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
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2011, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                           :: jk, jt            ! dummy loop index
  INTEGER(KIND=4)                           :: ierr              ! working integer
  INTEGER(KIND=4)                           :: narg, iargc       ! browse line
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
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: tim               ! time counter

  CHARACTER(LEN=256)                        :: cf_obs            ! observation-file name
  CHARACTER(LEN=256)                        :: cf_ref            ! reference-file name
  CHARACTER(LEN=256)                        :: cf_tst            ! test-file name
  CHARACTER(LEN=256)                        :: cv_in             ! cdf variable name
  CHARACTER(LEN=256)                        :: cf_out='chk.nc'   ! output filename

  TYPE (variable), DIMENSION(1)             :: stypvar           ! structure for attributes

  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfimprovechk IN-var OBS-file REF-file TST-file'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'        Estimate the improvement/deterioration of a test run,'
     PRINT *,'        compared with a reference run relative to some observations'
     PRINT *,'        This program computes the quantity zchk= ( REF - TEST )/(REF - OBS)'
     PRINT *,'        Where 0 < zchk <= 1, the TST is better than the reference'
     PRINT *,'        Where 1 < zchk, the TST  was corrected in the right sense but too much'
     PRINT *,'        Where  zchk < 0, the TST  was corrected was corrected in the wrong way.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'        IN-var    : netcdf input variable'
     PRINT *,'        OBS-file  : netcdf observation file'
     PRINT *,'        REF-file  : netcdf reference file'
     PRINT *,'        TST-file  : netcdf test file'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       none' 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : same as input variable.'
     STOP
  ENDIF

  CALL getarg (1, cv_in )
  CALL getarg (2, cf_obs)
  CALL getarg (3, cf_ref)
  CALL getarg (4, cf_tst)

  IF ( chkfile(cf_obs) .OR. chkfile(cf_ref) .OR. chkfile(cf_tst) ) STOP 99 ! missing files

  npiglo = getdim(cf_ref, cn_x)
  npjglo = getdim(cf_ref, cn_y)
  npk    = getdim(cf_ref, cn_z)
  npt    = getdim(cf_ref, cn_t)

  nvpk   = getvdim(cf_ref, cv_in)
  IF (nvpk == 2 ) nvpk = 1
  IF (nvpk == 3 ) nvpk = npk

  ipk(:)                       = nvpk  ! all variables 
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


  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  ALLOCATE (zobs(npiglo,npjglo), zref(npiglo,npjglo), ztst(npiglo,npjglo), zmask(npiglo,npjglo))
  ALLOCATE (zchk(npiglo,npjglo), tim(npt) )

  ! create output fileset

  ncout = create      (cf_out, cf_ref,  npiglo, npjglo, npk       )
  ierr  = createvar   (ncout,  stypvar, 1,      ipk,    id_varout )
  ierr  = putheadervar(ncout,  cf_ref,  npiglo, npjglo, npk       )

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

  tim  = getvar1d(cf_ref, cn_vtimec, npt     )
  ierr = putvar1d(ncout,  tim,       npt, 'T')

  ierr = closeout(ncout)

END PROGRAM cdfimprovechk
