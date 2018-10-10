PROGRAM cdfsig0
  !!======================================================================
  !!                     ***  PROGRAM  cdfsig0  ***
  !!=====================================================================
  !!  ** Purpose : Compute sigma0 3D field from gridT file
  !!               Store the results on a 'similar' cdf file.
  !!
  !!  ** Method  : Use NEMO equation of state
  !!
  !! History : 2.1  : 11/2006  : J.M. Molines : Original code
  !!           3.0  : 01/2011  : J.M. Molines : Doctor norm + Lic.
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
  !! @class Equation_of_state
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                           :: jk, jt             ! dummy loop index
  INTEGER(KIND=4)                           :: ierr               ! error status
  INTEGER(KIND=4)                           :: narg, ijarg,iargc        ! browse command line
  INTEGER(KIND=4)                           :: npiglo, npjglo     ! size of the domain
  INTEGER(KIND=4)                           :: npk, npkk, npt     ! size of the domain
  INTEGER(KIND=4)                           :: ncout              ! ncid of output file
  INTEGER(KIND=4), DIMENSION(1)             :: ipk, id_varout     ! level and  varid's

  REAL(KIND=4)                              :: zsps               ! missing value for salinity
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ztemp              ! temperature
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zsal               ! salinity
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zsig0              ! sigma-0
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zmask              ! 2D mask at current level

  REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: dtim               ! time counter

  CHARACTER(LEN=256)                        :: cf_tfil            ! input filename
  CHARACTER(LEN=256)                        :: cf_sfil            ! salinity file (option)
  CHARACTER(LEN=256)                        :: cf_out='sig0.nc'   ! output file name
  CHARACTER(LEN=256)                        :: cldum              ! dummy characte variable variable
  CHARACTER(LEN=256)                        :: cv_sal             ! salinity name in netcdf
  CHARACTER(LEN=256)                        :: cv_tem             ! temperature name in netcdf

  TYPE (variable), DIMENSION(1)             :: stypvar            ! structure for attributes
  LOGICAL                                   :: lnc4 = .FALSE.     ! flag for missing files
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfsig0 -t T-file [-s S-file] [-sal SAL-name] [-tem TEM-name] ...'
     PRINT *,'            ... [-nc4] [-o OUT-file]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute potential density (sigma-0) refered to the surface.' 
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -t T-file  : netcdf file with temperature and salinity.' 
     PRINT *,'           If salinity not in T-file, use -s option.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-s S-file]    : specify salinity file if not T-file.'
     PRINT *,'       [-sal SAL-name]  : name of salinity variable'
     PRINT *,'       [-tem TEM-name]  : name of temperature variable'
     PRINT *,'       [-nc4]  : enable chunking and compression'
     PRINT *,'       [-o OUT-file]    : specify output filename instead of ',TRIM(cf_out)
     PRINT *,'      '
     PRINT *,'     OPENMP SUPPORT : yes'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        none'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : ', TRIM(cn_vosigma0), ' ( kg/m3 - 1000 )'
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'       cdfsigi, cdfsigintu, signtr'
     STOP 
  ENDIF

  cv_sal=cn_vosaline
  cv_tem=cn_votemper
  cf_sfil='none'
  ijarg=1
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-t'   ) ; CALL getarg(ijarg, cf_tfil) ; ijarg=ijarg+1
!    options
     CASE ( '-s'   ) ; CALL getarg(ijarg, cf_sfil) ; ijarg=ijarg+1
     CASE ( '-o'   ) ; CALL getarg(ijarg, cf_out ) ; ijarg=ijarg+1
     CASE ( '-sal' ) ; CALL getarg(ijarg, cv_sal ) ; ijarg=ijarg+1
     CASE ( '-tem' ) ; CALL getarg(ijarg, cv_tem ) ; ijarg=ijarg+1
     CASE ( '-nc4' ) ; lnc4 = .TRUE.
     CASE DEFAULT    ; PRINT *,' ERROR : ', TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  ENDDO
  IF ( cf_sfil == 'none' ) cf_sfil=cf_tfil

  IF (chkfile(cf_tfil) .OR. chkfile(cf_sfil) ) STOP 99 ! missing file

  ! Look for missing value for salinity
  zsps = getspval(cf_sfil, cv_sal)

  npiglo = getdim (cf_tfil, cn_x)
  npjglo = getdim (cf_tfil, cn_y)
  npk    = getdim (cf_tfil, cn_z)
  npt    = getdim (cf_tfil, cn_t)

  npkk=npk
  IF ( npk == 0 ) npkk=1

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  ALLOCATE (ztemp(npiglo,npjglo), zsal (npiglo,npjglo) )
  ALLOCATE (zsig0(npiglo,npjglo), zmask(npiglo,npjglo) )
  ALLOCATE (dtim(npt) )

  CALL CreateOutput
  zsps = getspval( cf_sfil, cv_sal )

  DO jt=1,npt
     PRINT *,' TIME = ', jt, dtim(jt)/86400.,' days'
     DO jk = 1, npkk
        zmask(:,:)=1.

        ztemp(:,:)= getvar(cf_tfil, cv_tem, jk, npiglo, npjglo, ktime=jt)
        zsal(:,:) = getvar(cf_sfil, cv_sal, jk, npiglo, npjglo, ktime=jt)

        ! assuming spval is 0
        WHERE( zsal == zsps ) zmask = 0

        zsig0(:,:) = sigma0 (ztemp, zsal, npiglo, npjglo )* zmask(:,:)

        ierr = putvar(ncout, id_varout(1), zsig0, jk, npiglo, npjglo, ktime=jt)

     END DO  ! loop to next level
  END DO  ! next time frame

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

    ipk(:)                       = npkk  ! all variables (input and output are 3D)
    stypvar(1)%cname             = cn_vosigma0
    stypvar(1)%cunits            = 'kg/m3'
    stypvar(1)%rmissing_value    = 0.
    stypvar(1)%valid_min         = 0.001
    stypvar(1)%valid_max         = 40.
    stypvar(1)%clong_name        = 'Potential_density:sigma-0'
    stypvar(1)%cshort_name       = cn_vosigma0
    stypvar(1)%conline_operation = 'N/A'
    stypvar(1)%caxis             = 'TZYX'
    stypvar(1)%ichunk            = (/npiglo, MAX(1,npjglo/30), 1, 1 /)

    ! create output fileset
    ncout = create      (cf_out, cf_tfil, npiglo, npjglo, npk,       ld_nc4=lnc4  )
    ierr  = createvar   (ncout,  stypvar, 1,      ipk,    id_varout, ld_nc4=lnc4  )
    ierr  = putheadervar(ncout,  cf_tfil, npiglo, npjglo, npk       )

    dtim=getvar1d(cf_tfil, cn_vtimec, npt     )
    ierr=putvar1d(ncout,  dtim,       npt, 'T')

  END SUBROUTINE CreateOutput


END PROGRAM cdfsig0
