PROGRAM cdfsiginsitu
  !!======================================================================
  !!                     ***  PROGRAM  cdfsiginsitu  ***
  !!=====================================================================
  !!  ** Purpose : Compute sigma insitu 3D field from gridT file
  !!                Store the results on a 'similar' cdf file.
  !!
  !!  **  Method: read temp and salinity, compute sigma insitu
  !!              using depth taken from input T file

  !! History : 2.0  : 11/2004  : J.M. Molines : Original code
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
  INTEGER(KIND=4)                           :: narg, iargc, ijarg ! browse command line
  INTEGER(KIND=4)                           :: npiglo, npjglo     ! size of the domain
  INTEGER(KIND=4)                           :: npk, npkk, npt     ! size of the domain
  INTEGER(KIND=4)                           :: ncout              ! ncid of output file
  INTEGER(KIND=4), DIMENSION(1)             :: ipk, id_varout     ! level and  varid's

  REAL(KIND=4)                              :: zspval             ! missing value
  REAL(KIND=4)                              :: dep=0.0            ! depth to be used
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ztemp              ! temperature
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zsal               ! salinity
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zsigi              ! sigma-insitu
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zmask              ! 2D mask at current level
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: gdept              ! depth of T points

  REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: dtim               ! time counter

  CHARACTER(LEN=256)                        :: cf_tfil             ! input filename
  CHARACTER(LEN=256)                        :: cf_out='siginsitu.nc' ! output file name
  CHARACTER(LEN=256)                        :: cldum              ! dummy characte variable variable
  CHARACTER(LEN=256)                        :: cv_sal             ! salinity name in netcdf
  CHARACTER(LEN=256)                        :: cv_tem             ! temperature name in netcdf

  TYPE (variable), DIMENSION(1)             :: stypvar            ! structure for attributes
  LOGICAL                                   :: lnc4 = .FALSE.     ! flag for missing files
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfsiginsitu -t T-file [-sal SAL-name] [-tem TEM-name ] ...'
     PRINT *,'                [-dep depth] [-o OUT-file ] [-nc4 ] '
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute in situ density from temperature and salinity. Depths are taken' 
     PRINT *,'       from input file.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -t T-file : netcdf file with temperature and salinity.' 
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-sal SAL-name] : name of salinity variable'
     PRINT *,'       [-tem TEM-name] : name of temperature variable'
     PRINT *,'       [-dep depth ]   : depth to be used in case of 2D input file (only)'
     PRINT *,'       [-nc4]          : enable chunking and compression'
     PRINT *,'       [-o OUT-file]   : specify output filename instead of ',TRIM(cf_out)
     PRINT *,'      '
     PRINT *,'     OPENMP SUPPORT : yes'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        none'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) ,' or the file name specified '
     PRINT *,'                   with -o option'
     PRINT *,'         variables : vosigmainsitu (kg/m3 -1000 )'
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'      cdfsig0, cdfsigi, cdfsigntr '
     PRINT *,'      '
     STOP 
  ENDIF

  cv_sal=cn_vosaline
  cv_tem=cn_votemper
  ijarg=1
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-t'   ) ; CALL getarg(ijarg, cf_tfil) ; ijarg=ijarg+1
     CASE ( '-o'   ) ; CALL getarg(ijarg, cf_out ) ; ijarg=ijarg+1
     CASE ( '-sal' ) ; CALL getarg(ijarg, cv_sal ) ; ijarg=ijarg+1
     CASE ( '-tem' ) ; CALL getarg(ijarg, cv_tem ) ; ijarg=ijarg+1
     CASE ( '-dep' ) ; CALL getarg(ijarg, cldum  ) ; ijarg=ijarg+1 ; READ(cldum,*) dep 
     CASE ( '-nc4' ) ; lnc4 = .TRUE.
     CASE DEFAULT    ; PRINT *,' ERROR : ', TRIM(cldum) ,' : unknown option.' ; STOP 99
     END SELECT
  ENDDO

  IF ( chkfile(cf_tfil) ) STOP 99 ! missing file

  npiglo = getdim (cf_tfil,cn_x)
  npjglo = getdim (cf_tfil,cn_y)
  npk    = getdim (cf_tfil,cn_z)
  npt    = getdim (cf_tfil,cn_t)

  IF ( npk == 0 ) THEN ; npkk = 1
  ELSE                 ; npkk = npk  ! all variables (input and output are 3D)
  ENDIF

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  ALLOCATE (ztemp(npiglo,npjglo), zsal (npiglo,npjglo) )
  ALLOCATE (zsigi(npiglo,npjglo), zmask(npiglo,npjglo) )
  ALLOCATE (gdept(npkk), dtim(npt)                     )

  CALL CreateOutput
  zspval = getatt(cf_tfil, cv_sal, 'missing_value')

  IF ( npk == 0 ) THEN ; gdept(:)= dep
  ELSE                 ; gdept = getvar1d(cf_tfil, cn_vdeptht, npk    )
  ENDIF

  DO jt = 1, npt
     PRINT *,'time: ',jt
     DO jk = 1, npkk
        zmask(:,:) = 1.

        ztemp(:,:) = getvar(cf_tfil, cv_tem, jk, npiglo, npjglo, ktime=jt)
        zsal( :,:) = getvar(cf_tfil, cv_sal, jk, npiglo, npjglo, ktime=jt)

        WHERE( zsal == zspval ) zmask = 0

        zsigi(:,:) = sigmai(ztemp, zsal, gdept(jk), npiglo, npjglo )* zmask(:,:)

        ierr = putvar(ncout, id_varout(1), zsigi, jk, npiglo, npjglo, ktime=jt)

     END DO  ! loop to next level
  END DO  ! loop on time

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
    ipk(:) = npkk
    stypvar(1)%ichunk            = (/npiglo, MAX(1,npjglo/30), 1, 1 /)
    stypvar(1)%cname             = 'vosigmainsitu'
    stypvar(1)%cunits            = 'kg/m3'
    stypvar(1)%rmissing_value    = 0.
    stypvar(1)%valid_min         = 0.001
    stypvar(1)%valid_max         = 45.
    stypvar(1)%clong_name        = 'in situ density'
    stypvar(1)%cshort_name       = 'vosigmainsitu'
    stypvar(1)%conline_operation = 'N/A'
    stypvar(1)%caxis             = 'TZYX'

    ! create output fileset
    ncout = create      (cf_out, cf_tfil,  npiglo, npjglo, npk,       ld_nc4=lnc4 )
    ierr  = createvar   (ncout,  stypvar,  1,      ipk,    id_varout, ld_nc4=lnc4 )
    ierr  = putheadervar(ncout,  cf_tfil,  npiglo, npjglo, npk       )

    dtim  = getvar1d(cf_tfil, cn_vtimec, npt     )
    ierr  = putvar1d(ncout,  dtim,       npt, 'T')

  END SUBROUTINE CreateOutput

END PROGRAM cdfsiginsitu
