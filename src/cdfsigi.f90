PROGRAM cdfsigi
  !!======================================================================
  !!                     ***  PROGRAM  cdfsigi  ***
  !!=====================================================================
  !!  ** Purpose : Compute sigmai 3D field from gridT file
  !!                Store the results on a 'similar' cdf file.
  !!
  !!  **  Method: read temp and salinity, compute sigma-i
  !!              using depth given in argument (meters or dbar)
  !!
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
  INTEGER(KIND=4)                           :: npk, npt           ! size of the domain
  INTEGER(KIND=4)                           :: ncout              ! ncid of output file
  INTEGER(KIND=4), DIMENSION(1)             :: ipk, id_varout     ! level and  varid's

  REAL(KIND=4)                              :: ref_dep            ! reference depth in meters
  REAL(KIND=4)                              :: zspval             ! missing value
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ztemp              ! temperature
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zsal               ! salinity
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zsigi              ! sigma-i
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zmask              ! 2D mask at current level

  REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: dtim               ! time counter

  CHARACTER(LEN=256)                        :: cf_tfil            ! input filename
  CHARACTER(LEN=256)                        :: cf_sfil            ! salinity file (option)
  CHARACTER(LEN=256)                        :: cf_out='sigi.nc'   ! output file name
  CHARACTER(LEN=256)                        :: cldum              ! dummy string

  TYPE (variable), DIMENSION(1)             :: stypvar            ! structure for attributes

  LOGICAL                                   :: lnc4 = .FALSE.     ! Use nc4 with chunking and deflation
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfsigi -t T-file -r REF-dep(m) [-s S-file] [-o OUT-file] [-nc4] '
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute potential density referred to the depth given in arguments.' 
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -t T-file  : netcdf file with temperature and salinity.' 
     PRINT *,'          If salinity not in T-file, use -s option.'
     PRINT *,'       -r REF-dep : reference depth in meter.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'      [-s S-file   ] : Specify salinity file if not T-file.'
     PRINT *,'      [-o OUT-file ] : Specify output file name instead of ',TRIM(cf_out)
     PRINT *,'      [-nc4 ] : Use netcdf4 output with chunking and deflation level 1.'
     PRINT *,'           This option is effective only if cdftools are compiled with'
     PRINT *,'           a netcdf library supporting chunking and deflation.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        none'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : ', TRIM(cn_vosigmai),' (kg/m3 -1000 )'
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'      cdfsig0, cdfsiginintu' 
     PRINT *,'      '
     STOP 
  ENDIF

  ijarg=1
  cf_sfil='none'
  DO WHILE ( ijarg <= narg ) 
     CALL getarg (ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-t'   ) ; CALL getarg (ijarg, cf_tfil) ; ijarg=ijarg+1
     CASE ( '-r'   ) ; CALL getarg (ijarg, cldum  ) ; ijarg=ijarg+1 ; READ(cldum,*) ref_dep
        ! options
     CASE ( '-s'   ) ; CALL getarg (ijarg, cf_sfil) ; ijarg=ijarg+1
     CASE ( '-o'   ) ; CALL getarg (ijarg, cf_out ) ; ijarg=ijarg+1
     CASE ( '-nc4' ) ; lnc4 = .TRUE.
     CASE DEFAULT    ; PRINT *,' ERROR : ', TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  ENDDO

  IF ( cf_sfil == 'none' ) cf_sfil=cf_tfil

  IF ( chkfile(cf_tfil) .OR. chkfile(cf_sfil) ) STOP 99 ! missing file

  npiglo = getdim (cf_tfil, cn_x)
  npjglo = getdim (cf_tfil, cn_y)
  npk    = getdim (cf_tfil, cn_z)
  npt    = getdim (cf_tfil, cn_t)

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  ALLOCATE (ztemp(npiglo,npjglo), zsal (npiglo,npjglo) )
  ALLOCATE (zsigi(npiglo,npjglo), zmask(npiglo,npjglo) )
  ALLOCATE (dtim(npt) )

  CALL CreateOutput
  zspval= getatt(cf_sfil, cn_vosaline, cn_missing_value)
  DO jt = 1, npt
     PRINT *,'time: ',jt
     DO jk = 1, npk
        zmask(:,:) = 1.

        ztemp(:,:) = getvar(cf_tfil, cn_votemper,  jk, npiglo, npjglo, ktime=jt)
        zsal( :,:) = getvar(cf_sfil, cn_vosaline,  jk, npiglo, npjglo, ktime=jt)

        WHERE( zsal == zspval ) zmask = 0

        zsigi(:,:) = sigmai(ztemp, zsal, ref_dep, npiglo, npjglo )* zmask(:,:)

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
    ipk(:)= npk  ! all variables (input and output are 3D)
    stypvar(1)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(1)%cname             = cn_vosigmai
    stypvar(1)%cunits            = 'kg/m3'
    stypvar(1)%rmissing_value    = 0.
    stypvar(1)%valid_min         = 0.001
    stypvar(1)%valid_max         = 45.
    stypvar(1)%clong_name        = 'Potential_density:refered to '//TRIM(cldum)//' m'
    stypvar(1)%cshort_name       = cn_vosigmai
    stypvar(1)%conline_operation = 'N/A'
    stypvar(1)%caxis             = 'TZYX'

    ! create output fileset
    ncout = create      (cf_out, cf_tfil, npiglo, npjglo, npk       , ld_nc4=lnc4 )
    ierr  = createvar   (ncout,  stypvar, 1,      ipk,    id_varout , ld_nc4=lnc4 )
    ierr  = putheadervar(ncout,  cf_tfil, npiglo, npjglo, npk                     )

    dtim  = getvar1d(cf_tfil, cn_vtimec, npt     )
    ierr  = putvar1d(ncout,   dtim,      npt, 'T')

  END SUBROUTINE CreateOutput

END PROGRAM cdfsigi
