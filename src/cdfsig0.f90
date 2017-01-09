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
  !!----------------------------------------------------------------------
  USE cdfio
  USE eos
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2011, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                           :: jk, jt             ! dummy loop index
  INTEGER(KIND=4)                           :: ierr               ! error status
  INTEGER(KIND=4)                           :: narg, ijarg,iargc        ! browse command line
  INTEGER(KIND=4)                           :: npiglo, npjglo     ! size of the domain
  INTEGER(KIND=4)                           :: npk, npkk, npt     ! size of the domain
  INTEGER(KIND=4)                           :: ncout              ! ncid of output file
  INTEGER(KIND=4), DIMENSION(1)             :: ipk, id_varout     ! level and  varid's

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ztemp              ! temperature
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zsal               ! salinity
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zsig0              ! sigma-0
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zmask              ! 2D mask at current level
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: tim                ! time counter

  CHARACTER(LEN=256)                        :: cf_tfil            ! input filename
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
     PRINT *,' usage : cdfsig0 -t T-file [-sal SAL-name] [-tem TEM-name] [-nc4] [-o OUT-file]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute potential density (sigma-0) refered to the surface.' 
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -t T-file  : netcdf file with temperature and salinity.' 
     PRINT *,'         (for backward compatibility, -t can be ommited when T-file is the only'
     PRINT *,'         argument.)'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-sal SAL-name]  : name of salinity variable'
     PRINT *,'       [-tem TEM-name]  : name of temperature variable'
     PRINT *,'       [-nc4]  : enable chunking and compression'
     PRINT *,'       [-o OUT-file]    : specify output filename instead of ',TRIM(cf_out)
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        none'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : ', TRIM(cn_vosigma0), ' ( kg/m3 - 1000 )'
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'       cdfsigi'
     STOP
  ENDIF

  ijarg=1
  cv_sal=cn_vosaline
  cv_tem=cn_votemper
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-t'   ) ; CALL getarg(ijarg, cf_tfil) ; ijarg=ijarg+1
     CASE ( '-o'   ) ; CALL getarg(ijarg, cf_out ) ; ijarg=ijarg+1
     CASE ( '-sal' ) ; CALL getarg(ijarg, cv_sal ) ; ijarg=ijarg+1
     CASE ( '-tem' ) ; CALL getarg(ijarg, cv_tem ) ; ijarg=ijarg+1
     CASE ( '-nc4' ) ; lnc4 = .TRUE.
     CASE DEFAULT
      IF ( narg == 1 ) THEN
         cf_tfil = cldum
      ELSE
         PRINT *,' option ',TRIM(cldum),' not understood'
         STOP
      ENDIF
     END SELECT
  ENDDO
 
  IF (chkfile(cf_tfil) ) STOP ! missing file

  npiglo = getdim (cf_tfil, cn_x)
  npjglo = getdim (cf_tfil, cn_y)
  npk    = getdim (cf_tfil, cn_z)
  npt    = getdim (cf_tfil, cn_t)

  npkk=npk
  IF ( npk == 0 ) npkk=1

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


  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  ALLOCATE (ztemp(npiglo,npjglo), zsal (npiglo,npjglo) )
  ALLOCATE (zsig0(npiglo,npjglo), zmask(npiglo,npjglo) )
  ALLOCATE (tim(npt) )

  ! create output fileset
  ncout = create      (cf_out, cf_tfil, npiglo, npjglo, npk,       ld_nc4=lnc4  )
  ierr  = createvar   (ncout,  stypvar, 1,      ipk,    id_varout, ld_nc4=lnc4  )
  ierr  = putheadervar(ncout,  cf_tfil, npiglo, npjglo, npk       )

  tim=getvar1d(cf_tfil, cn_vtimec, npt     )
  ierr=putvar1d(ncout,  tim,       npt, 'T')

  DO jt=1,npt
     PRINT *,' TIME = ', jt, tim(jt)/86400.,' days'
     DO jk = 1, npkk
        zmask(:,:)=1.

        ztemp(:,:)= getvar(cf_tfil, cv_tem, jk, npiglo, npjglo, ktime=jt)
        zsal(:,:) = getvar(cf_tfil, cv_sal, jk, npiglo, npjglo, ktime=jt)

        ! assuming spval is 0
        WHERE( zsal <= 0 ) zmask = 0

        zsig0(:,:) = sigma0 (ztemp, zsal, npiglo, npjglo )* zmask(:,:)

        ierr = putvar(ncout, id_varout(1), zsig0, jk, npiglo, npjglo, ktime=jt)

     END DO  ! loop to next level
  END DO  ! next time frame

  ierr = closeout(ncout)

END PROGRAM cdfsig0
