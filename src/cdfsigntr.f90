PROGRAM cdfsigntr
  !!======================================================================
  !!                     ***  PROGRAM  cdfsigntr  ***
  !!=====================================================================
  !!  ** Purpose : Compute ntr volumic mass 3D field from gridT file
  !!               Store the results on a 'similar' cdf file.
  !!
  !!  ** Method  : Use Mc Dougall Jackett EOS (2005)
  !!
  !! History : 3.0  : 06/2013  : J.M. Molines from G. Madec idl code
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
  INTEGER(KIND=4)                           :: narg, iargc        ! browse command line
  INTEGER(KIND=4)                           :: npiglo, npjglo     ! size of the domain
  INTEGER(KIND=4)                           :: npk, npt           ! size of the domain
  INTEGER(KIND=4)                           :: ncout              ! ncid of output file
  INTEGER(KIND=4), DIMENSION(1)             :: ipk, id_varout     ! level and  varid's

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ztemp              ! temperature
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zsal               ! salinity
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zsigntr            ! sigma-0
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zmask              ! 2D mask at current level
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: tim                ! time counter

  CHARACTER(LEN=256)                        :: cf_tfil            ! input filename
  CHARACTER(LEN=256)                        :: cf_out='signtr.nc'   ! output file name

  TYPE (variable), DIMENSION(1)             :: stypvar            ! structure for attributes
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfsigntr T-file'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute neutral volumic mass (kg/m3) from temperature and salinity.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       T-file  : netcdf file with temperature and salinity.' 
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        none'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : ', TRIM(cn_vosigntr), ' ( kg/m3 )'
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'       cdfsig0, cdfsigi, cdfsiginsitu'
     STOP
  ENDIF

  CALL getarg (1, cf_tfil)
  IF (chkfile(cf_tfil) ) STOP ! missing file

  npiglo = getdim (cf_tfil, cn_x)
  npjglo = getdim (cf_tfil, cn_y)
  npk    = getdim (cf_tfil, cn_z)
  npt    = getdim (cf_tfil, cn_t)

  ipk(:)                       = npk  ! all variables (input and output are 3D)
  stypvar(1)%cname             = cn_vosigntr
  stypvar(1)%cunits            = 'kg/m3'
  stypvar(1)%rmissing_value    = 0.
  stypvar(1)%valid_min         = 0.001
  stypvar(1)%valid_max         = 1040.
  stypvar(1)%clong_name        = 'Neutral volumic mass'
  stypvar(1)%cshort_name       = cn_vosigntr
  stypvar(1)%conline_operation = 'N/A'
  stypvar(1)%caxis             = 'TZYX'

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  ALLOCATE (ztemp(npiglo,npjglo), zsal (npiglo,npjglo) )
  ALLOCATE (zsigntr(npiglo,npjglo), zmask(npiglo,npjglo) )
  ALLOCATE (tim(npt) )

  ! create output fileset
  ncout = create      (cf_out, cf_tfil, npiglo, npjglo, npk       )
  ierr  = createvar   (ncout,  stypvar, 1,      ipk,    id_varout )
  ierr  = putheadervar(ncout,  cf_tfil, npiglo, npjglo, npk       )

  tim=getvar1d(cf_tfil, cn_vtimec, npt     )
  ierr=putvar1d(ncout,  tim,       npt, 'T')

  DO jt=1,npt
     PRINT *,' TIME = ', jt, tim(jt)/86400.,' days'
     DO jk = 1, npk
        PRINT *, 'level : ', jk
        zmask(:,:)=1.

        ztemp(:,:)= getvar(cf_tfil, cn_votemper, jk, npiglo, npjglo, ktime=jt)
        zsal(:,:) = getvar(cf_tfil, cn_vosaline, jk, npiglo, npjglo, ktime=jt)

        ! assuming spval is 0
        WHERE( zsal == 0 ) zmask = 0

        zsigntr(:,:) = sigmantr (ztemp, zsal, npiglo, npjglo )* zmask(:,:)

        ierr = putvar(ncout, id_varout(1), zsigntr, jk, npiglo, npjglo, ktime=jt)

     END DO  ! loop to next level
  END DO  ! next time frame

  ierr = closeout(ncout)

END PROGRAM cdfsigntr
