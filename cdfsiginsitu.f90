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

  REAL(KIND=4)                              :: zspval             ! missing value
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ztemp              ! temperature
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zsal               ! salinity
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zsigi              ! sigma-insitu
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zmask              ! 2D mask at current level
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: tim                ! time counter
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: gdept              ! depth of T points

  CHARACTER(LEN=256)                        :: cf_tfil             ! input filename
  CHARACTER(LEN=256)                        :: cf_out='siginsitu.nc' ! output file name

  TYPE (variable), DIMENSION(1)             :: stypvar            ! structure for attributes
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfsiginsitu T-file'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute in situ density from temperature and salinity.' 
     PRINT *,'       Depths are taken from input file.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       T-file : netcdf file with temperature and salinity.' 
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        none'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : vosigmainsitu (kg/m3 -1000 )'
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'      cdfsig0, cdfsigi '
     PRINT *,'      '
     STOP
  ENDIF

  CALL getarg (1, cf_tfil)

  IF ( chkfile(cf_tfil) ) STOP ! missing file

  npiglo = getdim (cf_tfil,cn_x)
  npjglo = getdim (cf_tfil,cn_y)
  npk    = getdim (cf_tfil,cn_z)
  npt    = getdim (cf_tfil,cn_t)

  ipk(:)= npk  ! all variables (input and output are 3D)
  stypvar(1)%cname             = 'vosigmainsitu'
  stypvar(1)%cunits            = 'kg/m3'
  stypvar(1)%rmissing_value    = 0.
  stypvar(1)%valid_min         = 0.001
  stypvar(1)%valid_max         = 45.
  stypvar(1)%clong_name        = 'in situ density'
  stypvar(1)%cshort_name       = 'vosigmainsitu'
  stypvar(1)%conline_operation = 'N/A'
  stypvar(1)%caxis             = 'TZYX'

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  ALLOCATE (ztemp(npiglo,npjglo), zsal (npiglo,npjglo) )
  ALLOCATE (zsigi(npiglo,npjglo), zmask(npiglo,npjglo) )
  ALLOCATE (gdept(npk), tim(npt)                       )

  ! create output fileset
  ncout = create      (cf_out, cf_tfil,  npiglo, npjglo, npk       )
  ierr  = createvar   (ncout,  stypvar,  1,      ipk,    id_varout )
  ierr  = putheadervar(ncout,  cf_tfil,  npiglo, npjglo, npk       )

  zspval = getatt(cf_tfil, cn_vosaline, 'missing_value')

  gdept = getvar1d(cf_tfil, cn_vdeptht, npk    )
  tim   = getvar1d(cf_tfil, cn_vtimec, npt     )
  ierr  = putvar1d(ncout,  tim,        npt, 'T')

  DO jt = 1, npt
     PRINT *,'time: ',jt
     DO jk = 1, npk
        zmask(:,:) = 1.

        ztemp(:,:) = getvar(cf_tfil, cn_votemper, jk, npiglo, npjglo, ktime=jt)
        zsal( :,:) = getvar(cf_tfil, cn_vosaline, jk, npiglo, npjglo, ktime=jt)

        WHERE( zsal == zspval ) zmask = 0

        zsigi(:,:) = sigmai(ztemp, zsal, gdept(jk), npiglo, npjglo )* zmask(:,:)

        ierr = putvar(ncout, id_varout(1), zsigi, jk, npiglo, npjglo, ktime=jt)

     END DO  ! loop to next level
  END DO  ! loop on time

  ierr = closeout(ncout)

END PROGRAM cdfsiginsitu
