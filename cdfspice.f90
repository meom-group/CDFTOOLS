PROGRAM cdfspice
  !!======================================================================
  !!                     ***  PROGRAM  cdfspice  ***
  !!=====================================================================
  !!  ** Purpose : Compute spiciness 3D field from gridT file
  !!               Store the results on a 'similar' cdf file.
  !!
  !!  ** Method  :  spiciness = sum(i=0,5)[sum(j=0,4)[b(i,j)*theta^i*(s-35)^j]]
  !!                   with:  b     -> coefficients
  !!                          theta -> potential temperature
  !!                          s     -> salinity
  !!
  !!  **  Example:
  !!       spice(15,33)=   0.5445863      0.544586321373410  calcul en double
  !!       spice(15,33)=   0.5445864      (calcul en simple precision)
  !!
  !!  ** References : Flament (2002) "A state variable for characterizing 
  !!              water masses and their diffusive stability: spiciness."
  !!              Progress in Oceanography Volume 54, 2002, Pages 493-501.
  !!
  !! History : 2.1  : 03/2010  : C.O. Dufour  : Original code
  !!           3.0  : 01/2011  : J.M. Molines : Doctor norm + Lic.
  !!           3.0  : 06/2013  : J.M. Molines : Transfert spice in eos
  !!                                            Add OMP directive
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

  INTEGER(KIND=4)                           :: ji, jj, jk, jt     ! dummy loop index
  INTEGER(KIND=4)                           :: ierr               ! error status
  INTEGER(KIND=4)                           :: narg, iargc        ! browse command line
  INTEGER(KIND=4)                           :: npiglo, npjglo     ! size of the domain
  INTEGER(KIND=4)                           :: npk, npt           ! size of the domain
  INTEGER(KIND=4)                           :: ncout              ! ncid of output file
  INTEGER(KIND=4), DIMENSION(1)             :: ipk, id_varout     ! level and  varid's

  REAL(KIND=4)                              :: zspval             ! missing value
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: tim                ! time counter
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ztemp              ! temperature
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zsal               ! salinity
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zmask              ! 2D mask at current level

  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dtempt             ! temperature
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dsalt              ! salinity
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dsalref            ! reference salinity
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dspi               ! spiceness

  CHARACTER(LEN=256)                        :: cf_tfil            ! input filename
  CHARACTER(LEN=256)                        :: cf_out='spice.nc'  ! output file name

  TYPE (variable), DIMENSION(1)             :: stypvar            ! structure for attributes
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfspice T-file '
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the spiceness corresponding to temperatures and salinities'
     PRINT *,'       given in the input file.' 
     PRINT *,'      '
     PRINT *,'       spiciness = sum(i=0,5)[sum(j=0,4)[b(i,j)*theta^i*(s-35)^j]]'
     PRINT *,'                 with:  b     -> coefficients'
     PRINT *,'                        theta -> potential temperature'
     PRINT *,'                        s     -> salinity'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       T-file : netcdf file with temperature and salinity (gridT)' 
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       none'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : vospice'
     PRINT *,'      '
     PRINT *,'     REFERENCE :'
     PRINT *,'       Flament (2002) "A state variable for characterizing '
     PRINT *,'             water masses and their diffusive stability: spiciness."'
     PRINT *,'             Progress in Oceanography Volume 54, 2002, Pages 493-501.'
     STOP
  ENDIF
  IF ( narg == 0 ) THEN
     PRINT *,'usage : cdfspice  gridT '
     PRINT *,'    Output on spice.nc, variable vospice'
     STOP
  ENDIF

  CALL getarg (1, cf_tfil)

  IF ( chkfile(cf_tfil) ) STOP ! missing files

  npiglo = getdim (cf_tfil,cn_x)
  npjglo = getdim (cf_tfil,cn_y)
  npk    = getdim (cf_tfil,cn_z)
  npt    = getdim (cf_tfil,cn_t)

  ipk(:)                       = npk 
  stypvar(1)%cname             = 'vospice'
  stypvar(1)%cunits            = 'kg/m3'
  stypvar(1)%rmissing_value    = 0.
  stypvar(1)%valid_min         = -300.
  stypvar(1)%valid_max         = 300.
  stypvar(1)%clong_name        = 'spiciness'
  stypvar(1)%cshort_name       = 'vospice'
  stypvar(1)%conline_operation = 'N/A'
  stypvar(1)%caxis             = 'TZYX'

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  ALLOCATE (ztemp(npiglo,npjglo), zsal (npiglo,npjglo) )
  ALLOCATE (dspi( npiglo,npjglo), zmask(npiglo,npjglo) )
  ALLOCATE (dtempt(npiglo,npjglo), dsalt(npiglo,npjglo))
  ALLOCATE (dsalref(npiglo,npjglo))
  ALLOCATE (tim(npt))

  ! create output fileset
  ncout = create      (cf_out, cf_tfil, npiglo, npjglo, npk       )
  ierr  = createvar   (ncout,  stypvar, 1,      ipk,    id_varout )
  ierr  = putheadervar(ncout,  cf_tfil, npiglo, npjglo, npk       )

  tim  = getvar1d(cf_tfil, cn_vtimec, npt     )
  ierr = putvar1d(ncout,   tim,       npt, 'T')

  zspval = getspval( cf_tfil, cn_vosaline )

  ! Compute spiciness
  DO jt=1,npt
     PRINT *,' TIME = ', jt, tim(jt)/86400.,' days'
     DO jk = 1, npk
        PRINT *, 'Level ', jk
        zmask(:,:) = 1.e0

        ztemp(:,:) = getvar(cf_tfil, cn_votemper, jk, npiglo, npjglo, ktime=jt)
        zsal( :,:) = getvar(cf_tfil, cn_vosaline, jk, npiglo, npjglo, ktime=jt)

        WHERE(zsal == zspval ) zmask = 0.e0
     
        dspi(:,:) = spice ( ztemp, zsal, npiglo, npjglo )
        ierr      = putvar( ncout, id_varout(1), REAL(dspi*zmask), jk, npiglo, npjglo, ktime=jt)

     END DO  ! loop to next level
  END DO  ! next time frame

  ierr = closeout(ncout)

END PROGRAM cdfspice
