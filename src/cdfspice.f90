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

  INTEGER(KIND=4)                           :: ji, jj, jk, jt     ! dummy loop index
  INTEGER(KIND=4)                           :: ierr               ! error status
  INTEGER(KIND=4)                           :: narg, ijarg, iargc ! browse command line
  INTEGER(KIND=4)                           :: npiglo, npjglo     ! size of the domain
  INTEGER(KIND=4)                           :: npk, npkk, npt     ! size of the domain
  INTEGER(KIND=4)                           :: ncout              ! ncid of output file
  INTEGER(KIND=4), DIMENSION(1)             :: ipk, id_varout     ! level and  varid's

  REAL(KIND=4)                              :: zspval             ! missing value
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ztemp              ! temperature
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zsal               ! salinity
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zmask              ! 2D mask at current level

  REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: dtim               ! time counter
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dtempt             ! temperature
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dsalt              ! salinity
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dsalref            ! reference salinity
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dspi               ! spiceness

  CHARACTER(LEN=256)                        :: cf_tfil            ! input filename
  CHARACTER(LEN=256)                        :: cf_out='spice.nc'  ! output file name
  CHARACTER(LEN=256)                        :: cldum              ! dummy characte variable variable
  CHARACTER(LEN=256)                        :: cv_sal             ! salinity name in netcdf
  CHARACTER(LEN=256)                        :: cv_tem             ! temperature name in netcdf

  TYPE (variable), DIMENSION(1)             :: stypvar            ! structure for attributes
  LOGICAL                                   :: lnc4 = .FALSE.     ! flag for missing files

  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfspice -t T-file [-sal SAL-name] [-tem TEM-name] ...'
     PRINT *,'        ... [-o OUT-file] [-nc4]'
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
     PRINT *,'       -t T-file : netcdf file with temperature and salinity (gridT)' 
     PRINT *,'     '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-sal SAL-name]  : name of salinity variable'
     PRINT *,'       [-tem TEM-name]  : name of temperature variable'
     PRINT *,'       [-o OUT-file]    : specify output filename instead of ',TRIM(cf_out)
     PRINT *,'       [-nc4]  : enable chunking and compression'
     PRINT *,'      '
     PRINT *,'     OPENMP SUPPORT : yes'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       none'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) ,' unless option -o is used.'
     PRINT *,'         variables : vospice'
     PRINT *,'      '
     PRINT *,'     REFERENCE :'
     PRINT *,'       Flament (2002) "A state variable for characterizing '
     PRINT *,'             water masses and their diffusive stability: spiciness."'
     PRINT *,'             Progress in Oceanography Volume 54, 2002, Pages 493-501.'
     PRINT *,'      '
     STOP 
  ENDIF

  ! default name for salinity and temperature
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
     CASE ( '-nc4' ) ; lnc4 = .TRUE.
     CASE DEFAULT    ; PRINT *,'ERROR : ',TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  ENDDO

  IF ( chkfile(cf_tfil) ) STOP 99 ! missing files

  npiglo = getdim (cf_tfil,cn_x)
  npjglo = getdim (cf_tfil,cn_y)
  npk    = getdim (cf_tfil,cn_z)
  npt    = getdim (cf_tfil,cn_t)

  npkk=npk
  IF ( npk == 0 ) npkk=1

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  ALLOCATE (ztemp(npiglo,npjglo), zsal (npiglo,npjglo) )
  ALLOCATE (dspi( npiglo,npjglo), zmask(npiglo,npjglo) )
  ALLOCATE (dtempt(npiglo,npjglo), dsalt(npiglo,npjglo))
  ALLOCATE (dsalref(npiglo,npjglo))
  ALLOCATE (dtim(npt))

  CALL CreateOutput
  zspval = getspval( cf_tfil, cn_vosaline )

  ! Compute spiciness
  DO jt=1,npt
     PRINT *,' TIME = ', jt, dtim(jt)/86400.,' days'
     DO jk = 1, npkk
        PRINT *, 'Level ', jk
        zmask(:,:) = 1.e0

        ztemp(:,:) = getvar(cf_tfil, cv_tem, jk, npiglo, npjglo, ktime=jt)
        zsal( :,:) = getvar(cf_tfil, cv_sal, jk, npiglo, npjglo, ktime=jt)

        WHERE(zsal == zspval ) zmask = 0.e0

        dspi(:,:) = spice ( ztemp, zsal, npiglo, npjglo )
        ierr      = putvar( ncout, id_varout(1), REAL(dspi*zmask), jk, npiglo, npjglo, ktime=jt)

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
    ipk(:)                       = npkk
    stypvar(1)%ichunk            = (/npiglo, MAX(1,npjglo/30), 1, 1 /)
    stypvar(1)%cname             = 'vospice'
    stypvar(1)%cunits            = 'kg/m3'
    stypvar(1)%rmissing_value    = 0.
    stypvar(1)%valid_min         = -300.
    stypvar(1)%valid_max         = 300.
    stypvar(1)%clong_name        = 'spiciness'
    stypvar(1)%cshort_name       = 'vospice'
    stypvar(1)%conline_operation = 'N/A'
    stypvar(1)%caxis             = 'TZYX'

    ! create output fileset
    ncout = create      (cf_out, cf_tfil, npiglo, npjglo, npk,       ld_nc4=lnc4     )
    ierr  = createvar   (ncout,  stypvar, 1,      ipk,    id_varout, ld_nc4=lnc4     )
    ierr  = putheadervar(ncout,  cf_tfil, npiglo, npjglo, npk       )

    dtim = getvar1d(cf_tfil, cn_vtimec, npt     )
    ierr = putvar1d(ncout,   dtim,      npt, 'T')

  END SUBROUTINE CreateOutput


END PROGRAM cdfspice
