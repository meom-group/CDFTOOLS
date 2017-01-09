PROGRAM cdffracinv
  !!======================================================================
  !!                     ***  PROGRAM  cdffracinv  ***
  !!=====================================================================
  !!  ** Purpose : Computes fraction of inventory for passive tracers
  !!               output. This is the ratio between inventory at a
  !!               grid point and total inventory
  !!
  !! History : 2.1  : 07/2010  : C.O. Dufour  : Original code
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

  INTEGER(KIND=4)                           :: jt                    ! dummy loop index
  INTEGER(KIND=4)                           :: narg, iargc, ijarg    ! browse line
  INTEGER(KIND=4)                           :: npiglo, npjglo        ! size of the domain
  INTEGER(KIND=4)                           :: npk, npt              ! size of the domain
  INTEGER(KIND=4)                           :: ncout                 ! ncid of output file
  INTEGER(KIND=4)                           :: ierr                  ! error status
  INTEGER(KIND=4), DIMENSION(1)             :: ipk, id_varout        ! level and varid's of output vars

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: trcinvij              ! tracer inventory
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: fracinv               ! fraction of inventory
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: tim                   ! time counter

  CHARACTER(LEN=256)                        :: cf_trc                ! tracer file (for inventory)
  CHARACTER(LEN=256)                        :: cf_out='fracinv.nc'   ! output file name
  CHARACTER(LEN=256)                        :: cv_inv='invcfc'       ! inventory name
  CHARACTER(LEN=256)                        :: cv_out='fracinv'      ! output variable name
  CHARACTER(LEN=256)                        :: cglobal               ! global attribute
  CHARACTER(LEN=256)                        :: cldum                 ! dummy string

  TYPE(variable), DIMENSION(1)              :: stypvar               ! structure for attributes
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdffracinv TRC-file [-inv INV-name]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the fraction of inventory for passive tracers, which is '
     PRINT *,'       the ratio between inventory at a grid point and the total inventory.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       TRC-file : netcdf file with tracer inventory.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       -inv INV-name  : name of the netcdf name for inventory [ ',TRIM(cv_inv),' ]'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       none ... but : horizontal weight to be coded ?' 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : ', TRIM(cv_out)
     STOP
  ENDIF

  ijarg = 1 
  CALL getarg (ijarg, cf_trc) ; ijarg = ijarg + 1

  IF ( chkfile(cf_trc) ) STOP ! missing file

  DO WHILE (ijarg <= narg )
     CALL getarg(ijarg,cldum  ) ; ijarg = ijarg + 1
     SELECT CASE ( cldum )
     CASE ( '-inv' ) ; CALL getarg(ijarg, cv_inv) ; ijarg =ijarg + 1
     CASE DEFAULT ; PRINT *, 'option ', TRIM(cldum),' not understood' ; STOP
     END SELECT
  END DO

  npiglo = getdim (cf_trc,cn_x)
  npjglo = getdim (cf_trc,cn_y)
  npk    = getdim (cf_trc,cn_z)
  npt    = getdim (cf_trc,cn_t)

  ipk(1)                       = 1
  stypvar(1)%cname             = cv_out
  stypvar(1)%cunits            = ''
  stypvar(1)%rmissing_value    = 0.
  stypvar(1)%valid_min         = 0.
  stypvar(1)%valid_max         = 10000.
  stypvar(1)%clong_name        = 'Fraction of inventory'
  stypvar(1)%cshort_name       = cv_out
  stypvar(1)%conline_operation = 'N/A'
  stypvar(1)%caxis             = 'TYX'

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk

  ALLOCATE( trcinvij(npiglo,npjglo), fracinv(npiglo,npjglo) )
  ALLOCATE( tim(npt) )

  WRITE(cglobal,9000) TRIM(cf_trc), TRIM(cv_inv)
9000 FORMAT('cdffracinv ',a,' -inv ',a )

  ncout = create      (cf_out, cf_trc,  npiglo, npjglo, 1                           )
  ierr  = createvar   (ncout,  stypvar, 1,      ipk,    id_varout, cdglobal=cglobal )
  ierr  = putheadervar(ncout,  cf_trc,  npiglo, npjglo, 1                           )

  DO jt=1,npt
     fracinv( :,:) = 0.
     trcinvij(:,:) = getvar(cf_trc, cv_inv, 1, npiglo, npjglo, ktime=jt)
     ! JMM bug ?? : SUM(trcinij) is not the 'total inventory', should be weighted by model metrics ???
     !              also assume spval is 0
     fracinv( :,:) = trcinvij(:,:) / SUM(trcinvij(:,:))
     ierr = putvar(ncout, id_varout(1), fracinv, 1, npiglo, npjglo, ktime=jt)
  END DO

  tim  = getvar1d(cf_trc, cn_vtimec, npt     )
  ierr = putvar1d(ncout,  tim,       npt, 'T')
  ierr = closeout(ncout)

END PROGRAM cdffracinv
