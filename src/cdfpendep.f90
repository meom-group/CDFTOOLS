PROGRAM cdfpendep
  !!======================================================================
  !!                     ***  PROGRAM  cdfpendep  ***
  !!=====================================================================
  !!  ** Purpose : Computes penetration depth for passive tracer output. 
  !!               This is the ratio between inventory and surface 
  !!               concentration.
  !!
  !!  ** Method  : takes TRC files as input
  !!
  !! History : 2.1  : 02/2008  : J.M. Molines : Original code
  !!         : 2.1  : 09/2010  : C. Dufour    : Adapation to TOP evolution
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

  INTEGER(KIND=4)                           :: jt                   ! dummy loop index
  INTEGER(KIND=4)                           :: narg, iargc, ijarg   ! browse command line
  INTEGER(KIND=4)                           :: npiglo, npjglo       ! size of the domain
  INTEGER(KIND=4)                           :: npk, npt             ! size of the domain
  INTEGER(KIND=4)                           :: ncout                ! ncid of output file
  INTEGER(KIND=4)                           :: ierr                 ! error status
  INTEGER(KIND=4), DIMENSION(1)             :: ipk, id_varout       ! levels and varid's of output vats

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: trcinv, trcsurf      ! inventory, surface concentration
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: rpendep              ! penetration depth
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: tim                  ! time counter

  CHARACTER(LEN=256)                        :: cf_trcfil            ! tracer file name
  CHARACTER(LEN=256)                        :: cf_inv               ! inventory file name
  CHARACTER(LEN=256)                        :: cf_out='pendep.nc'   ! output file
  CHARACTER(LEN=256)                        :: cv_inv               ! inventory variable name
  CHARACTER(LEN=256)                        :: cv_trc               ! tracer variable name
  CHARACTER(LEN=256)                        :: cglobal              ! global attribute
  CHARACTER(LEN=256)                        :: cldum                ! dummy string

  TYPE(variable), DIMENSION(1)              :: typvar               ! structure for attributes
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  cv_inv = cn_invcfc
  cv_trc = cn_cfc11
  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage :  cdfpendep TRC-file INV-file  ... '
     PRINT *,'                    ... [-inv inventory_name -trc trc_name ]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'        Compute the penetration depth for passive tracers. It is the'
     PRINT *,'        ratio between the inventory and the surface concentration of'
     PRINT *,'        the tracer.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       TRC-file : netcdf file with tracer concentration.' 
     PRINT *,'       INV-file : netcdf file with inventory of the tracer.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-inv inventory_name ] : specify netcdf variable name for inventory.' 
     PRINT *,'                                Default is ', TRIM(cv_inv)
     PRINT *,'       [-trc tracer_name ]    : specify netcdf variable name for tracer.' 
     PRINT *,'                                Default is ', TRIM(cv_trc)
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        none'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : pendep (m)'
     STOP
  ENDIF

  ijarg = 1
  CALL getarg (ijarg, cf_trcfil) ; ijarg = ijarg + 1 
  CALL getarg (ijarg, cf_inv   ) ; ijarg = ijarg + 1

  IF ( chkfile(cf_trcfil) .OR. chkfile(cf_inv) ) STOP 99 ! missing file

  DO WHILE ( ijarg <= narg)
     CALL getarg(ijarg, cldum ) ; ijarg = ijarg + 1 
     SELECT CASE ( cldum )
     CASE ('-inv') ;  CALL getarg(ijarg, cv_inv) ; ijarg=ijarg+1
     CASE ('-trc') ;  CALL getarg(ijarg, cv_trc) ; ijarg=ijarg+1
     CASE DEFAULT  ; PRINT *, 'option ', TRIM(cldum),' not understood' ; STOP 99
     END SELECT
  END DO

  npiglo = getdim (cf_trcfil,cn_x)
  npjglo = getdim (cf_trcfil,cn_y)
  npk    = getdim (cf_trcfil,cn_z)
  npt    = getdim (cf_trcfil,cn_t)

  ipk(1)                      = 1
  typvar(1)%cname             = cn_pendep
  typvar(1)%cunits            = 'm'
  typvar(1)%rmissing_value    = 0.
  typvar(1)%valid_min         = 0.
  typvar(1)%valid_max         = 10000.
  typvar(1)%clong_name        = 'Penetration depth'
  typvar(1)%cshort_name       = cn_pendep
  typvar(1)%conline_operation = 'N/A'
  typvar(1)%caxis             = 'TYX'


  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  ALLOCATE( trcinv(npiglo,npjglo), trcsurf(npiglo,npjglo), rpendep(npiglo,npjglo) )
  ALLOCATE( tim(npt) )

  WRITE(cglobal,9000) TRIM(cf_trcfil), TRIM(cf_inv), TRIM(cv_inv), TRIM(cv_trc)
9000 FORMAT('cdfpendep ',a,' ',a,' -inv ',a,' -trc ',a )

  ncout = create      (cf_out, cf_trcfil, npiglo, npjglo, 1)
  ierr  = createvar   (ncout,  typvar,    1,      ipk,    id_varout, cdglobal=cglobal )
  ierr  = putheadervar(ncout,  cf_trcfil, npiglo, npjglo, 1)

  DO jt = 1,npt
     rpendep(:,:) = 0.
     trcinv( :,:) = getvar(cf_inv,    cv_inv, 1, npiglo, npjglo, ktime=jt)
     trcsurf(:,:) = getvar(cf_trcfil, cv_trc, 1, npiglo, npjglo, ktime=jt)

     WHERE( trcsurf /= 0 ) rpendep = trcinv/trcsurf
     ierr=putvar(ncout, id_varout(1), rpendep, 1, npiglo, npjglo, ktime=jt)
  END DO

  tim  = getvar1d(cf_trcfil, cn_vtimec, npt     )
  ierr = putvar1d(ncout,     tim,       npt, 'T')
  ierr = closeout(ncout)

END PROGRAM cdfpendep
