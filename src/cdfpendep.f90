PROGRAM cdfpendep
  !!======================================================================
  !!                     ***  PROGRAM  cdfpendep  ***
  !!=====================================================================
  !!  ** Purpose : Compute penetration depth for passive tracer output. 
  !!               This is the ratio between inventory and surface 
  !!               concentration.
  !!
  !!  ** Method  : takes TRC files as input
  !!
  !! History : 2.1  : 02/2008  : J.M. Molines : Original code
  !!         : 2.1  : 09/2010  : C. Dufour    : Adapation to TOP evolution
  !!           3.0  : 01/2011  : J.M. Molines : Doctor norm + Lic.
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  USE cdfio 
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class passive_tracer
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

  REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: dtim                 ! time counter

  CHARACTER(LEN=256)                        :: cf_trcfil            ! tracer file name
  CHARACTER(LEN=256)                        :: cf_inv               ! inventory file name
  CHARACTER(LEN=256)                        :: cf_out='pendep.nc'   ! output file
  CHARACTER(LEN=256)                        :: cv_inv               ! inventory variable name
  CHARACTER(LEN=256)                        :: cv_trc               ! tracer variable name
  CHARACTER(LEN=256)                        :: cglobal              ! global attribute
  CHARACTER(LEN=256)                        :: cldum                ! dummy string

  TYPE(variable), DIMENSION(1)              :: stypvar              ! structure for attributes

  LOGICAL                                   :: lnc4 = .FALSE.       ! Use nc4 with chunking and deflation
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  cv_inv = cn_invcfc
  cv_trc = cn_cfc11
  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage :  cdfpendep -trc TRC-file -i INV-file [-o OUT-file] [-nc4] ...'
     PRINT *,'              ... [-vinv inventory_name -vtrc trc_name ]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'        Compute the penetration depth for passive tracers. It is the ratio'
     PRINT *,'        between the inventory and the surface concentration of the tracer.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -trc TRC-file : netcdf file with tracer concentration.' 
     PRINT *,'       -i INV-file   : netcdf file with inventory of the tracer.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-vinv inventory_name ] : specify netcdf variable name for inventory.' 
     PRINT *,'                                Default is ', TRIM(cv_inv)
     PRINT *,'       [-vtrc trc_name   ]    : specify netcdf variable name for tracer.' 
     PRINT *,'                                Default is ', TRIM(cv_trc)
     PRINT *,'       [-o OUT-file ] : specify output file name instead of ', TRIM(cf_out)
     PRINT *,'       [ -nc4 ]     : Use netcdf4 output with chunking and deflation level 1.'
     PRINT *,'                 This option is effective only if cdftools are compiled with'
     PRINT *,'                 a netcdf library supporting chunking and deflation.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        none'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : pendep (m)'
     PRINT *,'      '
     STOP 
  ENDIF

  ijarg = 1
  DO WHILE ( ijarg <= narg)
     CALL getarg(ijarg, cldum ) ; ijarg = ijarg + 1 
     SELECT CASE ( cldum )
     CASE ('-trc' ) ;  CALL getarg(ijarg, cf_trcfil) ; ijarg=ijarg+1
     CASE ('-i'   ) ;  CALL getarg(ijarg, cf_inv   ) ; ijarg=ijarg+1
        ! options
     CASE ('-o'   ) ;  CALL getarg(ijarg, cf_out   ) ; ijarg=ijarg+1
     CASE ('-vinv') ;  CALL getarg(ijarg, cv_inv   ) ; ijarg=ijarg+1
     CASE ('-vtrc') ;  CALL getarg(ijarg, cv_trc   ) ; ijarg=ijarg+1
     CASE ('-nc4' ) ;  lnc4 = .TRUE.
     CASE DEFAULT  ; PRINT *, ' ERROR : ', TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  END DO

  IF ( chkfile(cf_trcfil) .OR. chkfile(cf_inv) ) STOP 99 ! missing file

  npiglo = getdim (cf_trcfil,cn_x)
  npjglo = getdim (cf_trcfil,cn_y)
  npk    = getdim (cf_trcfil,cn_z)
  npt    = getdim (cf_trcfil,cn_t)

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  ALLOCATE( trcinv(npiglo,npjglo), trcsurf(npiglo,npjglo), rpendep(npiglo,npjglo) )
  ALLOCATE( dtim(npt) )

  WRITE(cglobal,9000) TRIM(cf_trcfil), TRIM(cf_inv), TRIM(cv_inv), TRIM(cv_trc)
9000 FORMAT('cdfpendep ',a,' ',a,' -inv ',a,' -trc ',a )

  CALL CreateOutput

  DO jt = 1,npt
     rpendep(:,:) = 0.
     trcinv( :,:) = getvar(cf_inv,    cv_inv, 1, npiglo, npjglo, ktime=jt)
     trcsurf(:,:) = getvar(cf_trcfil, cv_trc, 1, npiglo, npjglo, ktime=jt)

     WHERE( trcsurf /= 0 ) rpendep = trcinv/trcsurf
     ierr=putvar(ncout, id_varout(1), rpendep, 1, npiglo, npjglo, ktime=jt)
  END DO

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
    ipk(1)                       = 1
    stypvar(1)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(1)%cname             = cn_pendep
    stypvar(1)%cunits            = 'm'
    stypvar(1)%rmissing_value    = 0.
    stypvar(1)%valid_min         = 0.
    stypvar(1)%valid_max         = 10000.
    stypvar(1)%clong_name        = 'Penetration depth'
    stypvar(1)%cshort_name       = cn_pendep
    stypvar(1)%conline_operation = 'N/A'
    stypvar(1)%caxis             = 'TYX'

    ncout = create      (cf_out, cf_trcfil, npiglo, npjglo, 1                          , ld_nc4=lnc4 )
    ierr  = createvar   (ncout,  stypvar,   1,      ipk,    id_varout, cdglobal=cglobal, ld_nc4=lnc4 )
    ierr  = putheadervar(ncout,  cf_trcfil, npiglo, npjglo, 1                                        )

    dtim = getvar1d(cf_trcfil, cn_vtimec, npt     )
    ierr = putvar1d(ncout,    dtim,       npt, 'T')

  END SUBROUTINE CreateOutput

END PROGRAM cdfpendep
