PROGRAM cdficb_clim
  !!======================================================================
  !!                     ***  PROGRAM  cdficb_clim  ***
  !!=====================================================================
  !!  ** Purpose : Concatenate ceberg mass and melt from 12 monthly files
  !!               into a single 12 time frame file.
  !!  ** Method  : read and write ...
  !!
  !! History :  3.0  : 06/2016  : N. Merino         : Original code
  !!         :  4.0  : 03/2017  : J.-M. Molines
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017
  !! $Id$
  !! Copyright (c) 2012, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class iceberg_processing
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                            :: jf, jt               ! dummy loop index
  INTEGER(KIND=4)                            :: ierr                 ! working integer
  INTEGER(KIND=4)                            :: nfiles               ! number of files
  INTEGER(KIND=4)                            :: narg, iargc, ijarg   ! command line 
  INTEGER(KIND=4)                            :: npiglo, npjglo, npt  ! size of the domain
  INTEGER(KIND=4)                            :: nvpk                 ! vertical levels in working variable
  INTEGER(KIND=4)                            :: nboutput=2           ! number of values to write in cdf output
  INTEGER(KIND=4)                            :: ncout                ! for netcdf output
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: ipk, id_varout

  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: ricbmass, ricbmelt   ! icbmass icbmelt

  REAL(KIND=8), DIMENSION(:),    ALLOCATABLE :: dtim                 ! time counter

  TYPE(variable), DIMENSION(:),  ALLOCATABLE :: stypvar              ! structure of output
  !
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cf_lst            ! input icb file
  CHARACTER(LEN=256)                         :: cf_out='icbdiags.nc' ! output file
  CHARACTER(LEN=256)                         :: cldum                ! dummy string
  !
  LOGICAL                                    :: lchk  = .FALSE.      ! missing file flag
  LOGICAL                                    :: lnc4  = .FALSE.      ! Use nc4 with chunking and deflation
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdficb_clim -l LST-ICB-monthly-files [-o OUT-file] [-nc4]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'        Concatenates 12 monthly input files, into a 12 frames output file.'
     PRINT *,'        This is done for the 2 variables corresponding to mass and melt.'
     PRINT *,'        No process done in this tool.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -l  LST-ICB-monthly-files : A list of 12 monthly-mean ICB files.'
     PRINT *,'            These files are likely produced by cdfmoy.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-o OUT-file] : specify output file name instead of ',TRIM(cf_out)
     PRINT *,'       [-nc4 ] : Use netcdf4 output with chunking and deflation level 1.'
     PRINT *,'            This option is effective only if cdftools are compiled with'
     PRINT *,'            a netcdf library supporting chunking and deflation.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'         none '
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : Mass  (Kg/m2 )'
     PRINT *,'                     Melt  (Kg/m2/s )'
     STOP 
  ENDIF

  ijarg = 1
  DO WHILE (ijarg <= narg )
     CALL getarg(ijarg, cldum) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-l'  ) ; CALL GetFileList
        ! options
     CASE ( '-o'  ) ; CALL getarg(ijarg, cf_out) ; ijarg=ijarg+1
     CASE ( '-nc4') ; lnc4 = .TRUE.
     CASE DEFAULT   ; PRINT *,' ERROR : ', TRIM(cldum),' unknown option.' ; STOP 99
     END SELECT
  ENDDO

  IF ( nfiles /= 12 ) THEN 
     PRINT *,' +++ ERROR : This program needs 12 monthly files in input.' ; STOP 99
  ENDIF

  DO jf=1, nfiles
     lchk = lchk .OR. chkfile(cf_lst(jf))
  ENDDO

  IF ( lchk ) STOP 99 ! missing file

  npiglo = getdim (cf_lst(1),cn_x)
  npjglo = getdim (cf_lst(1),cn_y)
  npt    = 12

  ALLOCATE ( ricbmass(npiglo,npjglo) )
  ALLOCATE ( ricbmelt(npiglo,npjglo) )
  ALLOCATE ( dtim(npt) )
  ALLOCATE ( stypvar(nboutput), ipk(nboutput), id_varout(nboutput) )

  CALL CreateOutput

  ! Check variable on first file ?
  IF (chkvar(cf_lst(1), cn_iicbmass)) THEN
     cn_iicbmass='missing'
     PRINT *,'' 
     PRINT *,' WARNING, ICEBERG MASS IS SET TO 0. '
     PRINT *,' '
  END IF
  !
  DO jt = 1, npt
     IF (TRIM(cn_iicbmass) /= 'missing') THEN ; ricbmass(:,:) = getvar(cf_lst(jt), cn_iicbmass, 1, npiglo, npjglo, ktime=1)
     ELSE                                     ; ricbmass(:,:) = 0.
     ENDIF
     ricbmelt(:,:) = getvar(cf_lst(jt), cn_iicbmelt, 1, npiglo, npjglo, ktime=1)

     ! netcdf output
     ierr = putvar(ncout,id_varout(1),REAL(ricbmass(:,:)),1,npiglo,npjglo, ktime=jt)
     ierr = putvar(ncout,id_varout(2),REAL(ricbmelt(:,:)),1,npiglo,npjglo, ktime=jt)
  END DO ! time loop
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

    ipk(:) = 1
    ! define new variables for output 
    stypvar%scale_factor      = 1.
    stypvar%add_offset        = 0.
    stypvar%savelog10         = 0.
    stypvar%conline_operation = 'N/A'
    stypvar%caxis             = 'T'

    stypvar(1)%ichunk         = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(1)%cname          = 'Mass'
    stypvar(1)%cunits         = 'Kg/m2'
    stypvar(1)%clong_name     = 'Icb mass per unit of area'
    stypvar(1)%cshort_name    = 'Mass'

    stypvar(2)%ichunk         = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(2)%cname          = 'Melt'
    stypvar(2)%cunits         = 'Kg/m2/s'
    stypvar(2)%clong_name     = 'Icb melt flux'
    stypvar(2)%cshort_name    = 'Melt'


    ! create output fileset
    ncout = create      (cf_out, 'none',  npiglo,      npjglo, 1,  cdep='depthw', ld_nc4=lnc4 )
    ierr  = createvar   (ncout,  stypvar, nboutput, ipk, id_varout              , ld_nc4=lnc4 )
    ierr  = putheadervar(ncout,  cf_lst(1), npiglo,      npjglo, 1)

    dtim = (/(jt,jt=1,12)/)
    ierr = putvar1d(ncout, dtim, npt, 'T')


  END SUBROUTINE CreateOutput

  SUBROUTINE GetFileList
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE GetFileList  ***
    !!
    !! ** Purpose :  Set up a file list given on the command line as 
    !!               blank separated list
    !!
    !! ** Method  :  Scan the command line until a '-' is found
    !!----------------------------------------------------------------------
    INTEGER (KIND=4)  :: ji
    INTEGER (KIND=4)  :: icur
    !!----------------------------------------------------------------------
    !!
    nfiles=0
    ! need to read a list of file ( number unknow ) 
    ! loop on argument till a '-' is found as first char
    icur=ijarg                          ! save current position of argument number
    DO ji = icur, narg                  ! scan arguments till - found
       CALL getarg ( ji, cldum )
       IF ( cldum(1:1) /= '-' ) THEN ; nfiles = nfiles+1
       ELSE                          ; EXIT
       ENDIF
    ENDDO
    ALLOCATE (cf_lst(nfiles) )
    DO ji = icur, icur + nfiles -1
       CALL getarg(ji, cf_lst( ji -icur +1 ) ) ; ijarg=ijarg+1
    END DO
  END SUBROUTINE GetFileList

END PROGRAM cdficb_clim
