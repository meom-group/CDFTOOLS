PROGRAM cdfcoast
  !!======================================================================
  !!                     ***  PROGRAM  cdfcoast  ***
  !!=====================================================================
  !!  ** Purpose : create a coast file from surface tmask
  !!
  !! History :  4.0  : 06/2022  : J.M. Molines  
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2022
  !! $Id$
  !! Copyright (c) 2022, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class mask
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                              :: ierr               ! error status
  INTEGER(KIND=4)                              :: ji, jj, jw         ! dummy loop index
  INTEGER(KIND=4)                              :: narg, iargc, ijarg ! command line 
  INTEGER(KIND=4)                              :: npiglo, npjglo     ! size of the domain
  INTEGER(KIND=4)                              :: isum               ! sum of the 9 tmask point
  INTEGER(KIND=4)                              :: ncout              ! ncid of output file
  INTEGER(KIND=4)                              :: iwidth=1           ! width of coastal band
  INTEGER(KIND=4), DIMENSION(1)                :: ipk, id_varout     ! outptut variables : number of levels,

  INTEGER(KIND=2), DIMENSION(:,:), ALLOCATABLE :: itmask             ! mask at jk level 
  INTEGER(KIND=2), DIMENSION(:,:), ALLOCATABLE :: icoast             ! masked cv_in at jk level

  CHARACTER(LEN=256)                           :: cf_out='coast.nc'  ! output file name
  CHARACTER(LEN=256)                           :: cf_msk             ! input mask file name
  CHARACTER(LEN=256)                           :: cldum              ! dummy char var

  TYPE (variable), DIMENSION(1)                :: stypvar                  ! output attribute

  LOGICAL                                      :: lnc4 =.FALSE.      ! use Netcdf4 chunking and deflation
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfcoast -f MASK-file [-v MASK-var] [-nc4] [-o OUT-file] [-w width ] '
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       This program is used to create a coastal mask file, in wich ocean'
     PRINT *,'       points next to the coast are set to 1, 0 elsewhere. These points'
     PRINT *,'       typically corresponds to runoff points.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f MASK-file  : name of the mask file ' 
     PRINT *,'      '
     PRINT *,'      OPTIONS : '
     PRINT *,'       -v MASK-var : input netcdf mask variable.' 
     PRINT *,'       -w width : width (pixel) of the coastal band'
     PRINT *,'       -nc4 : use netcdf4/Hdf5 chunking and deflation.'
     PRINT *,'       -o OUT-file     : name of coastal_mask file.[',TRIM(cf_out),']'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        none, all are given as arguments.'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       coastal.nc file unless -o option is used.'
     PRINT *,'       Variable : coast_mask'
     STOP 
  ENDIF

  ijarg = 1
  DO WHILE (ijarg <= narg)
    CALL getarg (ijarg, cldum ) ; ijarg = ijarg + 1
    SELECT CASE ( cldum)
    CASE ( '-f'   ) ; CALL getarg(ijarg, cf_msk   )  ; ijarg = ijarg + 1
    CASE ( '-v'   ) ; CALL getarg(ijarg, cn_tmask )  ; ijarg = ijarg + 1
    CASE ( '-w'   ) ; CALL getarg(ijarg, cldum    )  ; ijarg = ijarg + 1 ; READ(cldum,*) iwidth
    CASE ( '-nc4' ) ; lnc4 = .true. 
    CASE ( '-o'   ) ; CALL getarg(ijarg, cf_out   )  ; ijarg = ijarg + 1 
    CASE DEFAULT    ; PRINT *,' ERROR : ', TRIM(cldum),' : unknown option.' ; STOP 99
    END SELECT
  ENDDO

  IF (  chkfile(cf_msk) ) STOP 99 ! missing files

  ! get NEMO grid information
  npiglo = getdim (cf_msk,cn_x)
  npjglo = getdim (cf_msk,cn_y)

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo

  ! Allocate arrays
  ALLOCATE( itmask(npiglo,npjglo) )
  ALLOCATE( icoast(npiglo,npjglo) )
  itmask(:,:) = getvar(cf_msk, cn_tmask, 1, npiglo, npjglo)
  icoast(:,:) = 0
  DO jw = 1, iwidth
  DO jj = 2, npjglo-1
    DO ji = 2, npiglo-1
      IF ( itmask(ji,jj)  == 1 ) THEN
         isum=SUM(itmask(ji-1:ji+1, jj)) +  SUM(itmask(ji, jj-1:jj+1))
         ! check against 6 instead of 5 because central point (i,j) is
         ! counted twice
         IF ( isum /= 6 ) icoast(ji,jj) = icoast(ji,jj) + 1
      ENDIF
    ENDDO
  ENDDO
  WHERE (icoast == 1 ) itmask = 0
  ENDDO

  CALL CreateOutput
  ierr = putvar(ncout, id_varout(1), icoast, 1, npiglo, npjglo)

  ierr =closeout(ncout)

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
    REAL(KIND=8), DIMENSION(1) :: dltim
    !!----------------------------------------------------------------------
    ipk                          = 1

    stypvar(1)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(1)%cname             = 'coast_mask'

    stypvar(1)%cunits            = '1/0'
    stypvar(1)%rmissing_value    = 9999.
    stypvar(1)%valid_min         = 0.
    stypvar(1)%valid_max         = 1.

    stypvar(1)%clong_name        = 'Ocean points near the coast'

    stypvar(1)%cshort_name       = 'coast_mask'

    stypvar(1)%conline_operation = 'N/A'
    stypvar(1)%caxis             = 'TYX'
    stypvar(1)%cprecision        = 'i2'

    ncout = create      (cf_out, cf_msk,  npiglo, npjglo, 0, ld_nc4=lnc4 )
    ierr  = createvar   (ncout, stypvar, 1, ipk, id_varout,  ld_nc4=lnc4 )

    ierr  = putheadervar(ncout,  cf_msk,  npiglo, npjglo, 0)
    dltim = 0.d0
    ierr = putvar1d(ncout, dltim, 1,'T')


  END SUBROUTINE CreateOutput

END PROGRAM cdfcoast
