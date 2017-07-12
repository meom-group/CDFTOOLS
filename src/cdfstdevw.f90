PROGRAM cdfstdevw
  !!======================================================================
  !!                     ***  PROGRAM  cdfstdevw  ***
  !!=====================================================================
  !!  ** Purpose : Compute the RMS of W, from the mean squared value.
  !!
  !!  ** Method  : Read gridW and gridW2 and compute rms
  !!
  !! History : 2.1  : 11/2004  : J.M. Molines : Original code
  !!           3.0  : 05/2011  : J.M. Molines : Doctor norm + Lic.
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class statistics
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                            :: jk, jt            ! dummy loop index
  INTEGER(KIND=4)                            :: narg, iargc       ! command line
  INTEGER(KIND=4)                            :: ijarg             ! command line
  INTEGER(KIND=4)                            :: npiglo, npjglo    ! size of the domain
  INTEGER(KIND=4)                            :: npk, npt          ! size of the domain
  INTEGER(KIND=4)                            :: ncout             ! ncid of output variable
  INTEGER(KIND=4)                            :: ierr              ! error status
  INTEGER(KIND=4), DIMENSION(1)              :: ipko, id_varout   ! output variable

  REAL(KIND=4)                               :: rmiss             ! missing value attribute
  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: zvbar, zvba2      ! mean and mean2 variable

  REAL(KIND=8), DIMENSION(:),    ALLOCATABLE :: dtim              ! time counter
  REAL(KIND=8), DIMENSION(:,:),  ALLOCATABLE :: dsdev             ! standard deviation

  CHARACTER(LEN=256)                         :: cf_in             ! input mean file name
  CHARACTER(LEN=256)                         :: cf_in2            ! input mean2 file name
  CHARACTER(LEN=256)                         :: cf_out = 'rmsw.nc'! output file name
  CHARACTER(LEN=256)                         :: cv_in, cv_in2     ! input variable names
  CHARACTER(LEN=256)                         :: cldum             ! dummy character variable
  CHARACTER(LEN=256)                         :: cl_units, cl_longname, cl_shortname

  TYPE(variable), DIMENSION(1)               :: stypvaro          ! output data structure

  LOGICAL                                    :: lchk = .FALSE.    ! flag for missing files
  LOGICAL                                    :: lnc4 = .FALSE.    ! flag for netcdf4 with chunking and deflation
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  cv_in  = cn_vovecrtz
  cv_in2 = TRIM(cn_vovecrtz)//'_sqd'

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfstdevw -w W-file -w2 W2-file  [-o OUT-file] [-nc4] ...'
     PRINT *,'               .... [-var IN-var IN-var2]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the standard deviation of the vertical velocity from its mean'
     PRINT *,'       value and its mean square value. If a variable name is given ,then '
     PRINT *,'       the standard deviation of this variable instead of the vertical '
     PRINT *,'       velocity.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -w W-file  : netcdf file with mean values for w or <IN-var>' 
     PRINT *,'       -w2 W2-file : netcdf file with mean squared values for w or <IN-var>' 
     PRINT *,'      '
     PRINT *,'     OPTIONS: '
     PRINT *,'        [-o OUT-file] : specify the name of the output file instead of ',TRIM(cf_out)
     PRINT *,'        [-nc4 ]   : Use netcdf4 output with chunking and deflation level 1'
     PRINT *,'              This option is effective only if cdftools are compiled with'
     PRINT *,'              a netcdf library supporting chunking and deflation.'
     PRINT *,'        [-var IN-var IN-var2] : give name of mean variable if not ', TRIM(cn_vovecrtz)
     PRINT *,'              and ',TRIM(cn_vovecrtz)//'_sqd'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       none' 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) ,' (if IN-var specified, output file is rms_var.nc)'
     PRINT *,'         variables : ', TRIM(cv_in)//'_rms, (or IN-var_rms)  same unit than the input.'
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'       cdfstd, cdfrmsssh, cdfstdevts, cdfstats.'
     PRINT *,'      '
     STOP 
  ENDIF

  ijarg = 1  
  DO WHILE ( ijarg <= narg) 
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-w'   ) ; CALL getarg(ijarg, cf_in ) ; ijarg=ijarg+1
     CASE ( '-w2'  ) ; CALL getarg(ijarg, cf_in2) ; ijarg=ijarg+1
        ! options
     CASE ( '-o'   ) ; CALL getarg(ijarg, cf_out) ; ijarg=ijarg+1
     CASE ( '-nc4' ) ; lnc4 = .TRUE.
     CASE ( '-var' ) ; CALL getarg(ijarg, cv_in ) ; ijarg=ijarg+1 ; cf_out='rms_'//TRIM(cv_in)//'.nc'
       ;             ; CALL getarg(ijarg, cv_in2) ; ijarg=ijarg+1
     CASE DEFAULT    ; PRINT *,' ERROR : ',TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  ENDDO

  ! check existence of files
  lchk = lchk .OR. chkfile(cf_in )
  lchk = lchk .OR. chkfile(cf_in2)
  IF (lchk ) STOP 99 ! missing file

  npiglo = getdim (cf_in, cn_x)
  npjglo = getdim (cf_in, cn_y)
  npk    = getdim (cf_in, cn_z)
  npt    = getdim (cf_in, cn_t)

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  CALL CreateOutput

  ALLOCATE( zvbar(npiglo,npjglo), zvba2(npiglo,npjglo) )
  ALLOCATE( dsdev(npiglo,npjglo), dtim(npt)            )

  ierr = getvaratt(cf_in, cv_in, cl_units, rmiss, cl_longname, cl_shortname )

  DO jt = 1, npt
     DO jk = 1, npk
        zvbar(:,:) = getvar(cf_in,  cv_in,  jk, npiglo, npjglo, ktime=jt)
        zvba2(:,:) = getvar(cf_in2, cv_in2, jk, npiglo, npjglo, ktime=jt)

        dsdev(:,:) = SQRT ( DBLE(zvba2(:,:) - zvbar(:,:)*zvbar(:,:)) )

        ierr = putvar(ncout, id_varout(1), REAL(dsdev), jk, npiglo, npjglo, ktime=jt)
     END DO
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
    ipko(1) = npk
    stypvaro(1)%ichunk            = (/ npiglo, MAX(1,npjglo/30), 1,1 /)
    stypvaro(1)%cname             = TRIM(cv_in)//'_rms'
    stypvaro(1)%cunits            = TRIM(cl_units)
    stypvaro(1)%rmissing_value    = 0.
    stypvaro(1)%valid_min         = 0.
    stypvaro(1)%valid_max         = 10.
    stypvaro(1)%clong_name        = 'RMS_'//TRIM(cl_longname)
    stypvaro(1)%cshort_name       = TRIM(cv_in)//'_rms'
    stypvaro(1)%conline_operation = 'N/A'
    stypvaro(1)%caxis             = 'TZYX'

    ncout = create      (cf_out, cf_in,    npiglo, npjglo, npk       , ld_nc4=lnc4 )
    ierr  = createvar   (ncout,  stypvaro, 1,      ipko,   id_varout , ld_nc4=lnc4 )
    ierr  = putheadervar(ncout,  cf_in,    npiglo, npjglo, npk       )

    dtim = getvar1d(cf_in, cn_vtimec, npt     )
    ierr = putvar1d(ncout, dtim,      npt, 'T')

  END SUBROUTINE CreateOutput

END PROGRAM cdfstdevw
