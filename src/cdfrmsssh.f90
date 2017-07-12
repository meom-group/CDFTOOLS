PROGRAM cdfrmsssh
  !!======================================================================
  !!                     ***  PROGRAM  cdfrmsssh  ***
  !!=====================================================================
  !!  ** Purpose : Compute the RMS of SSH, from the mean squared value.
  !!
  !!  ** Method  : Read gridT and gridT2 and compute rms
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
  INTEGER(KIND=4)                            :: ijarg, ixtra      ! command line
  INTEGER(KIND=4)                            :: npiglo, npjglo    ! size of the domain
  INTEGER(KIND=4)                            :: npk, npt          ! size of the domain
  INTEGER(KIND=4)                            :: ncout             ! ncid of output variable
  INTEGER(KIND=4)                            :: ierr              ! error status
  INTEGER(KIND=4), DIMENSION(1)              :: ipko, id_varout   ! output variable

  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: zvbar, zvba2      ! mean and mean2 variable

  REAL(KIND=8), DIMENSION(:),    ALLOCATABLE :: dtim              ! time counter
  REAL(KIND=8), DIMENSION(:,:),  ALLOCATABLE :: dsdev             ! standard deviation

  CHARACTER(LEN=256)                         :: cf_in             ! input mean file name
  CHARACTER(LEN=256)                         :: cf_in2            ! input mean2 file name
  CHARACTER(LEN=256)                         :: cf_out = 'rms.nc' ! output file name
  CHARACTER(LEN=256)                         :: cv_in, cv_in2     ! input variable names
  CHARACTER(LEN=256)                         :: cldum             ! dummy character variable
  CHARACTER(LEN=256)                         :: cl_sossheig       ! local name for ssh
  CHARACTER(LEN=256)                         :: cl_sossheig2      ! local name for ssh2

  TYPE(variable), DIMENSION(1)               :: stypvaro          ! output data structure

  LOGICAL                                    :: lchk = .FALSE.    ! flag for missing files
  LOGICAL                                    :: lnc4 = .FALSE.    ! flag for missing files
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfrmsssh -t T-file -t2 T2-file [-o OUT-file] [-nc4] ... '
     PRINT *,'            ...  [-var VAR-ssh VAR-ssh2]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the standard deviation of the SSH from its mean value'
     PRINT *,'       its mean square value. '
     PRINT *,'      '
     PRINT *,'       Note that what is computed in this program is stictly the'
     PRINT *,'       standard deviation. It is very often called RMS, which is'
     PRINT *,'       an abuse. It is the same only in the case of zero mean value.'
     PRINT *,'       However, for historical reason, the name of this tool, remains'
     PRINT *,'       unchanged: cdfrmsssh'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -t T-file   : netcdf file with mean values for SSH' 
     PRINT *,'       -t2 T2-file : netcdf file with mean squared values for SSH' 
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-o OUT-file] : specify the name of the output file instead'
     PRINT *,'            of default name ', TRIM(cf_out)
     PRINT *,'       [-nc4] : use netcdf4 with chunking and deflation '
     PRINT *,'       [-var VAR-ssh VAR-ssh2] : specify the variable name for mean and '
     PRINT *,'             mean-squared ssh, if they differ from the standard names.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       none' 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) ,' unless option -o is used.'
     PRINT *,'         variables : ', TRIM(cn_sossheig)//'_rms, same unit than the input.'
     PRINT *,'      '
     PRINT *,'     SEA ALSO :'
     PRINT *,'       cdfstd, cdfstdevw, cdfstdevts.'
     STOP 
  ENDIF

  cl_sossheig  = cn_sossheig
  cl_sossheig2 = TRIM(cn_sossheig)//'_sqd'

  ijarg = 1  
  DO WHILE ( ijarg <= narg) 
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-t'   ) ; CALL getarg(ijarg, cf_in        ) ; ijarg=ijarg+1
     CASE ( '-t2'  ) ; CALL getarg(ijarg, cf_in2       ) ; ijarg=ijarg+1
        ! options
     CASE ( '-nc4' ) ; lnc4 = .TRUE.
     CASE ( '-o'   ) ; CALL getarg(ijarg, cf_out       ) ; ijarg=ijarg+1
     CASE ( '-var' ) ; CALL getarg(ijarg, cl_sossheig  ) ; ijarg=ijarg+1
       ;             ; CALL getarg(ijarg, cl_sossheig2 ) ; ijarg=ijarg+1
     CASE DEFAULT    ; PRINT *,' ERROR : ', TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  ENDDO

  cv_in  = cl_sossheig
  cv_in2 = cl_sossheig2

  ! check existence of files
  lchk = lchk .OR. chkfile(cf_in  )
  lchk = lchk .OR. chkfile(cf_in2 )
  IF (lchk ) STOP 99 ! missing file

  npiglo = getdim (cf_in, cn_x)
  npjglo = getdim (cf_in, cn_y)
  npk    = getdim (cf_in, cn_z)
  npt    = getdim (cf_in, cn_t)

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  ALLOCATE( zvbar(npiglo,npjglo), zvba2(npiglo,npjglo) )
  ALLOCATE( dsdev(npiglo,npjglo), dtim(npt)            )

  CALL CreateOutput

  DO jt = 1, npt
     zvbar(:,:) = getvar(cf_in,  cv_in,  1, npiglo, npjglo, ktime=jt)
     zvba2(:,:) = getvar(cf_in2, cv_in2, 1, npiglo, npjglo, ktime=jt)

     dsdev(:,:) = SQRT ( DBLE(zvba2(:,:) - zvbar(:,:)*zvbar(:,:)) )

     ierr = putvar(ncout, id_varout(1), REAL(dsdev), 1, npiglo, npjglo, ktime=jt)
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
    ipko(1) = 1
    stypvaro(1)%ichunk            = (/npiglo, MAX(1,npjglo/30), 1, 1 /)
    stypvaro(1)%cname             = TRIM(cv_in)//'_rms'
    stypvaro(1)%cunits            = 'm'
    stypvaro(1)%rmissing_value    = 0.
    stypvaro(1)%valid_min         = 0.
    stypvaro(1)%valid_max         = 100.
    stypvaro(1)%clong_name        = 'RMS_Sea_Surface_height'
    stypvaro(1)%cshort_name       = TRIM(cv_in)//'_rms'
    stypvaro(1)%conline_operation = 'N/A'
    stypvaro(1)%caxis             = 'TYX'

    ncout = create      (cf_out, cf_in,    npiglo, npjglo, npk       , ld_nc4=lnc4 )
    ierr  = createvar   (ncout,  stypvaro, 1,      ipko,   id_varout , ld_nc4=lnc4 )
    ierr  = putheadervar(ncout,  cf_in,    npiglo, npjglo, npk                     )

    dtim = getvar1d(cf_in, cn_vtimec, npt     )
    ierr = putvar1d(ncout, dtim,      npt, 'T')

  END SUBROUTINE CreateOutput

END PROGRAM cdfrmsssh
