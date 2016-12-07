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

  INTEGER(KIND=4)                            :: jk, jt            ! dummy loop index
  INTEGER(KIND=4)                            :: narg, iargc       ! command line
  INTEGER(KIND=4)                            :: ijarg, ixtra      ! command line
  INTEGER(KIND=4)                            :: npiglo, npjglo    ! size of the domain
  INTEGER(KIND=4)                            :: npk, npt          ! size of the domain
  INTEGER(KIND=4)                            :: ncout             ! ncid of output variable
  INTEGER(KIND=4)                            :: ierr              ! error status
  INTEGER(KIND=4), DIMENSION(1)              :: ipko, id_varout   ! output variable

  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: zvbar, zvba2      ! mean and mean2 variable
  REAL(KIND=4), DIMENSION(:),    ALLOCATABLE :: tim               ! time counter

  REAL(KIND=8), DIMENSION(:,:),  ALLOCATABLE :: dsdev             ! standard deviation

  CHARACTER(LEN=256)                         :: cf_in             ! input mean file name
  CHARACTER(LEN=256)                         :: cf_in2            ! input mean2 file name
  CHARACTER(LEN=256)                         :: cf_out = 'rms.nc' ! output file name
  CHARACTER(LEN=256)                         :: cv_in, cv_in2     ! input variable names
  CHARACTER(LEN=256)                         :: cldum             ! dummy character variable

  TYPE(variable), DIMENSION(1)               :: stypvaro          ! output data structure

  LOGICAL                                    :: lchk = .FALSE.    ! flag for missing files
  LOGICAL                                    :: lnc4 = .FALSE.    ! flag for missing files
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  cv_in = cn_sossheig

  narg= iargc()
  IF ( narg < 2 ) THEN
     PRINT *,' usage : cdfrmsssh T-file T2-file [-nc4] [-o outputfile]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the standard deviation of the SSH from its'
     PRINT *,'       mean value and its mean square value. '
     PRINT *,'      '
     PRINT *,'       Note that what is computed in this program is stictly the'
     PRINT *,'       standard deviation. It is very often called RMS, which is'
     PRINT *,'       an abuse. It is the same only in the case of zero mean value.'
     PRINT *,'       However, for historical reason, the name of this tool, remains'
     PRINT *,'       unchanged: cdfrmsssh'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       T-file  : netcdf file with mean values for SSH' 
     PRINT *,'       T2-file : netcdf file with mean squared values for SSH' 
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-nc4] : use netcdf4 with chunking and deflation '
     PRINT *,'       [-o output file ] : specify the name of the output file instead'
     PRINT *,'                          of default name ', TRIM(cf_out)
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       none' 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : ', TRIM(cv_in)//'_rms, same unit than the input.'
     PRINT *,'      '
     PRINT *,'     SEA ALSO :'
     PRINT *,'       cdfstd, cdfstdevw, cdfstdevts.'
     STOP
  ENDIF

  ijarg = 1  ; ixtra = 0
  DO WHILE ( ijarg <= narg) 
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-nc4' ) ; lnc4 = .TRUE.
     CASE ( '-o'   ) ; CALL getarg(ijarg, cf_out ) ; ijarg=ijarg+1
     CASE DEFAULT
        ixtra = ixtra + 1
        SELECT CASE ( ixtra ) 
        CASE ( 1 ) ; cf_in  = cldum
        CASE ( 2 ) ; cf_in2 = cldum
        CASE DEFAULT
           PRINT *, ' Too many variables ' ; STOP
        END SELECT
     END SELECT
  ENDDO

  ! check existence of files
  lchk = lchk .OR. chkfile(cf_in  )
  lchk = lchk .OR. chkfile(cf_in2 )
  IF (lchk ) STOP ! missing file

  npiglo = getdim (cf_in, cn_x)
  npjglo = getdim (cf_in, cn_y)
  npk    = getdim (cf_in, cn_z)
  npt    = getdim (cf_in, cn_t)

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

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  ALLOCATE( zvbar(npiglo,npjglo), zvba2(npiglo,npjglo) )
  ALLOCATE( dsdev(npiglo,npjglo), tim(npt)             )

  ncout = create      (cf_out, cf_in,    npiglo, npjglo, npk       , ld_nc4=lnc4 )
  ierr  = createvar   (ncout,  stypvaro, 1,      ipko,   id_varout , ld_nc4=lnc4 )
  ierr  = putheadervar(ncout,  cf_in,    npiglo, npjglo, npk                     )

  cv_in2 = TRIM(cv_in)//'_sqd'
  DO jt = 1, npt
     zvbar(:,:) = getvar(cf_in,  cv_in,  1, npiglo, npjglo, ktime=jt)
     zvba2(:,:) = getvar(cf_in2, cv_in2, 1, npiglo, npjglo, ktime=jt)

     dsdev(:,:) = SQRT ( DBLE(zvba2(:,:) - zvbar(:,:)*zvbar(:,:)) )

     ierr = putvar(ncout, id_varout(1), REAL(dsdev), 1, npiglo, npjglo, ktime=jt)
  END DO

  tim  = getvar1d(cf_in, cn_vtimec, npt     )
  ierr = putvar1d(ncout, tim,       npt, 'T')

  ierr = closeout(ncout)

END PROGRAM cdfrmsssh
