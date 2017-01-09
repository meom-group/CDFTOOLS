PROGRAM cdfstdevts
  !!======================================================================
  !!                     ***  PROGRAM  cdfstdevts  ***
  !!=====================================================================
  !!  ** Purpose : Compute the RMS of T and S, from the mean squared value.
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

  INTEGER(KIND=4)                            :: jk, jt, jvar      ! dummy loop index
  INTEGER(KIND=4)                            :: narg, iargc       ! command line
  INTEGER(KIND=4)                            :: ijarg, ireq       ! command line
  INTEGER(KIND=4)                            :: npiglo, npjglo    ! size of the domain
  INTEGER(KIND=4)                            :: npk, npt          ! size of the domain
  INTEGER(KIND=4)                            :: ncout             ! ncid of output variable
  INTEGER(KIND=4)                            :: ierr              ! error status
  INTEGER(KIND=4), DIMENSION(2)              :: ipko, id_varout   ! output variable

  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: zvbar, zvba2      ! mean and mean2 variable
  REAL(KIND=4), DIMENSION(:),    ALLOCATABLE :: tim               ! time counter

  REAL(KIND=8), DIMENSION(:,:),  ALLOCATABLE :: dsdev             ! standard deviation

  CHARACTER(LEN=256)                         :: cf_in             ! input mean file name
  CHARACTER(LEN=256)                         :: cf_in2            ! input mean2 file name
  CHARACTER(LEN=256)                         :: cf_out = 'stdevts.nc'! output file name
  CHARACTER(LEN=256)                         :: cv_in, cv_in2     ! input variable names
  CHARACTER(LEN=256)                         :: cldum             ! dummy character variable
  CHARACTER(LEN=256), DIMENSION(2)           :: cv_namesi         ! input variable names

  TYPE(variable), DIMENSION(2)               :: stypvaro          ! output data structure

  LOGICAL                                    :: lchk = .FALSE.    ! flag for missing files
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  cv_namesi(1) = cn_votemper
  cv_namesi(2) = cn_vosaline

  narg= iargc()
  IF ( narg /= 2 ) THEN
     PRINT *,' usage : cdfstdevts T-file T2-file '
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the standard deviation of the temperature'
     PRINT *,'       and salinity from their mean and  mean square values. '
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       T-file  : netcdf file with mean values for T, S' 
     PRINT *,'       T2-file : netcdf file with mean squared values for T,S' 
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       none' 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : ', TRIM(cn_votemper)//'_stdev, same unit than the input.'
     PRINT *,'                     ', TRIM(cn_vosaline)//'_stdev, same unit than the input.'
     PRINT *,'      '
     PRINT *,'     SEA ALSO :'
     PRINT *,'       cdfstd, cdfrmsssh, cdfstdevw.'
     STOP
  ENDIF

  ijarg = 1  ; ireq = 0
  DO WHILE ( ijarg <= narg) 
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE DEFAULT
        ireq = ireq + 1
        SELECT CASE ( ireq ) 
        CASE ( 1 ) ; cf_in  = cldum
        CASE ( 2 ) ; cf_in2 = cldum
        CASE DEFAULT
           PRINT *, ' Too many variables ' ; STOP
        END SELECT
     END SELECT
  ENDDO

  ! check existence of files
  lchk = lchk .OR. chkfile(cf_in )
  lchk = lchk .OR. chkfile(cf_in2)
  IF (lchk ) STOP ! missing file

  npiglo = getdim (cf_in, cn_x)
  npjglo = getdim (cf_in, cn_y)
  npk    = getdim (cf_in, cn_z)
  npt    = getdim (cf_in, cn_t)

  ipko(1) = npk
  stypvaro(1)%cname             = TRIM(cn_votemper)//'_stdev'
  stypvaro(1)%cunits            = 'DegC'
  stypvaro(1)%rmissing_value    = 0.
  stypvaro(1)%valid_min         = 0.
  stypvaro(1)%valid_max         = 20
  stypvaro(1)%clong_name        = 'STDEV_Temperature'
  stypvaro(1)%cshort_name       = TRIM(cn_votemper)//'_stdev'
  stypvaro(1)%conline_operation = 'N/A'
  stypvaro(1)%caxis             = 'TZYX'

  ipko(2) = npk
  stypvaro(2)%cname             = TRIM(cn_vosaline)//'_stdev'
  stypvaro(2)%cunits            = 'PSU'
  stypvaro(2)%rmissing_value    = 0.
  stypvaro(2)%valid_min         = 0.
  stypvaro(2)%valid_max         = 10
  stypvaro(2)%clong_name        = 'STDEV_Salinity'
  stypvaro(2)%cshort_name       = TRIM(cn_vosaline)//'_stdev'
  stypvaro(2)%conline_operation = 'N/A'
  stypvaro(2)%caxis             = 'TZYX'

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  ALLOCATE( zvbar(npiglo,npjglo), zvba2(npiglo,npjglo) )
  ALLOCATE( dsdev(npiglo,npjglo), tim(npt)             )

  ncout = create      (cf_out, cf_in,    npiglo, npjglo, npk       )
  ierr  = createvar   (ncout,  stypvaro, 2,      ipko,   id_varout )
  ierr  = putheadervar(ncout,  cf_in,    npiglo, npjglo, npk       )

  DO jvar = 1, 2
  cv_in  = cv_namesi(jvar)
  cv_in2 = TRIM(cv_in)//'_sqd'
  DO jt = 1, npt
     DO jk = 1, npk
        zvbar(:,:) = getvar(cf_in,  cv_in,  jk, npiglo, npjglo, ktime=jt)
        zvba2(:,:) = getvar(cf_in2, cv_in2, jk, npiglo, npjglo, ktime=jt)

        dsdev(:,:) = SQRT ( DBLE(zvba2(:,:) - zvbar(:,:)*zvbar(:,:)) )

        ierr = putvar(ncout, id_varout(jvar), REAL(dsdev), jk, npiglo, npjglo, ktime=jt)
     END DO
  END DO
  END DO

  tim  = getvar1d(cf_in, cn_vtimec, npt     )
  ierr = putvar1d(ncout, tim,       npt, 'T')

  ierr = closeout(ncout)

END PROGRAM cdfstdevts
