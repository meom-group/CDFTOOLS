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

  INTEGER(KIND=4)                            :: jk, jt, jvar      ! dummy loop index
  INTEGER(KIND=4)                            :: narg, iargc       ! command line
  INTEGER(KIND=4)                            :: ijarg             ! command line
  INTEGER(KIND=4)                            :: npiglo, npjglo    ! size of the domain
  INTEGER(KIND=4)                            :: npk, npt          ! size of the domain
  INTEGER(KIND=4)                            :: ncout             ! ncid of output variable
  INTEGER(KIND=4)                            :: ierr              ! error status
  INTEGER(KIND=4), DIMENSION(2)              :: ipko, id_varout   ! output variable

  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: zvbar, zvba2      ! mean and mean2 variable

  REAL(KIND=8), DIMENSION(:),    ALLOCATABLE :: dtim              ! time counter
  REAL(KIND=8), DIMENSION(:,:),  ALLOCATABLE :: dsdev             ! standard deviation

  CHARACTER(LEN=256)                         :: cf_tfil           ! input mean file name
  CHARACTER(LEN=256)                         :: cf_tfil2          ! input mean2 file name
  CHARACTER(LEN=256)                         :: cf_sfil           ! input mean file name
  CHARACTER(LEN=256)                         :: cf_sfil2          ! input mean2 file name
  CHARACTER(LEN=256)                         :: cf_out = 'stdevts.nc'! output file name
  CHARACTER(LEN=256)                         :: cv_in, cv_in2     ! input variable names
  CHARACTER(LEN=256)                         :: cldum             ! dummy character variable
  CHARACTER(LEN=256)                         :: cl_votemper       ! local name for T
  CHARACTER(LEN=256)                         :: cl_vosaline       ! local name for S
  CHARACTER(LEN=256)                         :: cl_votemper2      ! local name for T2
  CHARACTER(LEN=256)                         :: cl_vosaline2      ! local name for S2
  CHARACTER(LEN=256), DIMENSION(4)           :: cv_namesi         ! input variable names
  CHARACTER(LEN=256), DIMENSION(4)           :: cl_fil            ! input variable names

  TYPE(variable), DIMENSION(2)               :: stypvaro          ! output data structure

  LOGICAL                                    :: lchk = .FALSE.    ! flag for missing files
  LOGICAL                                    :: lnc4 = .FALSE.    ! Use nc4 with chunking and deflation
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfstdevts -t T-file -t2 T2-file [-o OUT-file] [-nc4] ...'
     PRINT *,'                   [-s S-file] [-s2 S2-file ] ...'
     PRINT *,'           [-var VAR-temp VAR-sal VAR-temp2 VAR-sal2 ]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the standard deviation of the temperature and salinity from'
     PRINT *,'       their mean and  mean squared values.  The mean squared values for T'
     PRINT *,'       and S are not automatically computed with cdfmoy (need some namelist'
     PRINT *,'       modification). Consider using the more generic program ''cdfstd'' if'
     PRINT *,'       you do not already have the mean T2 and mean S2.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -t T-file  : netcdf file with mean values for T, S' 
     PRINT *,'       -t T2-file : netcdf file with mean squared values for T,S' 
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-t S-file ]  : salinity file if not in T-file.'
     PRINT *,'       [-t S2-file]  : salinity sqd file if not in T2-file.'
     PRINT *,'       [-o OUT-file] : specify output file name instead of ',TRIM(cf_out)
     PRINT *,'       [-nc4 ] : Use netcdf4 output with chunking and deflation level 1.'
     PRINT *,'            This option is effective only if cdftools are compiled with'
     PRINT *,'            a netcdf library supporting chunking and deflation.'
     PRINT *,'       [-var VAR-temp VAR-sal VAR-temp2 VAR-sal2 ] : Specify variable names'
     PRINT *,'            for mean and mean-squared temperature and salinity.'
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
     PRINT *,'       cdfstd, cdfrmsssh, cdfstdevw, cdfstats.'
     PRINT *,'      '
     STOP 
  ENDIF

  cl_votemper  = cn_votemper
  cl_vosaline  = cn_vosaline
  cl_votemper2 = TRIM(cn_votemper)//'_sqd'
  cl_vosaline2 = TRIM(cn_vosaline)//'_sqd'

  cf_sfil ='none'
  cf_sfil2='none'

  ijarg = 1  
  DO WHILE ( ijarg <= narg) 
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-t' ) ; CALL getarg(ijarg, cf_tfil ) ; ijarg=ijarg+1
     CASE ( '-t2') ; CALL getarg(ijarg, cf_tfil2) ; ijarg=ijarg+1
        ! option
     CASE ( '-s' ) ; CALL getarg(ijarg, cf_sfil ) ; ijarg=ijarg+1
     CASE ( '-s2') ; CALL getarg(ijarg, cf_sfil2) ; ijarg=ijarg+1
     CASE ( '-o' ) ; CALL getarg(ijarg, cf_out  ) ; ijarg=ijarg+1
     CASE ('-nc4') ; lnc4 = .TRUE.
     CASE ('-var') ; CALL getarg(ijarg, cl_votemper ) ; ijarg=ijarg+1
       ;           ; CALL getarg(ijarg, cl_vosaline ) ; ijarg=ijarg+1
       ;           ; CALL getarg(ijarg, cl_votemper2) ; ijarg=ijarg+1
       ;           ; CALL getarg(ijarg, cl_vosaline2) ; ijarg=ijarg+1
     CASE DEFAULT  ; PRINT *,' ERROR : ', TRIM(cldum),' : unknown option.'
     END SELECT
  ENDDO

  IF ( cf_sfil  == 'none' ) cf_sfil  = cf_tfil
  IF ( cf_sfil2 == 'none' ) cf_sfil2 = cf_tfil2

  cv_namesi(1) = cl_votemper  ; cl_fil(1)=cf_tfil
  cv_namesi(2) = cl_vosaline  ; cl_fil(2)=cf_sfil
  cv_namesi(3) = cl_votemper2 ; cl_fil(3)=cf_tfil2
  cv_namesi(4) = cl_vosaline2 ; cl_fil(4)=cf_sfil2

  ! check existence of files
  lchk = lchk .OR. chkfile(cf_tfil )
  lchk = lchk .OR. chkfile(cf_tfil2)
  lchk = lchk .OR. chkfile(cf_sfil )
  lchk = lchk .OR. chkfile(cf_sfil2)
  IF (lchk ) STOP 99 ! missing file

  npiglo = getdim (cf_tfil, cn_x)
  npjglo = getdim (cf_tfil, cn_y)
  npk    = getdim (cf_tfil, cn_z)
  npt    = getdim (cf_tfil, cn_t)

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  ALLOCATE( zvbar(npiglo,npjglo), zvba2(npiglo,npjglo) )
  ALLOCATE( dsdev(npiglo,npjglo), dtim(npt)            )

  CALL CreateOutput

  DO jvar = 1, 2   ! temperature and salinity
     cv_in  = cv_namesi(jvar  )
     cv_in2 = cv_namesi(jvar+2)
     DO jt = 1, npt
        DO jk = 1, npk
           zvbar(:,:) = getvar(cl_fil(jvar)  , cv_in,  jk, npiglo, npjglo, ktime=jt)
           zvba2(:,:) = getvar(cl_fil(jvar+2), cv_in2, jk, npiglo, npjglo, ktime=jt)

           dsdev(:,:) = SQRT ( DBLE(zvba2(:,:) - zvbar(:,:)*zvbar(:,:)) )

           ierr = putvar(ncout, id_varout(jvar), REAL(dsdev), jk, npiglo, npjglo, ktime=jt)
        END DO
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
    stypvaro(1)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvaro(1)%cname             = TRIM(cl_votemper)//'_stdev'
    stypvaro(1)%cunits            = 'DegC'
    stypvaro(1)%rmissing_value    = 0.
    stypvaro(1)%valid_min         = 0.
    stypvaro(1)%valid_max         = 20
    stypvaro(1)%clong_name        = 'STDEV_Temperature'
    stypvaro(1)%cshort_name       = TRIM(cl_votemper)//'_stdev'
    stypvaro(1)%conline_operation = 'N/A'
    stypvaro(1)%caxis             = 'TZYX'

    ipko(2) = npk
    stypvaro(2)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvaro(2)%cname             = TRIM(cl_vosaline)//'_stdev'
    stypvaro(2)%cunits            = 'PSU'
    stypvaro(2)%rmissing_value    = 0.
    stypvaro(2)%valid_min         = 0.
    stypvaro(2)%valid_max         = 10
    stypvaro(2)%clong_name        = 'STDEV_Salinity'
    stypvaro(2)%cshort_name       = TRIM(cl_vosaline)//'_stdev'
    stypvaro(2)%conline_operation = 'N/A'
    stypvaro(2)%caxis             = 'TZYX'

    ncout = create      (cf_out, cf_tfil,    npiglo, npjglo, npk       , ld_nc4=lnc4)
    ierr  = createvar   (ncout,  stypvaro, 2,      ipko,   id_varout , ld_nc4=lnc4)
    ierr  = putheadervar(ncout,  cf_tfil,    npiglo, npjglo, npk       )

    dtim = getvar1d(cf_tfil, cn_vtimec, npt     )
    ierr = putvar1d(ncout, dtim,      npt, 'T')

  END SUBROUTINE CreateOutput

END PROGRAM cdfstdevts
