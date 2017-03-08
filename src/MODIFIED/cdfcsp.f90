PROGRAM cdfcsp
  !!======================================================================
  !!                     ***  PROGRAM  cdfcsp  ***
  !!=====================================================================
  !!  ** Purpose : Replace the masked part of the arrays (marked with
  !!               special values) with spval zero. Replace consistently
  !!               the definition of the spval in the variable attribut.
  !!
  !! History : 2.1  : 10/2006  : F. Castruccio : Original code
  !!         :  4.0  : 03/2017  : J.M. Molines  
  !!           3.0  : 12/2010  : J.M. Molines  : Doctor norm + Lic.
  !!----------------------------------------------------------------------
  USE cdfio 
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class file_operations
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                               :: jf, jk, jvar, jt ! dummy loop index
  INTEGER(KIND=4)                               :: narg, iargc      ! 
  INTEGER(KIND=4)                               :: ijarg            !
  INTEGER(KIND=4)                               :: nfiles           ! number of files in the list
  INTEGER(KIND=4)                               :: npiglo, npjglo   ! size of the domain
  INTEGER(KIND=4)                               :: npk , npt        ! size of the domain
  INTEGER(KIND=4)                               :: ncid, ierr       ! ncdf related integer
  INTEGER(KIND=4)                               :: nvars            ! Number of variables in a file
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: ipk              ! arrays of vertical level for each var
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: id_var           ! arrays of var id

  REAL(KIND=4), DIMENSION (:,:),    ALLOCATABLE :: tab              ! working array
  REAL(KIND=4)                                  :: zspval           ! special value read in file
  REAL(KIND=4)                                  :: spval=0.         ! special value written in file

  CHARACTER(LEN=256)                            :: cf_in            ! input file name
  CHARACTER(LEN=256)                            :: cunits           ! units attribute
  CHARACTER(LEN=256)                            :: clname           ! long name attribute
  CHARACTER(LEN=256)                            :: csname           ! short name attribute
  CHARACTER(LEN=256)                            :: cldum            ! working char variable
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_names         ! array of var name
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cf_lst           ! list of input files

  TYPE(variable), DIMENSION(:),     ALLOCATABLE :: stypvar          ! type for attributes
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()

  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfcsp -l LIST-files [-v value ]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Replace missing_values by 0 and update attribute.' 
     PRINT *,'       This program is not working properly with NETCDF4/HDF5 files!'
     PRINT *,'      '
     PRINT *,'     CAUTION :'
     PRINT *,'      ################################'
     PRINT *,'      # INPUT FILES ARE OVER-WRITTEN #'
     PRINT *,'      ################################'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -l LIST-files : The list of cdf file to process, all variables will '
     PRINT *,'              be processed.' 
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-v value ] : use value instead of 0 as the new missing_value'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        none'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : same as input file (modified)'
     PRINT *,'         variables : same as input file'
     STOP
  ENDIF

  !! Initialisation from 1st file (all file are assume to have the same geometry)
  ijarg=1
  DO WHILE ( ijarg <= narg )
     CALL getarg( ijarg,cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-f' ) ; CALL GetFileList
        ! options
     CASE ( '-v' ) ; CALL getarg (ijarg, cldum ) ; ijarg=ijarg+1  ; READ(cldum,*) spval
     CASE DEFAULT  ; PRINT *,' ERROR : ', TRIM(cldum) , ' : unknown option.' ; STOP
     END SELECT
  ENDDO
  cf_in = cf_lst(1)

  IF ( chkfile (cf_in) ) STOP ! missing file

  npiglo = getdim (cf_in, cn_x)
  npjglo = getdim (cf_in, cn_y)
  npk    = getdim (cf_in, cn_z, kstatus=ierr)
  npt    = getdim (cf_in, cn_t)

  IF (ierr /= 0 ) THEN
     npk   = getdim (cf_in,'z',kstatus=ierr)
     IF (ierr /= 0 ) THEN
        PRINT *, "ASSUME NO VERTICAL DIMENSIONS !"
        npk=0
     ENDIF
  ENDIF

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  ALLOCATE( tab(npiglo,npjglo) )

  nvars = getnvar(cf_in)

  ALLOCATE (cv_names(nvars), id_var(nvars),ipk(nvars), stypvar(nvars))

  cv_names(:) = getvarname(cf_in, nvars, stypvar)
  ipk(:)      = getipk    (cf_in, nvars         )
  id_var(:)   = getvarid  (cf_in, nvars         )

  DO jf = 1, nfiles
     cf_in = cf_lst(jf) 
     IF ( chkfile (cf_in) ) STOP ! missing file
     PRINT *, 'Change spval on file ', cf_in
     ncid = ncopen(cf_in)
     npt  = getdim (cf_in, cn_t)
     DO jvar = 1,nvars
        IF ( cv_names(jvar) == cn_vlon2d  .OR. &
             & cv_names(jvar) == cn_vlat2d  .OR. &
             & cv_names(jvar) == cn_vtimec  .OR. &
             & cv_names(jvar) == cn_vdeptht .OR. &
             & cv_names(jvar) == cn_vdepthu .OR. &
             & cv_names(jvar) == cn_vdepthv      )  THEN
           ! skip these variable
        ELSE
           ierr = getvaratt (cf_in, cv_names(jvar), cunits, zspval, clname, csname)
           ierr = cvaratt   (cf_in, cv_names(jvar), cunits,  spval, clname, csname)
           DO jt=1,npt
              DO jk = 1, ipk(jvar) 
                 tab(:,:) = getvar(cf_in, cv_names(jvar),   jk, npiglo, npjglo, ktime=jt )
                 WHERE( tab(:,:) == zspval ) tab(:,:) = spval
                 ierr     = putvar(ncid, id_var(jvar), tab, jk, npiglo, npjglo, ktime=jt )
              ENDDO
           END DO
        ENDIF
     ENDDO
  ENDDO

  ierr = closeout(ncid)

CONTAINS

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

END PROGRAM cdfcsp
