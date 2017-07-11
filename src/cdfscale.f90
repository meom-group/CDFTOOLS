PROGRAM cdfscale
  !!======================================================================
  !!                     ***  PROGRAM  cdfscale  ***
  !!=====================================================================
  !!  ** Purpose : Replace a variable in the file by its value x scale
  !!               given in argument
  !!
  !! History : 3.0  : 12/2011  : J.M. Molines : Original code
  !!         : 4.0  : 03/2017  : J.M. Molines  
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

  INTEGER(KIND=4)                               :: jk, jvar, jt           ! dummy loop index
  INTEGER(KIND=4)                               :: ivar                   ! index of working variable
  INTEGER(KIND=4)                               :: ierr, ireq             ! working integer
  INTEGER(KIND=4)                               :: narg, iargc, ijarg     ! browse line
  INTEGER(KIND=4)                               :: ncid                   ! ncid of input file for rewrite
  INTEGER(KIND=4)                               :: npiglo, npjglo         ! size of the domain
  INTEGER(KIND=4)                               :: npk, npt               ! size of the domain
  INTEGER(KIND=4)                               :: nvars                  ! Number of variables in a file
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: ipk                    ! arrays of vertical level for each var
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: id_var                 ! arrays of var id

  REAL(KIND=4)                                  :: vscale                 ! spval, replace value
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: tab                    ! Arrays for data

  CHARACTER(LEN=256)                            :: cldum                  ! dummy string for getarg
  CHARACTER(LEN=256)                            :: cf_inout               ! file name
  CHARACTER(LEN=256)                            :: cunits, clname, csname ! attributes
  CHARACTER(LEN=256)                            :: cv_inout               ! variable name
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_names               ! array of var name

  TYPE(variable), DIMENSION(:),     ALLOCATABLE :: stypvar                ! type for attributes

  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfscale -f INOUT-file -v IN-var -s SCAL-factor '
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Replace IN-var in INOUT-file by its values x SCAL-factor.'
     PRINT *,'      '
     PRINT *,'     CAUTION :'
     PRINT *,'      #############################'
     PRINT *,'      # INPUT FILE IS OVERWRITTEN #'
     PRINT *,'      #############################'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f INOUT-file : netcdf input file (!overwritten!).'
     PRINT *,'       -v IN-var : netcdf variable to be scaled.'
     PRINT *,'       -s SCAL-factor : Scale value to be used (multiplication factor).'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : input file is rewritten ' 
     PRINT *,'         variables : same name as input.' 
     STOP 
  ENDIF

  ijarg=1 ; ireq=0
  DO WHILE ( ijarg <= narg )
      CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
      SELECT CASE ( cldum )
      CASE ( '-f'   ) ; CALL getarg(ijarg, cf_inout ) ; ijarg=ijarg+1  ; ireq=ireq+1
      CASE ( '-v'   ) ; CALL getarg(ijarg, cv_inout ) ; ijarg=ijarg+1  ; ireq=ireq+1
      CASE ( '-s'   ) ; CALL getarg(ijarg, cldum    ) ; ijarg=ijarg+1  ; ireq=ireq+1 ;  READ(cldum,*) vscale
      CASE DEFAULT    ; PRINT *,' ERROR :', TRIM(cldum),' : unknown option.' ; STOP 99
      END SELECT
  ENDDO
  ! 3 arguments are mandatory : here ijarg must be 4
  IF ( ireq /= 3 ) THEN ; PRINT *,' ERROR : need at least 3 arguments !' ; STOP 99 ; ENDIF

  IF ( chkfile (cf_inout) )  STOP 99 ! missing file

  npiglo = getdim (cf_inout, cn_x              )
  npjglo = getdim (cf_inout, cn_y              )
  npk    = getdim (cf_inout, cn_z, kstatus=ierr)

  IF (ierr /= 0 ) THEN
     npk   = getdim (cf_inout,'z',kstatus=ierr)
     IF (ierr /= 0 ) THEN
        PRINT *, 'ASSUME NO VERTICAL DIMENSIONS !'
        npk=0
     ENDIF
  ENDIF
  ncid   = ncopen ( cf_inout       )
  npt    = getdim ( cf_inout, cn_t )

  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk
  PRINT *, 'npt   =', npt

  ALLOCATE( tab(npiglo,npjglo) )

  nvars = getnvar(cf_inout)

  ALLOCATE (cv_names(nvars), id_var(nvars),ipk(nvars), stypvar(nvars))

  cv_names(:) = getvarname(cf_inout,nvars,stypvar)
  ipk(:)      = getipk(cf_inout,nvars)
  id_var(:)   = getvarid(cf_inout,nvars)

  ! look for cv_inout in the list of variables
  DO jvar = 1, nvars
    IF ( cv_inout == cv_names(jvar) ) ivar = jvar
  ENDDO

  PRINT *,' Working with ',TRIM(cv_inout),' variable number ', ivar
  PRINT *,'    IPK   : ', ipk(ivar)
  PRINT *,'    scale : ', vscale

  ! work only for ivar
  DO jt=1,npt
    DO jk = 1, ipk(ivar) 
      tab(:,:) = getvar(cf_inout, cv_names(ivar), jk, npiglo, npjglo, ktime=jt )
      tab(:,:) = tab(:,:) * vscale
      ierr = putvar(ncid, id_var(ivar), tab, jk, npiglo, npjglo, ktime=jt)
    ENDDO
  END DO

  ierr = closeout(ncid)

END PROGRAM cdfscale
