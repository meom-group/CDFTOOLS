PROGRAM cdfscale
  !!======================================================================
  !!                     ***  PROGRAM  cdfscale  ***
  !!=====================================================================
  !!  ** Purpose : Replace a variable in the file by its value x scale
  !!               given in argument
  !!
  !! History : 3.0  : 12/2011  : J.M. Molines : Original code
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

  INTEGER(KIND=4)                               :: jk, jvar, jt           ! dummy loop index
  INTEGER(KIND=4)                               :: ivar                   ! index of working variable
  INTEGER(KIND=4)                               :: ierr                   ! working integer
  INTEGER(KIND=4)                               :: narg, iargc            ! browse line
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
     PRINT *,' usage : cdfscale INOUT-file IN-var scale '
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Replace IN-var in INOUT-file by its values x scale.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       INOUT-file : netcdf input file (!overwritten!).'
     PRINT *,'       IN-var : netcdf variable to be scaled.'
     PRINT *,'       scale : Scale value to be used (multiplication factor).'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : input file is rewritten ' 
     PRINT *,'         variables : same name as input.' 
     STOP
  ENDIF

  CALL getarg(1, cf_inout) 
  CALL getarg(2, cv_inout) 
  CALL getarg(3, cldum)  ; READ(cldum,*) vscale

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
    IF ( cv_inout == cv_names(jvar) ) THEN
      ivar = jvar
    ENDIF
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
