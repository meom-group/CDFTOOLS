PROGRAM cdfcsp
  !!======================================================================
  !!                     ***  PROGRAM  cdfcsp  ***
  !!=====================================================================
  !!  ** Purpose : Replace the masked part of the arrays (marked with
  !!               special values) with spval zero. Replace consistently
  !!               the definition of the spval in the variable attribut.
  !!
  !! History : 2.1  : 10/2006  : F. Castruccio : Original code
  !!           3.0  : 12/2010  : J.M. Molines  : Doctor norm + Lic.
  !!----------------------------------------------------------------------
  USE cdfio 
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2010, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                               :: jf, jk, jvar, jt ! dummy loop index
  INTEGER(KIND=4)                               :: narg, iargc      ! 
  INTEGER(KIND=4)                               :: npiglo, npjglo   ! size of the domain
  INTEGER(KIND=4)                               :: npk , npt        ! size of the domain
  INTEGER(KIND=4)                               :: ncid, ierr       ! ncdf related integer
  INTEGER(KIND=4)                               :: nvars            ! Number of variables in a file
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: ipk              ! arrays of vertical level for each var
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: id_var           ! arrays of var id

  REAL(KIND=4), DIMENSION (:,:),    ALLOCATABLE :: tab              ! working array
  REAL(KIND=4)                                  :: zspval           ! special value read in file

  CHARACTER(LEN=256)                            :: cf_in            ! input file name
  CHARACTER(LEN=256)                            :: cunits           ! units attribute
  CHARACTER(LEN=256)                            :: clname           ! long name attribute
  CHARACTER(LEN=256)                            :: csname           ! short name attribute
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_names         ! array of var name

  TYPE(variable), DIMENSION(:),     ALLOCATABLE :: stypvar          ! type for attributes
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()

  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfcsp list_of_files '
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Replace missing_values by 0 and update attribute' 
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       The list of cdf file to process, all variables will be processed' 
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
  CALL getarg (1, cf_in)
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

  DO jf = 1, narg
     CALL getarg (jf, cf_in)
     IF ( chkfile (cf_in) ) STOP ! missing file
     PRINT *, 'Change spval on file ', cf_in
     ncid = ncopen(cf_in)
     npt  = getdim (cf_in,cn_t)
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
           ierr = cvaratt   (cf_in, cv_names(jvar), cunits, 0.,     clname, csname)
          DO jt=1,npt
           DO jk = 1, ipk(jvar) 
              tab(:,:) = getvar(cf_in, cv_names(jvar),   jk, npiglo, npjglo, ktime=jt )
              WHERE( tab(:,:) == zspval ) tab(:,:) = 0.
              ierr     = putvar(ncid, id_var(jvar), tab, jk, npiglo, npjglo, ktime=jt )
           ENDDO
          END DO
        ENDIF 
     ENDDO
  ENDDO

  ierr = closeout(ncid)

END PROGRAM cdfcsp
