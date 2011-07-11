PROGRAM cdfmltmask
  !!======================================================================
  !!                     ***  PROGRAM  cdfmltmask  ***
  !!=====================================================================
  !!  ** Purpose : multiplication of file by a mask (0,1)
  !!
  !! History : 2.1  : 06/2007  : M. Juza      : Original code
  !!         : 2.1  : 06/2007  : P. Mathiot   : add forcing capabilities
  !!           3.0  : 01/2011  : J.M. Molines : Doctor norm + Lic.
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

  INTEGER(KIND=4)                             :: jk, jt          ! dummy loop index
  INTEGER(KIND=4)                             :: ierr            ! error status
  INTEGER(KIND=4)                             :: narg, iargc     ! command line 
  INTEGER(KIND=4)                             :: npiglo, npjglo  ! size of the domain
  INTEGER(KIND=4)                             :: npk, npt        ! size of the domain
  INTEGER(KIND=4)                             :: nvpk            ! vertical levels in working variable
  INTEGER(KIND=4)                             :: npkmask         ! vertical levels in mask file

  REAL(KIND=4)                                :: zspval          ! missing value attribute
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: zv              ! cv_in at jk level 
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: zmask           ! mask at jk level 
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: zvmask          ! masked cv_in at jk level

  CHARACTER(LEN=256)                          :: cunits          ! units attribute
  CHARACTER(LEN=256)                          :: clname          ! long_name attribute
  CHARACTER(LEN=256)                          :: csname          ! short_name attribute
  CHARACTER(LEN=256)                          :: cf_in           ! input file name
  CHARACTER(LEN=256)                          :: cf_msk          ! input mask file name
  CHARACTER(LEN=256)                          :: cv_in           ! cdf variable name
  CHARACTER(LEN=256)                          :: cvartype        ! variable position on Cgrid
  CHARACTER(LEN=256)                          :: cv_dep          ! depth dim name
  CHARACTER(LEN=256)                          :: ctmp            ! dummy string
  CHARACTER(LEN=20)                           :: cv_msk          ! mask variable name
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfmltmask IN-file MSK-file IN-var T| U | V | F | W | P'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Multiply IN-var of IN-file by the mask corresponding to the' 
     PRINT *,'       C-grid point position given as last argument.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       IN-file  : input netcdf file.' 
     PRINT *,'       MSK-file : input netcdf mask file.' 
     PRINT *,'       IN-var   : input variable name.'
     PRINT *,'       T| U | V | F | W | P : C-grid position of IN-var'
     PRINT *,'                P indicate a polygon mask created by cdfpoly.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        none, all are given as arguments.'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       The output file is a copy of the input file with only'
     PRINT *,'       the requested variable masked.'
     PRINT *,'       netcdf file : IN-file_masked'
     PRINT *,'         variables : IN-var (same as input).'
     STOP
  ENDIF

  CALL getarg (1, cf_in    )
  CALL getarg (2, cf_msk   )
  CALL getarg (3, cv_in    )
  CALL getarg (4, cvartype )

  IF ( chkfile (cf_in) .OR. chkfile(cf_msk) ) STOP ! missing files

  ! append _masked to input file name and copy initial file to new file, which will be modified
  !  using dd more efficient than cp for big files
  ctmp   = TRIM(cf_in)//'_masked'
  CALL system(' dd bs=10000000 if='//TRIM(cf_in)//' of='//TRIM(ctmp) )
  cf_in = ctmp

  PRINT *,' Working on copy : ', TRIM(cf_in)

  npiglo = getdim (cf_in,cn_x)
  npjglo = getdim (cf_in,cn_y)
  npk    = getdim (cf_in,cn_z, cdtrue=cv_dep, kstatus=ierr)

  IF (ierr /= 0 ) THEN
     npk   = getdim (cf_in, 'z', cdtrue=cv_dep, kstatus=ierr)
     IF (ierr /= 0 ) THEN
       npk   = getdim (cf_in, 'sigma', cdtrue=cv_dep, kstatus=ierr)
        IF ( ierr /= 0 ) THEN
          npk = getdim (cf_in, 'nav_lev', cdtrue=cv_dep, kstatus=ierr)
            IF ( ierr /= 0 ) THEN
              PRINT *,' assume file with no depth'
              npk=0
            ENDIF
        ENDIF
     ENDIF
  ENDIF

  npkmask = getdim (cf_msk, cn_z, cdtrue=cv_dep, kstatus=ierr)
  IF (ierr /= 0 ) THEN
     npkmask  = getdim (cf_msk, 'z', cdtrue=cv_dep, kstatus=ierr)
       IF ( ierr /= 0 ) THEN
            npkmask = getdim (cf_msk, 'nav_lev', cdtrue=cv_dep, kstatus=ierr)
            IF ( ierr /= 0 ) THEN
              PRINT *,' assume file with no depth'
              npkmask=0
            ENDIF
       ENDIF
  ENDIF

  npt   = getdim (cf_in, cn_t )
  nvpk  = getvdim(cf_in, cv_in)

  IF (nvpk == 2 ) nvpk = 1
  IF (nvpk == 3 ) nvpk = npk

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt
  PRINT *, 'nvpk   = ', nvpk

  IF (npk==0) npk=1

  ! Allocate arrays
  ALLOCATE( zmask(npiglo,npjglo) )
  ALLOCATE( zv   (npiglo,npjglo) )
  ALLOCATE(zvmask(npiglo,npjglo) )

  SELECT CASE (TRIM(cvartype))
  CASE ( 'T' )
     cv_msk='tmask'
  CASE ( 'U' )
     cv_msk='umask'
  CASE ( 'V' )
     cv_msk='vmask'
  CASE ( 'F' )
     cv_msk='fmask'
  CASE ( 'W' )
     cv_msk='tmask'
  CASE ( 'P' )   ! for polymask 
     cv_msk='polymask'
  CASE DEFAULT
     PRINT *, 'this type of variable is not known :', TRIM(cvartype)
     STOP
  END SELECT

  IF ( npkmask <= 1 ) THEN 
        zmask(:,:) = getvar(cf_msk, cv_msk, 1, npiglo, npjglo)
  ENDIF

  DO jt = 1, npt
     IF (MOD(jt,100)==0) PRINT *, jt,'/', npt
     DO jk = 1,nvpk
        ! Read cv_in
        zv(:,:) = getvar(cf_in, cv_in, jk, npiglo, npjglo, ktime=jt)
        IF ( npkmask > 1 ) THEN
        ! Read mask
          zmask(:,:) = getvar(cf_msk, cv_msk, jk, npiglo, npjglo)
        ENDIF
        ! Multiplication of cv_in by mask at level jk
        zvmask = zv * zmask
        ! Writing  on the copy of original file                 
        ierr = putvar(cf_in, cv_in, jk, npiglo, npjglo, 1, 1, ktime=jt, ptab=zvmask)
     END DO
  END DO
  ! set missing value attribute for cv_in as 0.
  ierr = getvaratt (cf_in, cv_in, cunits, zspval, clname, csname)
  ierr = cvaratt   (cf_in, cv_in, cunits, 0.,     clname, csname)

END PROGRAM cdfmltmask 
