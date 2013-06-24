PROGRAM cdfmltmask
  !!======================================================================
  !!                     ***  PROGRAM  cdfmltmask  ***
  !!=====================================================================
  !!  ** Purpose : multiplication of file by a mask (0,1)
  !!
  !! History : 2.1  : 06/2007  : M. Juza      : Original code
  !!         : 2.1  : 06/2007  : P. Mathiot   : add forcing capabilities
  !!           3.0  : 01/2011  : J.M. Molines : Doctor norm + Lic.
  !!           3.0  : 06/2013  : J.M. Molines : add multi variable capability
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

  INTEGER(KIND=4)                             :: jk, jt, jvar    ! dummy loop index
  INTEGER(KIND=4)                             :: ierr            ! error status
  INTEGER(KIND=4)                             :: narg, iargc     ! command line 
  INTEGER(KIND=4)                             :: npiglo, npjglo  ! size of the domain
  INTEGER(KIND=4)                             :: npk, npt        ! size of the domain
  INTEGER(KIND=4)                             :: npkmask         ! vertical levels in mask file
  INTEGER(KIND=4)                             :: nvar=1          ! number of variable to process(default 1)
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE  :: nvpk            ! vertical levels in working variable

  REAL(KIND=4)                                :: zspval          ! missing value attribute
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: zv              ! cv_in at jk level 
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: zmask           ! mask at jk level 
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: zvmask          ! masked cv_in at jk level

  CHARACTER(LEN=256)                          :: cunits          ! units attribute
  CHARACTER(LEN=256)                          :: clname          ! long_name attribute
  CHARACTER(LEN=256)                          :: csname          ! short_name attribute
  CHARACTER(LEN=256)                          :: cf_in           ! input file name
  CHARACTER(LEN=256)                          :: cf_msk          ! input mask file name
  CHARACTER(LEN=256)                          :: cvartype        ! variable position on Cgrid
  CHARACTER(LEN=256)                          :: cv_dep          ! depth dim name
  CHARACTER(LEN=256)                          :: ctmp            ! dummy string
  CHARACTER(LEN=20)                           :: cv_msk          ! mask variable name
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_in         ! cdf variable names to process
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfmltmask IN-file MSK-file IN-var1,var2,...  T| U | V | F | W | P'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Multiply IN-var of IN-file by the mask corresponding to the' 
     PRINT *,'       C-grid point position given as last argument.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       IN-file  : input netcdf file.' 
     PRINT *,'       MSK-file : input netcdf mask file.' 
     PRINT *,'       IN-var1,var2,...   : Comma separated list of variable names to mask.'
     PRINT *,'       T| U | V | F | W | P : C-grid position of IN-var'
     PRINT *,'                P indicate a polygon mask created by cdfpoly.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        none, all are given as arguments.'
     PRINT *,'      '
     PRINT *,'     OUTPUT : (jvar)'
     PRINT *,'       The output file is a copy of the input file with only'
     PRINT *,'       the requested variable masked.'
     PRINT *,'       netcdf file : IN-file_masked'
     PRINT *,'         variables : IN-var (same as input).'
     STOP
  ENDIF

  CALL getarg (1, cf_in    )
  CALL getarg (2, cf_msk   )
  CALL getarg (3, ctmp     ) ; CALL ParseVars(ctmp)
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

  npt   = getdim (cf_in, cn_t    )
  ALLOCATE( nvpk(nvar) )

  DO jvar = 1, nvar
     nvpk(jvar)  = getvdim(cf_in, cv_in(jvar))
     IF (nvpk(jvar) == 2 ) nvpk(jvar) = 1
     IF (nvpk(jvar) == 3 ) nvpk(jvar) = npk
  END DO

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt
  PRINT *, 'npkmask= ', npkmask

  DO jvar = 1, nvar
    PRINT *, 'nvpk(',TRIM(cv_in(jvar)),')   = ', nvpk(jvar)
  END DO

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
     DO jvar = 1, nvar  ! loop on variables
        DO jk = 1, nvpk(jvar)
        PRINT *,' Processing level ', jk,' variable ',TRIM(cv_in(jvar))
        IF ( npkmask > 1 ) THEN
        ! Read mask
          zmask(:,:) = getvar(cf_msk, cv_msk, jk, npiglo, npjglo)
        ENDIF
          ! Read cv_in
          zv(:,:) = getvar(cf_in, cv_in(jvar), jk, npiglo, npjglo, ktime=jt)
          ! Multiplication of cv_in by mask at level jk
          zvmask = zv * zmask
          ! Writing  on the copy of original file                 
          ierr = putvar(cf_in, cv_in(jvar), jk, npiglo, npjglo, 1, 1, ktime=jt, ptab=zvmask)
        END DO
     END DO
  END DO

  ! set missing value attribute for cv_in as 0.
  DO jvar = 1, nvar 
     ierr = getvaratt (cf_in, cv_in(jvar), cunits, zspval, clname, csname)
     ierr = cvaratt   (cf_in, cv_in(jvar), cunits, 0.,     clname, csname)
  END DO

CONTAINS
   SUBROUTINE ParseVars (cdum)
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE ParseVars  ***
      !!
      !! ** Purpose :  Decode variables names to be used
      !!
      !! ** Method  :  look for , in the argument string and set the number of
      !!         variable (nvaro), allocate cv_fix array and fill it with the
      !!         decoded  names.
      !!
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*), INTENT(in) :: cdum

      CHARACTER(LEN=80), DIMENSION(100) :: cl_dum  ! 100 is arbitrary
      INTEGER  :: ji
      INTEGER  :: inchar,  i1=1
      !!----------------------------------------------------------------------

      inchar= LEN(TRIM(cdum))
      ! scan the input string and look for ',' as separator
      DO ji=1,inchar
         IF ( cdum(ji:ji) == ',' ) THEN
            cl_dum(nvar) = cdum(i1:ji-1)
            i1=ji+1
            nvar=nvar+1
         ENDIF
      ENDDO

      ! last name of the list does not have a ','
      cl_dum(nvar) = cdum(i1:inchar)

      ALLOCATE ( cv_in(nvar) )
      DO ji=1, nvar
         cv_in(ji) = cl_dum(ji)
      ENDDO
   END SUBROUTINE ParseVars

END PROGRAM cdfmltmask 
