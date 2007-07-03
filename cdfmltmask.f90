PROGRAM cdfmltmask
  !!-------------------------------------------------------------------
  !!               ***  PROGRAM cdfmltmask  ***
  !!
  !!  **  Purpose  : multiplication of file by a mask (0,1)
  !!                 
  !!                 
  !!  * history:
  !!        Original :   Melanie JUZA  (june 2007) 
  !!        Modified :   Pierre Mathiot(june 2007) update for forcing fields
  !!-------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !! * Modules used
  USE cdfio

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: jk, jt, jkk
  INTEGER   :: narg, iargc                         !: command line 
  INTEGER   :: npiglo,npjglo, npk,npt              !: size of the domain
  INTEGER   :: nvpk                                !: vertical levels in working variable

  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::  zv !: cvar at jk level 
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::  zmask !: mask at jk level 
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE :: zvmask !: masked cvar at jk level

  CHARACTER(LEN=80) :: cfilev , cfilemask, ctmp
  CHARACTER(LEN=80) :: cvar, cvartype, cdep
  CHARACTER(LEN=20) ::  cvmask
  
  LOGICAL :: lforcing= .FALSE.  ! PM

  INTEGER    :: istatus

  ! constants

  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfmltmask ncfile maskfile cdfvar T| U | V | F | W '
     PRINT *,' Mask the file  '
     PRINT *,' output on ncfile_masked'
     STOP
  ENDIF

  CALL getarg (1, cfilev)
  CALL getarg (2, cfilemask)
  CALL getarg (3, cvar)
  CALL getarg (4, cvartype)

  ! append _masked to input file name and copy initial file to new file, which will be modified
  !  using dd more efficient than cp for big files
  ctmp=TRIM(cfilev)//'_masked'
   CALL system(' dd bs=10000000 if='//TRIM(cfilev)//' of='//TRIM(ctmp) )
  cfilev=ctmp
  print *, TRIM(cfilev)

  npiglo= getdim (cfilev,'x')
  npjglo= getdim (cfilev,'y')
  npk   = getdim (cfilev,'depth',cdtrue=cdep, kstatus=istatus) ; !print *, istatus
  IF (istatus /= 0 ) THEN
     npk   = getdim (cfilev,'z',cdtrue=cdep,kstatus=istatus) ; !print *, istatus
     IF (istatus /= 0 ) THEN
       npk   = getdim (cfilev,'sigma',cdtrue=cdep,kstatus=istatus) ; !print *, istatus
        IF ( istatus /= 0 ) THEN
          PRINT *,' assume file with no depth'
          npk=0
        ENDIF
     ENDIF
  ENDIF

  npt   = getdim (cfilev,'time')
  nvpk  = getvdim(cfilev,cvar)

  IF (nvpk == 2 ) nvpk = 1
  IF (nvpk == 3 ) nvpk = npk

  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk
  PRINT *, 'npt   =', npt
  PRINT *, 'nvpk  =', nvpk

  IF (npk==0) THEN  ! PM
     npk=1          ! PM
     lforcing= .TRUE.  ! PM
     PRINT *, 'it is a forcing field, assume npk=1'  ! PM
  END IF  ! PM

  ! Allocate arrays
  ALLOCATE ( zmask(npiglo,npjglo) )
  ALLOCATE ( zv(npiglo,npjglo) )
  ALLOCATE(zvmask(npiglo,npjglo))

  SELECT CASE (TRIM(cvartype))
  CASE ( 'T' )
     cvmask='tmask'
  CASE ( 'U' )
     cvmask='umask'
  CASE ( 'V' )
     cvmask='vmask'
  CASE ( 'F' )
     cvmask='fmask'
  CASE ( 'W' )
     cvmask='tmask'
  CASE DEFAULT
     PRINT *, 'this type of variable is not known :', TRIM(cvartype)
     STOP
  END SELECT

  DO jt = 1, npt
     IF (MOD(jt,100)==0) PRINT *, jt,'/', npt
     DO jkk = 1,nvpk
        ! Read cvar
        IF (.NOT. lforcing)jk=jkk !PM 
        IF (lforcing) jk=jt  !PM
        zv(:,:)= getvar(cfilev, cvar, jk ,npiglo,npjglo, ktime=jt)
        ! Read mask
        zmask(:,:)=getvar(cfilemask,cvmask,jkk,npiglo,npjglo)
        ! Multiplication of cvar by mask à level jk
        zvmask=zv*zmask
        ! Writing  on the original file                 
        istatus=putvar(cfilev,cvar,jk,npiglo,npjglo,1,1,ktime=jt, ptab=zvmask)
     END DO
  END DO

END PROGRAM cdfmltmask 
