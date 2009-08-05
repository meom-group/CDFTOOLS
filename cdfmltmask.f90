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
  INTEGER   :: npkmask                             !: vertical levels in mask file

  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::  zv !: cvar at jk level 
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::  zmask !: mask at jk level 
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE :: zvmask !: masked cvar at jk level

  CHARACTER(LEN=256) :: cfilev , cfilemask, ctmp
  CHARACTER(LEN=256) :: cvar, cvartype, cdep
  CHARACTER(LEN=20) ::  cvmask

  INTEGER    :: istatus

  ! constants

  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfmltmask ncfile maskfile cdfvar T| U | V | F | W | P'
     PRINT *,' Mask the file  '
     PRINT *,' output on ncfile_masked'
     PRINT *,'  Point type P correspond to polymask'
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
  npkmask = getdim (cfilemask,'depth',cdtrue=cdep, kstatus=istatus)
  IF (istatus /= 0 ) THEN
     npkmask  = getdim (cfilemask,'z',cdtrue=cdep,kstatus=istatus) ; !print *, istatus
        IF ( istatus /= 0 ) THEN
          PRINT *,' assume mask file with no depth'
          npkmask=0
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

  IF (npk==0) npk=1

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
  CASE ( 'P' )   ! for polymask 
     cvmask='polymask'
  CASE DEFAULT
     PRINT *, 'this type of variable is not known :', TRIM(cvartype)
     STOP
  END SELECT

  IF ( npkmask <= 1 ) THEN 
        zmask(:,:)=getvar(cfilemask,cvmask,1,npiglo,npjglo)
  ENDIF
  DO jt = 1, npt
     IF (MOD(jt,100)==0) PRINT *, jt,'/', npt
     DO jk = 1,nvpk
        ! Read cvar
        zv(:,:)= getvar(cfilev, cvar, jk ,npiglo,npjglo, ktime=jt)
        IF ( npkmask > 1 ) THEN
        ! Read mask
          zmask(:,:)=getvar(cfilemask,cvmask,jk,npiglo,npjglo)
        ENDIF
        ! Multiplication of cvar by mask at level jk
        zvmask=zv*zmask
        ! Writing  on the original file                 
        istatus=putvar(cfilev,cvar,jk,npiglo,npjglo,1,1,ktime=jt, ptab=zvmask)
     END DO
  END DO

END PROGRAM cdfmltmask 
