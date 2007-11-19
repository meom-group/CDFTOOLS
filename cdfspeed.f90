PROGRAM cdfspeed
  !!-------------------------------------------------------------------
  !!               ***  PROGRAM cdfmeanvar  ***
  !!
  !!  **  Purpose  :  Compute the Mean Value and variance  over the ocean
  !!                  PARTIAL STEPS
  !!  
  !!  **  Method   :  compute the sum ( V * e1 *e2 * e3 *mask )/ sum( e1 * e2 * e3 *mask )
  !!
  !!
  !! history ;
  !!  Original :  J.M. Molines (Oct. 2005) 
  !!              J.M. Molines  Add variance Nov. 2006
  !!-------------------------------------------------------------------
  !!  $Rev: 94 $
  !!  $Date: 2007-09-11 20:34:24 +0200 (Tue, 11 Sep 2007) $
  !!  $Id: cdfmeanvar.f90 94 2007-09-11 18:34:24Z molines $
  !!--------------------------------------------------------------
  !! * Modules used
  USE cdfio

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: jk, ik, jt
  INTEGER   :: narg, iargc, ncout, ierr            !: command line 
  INTEGER   :: npiglo,npjglo,npk,nt                !: size of the domain
  INTEGER   :: nvpk                                !: vertical levels in working variable
  INTEGER, DIMENSION(1) ::  ipk, id_varout         ! 

  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE  :: zu, zv, U

  CHARACTER(LEN=80) :: cfilev, cfileu
  CHARACTER(LEN=80) :: cfileout='speed.nc'
  CHARACTER(LEN=80) :: cvaru, cvarv, cvartype
  CHARACTER(LEN=20) :: ce1, ce2, ce3, cvmask, cvtype, cdep

  LOGICAL    :: lforcing
  INTEGER    :: istatus

  TYPE (variable), DIMENSION(1) :: typvar         !: structure for attibutes
  ! constants

  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfmeanvar  ncfileU ncfileV' 
     PRINT *,' Computes the speed current or wind'
     PRINT *,' Output on speed.nc'
     STOP
  ENDIF

  CALL getarg (1, cfileu)
  CALL getarg (2, cfilev)
  CALL getarg (3, cvaru)
  CALL getarg (4, cvarv)

  IF (narg > 4 ) THEN
       PRINT *, ' ERROR : You must give just fileU and fileV and cvaru and cvarv'
       STOP
  ENDIF

  npiglo= getdim (cfilev,'x')
  npjglo= getdim (cfilev,'y')
  npk   = getdim (cfilev,'depth')
  nvpk  = getvdim(cfilev,cvarv)
  nt    = getdim (cfilev,'time_counter')

  IF (nvpk == 2 ) nvpk = 1
  IF (nvpk == 3 ) nvpk = npk

  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk
  PRINT *, 'nvpk  =', nvpk
  PRINT *, 'nt    =', nt

  lforcing=.FALSE.
  IF ((npk .EQ. 0) .AND. (nt .GT. 1)) THEN
     lforcing=.TRUE.
     npk=1
     PRINT *, 'W A R N I N G : you used a forcing field'
  END IF

  ipk(1) = 1  !  2D
  ! define new variables for output ( must update att.txt)
  typvar(1)%name='U'
  typvar(1)%units='m.s-1'
  typvar(1)%missing_value=0.
  typvar(1)%valid_min= -1000.
  typvar(1)%valid_max= 1000.
  typvar(1)%long_name='Current or wind speed'
  typvar(1)%short_name='U'
  typvar(1)%online_operation='N/A'
  typvar(1)%axis='TZYX'

  ! create output fileset
  ncout =create(cfileout, cfilev, npiglo,npjglo,0)
  ierr= createvar(ncout ,typvar,1, ipk,id_varout )
  ierr= putheadervar(ncout, cfilev, npiglo, npjglo,0)
  ! Allocate arrays
  ALLOCATE ( zv(npiglo,npjglo), zu(npiglo,npjglo), U(npiglo,npjglo))

  DO jt = 1,nt
     DO jk = 1,nvpk
        ! Get velocities v at ik
        IF ( lforcing ) THEN
           ik = jt
           zv(:,:)= getvar(cfilev, cvarv,ik,npiglo,npjglo,ktime=jt)
           zu(:,:)= getvar(cfileu, cvaru,ik,npiglo,npjglo,ktime=jt)
           ik=1
        ELSE
           zv(:,:)= getvar(cfilev,cvarv,ik,npiglo,npjglo,ktime=jt)
           zu(:,:)= getvar(cfileu,cvaru,ik,npiglo,npjglo,ktime=jt)
        END IF
        
        U=SQRT(zv*zv+zu*zu)

        IF (lforcing ) THEN
           PRINT *, jt
           ierr = putvar(ncout, id_varout(1) ,U, jt ,npiglo, npjglo, jt)
        ELSE
           ierr = putvar(ncout, id_varout(1) ,U, jk ,npiglo, npjglo, jt)
        END IF        
     END DO
  END DO
  ierr = closeout(ncout)

END PROGRAM cdfspeed
