PROGRAM cdfspeed
  !!-------------------------------------------------------------------
  !!               ***  PROGRAM cdfspeed  ***
  !!
  !!  **  Purpose  :  combine u and v to obtains the wind speed
  !!  
  !!  **  Method   :  sqrt(u**2 + v**2)
  !!
  !!
  !! history ;
  !!  Original :  P. Mathiot (Nov. 2007) from cdfmeanvar
  !!-------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !! * Modules used
  USE cdfio

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: jk, ik, jt, ji, jj
  INTEGER   :: narg, iargc, ncout, ierr            !: command line 
  INTEGER   :: npiglo,npjglo,npk,nt                !: size of the domain
  INTEGER   :: nvpk                                !: vertical levels in working variable
  INTEGER, DIMENSION(1) ::  ipk, id_varout         ! 

  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE  :: zu, zv, U

  REAL(kind=4), DIMENSION(:), ALLOCATABLE  :: tim

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
     PRINT *,' Usage : cdfspeed  ncfileU ncfileV cdfvarU cdfvarV' 
     PRINT *,' Computes the speed current or wind'
     PRINT *,' If the input files are 3D, the input is assumed to be '
     PRINT *,' a model output on native C-grid. Speed is computed on the A-grid.'
     PRINT *,' If the input file is 2D and have many time steps, then '
     PRINT *,' we assume that this is a forcing file already on the A-grid.'
     PRINT *,' Output on speed.nc, variable U'
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

  ! define new variables for output
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
  IF (lforcing ) THEN
     ipk(1) = 1 !  2D
     ncout =create(cfileout, cfilev, npiglo,npjglo,0)
     ierr= createvar(ncout ,typvar,1, ipk,id_varout )
     ierr= putheadervar(ncout, cfilev, npiglo, npjglo,0)
  ELSE
     ipk(1)=npk
     ncout =create(cfileout, cfilev, npiglo,npjglo,npk)
     ierr= createvar(ncout ,typvar,1, ipk,id_varout )
     ierr= putheadervar(ncout, cfilev, npiglo, npjglo,npk)
  END IF
  ! Allocate arrays
  ALLOCATE ( zv(npiglo,npjglo), zu(npiglo,npjglo), U(npiglo,npjglo), tim(nt))

  DO jt=1,nt
     tim(jt)=jt
  END DO
  ierr=putvar1d(ncout,tim,nt,'T')

  DO jt = 1,nt
     DO jk = 1,nvpk
        ! Get velocities v at ik
           zv(:,:)= getvar(cfilev, cvarv,jk,npiglo,npjglo,ktime=jt)
           zu(:,:)= getvar(cfileu, cvaru,jk,npiglo,npjglo,ktime=jt)
        IF ( lforcing ) THEN
          ! u and v are already on the T grid points
        ELSE
           ! in this case we are on the C-grid and the speed mus be computed on the A-grid
           DO ji=1,npiglo -1
             DO jj=1,npjglo
               zu(ji,jj)=0.5*(zu(ji,jj)+zu(ji+1,jj))
             ENDDO
           ENDDO
           DO ji=1,npiglo 
             DO jj=1,npjglo-1
               zv(ji,jj)=0.5*(zv(ji,jj)+zv(ji,jj+1))
             ENDDO
           ENDDO
        END IF
        
        U=SQRT(zv*zv+zu*zu)
        ierr = putvar(ncout, id_varout(1) ,U, jk ,npiglo, npjglo, ktime=jt)
     END DO
  END DO
  ierr = closeout(ncout)

END PROGRAM cdfspeed
