PROGRAM cdfmeanvar
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
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !! * Modules used
  USE cdfio

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: jk, ik, jt
  INTEGER   :: imin=0, imax=0, jmin=0, jmax=0      !: domain limitation for computation
  INTEGER   :: kmin=0, kmax=0                      !: domain limitation for computation
  INTEGER   :: ierr                                !: working integer
  INTEGER   :: narg, iargc                         !: command line 
  INTEGER   :: npiglo,npjglo,npk,nt                !: size of the domain
  INTEGER   :: nvpk                                !: vertical levels in working variable

  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::  e1, e2, e3,  zv   !:  metrics, velocity
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::  zmask             !:   npiglo x npjglo
  REAL(KIND=4), DIMENSION (:),     ALLOCATABLE :: gdep               !:  depth 

  REAL(KIND=8)      :: zvol, zsum, zvol2d, zsum2d, zsurf, zvar, zvar2d
  CHARACTER(LEN=256) :: cfilev , cdum
  CHARACTER(LEN=256) :: coordhgr='mesh_hgr.nc',  coordzgr='mesh_zgr.nc',cmask='mask.nc'
  CHARACTER(LEN=256) :: cvar, cvartype
  CHARACTER(LEN=20) :: ce1, ce2, ce3, cvmask, cvtype, cdep

  INTEGER    :: istatus

  ! constants

  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfmeanvar  ncfile cdfvar T| U | V | F | W [imin imax jmin jmax kmin kmax] '
     PRINT *,' Computes the mean value, and the spatial variance of the field (3D, weighted) '
     PRINT *,' imin imax jmin jmax kmin kmax can be given in option '
     PRINT *,'    if imin = 0 then ALL i are taken'
     PRINT *,'    if jmin = 0 then ALL j are taken'
     PRINT *,'    if kmin = 0 then ALL k are taken'
     PRINT *,' PARTIAL CELLS VERSION'
     PRINT *,' Files mesh_hgr.nc, mesh_zgr.nc ,mask.nc '
     PRINT *,'  must be in the current directory'
     PRINT *,' Output on standard output'
     STOP
  ENDIF

  CALL getarg (1, cfilev)
  CALL getarg (2, cvar)
  CALL getarg (3, cvartype)

  IF (narg > 3 ) THEN
    IF ( narg /= 9 ) THEN
       PRINT *, ' ERROR : You must give 6 optional values (imin imax jmin jmax kmin kmax)'
       STOP
    ELSE
    ! input optional imin imax jmin jmax
      CALL getarg ( 4,cdum) ; READ(cdum,*) imin
      CALL getarg ( 5,cdum) ; READ(cdum,*) imax
      CALL getarg ( 6,cdum) ; READ(cdum,*) jmin
      CALL getarg ( 7,cdum) ; READ(cdum,*) jmax
      CALL getarg ( 8,cdum) ; READ(cdum,*) kmin
      CALL getarg ( 9,cdum) ; READ(cdum,*) kmax
    ENDIF
  ENDIF

  npiglo= getdim (cfilev,'x')
  npjglo= getdim (cfilev,'y')
  npk   = getdim (cfilev,'depth')
  nvpk  = getvdim(cfilev,cvar)
  nt    = getdim (cfilev,'time_counter')

  IF (imin /= 0 ) THEN ; npiglo=imax -imin + 1;  ELSE ; imin=1 ; ENDIF
  IF (jmin /= 0 ) THEN ; npjglo=jmax -jmin + 1;  ELSE ; jmin=1 ; ENDIF
  IF (kmin /= 0 ) THEN ; npk   =kmax -kmin + 1;  ELSE ; kmin=1 ; ENDIF

  IF (nvpk == 2 ) nvpk = 1
  IF (nvpk == 3 ) nvpk = npk

  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk
  PRINT *, 'nvpk  =', nvpk
  PRINT *, 'nt    =', nt

  IF ((npk .EQ. 0) .AND. (nt .GT. 1)) THEN
     npk=1
     PRINT *, 'W A R N I N G : you used a forcing field'
  END IF
  ! Allocate arrays
  ALLOCATE ( zmask(npiglo,npjglo) )
  ALLOCATE ( zv(npiglo,npjglo) )
  ALLOCATE ( e1(npiglo,npjglo),e2(npiglo,npjglo), e3(npiglo,npjglo) )
  ALLOCATE ( gdep(npk) )
  SELECT CASE (TRIM(cvartype))
  CASE ( 'T' )
     ce1='e1t'
     ce2='e2t'
     ce3='e3t_ps'
     cvmask='tmask'
     cdep='gdept'
  CASE ( 'U' )
     ce1='e1u'
     ce2='e2u'
     ce3='e3t_ps'
     cvmask='umask'
     cdep='gdept'
  CASE ( 'V' )
     ce1='e1v'
     ce2='e2v'
     ce3='e3t_ps'
     cvmask='vmask'
     cdep='gdept'
  CASE ( 'F' )
     ce1='e1f'
     ce2='e2f'
     ce3='e3t_ps'
     cvmask='fmask'
     cdep='gdept'
  CASE ( 'W' )
     ce1='e1t'
     ce2='e2t'
     ce3='e3w_ps'
     cvmask='tmask'
     cdep='gdepw'
  CASE DEFAULT
      PRINT *, 'this type of variable is not known :', trim(cvartype)
      STOP
  END SELECT

  e1(:,:) = getvar(coordhgr, ce1, 1,npiglo,npjglo,kimin=imin,kjmin=jmin)
  e2(:,:) = getvar(coordhgr, ce2, 1,npiglo,npjglo,kimin=imin,kjmin=jmin)
  gdep(:) = getvare3(coordzgr,cdep,npk)

  zvol=0.d0
  zsum=0.d0
  DO jt = 1,nt
     DO jk = 1,nvpk
        ik = jk+kmin-1
        ! Get velocities v at ik
        zv(:,:)= getvar(cfilev, cvar,  ik ,npiglo,npjglo,kimin=imin,kjmin=jmin, ktime=jt)
        IF ( nvpk /= 1 .OR. jt == 1 ) THEN
          ! if there is only one level do not read mask and e3 every time step ...
          zmask(:,:)=getvar(cmask,cvmask,ik,npiglo,npjglo,kimin=imin,kjmin=jmin)
          ! get e3 at level ik ( ps...)
          e3(:,:) = getvar(coordzgr, ce3, ik,npiglo,npjglo,kimin=imin,kjmin=jmin, ldiom=.true.)
        END IF
        ! 
        zsurf=sum(e1 * e2 * zmask)
        zvol2d=sum(e1 * e2 * e3 * zmask)
        zvol=zvol+zvol2d
        zsum2d=sum(zv*e1*e2*e3*zmask)
        zvar2d=sum(zv*zv*e1*e2*e3*zmask)
        zsum=zsum+zsum2d
        zvar=zvar+zvar2d
        IF (zvol2d /= 0 )THEN
           PRINT *, ' Mean value at level ',ik,'(',gdep(ik),' m) ',zsum2d/zvol2d, 'surface = ',zsurf/1.e6,' km^2 jt=', jt 
           PRINT *, ' Mean value2 at level ',ik,'(',gdep(ik),' m) ',zvar2d/zvol2d, 'variance=', &
                &    zvar2d/zvol2d - (zsum2d/zvol2d)*(zsum2d/zvol2d)
        ELSE
           PRINT *, ' No points in the water at level ',ik,'(',gdep(ik),' m) '
        ENDIF
     END DO
  END DO

  PRINT * ,' Mean value over the ocean: ', zsum/zvol
  PRINT * ,' Global variance over the ocean: ', zvar/zvol - (zsum/zvol)*(zsum/zvol)
  PRINT * ,' Global std dev over the ocean: ', sqrt(zvar/zvol - (zsum/zvol)*(zsum/zvol))
  
END PROGRAM cdfmeanvar
