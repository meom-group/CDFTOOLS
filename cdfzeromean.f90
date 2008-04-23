PROGRAM cdfzeromean
  !!-------------------------------------------------------------------
  !!               ***  PROGRAM cdfzeromean  ***
  !!
  !!  **  Purpose  :  Compute the Mean Value over the ocean
  !!                  Produce a file with a 'zeromean' variable
  !!                  PARTIAL STEPS
  !!  
  !!  **  Method   :  compute the sum ( V * e1 *e2 * e3 *mask )/ sum( e1 * e2 * e3 *mask ) 
  !!                  The mean( 3D) value is rested from the initial field
  !!
  !!
  !! history ;
  !!  Original :  J.M. Molines (Oct. 2005) 
  !!-------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !! * Modules used
  USE cdfio

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: jk, ik, jt, ivar                    !: dummy loop index
  INTEGER   :: imin=0, imax=0, jmin=0, jmax=0      !: domain limitation for computation
  INTEGER   :: kmin=0, kmax=0                      !: domain limitation for computation
  INTEGER   :: ierr                                !: working integer
  INTEGER   :: narg, iargc                         !: command line 
  INTEGER   :: npiglo,npjglo, npk, nt              !: size of the domain
  INTEGER   :: npiglo_fi,npjglo_fi
  INTEGER   :: nvpk                                !: vertical levels in working variable

  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::  e1, e2, e3,  zv   !:  metrics, velocity
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::  zmask             !:   npiglo x npjglo
  REAL(KIND=4), DIMENSION (:),     ALLOCATABLE ::  gdep              !:  depth 

  REAL(KIND=8)      :: zvol, zsum, zvol2d, zsum2d, zsurf, zmean !: double precision cumul/mean
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: zmean2d            !: per level mean

  CHARACTER(LEN=80) :: cfilev , cdum, cfileout='zeromean.nc'    
  CHARACTER(LEN=80) :: coordhgr='mesh_hgr.nc',  coordzgr='mesh_zgr.nc',cmask='mask.nc'
  CHARACTER(LEN=80) :: cvar, cvartype
  CHARACTER(LEN=20) :: ce1, ce2, ce3, cvmask, cvtype, cdep

  ! output stuff variables
  INTEGER :: ncout, nvars  ! number of variables in the input file
  INTEGER, DIMENSION(1) ::  ipk, id_varout
  REAL(KIND=4), DIMENSION(:), ALLOCATABLE :: tim, dep
  TYPE (variable), DIMENSION(1) :: typvar         !: structure for attibutes
  TYPE (variable), DIMENSION(:),ALLOCATABLE :: typvarin         !: structure for attibutes
  CHARACTER(LEN=80), DIMENSION(:), ALLOCATABLE :: cvarname
  LOGICAL :: lnodep=.FALSE.

  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfzeromean  ncfile cdfvar T| U | V | F | W [imin imax jmin jmax kmin kmax] '
     PRINT *,' Computes the mean value of the field (3D, weighted) '
     PRINT *,' and return a ncdf file with the variable (field - mean) '
     PRINT *,' imin imax jmin jmax kmin kmax can be given in option '
     PRINT *,'    if imin = 0 then ALL i are taken'
     PRINT *,'    if jmin = 0 then ALL j are taken'
     PRINT *,'    if kmin = 0 then ALL k are taken'
     PRINT *,' PARTIAL CELLS VERSION'
     PRINT *,' Files mesh_hgr.nc, mesh_zgr.nc ,mask.nc '
     PRINT *,'  must be in the current directory'
     PRINT *,' Output on standard output and on zeromean.nc file, variable same as input'
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

  ! get dimensions from input file
  npiglo= getdim (cfilev,'x')
  npjglo= getdim (cfilev,'y')
  npk   = getdim (cfilev,'depth')
  nt    = getdim (cfilev,'time')
  nvpk  = getvdim(cfilev,cvar)
  IF (npk == 0 ) THEN ;  lnodep=.TRUE. ; npk = 1 ; ENDIF  ! no depth dimension ==> 1 level
  ! save original npiglo, npiglo
  npiglo_fi=npiglo
  npjglo_fi=npjglo

  IF (imin /= 0 ) THEN ; npiglo=imax -imin + 1;  ELSE ; imin=1 ; ENDIF
  IF (jmin /= 0 ) THEN ; npjglo=jmax -jmin + 1;  ELSE ; jmin=1 ; ENDIF
  IF (kmin /= 0 ) THEN ; npk   =kmax -kmin + 1;  ELSE ; kmin=1 ; ENDIF

  IF (nvpk == 2 ) nvpk = 1
  IF (nvpk == 3 ) nvpk = npk

  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk
  PRINT *, 'nt    =', nt
  PRINT *, 'nvpk  =', nvpk

  ! Allocate arrays
  ALLOCATE ( zmask(npiglo,npjglo) )
  ALLOCATE ( zv(npiglo,npjglo) )
  ALLOCATE ( e1(npiglo,npjglo),e2(npiglo,npjglo), e3(npiglo,npjglo) )
  ALLOCATE ( gdep (npk) ,zmean2d(nvpk) )
  ALLOCATE ( tim(nt) ,dep (nvpk))

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
     PRINT *, 'this type of variable is not known :', TRIM(cvartype)
     STOP
  END SELECT

  e1(:,:) = getvar(coordhgr, ce1, 1,npiglo,npjglo,kimin=imin,kjmin=jmin)
  e2(:,:) = getvar(coordhgr, ce2, 1,npiglo,npjglo,kimin=imin,kjmin=jmin)
  gdep(:) = getvare3(coordzgr,cdep,npk)

  DO jt=1,nt
     zvol=0.d0
     zsum=0.d0
     DO jk = 1,nvpk
        ik = jk+kmin-1
        dep(ik)=gdep(jk)
        ! Get velocities v at ik
        zv(:,:)= getvar(cfilev, cvar,  ik ,npiglo,npjglo,kimin=imin,kjmin=jmin)
        !     zv(:,:)= getvar(cfilev, cvar,  jt ,npiglo,npjglo,kimin=imin,kjmin=jmin,ktime=jt)
        zmask(:,:)=getvar(cmask,cvmask,ik,npiglo,npjglo,kimin=imin,kjmin=jmin)
        !    zmask(:,npjglo)=0.

        ! get e3 at level ik ( ps...)
        e3(:,:) = getvar(coordzgr, ce3, ik,npiglo,npjglo,kimin=imin,kjmin=jmin, ldiom=.TRUE.)

        ! 
        zsurf=SUM(e1 * e2 * zmask)
        zvol2d=SUM(e1 * e2 * e3 * zmask)
        zvol=zvol+zvol2d
        zsum2d=SUM(zv*e1*e2*e3*zmask)
        zsum=zsum+zsum2d
        IF (zvol2d /= 0 )THEN
           PRINT *, ' Mean value at level ',ik,'(',gdep(ik),' m) ',zsum2d/zvol2d, 'surface = ',zsurf/1.e6,' km^2'
           zmean2d(ik) = zsum2d/zvol2d
        ELSE
           PRINT *, ' No points in the water at level ',ik,'(',gdep(ik),' m) '
        ENDIF

     END DO
     zmean=zsum/zvol
     PRINT * ,' Mean value over the ocean: ', zmean, jt
  END DO
  DEALLOCATE ( zv, zmask)
  npiglo=npiglo_fi ; npjglo=npjglo_fi
  ALLOCATE (zv(npiglo,npjglo), zmask(npiglo,npjglo) )
  ! re-read file and rest mean value from the variable and store on file
  nvars = getnvar(cfilev)
  ALLOCATE ( typvarin(nvars), cvarname(nvars) )
  cvarname(:) = getvarname(cfilev,nvars,typvarin)
  ! look for the working variable
  DO ivar = 1, nvars
     IF ( TRIM(cvarname(ivar)) == TRIM(cvar) ) EXIT
  END DO

  typvar(1)%name= cvar
  typvar%units=typvarin(ivar)%units
  typvar%missing_value=typvarin(ivar)%missing_value
  typvar%valid_min=typvarin(ivar)%valid_min-zmean
  typvar%valid_max=typvarin(ivar)%valid_max-zmean
  typvar(1)%long_name=typvarin(ivar)%long_name//' zero mean '
  typvar(1)%short_name=cvar
  typvar%online_operation='N/A'
  typvar%axis=typvarin(ivar)%axis
  ipk(1)=nvpk

  ik=nvpk
  IF ( lnodep ) ik = 0  ! no depth variable in input file : the same in output file
  ncout = create(cfileout, cfilev, npiglo,npjglo,ik)
  ierr = createvar(ncout ,typvar ,1, ipk,id_varout )
  ierr= putheadervar(ncout, cfilev,npiglo, npjglo,ik,pdep=dep)
  tim=getvar1d(cfilev,'time_counter',nt)

  DO jt=1,nt
     DO jk = 1,nvpk
        ik = jk+kmin-1
        ! Get velocities v at ik
        zv(:,:)= getvar(cfilev, cvar,  ik ,npiglo,npjglo,kimin=imin,kjmin=jmin)
        !     zv(:,:)= getvar(cfilev, cvar,  jt ,npiglo,npjglo,kimin=imin,kjmin=jmin,ktime=jt)
        zmask(:,:)=getvar(cmask,cvmask,ik,npiglo,npjglo,kimin=imin,kjmin=jmin)
        !    zmask(:,npjglo)=0.
        WHERE (zmask /= 0 ) zv=zv - zmean
        ierr = putvar(ncout, id_varout(1) ,zv, ik,npiglo, npjglo )
     END DO
  END DO
  ierr=putvar1d(ncout,tim,nt,'T')
  ierr=closeout(ncout)

END PROGRAM cdfzeromean
