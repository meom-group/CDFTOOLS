PROGRAM cdftrp_bathy
  !!-------------------------------------------------------------------
  !!              PROGRAM CDFVITA
  !!              **************
  !!
  !!  **  Purpose: Compute surface velocity on t grid
  !!                 gridU ,  gridV   gridT (reference)
  !!  
  !!  **  Method: Try to avoid 3 d arrays 
  !!
  !! history:
  !!    Original:  J.M. Molines (Nov 2006 ) for ORCA025
  !!-------------------------------------------------------------------
  !!  $Rev: 175 $
  !!  $Date: 2008-03-25 10:23:47 +0100 (Tue, 25 Mar 2008) $
  !!  $Id: cdfvita.f90 175 2008-03-25 09:23:47Z molines $
  !!--------------------------------------------------------------
  !!
  !! * Modules used
  USE cdfio 

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: ji,jj,jk, jim1,jip1,jjm1,jjp1
  INTEGER   :: narg, iargc                                  !: 
  INTEGER   :: npiglo,npjglo, npk                                !: size of the domain
  INTEGER, DIMENSION(:),ALLOCATABLE ::  ipk, id_varout
  TYPE(variable), DIMENSION(:), ALLOCATABLE :: typvar
  REAL(KIND=8) , DIMENSION (:,:), ALLOCATABLE :: zu, zv, u, v, vmod, hdept, e1u, e2v, dhdx, dhdy, alpha, tmask
  REAL(KIND=4) ,DIMENSION(1)                  :: timean

  CHARACTER(LEN=80) :: cfile ,cfilev, cfilew,  cfilet, cfileout='trpiso.nc'            !: file name

  INTEGER    :: ncout
  INTEGER    :: istatus, ierr

  !!  Read command line
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdftrp_bathy trp.nc (build with cdfvtrp)'
     PRINT *,' need mask;nc, mesh_zgr.nc and mesh_hgr.nc in the current directory '
     STOP
  ENDIF
  !!
  !! Initialisation from 1st file (all file are assume to have the same geometry)
  CALL getarg (1, cfile)

  ALLOCATE ( ipk(2), id_varout(2), typvar(2) )

  npiglo = getdim (cfile,'x')
  npjglo = getdim (cfile,'y')
  npk    = getdim (cfile,'depth')
  !ipk(1)=1
  !ipk(2)=1
  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk

  ALLOCATE( zu(npiglo,npjglo),  zv(npiglo,npjglo)  )
  ALLOCATE( u(npiglo,npjglo), v(npiglo,npjglo) )
  ALLOCATE( e1u(npiglo,npjglo), e2v(npiglo,npjglo), hdept(npiglo,npjglo), tmask(npiglo,npjglo))
  ALLOCATE( dhdy(npiglo,npjglo), dhdx(npiglo,npjglo), alpha(npiglo,npjglo) )
  
  ipk(1)      = npk
  typvar(1)%name='soualz'
  typvar(1)%units='m/s'
  typvar(1)%missing_value=0.
  typvar(1)%valid_min= -10000000.
  typvar(1)%valid_max= 10000000.
  typvar(1)%long_name='Velocity T point along isodepth'
  typvar(1)%short_name='soualz'
  typvar(1)%online_operation='N/A'
  typvar(1)%axis='TYX'

  ipk(2)      = npk
  typvar(2)%name='sovacz'
  typvar(2)%units='m/s'
  typvar(2)%missing_value=0.
  typvar(2)%valid_min= -10000000.
  typvar(2)%valid_max= 10000000.
  typvar(2)%long_name='Velocity T point across isodepth'
  typvar(2)%short_name='sovacz'
  typvar(2)%online_operation='N/A'
  typvar(2)%axis='TYX'

  ncout =create(cfileout, cfile,npiglo,npjglo,1)
  ierr= createvar(ncout ,typvar,2, ipk,id_varout )
  ierr= putheadervar(ncout, cfile, npiglo, npjglo,1)

    zu(:,:) = getvar(cfile,'sozoutrp',1 ,npiglo, npjglo)
    zv(:,:) = getvar(cfile,'somevtrp',1 ,npiglo, npjglo)


    u = 0. ; v = 0. ; u(:,:) = 0. ; v(:,:)=0
    DO ji=1, npiglo
      DO jj=1,npjglo
         jjm1=jj-1
         jim1=ji-1
         IF (jj-1 == 0     ) jjm1=npjglo
         IF (ji-1 == 0     ) jim1=npiglo
        u(ji,jj) = 0.5* (zu(ji,jj)+ zu(jim1,jj))
        v(ji,jj) = 0.5* (zv(ji,jj)+ zv(ji,jjm1))
      END DO
    END DO
    tmask(:,:) = getvar('mask.nc','tmask',1 ,npiglo, npjglo)
    hdept(:,:) = getvar('hdept.nc','hdept',1 ,npiglo, npjglo)
    e1u(:,:) = getvar('mesh_hgr.nc','e1u',1 ,npiglo, npjglo)
    e2v(:,:) = getvar('mesh_hgr.nc','e2v',1 ,npiglo, npjglo)
    PRINT *, '',MAXVAL(hdept)
    hdept=hdept*tmask
    !DO jj = 1, npjglo
    !   DO ji = 1, npiglo
    !      hdept(ji,jj)=jj+ji
    !!   END DO
    !END DO
    !u=1
    !v=0
    !hdept=
     DO jj = 1, npjglo 
        DO ji = 1, npiglo   ! vector opt.
           jjm1=jj-1
           jim1=ji-1
           jjp1=jj+1
           jip1=ji+1
           IF (jj-1 == 0     ) jjm1=npjglo
           IF (ji-1 == 0     ) jim1=npiglo
           IF (jj+1 == npjglo) jjp1=1
           IF (ji+1 == npiglo) jip1=1
           dhdx(ji,jj)=(hdept(jip1,jj  )-hdept(jim1,jj  ))/(e1u(ji,jj) +e1u(jim1,jj  ))*tmask(ji,jj)
           dhdy(ji,jj)=(hdept(ji  ,jjp1)-hdept(ji  ,jjm1))/(e2v(ji,jj) +e2v(ji  ,jjm1))*tmask(ji,jj)
        END DO
     END DO

     !!!! calcul de l'angle entre la ligne de plus grande pente et le repere u,v
     zv=v*tmask
     zu=u*tmask
     !dhdy=-1
     !dhdx=-1
     alpha=atan2(dhdy,dhdx)*tmask!*180/3.14159*tmask
     !alpha=0.
     !!!! rotation
     PRINT *, 'atan2(1,1) = ',atan2(1.,1.)
     PRINT *, 'cos(60) = ',cos(3.14159/6.)
     u=(zu*cos(alpha)+zv*sin(alpha))*tmask    !!!vitesse across isoline (oriented from shelf to abyssal plain)
     v=-(-zu*sin(alpha)+zv*cos(alpha))*tmask    !!!vitesse along  isoline (oriented at right of u
     PRINT *, 'iso : ',MAXVAL(sqrt(u**2+v**2)), MAXVAL(u), MAXVAL(v)
     PRINT *, 'normal : ',MAXVAL(sqrt(zu**2+zv**2)), MAXVAL(zu), MAXVAL(zv)
    ierr=putvar(ncout,id_varout(1), REAL(v), 1, npiglo, npjglo)
    ierr=putvar(ncout,id_varout(2), REAL(u), 1, npiglo, npjglo)
    PRINT *, ' SUM DRAKE : ', SUM(sqrt(u(437,51:118)**2))
    PRINT *, ' SUM DRAKE : ', SUM(sqrt(v(437,51:118)**2))

    PRINT *, ' SUM DRAKE : ', SUM(sqrt(u(443,68:118)**2))
    PRINT *, ' SUM DRAKE : ', SUM(sqrt(v(443,68:118)**2))
! END DO

    timean=getvar1d(cfile,'time_counter',1)
    ierr=putvar1d(ncout,timean,1,'T')
    istatus = closeout(ncout)

END PROGRAM cdftrp_bathy
