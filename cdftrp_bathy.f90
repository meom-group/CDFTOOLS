PROGRAM cdftrp_bathy
  !!-------------------------------------------------------------------
  !!              PROGRAM cdftrp_bathy
  !!              ********************
  !!
  !!  **  Purpose: Compute  vertically integrated transport components,
  !!               along bathymetry( horizontal) and across bathy
  !!  
  !!  **  Method: Use output from cdfvtrp
  !!
  !! history:
  !!    Original:  P. Mathiot 2008.
  !!-------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !!
  !! * Modules used
  USE cdfio 

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: ji,jj, jim1,jip1,jjm1,jjp1
  INTEGER   :: narg, iargc                                  !: 
  INTEGER   :: npiglo,npjglo, npk                                !: size of the domain
  INTEGER, DIMENSION(:),ALLOCATABLE ::  ipk, id_varout
  TYPE(variable), DIMENSION(:), ALLOCATABLE :: typvar
  REAL(KIND=8) , DIMENSION (:,:), ALLOCATABLE :: zu, zv, u, v, vmod, hdept
  REAL(KIND=8) , DIMENSION (:,:), ALLOCATABLE :: e1u, e2v, dhdx, dhdy, alpha, tmask
  REAL(KIND=4) ,DIMENSION(1)                  :: timean

  CHARACTER(LEN=256) :: cfile ,cfilev, cfilew,  cfilet, cfileout='trpiso.nc'            !: file name

  INTEGER    :: ncout
  INTEGER    :: ierr

  !!  Read command line
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdftrp_bathy trp.nc (build with cdfvtrp)'
     PRINT *,' need mask.nc, mesh_hgr.nc and hdept.nc in the current directory '
     PRINT *,' (*) hdept.nc is a file with only hdept (2D) variable extracted from'
     PRINT *,'     mesh_zgr.nc. It can be a lonk to mesh_zgr'
     STOP
  ENDIF
  !!
  !! Initialisation from 1st file (all file are assume to have the same geometry)
  CALL getarg (1, cfile)

  ALLOCATE ( ipk(2), id_varout(2), typvar(2) )

  npiglo = getdim (cfile,'x')
  npjglo = getdim (cfile,'y')
  npk    = getdim (cfile,'depth')

  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk

  ALLOCATE( zu(npiglo,npjglo),  zv(npiglo,npjglo)  )
  ALLOCATE( u(npiglo,npjglo), v(npiglo,npjglo) )
  ALLOCATE( e1u(npiglo,npjglo), e2v(npiglo,npjglo), hdept(npiglo,npjglo), tmask(npiglo,npjglo))
  ALLOCATE( dhdy(npiglo,npjglo), dhdx(npiglo,npjglo), alpha(npiglo,npjglo) )

  ipk(1)      = npk
  typvar(1)%name='soualz'
  typvar(1)%units='m3/s'
  typvar(1)%missing_value=0.
  typvar(1)%valid_min= -10000000.
  typvar(1)%valid_max= 10000000.
  typvar(1)%long_name='Transport at T point along isodepth'
  typvar(1)%short_name='soualz'
  typvar(1)%online_operation='N/A'
  typvar(1)%axis='TYX'

  ipk(2)      = npk
  typvar(2)%name='sovacz'
  typvar(2)%units='m3/s'
  typvar(2)%missing_value=0.
  typvar(2)%valid_min= -10000000.
  typvar(2)%valid_max= 10000000.
  typvar(2)%long_name='Transport at T point across isodepth'
  typvar(2)%short_name='sovacz'
  typvar(2)%online_operation='N/A'
  typvar(2)%axis='TYX'

  ncout =create(cfileout, cfile,npiglo,npjglo,1)
  ierr= createvar(ncout ,typvar,2, ipk,id_varout )
  ierr= putheadervar(ncout, cfile, npiglo, npjglo,1)

  zu(:,:) = getvar(cfile,'sozoutrp',1 ,npiglo, npjglo)
  zv(:,:) = getvar(cfile,'somevtrp',1 ,npiglo, npjglo)

  ! put velocity components on T points
  u(:,:) = 0. ; v(:,:)=0
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

  ! compute bathymetric gradient
  DO jj = 1, npjglo 
     DO ji = 1, npiglo   
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

  ! Compute  the angle between the bathymetric slope and model coordinates.
  zv=v*tmask
  zu=u*tmask

  alpha=ATAN2(dhdy,dhdx)*tmask!*180/3.14159*tmask

  u=(zu*COS(alpha)+zv*SIN(alpha))*tmask      ! transport accross isoline (oriented from shelf to abyssal plain)
  v=-(-zu*SIN(alpha)+zv*COS(alpha))*tmask    ! transport along  isoline (oriented at right of u
  PRINT *, 'iso : ',MAXVAL(SQRT(u**2+v**2)), MAXVAL(u), MAXVAL(v)
  PRINT *, 'normal : ',MAXVAL(SQRT(zu**2+zv**2)), MAXVAL(zu), MAXVAL(zv)
  ierr=putvar(ncout,id_varout(1), REAL(v), 1, npiglo, npjglo)
  ierr=putvar(ncout,id_varout(2), REAL(u), 1, npiglo, npjglo)
  PRINT *, ' SUM DRAKE : ', SUM(SQRT(u(437,51:118)**2))
  PRINT *, ' SUM DRAKE : ', SUM(SQRT(v(437,51:118)**2))

  PRINT *, ' SUM DRAKE : ', SUM(SQRT(u(443,68:118)**2))
  PRINT *, ' SUM DRAKE : ', SUM(SQRT(v(443,68:118)**2))

  timean = getvar1d(cfile,'time_counter',1)
  ierr = putvar1d(ncout,timean,1,'T')
  ierr = closeout(ncout)

END PROGRAM cdftrp_bathy
