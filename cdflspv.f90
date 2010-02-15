PROGRAM cdflspv
  !! --------------------------------------------------------------
  !!               ***   PROGRAM CDFLSPV ***
  !!  ** Purpose:  This program is used to compute the 
  !!               large scale potential vorticity
  !!               from a set of T S  files.
  !!
  !!  ** Method:   pv = 1/rho0 *  f * d(rho)/d(z)
  !!                rho0 = 1020. kg/m3
  !!                f is the coriolis factor
  !!                zeta is the relative vorticity
  !!               Output is done for f (2D) (at f-points)
  !!                                 f/rho0 d(rho)/d(z) (3D) at W points
  !!
  !!  ** Usage :
  !!         cdfpv gridT gridU gridV files.
  !!             output is done on pv.nc, with variable name
  !!                volspv  (PV) 
  !!
  !!  * history:
  !!        Original : J.M. Molines for SPEM in Dynamo (1996)
  !!        Modif    : J-O. Beismann for OPA (1999)
  !!        Modif    : J.M. Molines for normalization Clipper (March 2000)
  !!                 : J.M. Molines in cdftools, f90 dor DRAKKAR (Nov. 2005)
  !! ---------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !! * Used modules
  USE cdfio
  USE eos

  !! * Local declaration
  IMPLICIT NONE

  INTEGER :: npiglo, npjglo, npk, npt  
  INTEGER :: narg, iargc
  INTEGER :: ji,jj,jk, jt
  INTEGER :: ncout, ierr
  INTEGER :: iup=1 , idown=2, itmp
  INTEGER, DIMENSION(1) :: ipk, id_varout !: for output variables
  !
  REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE ::  sigma
  REAL(KIND=4), DIMENSION(:,:)  , ALLOCATABLE ::  ztemp, zsal, dsig, zmask, fcorio, pv,&
       &                                          e3w, gphit
  REAL(KIND=4), DIMENSION(:)    , ALLOCATABLE ::  time_tag, h1d, gdepw
  REAL(KIND=4)                                ::  zrot, pi, rho0=1020.

  CHARACTER(LEN=256) ::   cfilet,  cfilout
  CHARACTER(LEN=256) ::  coordhgr='mesh_hgr.nc', coordzgr='mesh_zgr.nc'

  TYPE(variable) , DIMENSION(1)   :: typvar    !: structure for attributes
  !

  !! * Read command line
  narg=iargc()
  IF (narg == 0 ) THEN
     PRINT *, &
          &' >>>> usage: cdflspv gridT  files '
     PRINT *,'   Output is done on lspv.nc'
     PRINT *,'   variables  volspv '
     PRINT *,'  mesh_hgr.nc, mesh_zgr.nc are required'
     STOP
  ENDIF
  CALL getarg(1,cfilet)

  npiglo=getdim(cfilet,'x')
  npjglo=getdim(cfilet,'y')
  npk   =getdim(cfilet,'depth')
  npt   =getdim(cfilet,'time')

  ALLOCATE( sigma(npiglo,npjglo,2) )
  ALLOCATE( ztemp(npiglo,npjglo), zsal(npiglo,npjglo)  )
  ALLOCATE (fcorio(npiglo,npjglo),pv(npiglo,npjglo) )
  ALLOCATE( zmask(npiglo,npjglo), dsig(npiglo,npjglo)  )
  ALLOCATE( time_tag(npt), h1d(npk) ,gdepw(npk))
  ALLOCATE ( gphit(npiglo,npjglo))
  ALLOCATE ( e3w(npiglo,npjglo) )


  ! read mesh_mask/ time information
  time_tag(:)=getvar1d(cfilet,'time_counter', npt)
  h1d(:)=getvar1d(cfilet,'deptht',npk)
  gdepw(:) = getvare3(coordzgr, 'gdepw',npk)

  gphit(:,:) = getvar(coordhgr,'gphit',1,npiglo,npjglo)
  
  ! Compute coriolis factor
  pi=ACOS(-1.)
  fcorio(:,:)=4*pi/86400.*ABS(SIN(pi/180*gphit(:,:)))

  ! ... open output file and write header
  ipk(:)=npk
  typvar(1)%name= 'volspv'
  typvar(1)%units='kg.m-4.s-1 x 1e7'
  typvar(1)%missing_value=0.
  typvar(1)%valid_min= -1000.
  typvar(1)%valid_max= 1000.
  typvar(1)%long_name='Large Scale Potential_vorticity'
  typvar(1)%short_name='volspv'
  typvar(1)%online_operation='N/A'
  typvar(1)%axis='TZYX'

  cfilout='lspv.nc'

  ncout = create(cfilout,'none' ,npiglo,npjglo,npk,cdep='depthw')
  ierr = createvar(ncout, typvar,1,ipk, id_varout )
  ierr = putheadervar(ncout , cfilet, npiglo, npjglo, npk,pdep=gdepw)
  ierr = putvar1d(ncout,time_tag,npt,'T')

  DO jt=1,npt
     PRINT *, 'time ',jt
  ! surface PV is unknown ...
  pv(:,:) = 0.
  ierr = putvar(ncout,id_varout(1), pv,1,npiglo,npjglo,ktime=jt)

  ! initialize first level
  ztemp(:,:) =   getvar(cfilet,'votemper',1,npiglo,npjglo,ktime=jt)
  zsal(:,:)  =   getvar(cfilet,'vosaline',1,npiglo,npjglo,ktime=jt)

  zmask = 1.0
  WHERE(zsal == 0 ) zmask = 0.0
  sigma(:,:,iup) =  sigma0 ( ztemp,zsal,npiglo,npjglo )* zmask(:,:)

  ! Main vertical loop
  DO jk=2,npk
     ztemp(:,:) =   getvar(cfilet,'votemper',jk,npiglo,npjglo,ktime=jt)
     zsal(:,:)  =   getvar(cfilet,'vosaline',jk,npiglo,npjglo,ktime=jt)
     e3w (:,:)  =   getvar(coordzgr,'e3w_ps', jk, npiglo,npjglo, ldiom=.true.)
     WHERE (e3w == 0 ) e3w = 1.

     zmask=1.0
     WHERE(zsal == 0 ) zmask = 0.0
     sigma(:,:,idown) =  sigma0 ( ztemp,zsal,npiglo,npjglo )* zmask(:,:)

     !  d(sigma0)/dz at W point ( masked if down level is masked )
     dsig(:,:)=(sigma(:,:,idown) - sigma(:,:,iup)) /e3w *zmask

     ! Full pv:
     DO ji=1,npiglo 
        DO jj = 1, npjglo
           pv(ji,jj) = (fcorio(ji,jj))*dsig(ji,jj)*1.e7
        END DO
     END DO
     ierr = putvar(ncout,id_varout(1), pv,jk,npiglo,npjglo,ktime=jt)
     
     ! swap index up and down
     itmp=iup
     iup=idown
     idown=itmp
  END DO   ! level loop
  END DO   ! time loop

  ierr = closeout(ncout)
  PRINT *,'cdflspv completed successfully'
END  PROGRAM cdflspv
