PROGRAM cdfpv
  !! --------------------------------------------------------------
  !!               ***   PROGRAM CDFPV ***
  !!  ** Purpose:  This program is used to compute the potential vorticity
  !!               from a set of T S U V files.
  !!
  !!  ** Method:   pv = 1/rho0 * ( f + zeta) d(rho)/d(z)
  !!                rho0 = 1020. kg/m3
  !!                f is the coriolis factor
  !!                zeta is the relative vorticity
  !!               Output is done for f (2D) (at f-points)
  !!                                 zeta (3D) at f-points
  !!                                 f/rho0 d(rho)/d(z) (3D) at W points
  !!                                 PV at T point.
  !!
  !!  ** Usage :
  !!         cdfpv gridT gridU gridV files.
  !!             output is done on pv.nc, with variable name
  !!                vopv  (PV) 
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
  INTEGER :: ji,jj,jk
  INTEGER :: ncout, ierr
  INTEGER :: iup=1 , idown=2, itmp
  INTEGER, DIMENSION(1) :: ipk, id_varout !: for output variables
  !
  REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE ::  sigma, rotn
  REAL(KIND=4), DIMENSION(:,:)  , ALLOCATABLE ::  ztemp, zsal,un, vn, dsig, rot, fmask, zmask, fcorio, pv,&
       &                                          e1u, e2f, e1f, e2v, e3w, gphit
  REAL(KIND=4), DIMENSION(:)    , ALLOCATABLE ::  time_tag, h1d, gdepw
  REAL(KIND=4)                                ::  zrot, pi, rho0=1020.

  CHARACTER(LEN=80) ::   cfilet,cfileu, cfilev,  cfilout
  CHARACTER(LEN=80) ::  coordhgr='mesh_hgr.nc', coordzgr='mesh_zgr.nc'

  TYPE(variable) , DIMENSION(1)   :: typvar    !: structure for attributes
  !

  !! * Read command line
  narg=iargc()
  IF (narg < 3 ) THEN
     PRINT *, &
          &' >>>> usage: cdfpv gridT gridU gridV files '
     PRINT *,'   Output is done on pv.nc'
     PRINT *,'   variables  vopv '
     PRINT *,'  mesh_hgr.nc, mesh_zgr.nc are required'
     STOP
  ENDIF
  CALL getarg(1,cfilet)
  CALL getarg(2,cfileu)
  CALL getarg(3,cfilev)

  npiglo=getdim(cfilet,'x')
  npjglo=getdim(cfilet,'y')
  npk   =getdim(cfilet,'depth')
  npt   =getdim(cfilet,'time')

  ALLOCATE( sigma(npiglo,npjglo,2) )
  ALLOCATE ( rotn(npiglo,npjglo,2)  )
  ALLOCATE( ztemp(npiglo,npjglo), zsal(npiglo,npjglo)  )
  ALLOCATE ( un(npiglo,npjglo)  , vn(npiglo,npjglo),fcorio(npiglo,npjglo),pv(npiglo,npjglo) )
  ALLOCATE( zmask(npiglo,npjglo), fmask(npiglo,npjglo) ,dsig(npiglo,npjglo), rot(npiglo,npjglo) )
  ALLOCATE( time_tag(npt), h1d(npk) ,gdepw(npk))
  ALLOCATE ( e1u(npiglo,npjglo) , e1f(npiglo,npjglo) ,gphit(npiglo,npjglo))
  ALLOCATE ( e2v(npiglo,npjglo) , e2f(npiglo,npjglo) ,e3w(npiglo,npjglo) )


  ! read mesh_mask/ time information
  time_tag(:)=getvar1d(cfilet,'time_counter', npt)
  h1d(:)=getvar1d(cfilet,'deptht',npk)
  gdepw(:) = getvare3(coordzgr, 'gdepw',npk)

  e1u=  getvar(coordhgr, 'e1u', 1,npiglo,npjglo)
  e1f=  getvar(coordhgr, 'e1f', 1,npiglo,npjglo)
  e2v=  getvar(coordhgr, 'e2v', 1,npiglo,npjglo)
  e2f=  getvar(coordhgr, 'e2f', 1,npiglo,npjglo)
  gphit(:,:) = getvar(coordhgr,'gphit',1,npiglo,npjglo)
  
  ! Compute coriolis factor
  pi=ACOS(-1.)
  fcorio(:,:)=4*pi/86400.*SIN(pi/180*gphit(:,:))

  ! ... open output file and write header
  ipk(:)=npk
  typvar(1)%name= 'vopv'
  typvar(1)%units='kg.m-4.s-1 x 1e7'
  typvar(1)%missing_value=0.
  typvar(1)%valid_min= -1000.
  typvar(1)%valid_max= 1000.
  typvar(1)%long_name='Full_Potential_vorticity'
  typvar(1)%short_name='vopv'
  typvar(1)%online_operation='N/A'
  typvar(1)%axis='TZYX'

  cfilout='pv.nc'

  ncout = create(cfilout,cfilet ,npiglo,npjglo,npk,cdep='depthw')
  ierr = createvar(ncout, typvar,1,ipk, id_varout )
  ierr = putheadervar(ncout , cfilet, npiglo, npjglo, npk,pdep=gdepw)
  ierr = putvar1d(ncout,time_tag,1,'T')
  pv(:,:) = 0.
  ierr = putvar(ncout,id_varout(1), pv,1,npiglo,npjglo)

  ! initialize first level
  ztemp(:,:) =   getvar(cfilet,'votemper',1,npiglo,npjglo)
  zsal(:,:)  =   getvar(cfilet,'vosaline',1,npiglo,npjglo)
  un  (:,:)  =   getvar(cfileu,'vozocrtx',1,npiglo,npjglo)
  vn  (:,:)  =   getvar(cfilev,'vomecrty',1,npiglo,npjglo)

  ! compute the mask
  DO jj = 1, npjglo - 1
     DO ji = 1, npiglo - 1
        fmask(ji,jj)=0.
        fmask(ji,jj)= un(ji,jj)*un(ji,jj+1) * vn(ji,jj)*vn(ji+1,jj)
        IF (fmask(ji,jj) /= 0.) fmask(ji,jj)=1.
     ENDDO
  ENDDO

  zmask = 1.0
  WHERE(zsal == 0 ) zmask = 0.0
  sigma(:,:,iup) =  sigma0 ( ztemp,zsal,npiglo,npjglo )* zmask(:,:)
  rotn(:,:,iup) = 0.
  DO jj = 1, npjglo -1
     DO ji = 1, npiglo -1   ! vector opt.
        rotn(ji,jj,iup) = (  e2v(ji+1,jj  ) * vn(ji+1,jj  ) - e2v(ji,jj) * vn(ji,jj)    &
             &              - e1u(ji  ,jj+1) * un(ji  ,jj+1) + e1u(ji,jj) * un(ji,jj)  ) &
             &           * fmask(ji,jj) / ( e1f(ji,jj) * e2f(ji,jj) )
     END DO
  END DO

  ! Main vertical loop
  DO jk=2,npk
     PRINT *, 'Level ',jk
     ztemp(:,:) =   getvar(cfilet,'votemper',jk,npiglo,npjglo)
     zsal(:,:)  =   getvar(cfilet,'vosaline',jk,npiglo,npjglo)
     un  (:,:)  =   getvar(cfileu,'vozocrtx',jk,npiglo,npjglo)
     vn  (:,:)  =   getvar(cfilev,'vomecrty',jk,npiglo,npjglo)
     e3w (:,:)  =   getvar(coordzgr,'e3w_ps', jk, npiglo,npjglo)
     WHERE (e3w == 0 ) e3w = 1.

     ! compute the mask at level jk
     DO jj = 1, npjglo - 1
        DO ji = 1, npiglo - 1
           fmask(ji,jj)=0.
           fmask(ji,jj)= un(ji,jj)*un(ji,jj+1) * vn(ji,jj)*vn(ji+1,jj)
           IF (fmask(ji,jj) /= 0.) fmask(ji,jj)=1.
        ENDDO
     ENDDO
     zmask=1.0
     WHERE(zsal == 0 ) zmask = 0.0
     sigma(:,:,idown) =  sigma0 ( ztemp,zsal,npiglo,npjglo )* zmask(:,:)

     !  d(sigma0)/dz at W point ( masked if down level is masked )
     dsig(:,:)=(sigma(:,:,idown) - sigma(:,:,iup)) /e3w *zmask

     rotn(:,:,idown) = 0.
     DO jj = 1, npjglo -1
        DO ji = 1, npiglo -1   ! vector opt.
           rotn(ji,jj,idown) = (  e2v(ji+1,jj  ) * vn(ji+1,jj  ) - e2v(ji,jj) * vn(ji,jj)    &
                &              - e1u(ji  ,jj+1) * un(ji  ,jj+1) + e1u(ji,jj) * un(ji,jj)  ) &
                &           * fmask(ji,jj) / ( e1f(ji,jj) * e2f(ji,jj) )
        END DO
     END DO

     ! curl at f point, w level
     rot(:,:)= 0.5*( rotn(:,:,idown) + rotn(:,:,iup) )

     ! Full pv:
     DO ji=2,npiglo 
        DO jj = 2, npjglo
           zrot=0.25*( rot(ji,jj) + rot(ji-1,jj) + rot(ji,jj-1) + rot(ji-1,jj-1) )
!          pv(ji,jj) = 1/rho0*(fcorio(ji,jj)+zrot)*dsig(ji,jj)*1.e11
           pv(ji,jj) = (fcorio(ji,jj)+zrot)*dsig(ji,jj)*1.e7
           !        pv(ji,jj) = dsig(ji,jj)*1000.
        END DO
     END DO
     ierr = putvar(ncout,id_varout(1), pv,jk,npiglo,npjglo)
     
     ! swap index up and down
     itmp=iup
     iup=idown
     idown=itmp
  END DO

  ierr = closeout(ncout)
  PRINT *,'cdfpv completed successfully'
END  PROGRAM cdfpv
