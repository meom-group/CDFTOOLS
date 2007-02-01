PROGRAM cdfvhst
  !!-------------------------------------------------------------------
  !!               ***  PROGRAM cdfvhst  ***
  !!
  !!  **  Purpose  :  Compute Verticaly integrated  Heat Salt Transport. 
  !!                  PARTIAL STEPS
  !!  
  !!  **  Method   :  Compute the 2D fields somevt, somevs and sozout, sozous
  !!                  as the integral on the vertical of ut, vt, us, vs
  !!                  Save on the nc file
  !!
  !!
  !! history ;
  !!  Original :  J.M. Molines (jan. 2005) (known then as cdfheattrp-save.f90 )
  !!              J.M. Molines : use module
  !!-------------------------------------------------------------------
  !! * Modules used
  USE cdfio

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: jk                            !: dummy loop index
  INTEGER   :: ierr                                !: working integer
  INTEGER   :: narg, iargc                         !: command line 
  INTEGER   :: npiglo,npjglo, npk                  !: size of the domain
  INTEGER   ::  ncout
  INTEGER, DIMENSION(4) ::  ipk, id_varout         !

  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::  zmask, e1v, e3v ,gphiv, zvt, zvs !: mask, metrics
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::         e2u, e3u ,gphiu, zut, zus !: mask, metrics
  REAL(KIND=4) ,DIMENSION(1)                    ::  tim
  REAL(KIND=4) ,DIMENSION(:) , ALLOCATABLE ::  gphimean_glo, gphimean_atl, gphimean_pac, &
       &                                       gphimean_ind, gphimean_aus, gphimean_med

  REAL(KIND=8),   DIMENSION (:,:), ALLOCATABLE :: zwk , zwks, zwkut, zwkus
  REAL(KIND=8),   DIMENSION (:,:), ALLOCATABLE :: ztrp, ztrps,ztrput, ztrpus
  REAL(KIND=8) ,DIMENSION(:) , ALLOCATABLE ::  zonal_heat_glo, zonal_heat_atl, zonal_heat_pac, &
       &                                       zonal_heat_ind, zonal_heat_aus, zonal_heat_med
  REAL(KIND=8) ,DIMENSION(:) , ALLOCATABLE ::  zonal_salt_glo, zonal_salt_atl, zonal_salt_pac, &
       &                                       zonal_salt_ind, zonal_salt_aus, zonal_salt_med

  CHARACTER(LEN=80) :: cfilet , cfileoutnc='trp.nc'
  CHARACTER(LEN=80) :: coordhgr='mesh_hgr.nc',  coordzgr='mesh_zgr.nc'
  TYPE (variable), DIMENSION(4)    :: typvar      !: structure for attribute


  INTEGER    :: istatus

  ! constants
  REAL(KIND=4), PARAMETER   ::  rau0=1000., rcp=4000.

  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfvhst  VTfile '
     PRINT *,' Computes the vertically integrated transports at each grid cell'
     PRINT *,' PARTIAL CELLS VERSION'
     PRINT *,' Files mesh_hgr.nc, mesh_zgr.nc ,mask.nc must be in te current directory'
     PRINT *,' Output on trp.nc, variables somevt somevs sozout sozous '
     STOP
  ENDIF

  CALL getarg (1, cfilet)
  npiglo= getdim (cfilet,'x')
  npjglo= getdim (cfilet,'y')
  npk   = getdim (cfilet,'depth')

 ! define new variables for output 
  typvar(1)%name= 'somevt'
  typvar(2)%name= 'somevs'
  typvar(3)%name= 'sozout'
  typvar(4)%name= 'sozous'

  typvar(1)%units='W'
  typvar(2)%units='kg.s-1'
  typvar(3)%units='W'
  typvar(4)%units='kg.s-1'

  typvar%missing_value=0.
  typvar%valid_min= -100.
  typvar%valid_max= 100.

  typvar(1)%long_name='Meridional_heat_transport'
  typvar(2)%long_name='Meridional_salt_transport'
  typvar(3)%long_name='Zonal_heat_transport'
  typvar(4)%long_name='Zonal_salt_transport'

  typvar(1)%short_name='somevt'
  typvar(2)%short_name='somevs'
  typvar(3)%short_name='sozout'
  typvar(4)%short_name='sozous'
  typvar%online_operation='N/A'
  typvar%axis='TYX'

  ipk(1) = 1  !  2D
  ipk(2) = 1  !  2D
  ipk(3) = 1  !  2D
  ipk(4) = 1  !  2D

  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk

  ! Allocate arrays
  ALLOCATE ( zmask(npiglo,npjglo) )
  ALLOCATE ( zwk(npiglo,npjglo)  ,zvt(npiglo,npjglo) )
  ALLOCATE ( zwks(npiglo,npjglo) ,zvs(npiglo,npjglo) )
  ALLOCATE ( zwkut(npiglo,npjglo)  ,zut(npiglo,npjglo) )
  ALLOCATE ( zwkus(npiglo,npjglo) ,zus(npiglo,npjglo) )
  ALLOCATE ( e1v(npiglo,npjglo),e3v(npiglo,npjglo), gphiv(npiglo,npjglo))
  ALLOCATE ( e2u(npiglo,npjglo),e3u(npiglo,npjglo), gphiu(npiglo,npjglo))
  ALLOCATE ( ztrp(npiglo,npjglo))
  ALLOCATE ( ztrps(npiglo,npjglo))
  ALLOCATE ( ztrput(npiglo,npjglo))
  ALLOCATE ( ztrpus(npiglo,npjglo))
  ALLOCATE ( zonal_heat_glo(npjglo), zonal_heat_atl(npjglo), zonal_heat_pac(npjglo))
  ALLOCATE ( zonal_heat_ind(npjglo), zonal_heat_aus(npjglo), zonal_heat_med(npjglo) )
  ALLOCATE ( zonal_salt_glo(npjglo), zonal_salt_atl(npjglo), zonal_salt_pac(npjglo))
  ALLOCATE ( zonal_salt_ind(npjglo), zonal_salt_aus(npjglo), zonal_salt_med(npjglo) )
  ALLOCATE ( gphimean_glo(npjglo) , gphimean_atl(npjglo), gphimean_pac(npjglo))
  ALLOCATE ( gphimean_ind(npjglo),gphimean_aus(npjglo),gphimean_med(npjglo))

  ! create output fileset
  ncout =create(cfileoutnc, cfilet, npiglo,npjglo,npk)
  ierr= createvar(ncout ,typvar,4, ipk,id_varout )
  ierr= putheadervar(ncout, cfilet,npiglo, npjglo,npk)
  tim=getvar1d(cfilet,'time_counter',1)
  ierr=putvar1d(ncout,tim,1,'T')


  e1v(:,:) = getvar(coordhgr, 'e1v', 1,npiglo,npjglo)
  e2u(:,:) = getvar(coordhgr, 'e2u', 1,npiglo,npjglo)
  gphiv(:,:) = getvar(coordhgr, 'gphiv', 1,npiglo,npjglo)
  gphiu(:,:) = getvar(coordhgr, 'gphiu', 1,npiglo,npjglo)

  ztrp(:,:)= 0
  ztrps(:,:)= 0
  ztrput(:,:)= 0
  ztrpus(:,:)= 0
  DO jk = 1,npk
     PRINT *,'level ',jk
     ! Get temperature and salinity at jk
     zvt(:,:)= getvar(cfilet, 'vomevt',  jk ,npiglo,npjglo)
     zvs(:,:)= getvar(cfilet, 'vomevs',  jk ,npiglo,npjglo)
     zut(:,:)= getvar(cfilet, 'vozout',  jk ,npiglo,npjglo)
     zus(:,:)= getvar(cfilet, 'vozous',  jk ,npiglo,npjglo)

     ! get e3v at level jk
     e3v(:,:) = getvar(coordzgr, 'e3v_ps', jk,npiglo,npjglo)
     e3u(:,:) = getvar(coordzgr, 'e3u_ps', jk,npiglo,npjglo)
     zwk(:,:) = zvt(:,:)*e1v(:,:)*e3v(:,:)
     zwks(:,:) = zvs(:,:)*e1v(:,:)*e3v(:,:)
     zwkut(:,:) = zut(:,:)*e2u(:,:)*e3u(:,:)
     zwkus(:,:) = zus(:,:)*e2u(:,:)*e3u(:,:)

     ! integrates vertically 
     ztrp(:,:) = ztrp(:,:) + zwk(:,:) * rau0*rcp
     ztrps(:,:) = ztrps(:,:) + zwks(:,:)  
     ztrput(:,:) = ztrput(:,:) + zwkut(:,:) * rau0*rcp
     ztrpus(:,:) = ztrpus(:,:) + zwkus(:,:)  

  END DO  ! loop to next level

     ierr = putvar(ncout, id_varout(1) ,SNGL(ztrp),   1, npiglo, npjglo)
     ierr = putvar(ncout, id_varout(2) ,SNGL(ztrps),  1, npiglo, npjglo)
     ierr = putvar(ncout, id_varout(3) ,SNGL(ztrput), 1, npiglo, npjglo)
     ierr = putvar(ncout, id_varout(4) ,SNGL(ztrpus), 1, npiglo, npjglo)

     istatus = closeout (ncout)

   END PROGRAM cdfvhst
