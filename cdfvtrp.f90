PROGRAM cdfvtrp
  !!-------------------------------------------------------------------
  !!               ***  PROGRAM cdfvtrp  ***
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
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
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

  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::  zmask, e1v, e3v !: mask, metrics
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::         e2u, e3u !: mask, metrics
  REAL(KIND=4) ,DIMENSION(1)                    ::  tim

  REAL(KIND=8),   DIMENSION (:,:), ALLOCATABLE :: zwku , zwkv, zu, zv
  REAL(KIND=8),   DIMENSION (:,:), ALLOCATABLE :: ztrpu, ztrpv

  CHARACTER(LEN=256) :: cfileu, cfilev , cfileoutnc='trp.nc'
  CHARACTER(LEN=256) :: coordhgr='mesh_hgr.nc',  coordzgr='mesh_zgr.nc'
  TYPE (variable), DIMENSION(4)    :: typvar      !: structure for attribute


  INTEGER    :: istatus

  ! constants

  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfvtrp  Ufile Vfile '
     PRINT *,' Computes the vertically integrated transports at each grid cell'
     PRINT *,' PARTIAL CELLS VERSION'
     PRINT *,' Files mesh_hgr.nc, mesh_zgr.nc ,mask.nc must be in te current directory'
     PRINT *,' Output on trp.nc, variables somevtrp sozoutrp '
     STOP
  ENDIF

  CALL getarg (1, cfileu)
  CALL getarg (2, cfilev)
  npiglo= getdim (cfileu,'x')
  npjglo= getdim (cfileu,'y')
  npk   = getdim (cfileu,'depth')

 ! define new variables for output 
  typvar(2)%name= 'somevtrp'
  typvar(1)%name= 'sozoutrp'

  typvar(1)%units='m3/s'
  typvar(2)%units='m3/s'

  typvar%missing_value=0.
  typvar%valid_min= -100.
  typvar%valid_max= 100.

  typvar(2)%long_name='Z_Integrated_Meridional_mass_transport'
  typvar(1)%long_name='Z_Integrated_Zonal_mass_transport'

  typvar(2)%short_name='somevtrp'
  typvar(1)%short_name='sozoutrp'

  typvar%online_operation='N/A'
  typvar%axis='TYX'

  ipk(1) = 1  !  2D
  ipk(2) = 1  !  2D
 
  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk

  ! Allocate arrays
  ALLOCATE ( zmask(npiglo,npjglo) )

  ALLOCATE ( e1v(npiglo,npjglo),e3v(npiglo,npjglo))
  ALLOCATE ( e2u(npiglo,npjglo),e3u(npiglo,npjglo))
  ALLOCATE ( zwku(npiglo,npjglo), zwkv(npiglo,npjglo))
  ALLOCATE ( ztrpu(npiglo,npjglo), ztrpv(npiglo,npjglo))
  ALLOCATE ( zu(npiglo,npjglo), zv(npiglo,npjglo))
 

  ! create output fileset
  ncout =create(cfileoutnc, cfileu, npiglo,npjglo,1)
  ierr= createvar(ncout ,typvar,2, ipk,id_varout )
  ierr= putheadervar(ncout, cfileu,npiglo, npjglo,1)
  tim=getvar1d(cfileu,'time_counter',1)
  ierr=putvar1d(ncout,tim,1,'T')


  e1v(:,:) = getvar(coordhgr, 'e1v', 1,npiglo,npjglo)
  e2u(:,:) = getvar(coordhgr, 'e2u', 1,npiglo,npjglo)

  ztrpu(:,:)= 0
  ztrpv(:,:)= 0
  DO jk = 1,npk
     PRINT *,'level ',jk
     ! Get temperature and salinity at jk
     zu(:,:)= getvar(cfileu, 'vozocrtx',  jk ,npiglo,npjglo)
     zv(:,:)= getvar(cfilev, 'vomecrty',  jk ,npiglo,npjglo)

     ! get e3v at level jk
     e3v(:,:) = getvar(coordzgr, 'e3v_ps', jk,npiglo,npjglo, ldiom=.true.)
     e3u(:,:) = getvar(coordzgr, 'e3u_ps', jk,npiglo,npjglo, ldiom=.true.)
     zwku(:,:) = zu(:,:)*e2u(:,:)*e3u(:,:)
     zwkv(:,:) = zv(:,:)*e1v(:,:)*e3v(:,:)
     ! integrates vertically 
     ztrpu(:,:) = ztrpu(:,:) + zwku(:,:)
     ztrpv(:,:) = ztrpv(:,:) + zwkv(:,:)
 
  END DO  ! loop to next level

  ierr = putvar(ncout, id_varout(1) ,REAL(ztrpu(:,:)), 1, npiglo, npjglo)
  ierr = putvar(ncout, id_varout(2) ,REAL(ztrpv(:,:)), 1, npiglo, npjglo)

  istatus = closeout (ncout)

   END PROGRAM cdfvtrp
