PROGRAM cdfpsi_full
  !!-------------------------------------------------------------------
  !!               ***  PROGRAM cdfpsi_full  ***
  !!
  !!  **  Purpose  :  Compute Barotropic Stream Function
  !!                  FULL STEPS
  !!  
  !!  **  Method   :  Compute the 2D fields ztrpu, ztrpv 
  !!                  as the integral on the vertical of u, v on their
  !!                  respective points. 
  !!                  Then integrate from south to north : ==> psiu
  !!                  Then integrate from West to East   : ==> psiv
  !!                  (should be almost the same (if no error )
  !!                  Then normalize the values setting psi (jpi,jpj) = 0
  !!                  Following Anne-Marie matlab program, we only take psiu.
  !!                  An alternative could be to average psiu and psiv ...
  !!
  !! history ;
  !!  Original :  J.M. Molines (May 2005 )
  !!-------------------------------------------------------------------
  !! * Modules used
  USE cdfio

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: ji,jj,jk                            !: dummy loop index
  INTEGER   :: ierr                                !: working integer
  INTEGER   :: narg, iargc                         !: command line 
  INTEGER   :: npiglo,npjglo, npk                  !: size of the domain
  INTEGER   :: ncout
  INTEGER, DIMENSION(1) ::  ipk, id_varout         !

  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::  zmask, e1v, zv !: mask, metrics
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::         e2u, zu !: mask, metrics
  REAL(KIND=4), DIMENSION (:),     ALLOCATABLE ::         e3t     !:  full step vertical metrics
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::         glamf, gphif
  REAL(KIND=4) ,DIMENSION(1)                    ::  tim

  REAL(KIND=8),   DIMENSION (:,:), ALLOCATABLE :: ztrpu, ztrpv, psiu, psiv

  CHARACTER(LEN=80) :: cfileu ,cfilev, cfileoutnc='psi.nc'
  CHARACTER(LEN=80) :: coordhgr='mesh_hgr.nc',  coordzgr='mesh_zgr.nc', cmask='mask.nc'

  TYPE (variable), DIMENSION(1)  :: typvar        !: structure for attributes

  INTEGER    :: istatus

  ! constants
  REAL(KIND=4), PARAMETER   ::  rau0=1000., rcp=4000.

  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfpsi  Ufile Vfile '
     PRINT *,' Computes the barotropic stream function as the integral of the transport'
     PRINT *,' FULL STEPS  VERSION'
     PRINT *,' Files mesh_hgr.nc, mesh_zgr.nc ,mask.nc must be in te current directory'
     PRINT *,' Output on psi.nc, variables sobarstf'
     STOP
  ENDIF

  CALL getarg (1, cfileu)
  CALL getarg (2, cfilev)

  npiglo= getdim (cfileu,'x')
  npjglo= getdim (cfileu,'y')
  npk   = getdim (cfileu,'depth')

 ! define new variables for output ( must update att.txt)
  typvar(1)%name= 'sobarstf'
  typvar(1)%units='m3/s'
  typvar(1)%missing_value=0.
  typvar(1)%valid_min= -300.e6
  typvar(1)%valid_max= 300.e6
  typvar(1)%long_name='Barotropic_Stream_Function'
  typvar(1)%short_name='sobarstf'
  typvar(1)%online_operation='N/A'
  typvar(1)%axis='TYX'

  ipk(1) = 1  !  2D ( X, Y , T )

  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk

  ! Allocate arrays
  ALLOCATE ( zmask(npiglo,npjglo) )
  ALLOCATE ( e1v(npiglo,npjglo) )
  ALLOCATE ( e2u(npiglo,npjglo))
  ALLOCATE ( e3t(npk) )
  ALLOCATE ( zu(npiglo,npjglo),ztrpu(npiglo,npjglo), psiu(npiglo,npjglo) )
  ALLOCATE ( zv(npiglo,npjglo),ztrpv(npiglo,npjglo), psiv(npiglo,npjglo))
  ALLOCATE ( glamf(npiglo,npjglo), gphif(npiglo,npjglo))

  glamf(:,:) = getvar(coordhgr, 'glamf',1,npiglo,npjglo)
  gphif(:,:) = getvar(coordhgr, 'gphif',1,npiglo,npjglo)

  ! create output fileset
  ncout =create(cfileoutnc, cfileu, npiglo,npjglo,1)
  ierr= createvar(ncout ,typvar,1, ipk,id_varout )
  ierr= putheadervar(ncout, cfileu,npiglo, npjglo,1,glamf, gphif)
  tim=getvar1d(cfileu,'time_counter',1)
  ierr=putvar1d(ncout,tim,1,'T')


  e1v(:,:) = getvar(coordhgr, 'e1v', 1,npiglo,npjglo)
  e2u(:,:) = getvar(coordhgr, 'e2u', 1,npiglo,npjglo)
  zmask(:,:) = getvar(cmask, 'fmask', 1,npiglo,npjglo)
  ! get rid of the free-slip/no-slip condition
  WHERE ( zmask == 2 ) zmask = 1

  ztrpu(:,:)= 0.d0
  ztrpv(:,:)= 0.d0
  e3t(:)   = getvare3(coordzgr,'e3t',npk)

  DO jk = 1,npk
     PRINT *,'level ',jk
     ! Get temperature and salinity at jk
     zu(:,:)= getvar(cfileu, 'vozocrtx',  jk ,npiglo,npjglo)
     zv(:,:)= getvar(cfilev, 'vomecrty',  jk ,npiglo,npjglo)

     ! integrates vertically 
     ztrpu(:,:) = ztrpu(:,:) + zu(:,:)*e2u(:,:)*e3t(jk)  ! zonal transport of each grid cell
     ztrpv(:,:) = ztrpv(:,:) + zv(:,:)*e1v(:,:)*e3t(jk)  ! meridional transport of each grid cell

  END DO  ! loop to next level

  ! integrate from the south to the north with zonal transport
  psiu(:,:) = 0.d0

  DO jj = 2, npjglo
    psiu(:,jj) = psiu(:,jj-1) - ztrpu(:,jj)   ! psi at f point
  END DO


  ! integrate zonally form west to east
  psiv(1,:)=psiu(1,:)

  DO ji=2, npiglo
     psiv(ji,:) = psiv(ji-1,:) + ztrpv(ji,:)  ! psi at f point
  END DO
  psiu(:,:) = (psiu(:,:) -psiu(npiglo,npjglo) ) * zmask(:,:)
  psiv(:,:) = (psiv(:,:) -psiv(npiglo,npjglo) ) * zmask(:,:)
  psiv=0.5 * (psiu+psiv)

  ierr = putvar(ncout, id_varout(1) ,SNGL(psiu),   1, npiglo, npjglo)
!  ierr = putvar(ncout, id_varout(2) ,SNGL(psiv),   1, npiglo, npjglo)

     istatus = closeout (ncout)

   END PROGRAM cdfpsi_full
