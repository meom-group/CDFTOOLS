PROGRAM cdfpsi_level
  !!-------------------------------------------------------------------
  !!               ***  PROGRAM cdfpsi_level  ***
  !!
  !!  **  Purpose  :  Compute Stream Function for each level
  !!                  PARTIAL STEPS
  !!  
  !!  **  Method   :  Compute the 3D fields ztrpu, ztrpv as the vertical
  !!                  integral of u, v on their respective points.
  !!                  Then integrate from south to north : ==> psiu
  !!                  Then integrate from West to East   : ==> psiv
  !!                  (should be almost the same (if no error ))
  !!   Default (appropriate for global model): output psiu;
  !!                    normalizes the values setting psi (jpi,jpj) = 0
  !!   If option "V" is given as last argument, output psiv, 
  !!                    normalizes values setting psi(jpi,1) = 0.
  !!                    This is appropriate for North Atlantic 
  !!
  !! history ;
  !!  Original :  J.M. Molines (May 2005 )
  !!-------------------------------------------------------------------
  !!  $Rev: 256 $
  !!  $Date: 2009-07-21 17:49:27 +0200 (mar. 21 juil. 2009) $
  !!  $Id$
  !!--------------------------------------------------------------
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

  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::  zmask, e1v, e3v , zv !: mask, metrics
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::         e2u, e3u , zu !: mask, metrics
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::         glamf, gphif
  REAL(KIND=4) ,DIMENSION(1)                    ::  tim

  REAL(KIND=8),   DIMENSION (:,:), ALLOCATABLE :: ztrpu, ztrpv, psiu, psiv

  CHARACTER(LEN=256) :: cfileu ,cfilev, cfileoutnc='psi_level.nc'
  CHARACTER(LEN=256) :: coordhgr='mesh_hgr.nc',  coordzgr='mesh_zgr.nc', cmask='mask.nc'
  CHARACTER(LEN=1)  :: coption
  CHARACTER(LEN=256) :: cdep

  TYPE(variable), DIMENSION(1) :: typvar         !: structure for attributes

  INTEGER    :: istatus

  ! constants

  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfpsi_level  Ufile Vfile <V> (optional argument)'
     PRINT *,' Compute the barotropic stream function as the integral of the transport'
     PRINT *,' PARTIAL CELLS VERSION'
     PRINT *,' Files mesh_hgr.nc, mesh_zgr.nc ,mask.nc must be in te current directory'
     PRINT *,' Output on psi_level.nc, variables sobarstf on f-points'
     PRINT *,' Default works well for a global ORCA grid. use V 3rdargument for North Atlantic'
     STOP
  ENDIF

  CALL getarg (1, cfileu  )
  CALL getarg (2, cfilev  )
  CALL getarg (3, coption )

  npiglo= getdim (cfileu,'x')
  npjglo= getdim (cfileu,'y')
  npk   = getdim (cfileu,'depth')

 ! define new variables for output ( must update att.txt)
  typvar(1)%cname= 'sobarstf'
  typvar(1)%cunits='m3/s'
  typvar(1)%rmissing_value=0.
  typvar(1)%valid_min= -300.e6
  typvar(1)%valid_max= 300.e6
  typvar(1)%clong_name='Barotropic_Stream_Function'
  typvar(1)%cshort_name='sobarstf'
  typvar(1)%conline_operation='N/A'
  typvar(1)%caxis='TZYX'
  ipk(1) = npk  !  3D ( X, Y , Z, T )

  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk
  IF ( coption == 'V') PRINT *, ' Use psiv (ex. North Atlantic case)'

  ! Allocate arrays
  ALLOCATE ( zmask(npiglo,npjglo) )
  ALLOCATE ( e1v(npiglo,npjglo),e3v(npiglo,npjglo))
  ALLOCATE ( e2u(npiglo,npjglo),e3u(npiglo,npjglo))
  ALLOCATE ( zu(npiglo,npjglo),ztrpu(npiglo,npjglo), psiu(npiglo,npjglo) )
  ALLOCATE ( zv(npiglo,npjglo),ztrpv(npiglo,npjglo), psiv(npiglo,npjglo))
  ALLOCATE ( glamf(npiglo,npjglo), gphif(npiglo,npjglo))

  glamf(:,:) = getvar(coordhgr, 'glamf',1,npiglo,npjglo)
  gphif(:,:) = getvar(coordhgr, 'gphif',1,npiglo,npjglo)

  ! create output fileset
  ncout =create(cfileoutnc, cfileu, npiglo,npjglo,npk)
  ierr= createvar(ncout ,typvar,1, ipk,id_varout )
  ierr= putheadervar(ncout , cfileu, npiglo, npjglo, npk)
!  ierr= putheadervar(ncout , cfileu, npiglo, npjglo, npk,cdep=cdep)
!  ierr= putheadervar(ncout, cfileu,npiglo, npjglo,1,glamf, gphif)
  tim=getvar1d(cfileu,'time_counter',1)
  ierr=putvar1d(ncout,tim,1,'T')

  e1v(:,:) = getvar(coordhgr, 'e1v', 1,npiglo,npjglo)
  e2u(:,:) = getvar(coordhgr, 'e2u', 1,npiglo,npjglo)

     ztrpu(:,:)= 0.d0
     ztrpv(:,:)= 0.d0

  DO jk = 1,npk
  
     zmask(:,:) = getvar(cmask, 'fmask', jk,npiglo,npjglo)
     ! get rid of the free-slip/no-slip condition
     WHERE ( zmask >= 2 ) zmask = 1

     PRINT *,'level ',jk
     IF ( coption == 'V' ) THEN
        zv(:,:)= getvar(cfilev, 'vomecrty',  jk ,npiglo,npjglo)
        e3v(:,:) = getvar(coordzgr, 'e3v_ps', jk,npiglo,npjglo, ldiom=.true.)
        ztrpv(:,:) = zv(:,:)*e1v(:,:)*e3v(:,:)*1.d0  ! meridional transport of each grid cell
     ELSE
     ! Get zonal velocity  at jk
        zu(:,:)= getvar(cfileu, 'vozocrtx',  jk ,npiglo,npjglo)
     ! get e3v at level jk
        e3u(:,:) = getvar(coordzgr, 'e3u_ps', jk,npiglo,npjglo, ldiom=.true.)
     ! integrates vertically 
        ztrpu(:,:) = zu(:,:)*e2u(:,:)*e3u(:,:)*1.d0  ! zonal transport of each grid cell
     ENDIF

  IF (coption == 'V' ) THEN
  ! integrate zonally from east to west 
     psiv(npiglo,:)= 0.0
     DO ji=npiglo-1,1,-1
        psiv(ji,:) = psiv(ji+1,:) - ztrpv(ji,:)  ! psi at f point
     END DO
     psiv(:,:) = psiv(:,:) *zmask(:,:)
     ierr = putvar(ncout, id_varout(1) ,REAL(psiv),   jk, npiglo, npjglo)
     !ierr = putvar(ncout, id_varout(1) ,REAL(ztrpv),   jk, npiglo, npjglo)

  ELSE
  ! integrate from the south to the north with zonal transport
     psiu(:,:) = 0.d0

     DO jj = 2, npjglo
       psiu(:,jj) = psiu(:,jj-1) - ztrpu(:,jj)   ! psi at f point
     END DO
     psiu(:,:) = (psiu(:,:) -psiu(npiglo,npjglo) ) * zmask(:,:)
     ierr = putvar(ncout, id_varout(1) ,REAL(psiu),   jk, npiglo, npjglo)
     !ierr = putvar(ncout, id_varout(1) ,REAL(ztrpu),   jk, npiglo, npjglo)
   ENDIF

  END DO  ! loop to next level

   istatus = closeout (ncout)

   END PROGRAM cdfpsi_level
