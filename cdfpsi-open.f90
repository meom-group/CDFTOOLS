PROGRAM cdfpsi_open
  !!-------------------------------------------------------------------
  !!               ***  PROGRAM cdfpsi_open  ***
  !!
  !!  **  Purpose  :  Compute Barotropic Stream Function
  !!                  PARTIAL STEPS
  !!  
  !!  **  Method   :  Compute the 2D fields ztrpu, ztrpv 
  !!                  as the integral on the vertical of u, v on their
  !!                  respective points. 
  !!                  Then, starting from the upper left  point, 
  !!                  initialize psi from W to E, on the northern line, using ztrpv.
  !!                  Then, from this first line (N to S) cumulate transport
  !!                  using ztrpu.
  !!                 This works in any case. (except perharps when northern line is a 
  !!                 folding line ( orca config)
  !!              REM: in this version PSI is not masked, which allows the exact calculation
  !!                 of barotropic transport just making a difference.
  !!
  !! history ;
  !!  Original :  J.M. Molines (May 2005 )
  !!             open version   J.M. Molines (March 2007 )
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

  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::  zmask, e1v, e3v , zv !: mask, metrics
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::         e2u, e3u , zu !: mask, metrics
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::         glamf, gphif
  REAL(KIND=4) ,DIMENSION(1)                    ::  tim

  REAL(KIND=8),   DIMENSION (:,:), ALLOCATABLE :: ztrpu, ztrpv, psi

  CHARACTER(LEN=80) :: cfileu ,cfilev, cfileoutnc='psi.nc'
  CHARACTER(LEN=80) :: coordhgr='mesh_hgr.nc',  coordzgr='mesh_zgr.nc', cmask='mask.nc'
  CHARACTER(LEN=1)  :: coption

  TYPE(variable), DIMENSION(1)  :: typvar         !: structure for attributes

  INTEGER    :: istatus

  ! constants

  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfpsi_open  Ufile Vfile '
     PRINT *,' Computes the barotropic stream function as the integral of the transport'
     PRINT *,' PARTIAL CELLS VERSION'
     PRINT *,' Files mesh_hgr.nc, mesh_zgr.nc ,mask.nc must be in te current directory'
     PRINT *,' Output on psi.nc, variables sobarstf on f-points'
     STOP
  ENDIF

  CALL getarg (1, cfileu  )
  CALL getarg (2, cfilev  )

  npiglo= getdim (cfileu,'x')
  npjglo= getdim (cfileu,'y')
  npk   = getdim (cfileu,'depth')

 ! define new variables for output
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
  ALLOCATE ( e1v(npiglo,npjglo),e3v(npiglo,npjglo))
  ALLOCATE ( e2u(npiglo,npjglo),e3u(npiglo,npjglo))
  ALLOCATE ( zu(npiglo,npjglo),ztrpu(npiglo,npjglo), psi(npiglo,npjglo) )
  ALLOCATE ( zv(npiglo,npjglo),ztrpv(npiglo,npjglo))
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
  DO jk = 1,npk
     ! Get  velocity  at jk
        zu(:,:)= getvar(cfileu, 'vozocrtx',  jk ,npiglo,npjglo)
        zv(:,:)= getvar(cfilev, 'vomecrty',  jk ,npiglo,npjglo)
     ! get e3 at level jk
        e3u(:,:) = getvar(coordzgr, 'e3u_ps', jk,npiglo,npjglo)
        e3v(:,:) = getvar(coordzgr, 'e3v_ps', jk,npiglo,npjglo)
     ! integrates vertically 
        ztrpu(:,:) = ztrpu(:,:) + zu(:,:)*e2u(:,:)*e3u(:,:)  ! zonal transport of each grid cell
        ztrpv(:,:) = ztrpv(:,:) + zv(:,:)*e1v(:,:)*e3v(:,:)  ! meridional transport of each grid cell
  END DO  ! loop to next level

  ! compute transport along line jj=jpj
  psi(1,npjglo) = ztrpv(1,npjglo)
    DO ji=2,npiglo 
      psi(ji,npjglo) = psi(ji-1,npjglo) + ztrpv(ji,npjglo)
    END DO
  ! Then compute from N to S the transport using zonal contribution
    DO jj=npjglo-1,1,-1
      DO ji=1,npiglo
        psi(ji,jj)=psi(ji,jj+1)  + ztrpu(ji,jj+1)
      END DO
    END DO
  
   ierr = putvar(ncout, id_varout(1) ,SNGL(psi),   1, npiglo, npjglo)
   istatus = closeout (ncout)

   END PROGRAM cdfpsi_open
