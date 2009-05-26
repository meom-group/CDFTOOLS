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
  !!  **  Usage : cdfpsi-open fileU fileV [-mask] [-moy]
  !!                when option -mask is used, output is masked (by fmask(k=1) )
  !!                when option -moy is used, the result is the mean value of the meridional and
  !!                     and zonal based calculation.
  !!                
  !!
  !! history ;
  !!  Original :  J.M. Molines (May 2005 )
  !!             open version   J.M. Molines (March 2007 )
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!-------------------------------------------------------------------
  !! * Modules used
  USE cdfio

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: ji,jj,jk, jarg                      !: dummy loop index
  INTEGER   :: ierr                                !: working integer
  INTEGER   :: narg, iargc                         !: command line 
  INTEGER   :: npiglo,npjglo, npk                  !: size of the domain
  INTEGER   :: ncout
  INTEGER, DIMENSION(3) ::  ipk, id_varout         !

  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::  zmask, e1v, e3v , zv !: mask, metrics
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::         e2u, e3u , zu !: mask, metrics
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::         glamf, gphif
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::         sshu, sshv, ztmp
  REAL(KIND=4) ,DIMENSION(1)                    ::  tim

  REAL(KIND=8),   DIMENSION (:,:), ALLOCATABLE :: ztrpu, ztrpv, psi1, psi2, psissh
  REAL(KIND=8),   DIMENSION (:,:), ALLOCATABLE :: ztrpussh, ztrpvssh

  CHARACTER(LEN=80) :: cfileu ,cfilev, cfileoutnc='psi.nc', cfilet
  CHARACTER(LEN=80) :: coordhgr='mesh_hgr.nc',  coordzgr='mesh_zgr.nc', cmask='mask.nc'
  CHARACTER(LEN=10)  :: coption

  TYPE(variable), DIMENSION(3)  :: typvar         !: structure for attributes

  INTEGER    :: istatus
  LOGICAL    :: lmask=.FALSE., lmoy=.FALSE.

  ! constants

  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfpsi_open  Ufile Vfile Tfile -mask -moy'
     PRINT *,' Computes the barotropic stream function as the integral of the transport'
     PRINT *,'    Option -mask : result are masked (default : no masked) '
     PRINT *,'    Option -moy : results is the mean between meridional and zonal calculation'
     PRINT *,'       that should be the same if all was perfect ...'
     PRINT *,' PARTIAL CELLS VERSION'
     PRINT *,' Files mesh_hgr.nc, mesh_zgr.nc ,mask.nc must be in te current directory'
     PRINT *,' Output on psi.nc, variables sobarstf on f-points'
     STOP
  ENDIF

  CALL getarg (1, cfileu  )
  CALL getarg (2, cfilev  )
  CALL getarg (3, cfilet  )
  IF (narg > 3 ) THEN
     DO  jarg = 4, narg 
        CALL getarg(jarg, coption)
        SELECT CASE ( coption )
        CASE ( '-mask' )
           lmask=.TRUE.
        CASE ( '-moy' )
           lmoy=.TRUE.
        CASE DEFAULT
           PRINT *,' Unknown option : ', TRIM(coption)
           STOP
        END SELECT
     END DO
  ENDIF

  npiglo= getdim (cfileu,'x')
  npjglo= getdim (cfileu,'y')
  npk   = getdim (cfileu,'depth')

  ! define new variables for output
  typvar(1)%name= 'sobarstf'
  typvar(2)%name= 'sobarstfssh'
  typvar(3)%name= 'sobarstftotal'
  typvar(:)%units='m3/s'
  typvar(:)%missing_value=0.
  typvar(:)%valid_min= -300.e6
  typvar(:)%valid_max= 300.e6
  typvar(1)%long_name='Barotropic_Stream_Function '
  typvar(2)%long_name='Barotropic_Stream_Function SSH contribution'
  typvar(3)%long_name='Barotropic_Stream_Function total'
  typvar(1)%short_name='sobarstf '
  typvar(2)%short_name='sobarstfssh '
  typvar(3)%short_name='sobarstftotal '
  typvar(:)%online_operation='N/A'
  typvar(:)%axis='TYX'
  ipk(:) = 1  !  2D ( X, Y , T )

  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk

  ! Allocate arrays
  ALLOCATE ( zmask(npiglo,npjglo) )
  ALLOCATE ( e1v(npiglo,npjglo),e3v(npiglo,npjglo))
  ALLOCATE ( e2u(npiglo,npjglo),e3u(npiglo,npjglo), sshu(npiglo,npjglo), sshv(npiglo,npjglo), ztmp(npiglo,npjglo))
  ALLOCATE ( zu(npiglo,npjglo),ztrpu(npiglo,npjglo), psi1(npiglo,npjglo), psi2(npiglo,npjglo) ,psissh(npiglo,npjglo))
  ALLOCATE ( zv(npiglo,npjglo),ztrpv(npiglo,npjglo),ztrpussh(npiglo,npjglo), ztrpvssh(npiglo,npjglo))
  ALLOCATE ( glamf(npiglo,npjglo), gphif(npiglo,npjglo))

  glamf(:,:) = getvar(coordhgr, 'glamf',1,npiglo,npjglo)
  gphif(:,:) = getvar(coordhgr, 'gphif',1,npiglo,npjglo)

  ! create output fileset
  ncout =create(cfileoutnc, cfileu, npiglo,npjglo,1)
  ierr= createvar(ncout ,typvar,3, ipk,id_varout )
  ierr= putheadervar(ncout, cfileu,npiglo, npjglo,1,glamf, gphif)
  tim=getvar1d(cfileu,'time_counter',1)
  ierr=putvar1d(ncout,tim,1,'T')


  e1v(:,:) = getvar(coordhgr, 'e1v', 1,npiglo,npjglo)
  e2u(:,:) = getvar(coordhgr, 'e2u', 1,npiglo,npjglo)
  zmask(:,:) = getvar(cmask, 'fmask', 1,npiglo,npjglo)
  ! get rid of the free-slip/no-slip condition
  WHERE ( zmask >= 2 ) zmask = 1

  ztrpu(:,:)= 0.d0
  ztrpv(:,:)= 0.d0
  ztmp(:,:) = getvar(cfilet, 'sossheig', 1, npiglo,npjglo)
  ! put it on sshu and sshv
  DO jj=1,npjglo
    DO ji=1,npiglo-1
     sshu(ji,jj)=0.5*(ztmp(ji+1,jj) + ztmp(ji,jj) )
    ENDDO
  ENDDO
  DO jj=1,npjglo-1
    DO ji=1,npiglo
     sshv(ji,jj)=0.5*(ztmp(ji,jj+1) + ztmp(ji,jj) )
    ENDDO
  ENDDO
  
     
  DO jk = 1,npk
     ! Get  velocity  at jk
     zu(:,:)= getvar(cfileu, 'vozocrtx',  jk ,npiglo,npjglo)
     zv(:,:)= getvar(cfilev, 'vomecrty',  jk ,npiglo,npjglo)
     ! get e3 at level jk
     e3u(:,:) = getvar(coordzgr, 'e3u_ps', jk,npiglo,npjglo, ldiom=.true.)
     e3v(:,:) = getvar(coordzgr, 'e3v_ps', jk,npiglo,npjglo, ldiom=.true.)
     ! integrates vertically 
     ztrpu(:,:) = ztrpu(:,:) + zu(:,:)*e2u(:,:)*e3u(:,:)  ! zonal transport of each grid cell
     ztrpv(:,:) = ztrpv(:,:) + zv(:,:)*e1v(:,:)*e3v(:,:)  ! meridional transport of each grid cell
     IF ( jk == 1 ) THEN
     ztrpussh(:,:)= zu(:,:)*e2u(:,:)*sshu(:,:)            ! contrib of SSH
     ztrpvssh(:,:)= zv(:,:)*e1v(:,:)*sshv(:,:)            ! contrib of SSH
     ENDIF
  END DO  ! loop to next level

  ! compute psissh
  ! compute transport along line jj=jpj
  psissh(1,npjglo-1) = ztrpvssh(1,npjglo-1)
  DO ji=2,npiglo
     psissh(ji,npjglo-1) = psissh(ji-1,npjglo-1) + ztrpvssh(ji,npjglo-1)
  END DO
  DO jj=npjglo-2,1,-1
     DO ji=1,npiglo
        psissh(ji,jj)=psissh(ji,jj+1)  + ztrpussh(ji,jj+1)
     END DO
  END DO

  ! compute transport along line jj=jpj
  psi1(1,npjglo) = ztrpv(1,npjglo)  
  DO ji=2,npiglo 
     psi1(ji,npjglo) = psi1(ji-1,npjglo) + ztrpv(ji,npjglo)
  END DO
  ! set transport to 0 on Africa ( ji=1241    ; jj= 364  )
  ! Then compute from N to S the transport using zonal contribution
     psi1(:,npjglo)=psi1(:,npjglo) - psi1(1241,npjglo)
  !
  DO jj=npjglo-1,1,-1
     DO ji=1,npiglo
        psi1(ji,jj)=psi1(ji,jj+1)  + ztrpu(ji,jj+1)
     END DO
  END DO


  IF ( lmoy ) THEN
     ! compute transport along line ji=jpi
     psi2(npiglo,npjglo) = psi1(npiglo,npjglo)
     DO jj=npjglo-1, 1, -1
        psi2(npiglo,jj) = psi2(npiglo,jj+1) + ztrpu(npiglo,jj+1)
     END DO
     ! Then compute from W to E the transport using meridional contribution
     DO jj=npjglo,1,-1
        DO ji=npiglo-1,1,-1
           psi2(ji,jj)=psi2(ji+1,jj)  - ztrpv(ji+1,jj)
        END DO
     END DO
     psi1=0.5*(psi1 +psi2)
  END IF

  IF ( lmask ) psi1=psi1*zmask
  IF ( lmask ) psissh=psissh*zmask

  ierr = putvar(ncout, id_varout(1) ,SNGL(psi1),          1, npiglo, npjglo)
  ierr = putvar(ncout, id_varout(2) ,SNGL(psissh),        1, npiglo, npjglo)
  ierr = putvar(ncout, id_varout(3) ,SNGL(psissh+psi1),   1, npiglo, npjglo)
!   ierr = putvar(ncout, id_varout(2) ,SNGL(ztrpu),   1, npiglo, npjglo)
!   ierr = putvar(ncout, id_varout(3) ,SNGL(ztrpv),   1, npiglo, npjglo)
  istatus = closeout (ncout)

END PROGRAM cdfpsi_open
