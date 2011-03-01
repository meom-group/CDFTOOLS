PROGRAM cdftransig_xy3d
  !!-------------------------------------------------------------------
  !!               ***  PROGRAM cdftransig_xy3d  ***
  !!
  !!  **  Purpose  :  calculates u and v transports  
  !!                  in rho coordinates.  produces a 3D field.
  !! allow two 3D arrays for more efficient reading 
  !!
  !! history ;
  !!  Original : A.M. Treguier (feb 2006)
  !!  Allow increased resolution in density in deeper layers (feb 2011)
  !!-------------------------------------------------------------------
  !! * Modules used
  USE cdfio
  USE eos

  !! * Local variables
  IMPLICIT NONE
!  FOR sigma 0 as the density coordinate 
!  REAL(KIND=4), PARAMETER :: pref =  0                       !:  reference for density 
!  INTEGER, PARAMETER      :: jpbin = 101                     !:   density  bins
!  REAL(KIND=4), PARAMETER :: s1min = 23.,s1scal=0.05         !:  reference for density 
!  CHARACTER (LEN=7)       :: clsigma = 'sigma_0'
!  FOR sigma 1 as the density coordinate 
  REAL(KIND=4), PARAMETER :: pref = 1000                      !:  reference for density 
  INTEGER, PARAMETER      :: jpbin = 93                       !:   density  bins
  REAL(KIND=4), PARAMETER :: s1min = 24.2,s1scal=0.1          !:  min sigma and delta_sigma
  REAL(KIND=4), PARAMETER :: s1zoom = 32.3,s1scalmin=0.05     !:  min sigma for increased resolution
  CHARACTER (LEN=7)       :: clsigma = 'sigma_1'
!  FOR sigma 1 as the density coordinate  for ACC region 
!  REAL(KIND=4), PARAMETER :: pref = 1000                    !:  reference for density 
!  INTEGER, PARAMETER      :: jpbin = 88                     !:   density  bins
!  REAL(KIND=4), PARAMETER :: s1min = 24.5,s1scal=0.1          !:  reference for density 
!  CHARACTER (LEN=7)       :: clsigma = 'sigma_1'
!  FOR sigma 2 as the density coordinate 
!  REAL(KIND=4), PARAMETER :: pref = 2000                     !:  reference for density 
!  INTEGER, PARAMETER      :: jpbin = 174                     !:   density  bins
!  REAL(KIND=4), PARAMETER :: s1min = 29,s1scal=0.05          !:  reference for density 
!  CHARACTER (LEN=7)       :: clsigma = 'sigma_2'
  
  INTEGER   :: jj, jk ,ji, jt , jib                !: dummy loop index
  INTEGER   :: ierr                                !: working integer
  INTEGER   :: narg, iargc                         !: command line 
  INTEGER   :: npiglo,npjglo, npk                  !: size of the domain
  INTEGER   :: ncout, ntags
  INTEGER, DIMENSION(2) ::   id_varout , ipk           !
  INTEGER, DIMENSION (:)          ,   ALLOCATABLE ::  itab                 !: look up table for density intervals
  INTEGER   :: jpsigmax , jitrans                                          !: dimension for itab, intermediate index
  REAL(KIND=4), DIMENSION (:,:)   ,   ALLOCATABLE ::  e1v,  gphiv          !:  2D x,y metrics, velocity
  REAL(KIND=4), DIMENSION (:,:)   ,   ALLOCATABLE ::  e2u                  !:  metrics, velocity
!!!                               
  REAL(KIND=4), DIMENSION (:,:)   ,   ALLOCATABLE ::  zt,zs, zv, e3v     !:  x,1,z arrays metrics, velocity
  REAL(KIND=4), DIMENSION (:,:)   ,   ALLOCATABLE ::  zu, e3u            !:  metrics, velocity
  REAL(KIND=4), DIMENSION (:,:)   ,   ALLOCATABLE ::  zmasku,zmaskv      !:  masks x,1,jpbin
  INTEGER,      DIMENSION (:,:)   ,   ALLOCATABLE ::  ibinu, ibinv   !: integer value corresponding to density for binning
  REAL(KIND=4), DIMENSION (:)     ,   ALLOCATABLE ::  gdept               !: array for depth of T points  
  REAL(KIND=4), DIMENSION (jpbin)                 ::  sigma               !: density coordinate, center of bins
  REAL(KIND=4), DIMENSION (jpbin+1)               ::  sig_edge            !: density coordinate, edge of bins.
  REAL(KIND=4),DIMENSION(1)                       ::  timean, tim
  REAL(KIND=4) ,DIMENSION(:,:)   , ALLOCATABLE   ::   zdensu, zdensv   !: density on u and v points 
!!!        3D arrays below are x,y,z
   REAL(KIND=8) ,DIMENSION(:,:,:) , ALLOCATABLE   ::  dusigsig,dvsigsig         !: cumulated transports,   
   REAL(KIND=8) ,DIMENSION(:,:)   , ALLOCATABLE   ::  dens2d 
   REAL(KIND=8)                                   :: total_time
   REAL(KIND=4)     :: sigtest
!!!  below 2D arrays npiglo,1
  CHARACTER(LEN=80) :: cfilev , cfilet,  cfileu, config , ctag 
  CHARACTER(LEN=80) :: cfileout='uvxysig.nc'
  CHARACTER(LEN=80) :: coordhgr='mesh_hgr.nc',  coordzgr='mesh_zgr.nc'
  CHARACTER(LEN=1)                  :: clanswer
  TYPE (variable), DIMENSION(2)     :: typvar     !: structure for attributes
    
  INTEGER    :: istatus 
  LOGICAL    :: lprint = .false.

  ! constants
!   lprint = .true. 
  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  ntags = narg-1
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdftransig_xyz CONFIG ''list_of_tags'' '
     PRINT *,' Computes the density transport in density space '
     PRINT *,' PARTIAL CELLS VERSION'
     PRINT *,' Files mesh_hgr.nc, mesh_zgr.nc, U, V, and T '
     PRINT *,'  must be in the current directory'
     PRINT *,'  Output on uvsigsig'
     PRINT *,'  variables vouxysig, vovxysig '
     STOP
  ENDIF
  !! Initialisation from 1st file (all file are assume to have the same geometry)
  CALL getarg (1, config)
  CALL getarg (2, ctag)
  WRITE(cfilev,'(a,"_",a,"_gridV.nc")') TRIM(config),TRIM(ctag)

  npiglo= getdim (cfilev,'x')
  npjglo= getdim (cfilev,'y')
  npk   = getdim (cfilev,'depth')

! define densities at middle of bins and edge
  jitrans = 0
  DO ji=1,jpbin
    sigtest  = s1min +(ji-0.5)*s1scal
    if ( sigtest > s1zoom ) THEN
       if ( jitrans == 0 ) jitrans = ji
       sigma(ji) = s1zoom + (ji-jitrans+0.5)*s1scalmin
    else
       sigma(ji) = sigtest
    endif
  ENDDO
  IF (lprint) print *, ' min density:',sigma(1), ' max density:', sigma(jpbin)
  IF (lprint) print *, ' verify sigma:', sigma
  sig_edge(1) = s1min
  DO ji=2,jpbin 
   sig_edge(ji) = 0.5* (sigma(ji)+sigma(ji-1))
  end do 
  sig_edge(jpbin+1) = sig_edge(jpbin)+s1scalmin
  IF (lprint) print *, ' sig_edge : ', sig_edge
 !
 !  define a lookup table array so that the density can be binned according to 
 !  the smallest interval s1scalmin
 jpsigmax = (sig_edge(jpbin+1)-sig_edge(1))/s1scalmin +1
 allocate ( itab(jpsigmax))
 itab(:) = 0
 DO ji=1,jpsigmax
    sigtest = s1min+ (ji-0.5)*s1scalmin
    DO jj=1,jpbin
      if ( sigtest > sig_edge(jj) .AND. sigtest <= sig_edge(jj+1) ) THEN
        itab(ji) = jj
      endif
    end do
 enddo
 IF (lprint) print *, ' jpsigmax=' , jpsigmax
 IF (lprint) print *, ' verify itab:', itab


 ! define new variables for output ( must update att.txt)
  ! define output variables
  typvar(1)%name= 'vouxysig'
  typvar(2)%name= 'vovxysig'

  typvar(1)%units='m/s'
  typvar(2)%units='m/s'
  typvar%missing_value=0.
  typvar%valid_min= -10.
  typvar%valid_max= 10.

  typvar(1)%long_name='Zonal_Velocity_sig_coord'
  typvar(2)%long_name='Meridional_Velocity_sig_coord'

  typvar(1)%short_name='vouxysig'
  typvar(2)%short_name='vovxysig'

  typvar%online_operation='N/A'
  typvar%axis='TZYX'

!                  output file has  jpbin sigma values
  ipk(:) = jpbin     

  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk, ' jpbin:', jpbin

  ! Allocate arrays
  ALLOCATE ( zv (npiglo,npjglo),  zu (npiglo,npjglo) )
  ALLOCATE ( zt (npiglo,npjglo),  zs (npiglo,npjglo) )
  ALLOCATE ( e3v(npiglo,npjglo),  e3u(npiglo,npjglo) )
  ALLOCATE ( ibinu(npiglo, npjglo), ibinv(npiglo, npjglo) )
  ALLOCATE ( e1v(npiglo,npjglo), gphiv(npiglo,npjglo) ,gdept(npk) )
  ALLOCATE ( e2u(npiglo,npjglo) )
  ALLOCATE ( dusigsig(npiglo,npjglo,jpbin), dvsigsig(npiglo,npjglo,jpbin))
  ALLOCATE ( dens2d(npiglo,npjglo) )
  ALLOCATE ( zdensu(npiglo,npjglo) ,zdensv(npiglo,npjglo) )
  ALLOCATE ( zmasku(npiglo,npjglo), zmaskv(npiglo,npjglo))


  e1v(:,:)   = getvar  (coordhgr, 'e1v', 1,npiglo,npjglo) 
  e2u(:,:)   = getvar  (coordhgr, 'e2u', 1,npiglo,npjglo) 
  gphiv(:,:) = getvar  (coordhgr, 'gphiv', 1,npiglo,npjglo)
  IF (lprint) PRINT *, '  read in hgr file OK'
  gdept(:)   = getvare3(coordzgr, 'gdept_0',npk)
 

  ! create output fileset
   IF (lprint) PRINT *, ' ready to create file:',trim( cfileout), ' from reference:',trim(cfilev )
   ncout =create(cfileout, cfilev, npiglo,npjglo,jpbin,cdep=clsigma)
   IF (lprint) print *, ' ncout=',ncout, ' ready to create variables:'
   ierr= createvar(ncout ,typvar,2, ipk  ,id_varout )
   IF (lprint) print *, ' ierr=',ierr, ' writing variables headers:'
   ierr= putheadervar(ncout, cfilev, npiglo, npjglo,jpbin,pdep=sigma)

   total_time=0

! initialize transport to 0
   dusigsig (:,:,:) = 0.; dvsigsig (:,:,:) =0;
!    loop on time and depth ---------------------------------------------------
! 
   
DO jk= 1, npk-1
   PRINT *, ' working on depth jk=',jk
   e3v(:,:) = getvar(coordzgr, 'e3v', jk, npiglo,npjglo )
   e3u(:,:) = getvar(coordzgr, 'e3u', jk, npiglo,npjglo )
 

    DO jt = 2, narg
 
      CALL getarg (jt, ctag)
      IF (lprint   ) PRINT *, ' working on  ctag=',trim(ctag)
 
      WRITE(cfilet,'(a,"_",a,"_gridT.nc")') TRIM(config),TRIM(ctag)
      WRITE(cfileu,'(a,"_",a,"_gridU.nc")') TRIM(config),TRIM(ctag)
      WRITE(cfilev,'(a,"_",a,"_gridV.nc")') TRIM(config),TRIM(ctag)

      IF (jk== 1 ) THEN
        tim=getvar1d(cfilet,'time_counter',1)
        total_time = total_time + tim(1)
      ENDIF

     ! Get velocities u, v  and mask   if first time slot only 
     zv(:,:)= getvar ( cfilev, 'vomecrty', jk ,npiglo,npjglo )
     zu(:,:)= getvar ( cfileu, 'vozocrtx', jk ,npiglo,npjglo )
     IF (jt == 2) THEN 
       zmasku(:,:)= 1; zmaskv(:,:)= 1;
       WHERE( zu == 0) zmasku(:,:)= 0.0;
       WHERE( zv == 0) zmaskv(:,:)= 0.0;
       IF (lprint  ) PRINT *, ' min,max u:',minval(zu),maxval(zu)
     ENDIF
!                     density  
     zt(:,:)= getvar ( cfilet, 'votemper', jk ,npiglo,npjglo )
     zs(:,:)= getvar ( cfilet, 'vosaline', jk ,npiglo,npjglo )
     
     IF ( pref == 0. ) THEN
        dens2d = sigma0(zt,zs,npiglo,npjglo)
     ELSE
        dens2d = sigmai(zt,zs,pref,npiglo,npjglo)
     ENDIF
!  density on u points masked by u  , single precision 
     zdensu(1:npiglo-1,:) = 0.5*( dens2d(1:npiglo-1,:) + dens2d(2:npiglo,:))
     zdensu(npiglo,:) = zdensu(2,:)
     zdensu(:,:) = zdensu(:,:) * zmasku(:,:)
!  density on v points masked by v  , single precision
     zdensv(:,1:npjglo-1) = 0.5*( dens2d(:,1:npjglo-1) + dens2d(:,2:npjglo) )
     zdensv(:,:) = zdensv(:,:) * zmaskv(:,:)

!  bins density - bins based on dens2d 
     DO jj=1,npjglo
        DO ji=1,npiglo
           jib   = ifix( (zdensu(ji,jj) - s1min)/s1scalmin )+1
           jib   = max( jib ,1   )
           jib   = min( jib,jpsigmax)
           ibinu(ji,jj) = itab (jib)
           jib   = ifix( (zdensv(ji,jj) - s1min)/s1scalmin )+1
           jib   = max( jib ,1   )
           jib   = min( jib,jpsigmax)
           ibinv(ji,jj) =  itab(jib)
        enddo
      enddo
       zu(:,:) = zu(:,:)*e3u(:,:)
       zv(:,:) = zv(:,:)*e3v(:,:)
       DO jj=1,npjglo
          DO ji=1,npiglo
             dusigsig(ji,jj,ibinu(ji,jj)) = dusigsig(ji,jj,ibinu(ji,jj))+ e2u(ji,jj)*zu(ji,jj) 
             dvsigsig(ji,jj,ibinv(ji,jj)) = dvsigsig(ji,jj,ibinv(ji,jj))+ e1v(ji,jj)*zv(ji,jj)            
          END DO
       END DO

!  -----------------------------------------end of loop on ctags
    END DO
!           
! -----------------  end of loop on jk
  END DO
   
   timean(1)= total_time/ntags
   ierr=putvar1d(ncout,timean,1,'T')
   DO jk=1, jpbin
     zt = dusigsig(:,:,jk) / ntags
     ierr = putvar (ncout, id_varout(1), zt, jk, npiglo, npjglo)
   ENDDO
   DO jk=1, jpbin
     zt = dvsigsig(:,:,jk) / ntags
     ierr = putvar (ncout, id_varout(2), zt, jk, npiglo, npjglo)
   ENDDO
 

  ierr = closeout(ncout)
 
END PROGRAM cdftransig_xy3d

   
