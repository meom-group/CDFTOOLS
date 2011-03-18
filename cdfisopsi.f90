PROGRAM cdfisopsi
  !!-------------------------------------------------------------------
  !!             ***  PROGRAM cdfisopsi  ***
  !!
  !!  **  Purpose: Compute a geostrophic streamfunction projected
  !!               on an isopycn (Ref: McDougall and ?, need reference)
  !!  
  !!  **  Method: read temp and salinity, compute sigmainsitu and sigma
  !!              at a reference level, projection of p,T,S on a given
  !!              isopycnal, compute specific volume anomaly and
  !!              integrates it.
  !!
  !! history: 
  !!     Original :   R. Dussin Dec 2010 (from various existing cdftools)
  !!
  !!-------------------------------------------------------------------
  !!  $Rev: 256 $
  !!  $Date: 2009-07-21 17:49:27 +0200 (mar. 21 juil. 2009) $
  !!  $Id: cdfsiginsitu.f90 256 2009-07-21 15:49:27Z molines $
  !!--------------------------------------------------------------
  !! * Modules used
  USE cdfio
  USE eos

  !! * Local variables
  IMPLICIT NONE
  INTEGER, PARAMETER  :: nvars=7
  INTEGER             :: jj,ji,jk,jt               !: dummy loop index
  INTEGER             :: ierr                      !: working integer
  INTEGER             :: narg, iargc               !: command line arguments
  INTEGER             :: npiglo,npjglo, npk ,npt   !: size of the domain
  INTEGER             :: k0                        !: 
  INTEGER, DIMENSION(nvars) ::  ipk, &             !: outptut variables : number of levels,
       &                        id_varout          !: ncdf varid's
  REAL(KIND=4) , DIMENSION (:,:), ALLOCATABLE :: ztemp, zsal , zssh, & !: Array to read a layer of data
       &                                         ztemp0, zsal0,      & !: Arrays for reference profile
       &                                         zsiginsitu , &        !: in-situ density
       &                                         zsig0, zsigsurf, &    !: potential density of ref profile and surface
       &                                         zmask, zdep           !: 2D mask at current level, level depths
  REAL(KIND=4),DIMENSION(:,:), ALLOCATABLE    :: v2d, ztempint, zsalint, zint, pint, alpha !: 2d working arrays
  REAL(KIND=4),DIMENSION(:,:,:), ALLOCATABLE  :: v3d, ztemp3, zsal3, sva3                  !: 3d array
  REAL(KIND=4),DIMENSION(:,:), ALLOCATABLE    :: e1t, e2t
  REAL(KIND=4),DIMENSION(:,:), ALLOCATABLE    :: deltapsi1, deltapsi2, psi0, psi, sva2
  REAL(KIND=4),DIMENSION(:), ALLOCATABLE      :: prof, tim              !: prof (m) and time (sec)
  REAL(KIND=4)                                :: P1, P2
  REAL(KIND=4)                                :: spval  !: missing value
  REAL(KIND=4)                                :: refdepth
  REAL(KIND=4)                                :: sigmaref
  REAL(KIND=4)                                :: tmean, smean, hmean, pmean

  CHARACTER(LEN=256) :: cfilet ,cfileout='isopsi.nc', coordhgr='mesh_hgr.nc', coordzgr='mesh_zgr.nc' 
  CHARACTER(LEN=256) :: cdum

  TYPE(variable) , DIMENSION(nvars) :: typvar         !: structure for attributes

  INTEGER    :: ncout
  INTEGER    :: istatus

  !--------------------------------------------------------------------
  !!  Read command line
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfisopsi ref_level sigma_ref gridT '
     PRINT *,' Output on isopsi.nc, variable soisopsi'
     PRINT *,' Depths are taken from input file '
     PRINT *,' requires mesh_hgr.nc and mesh_zgr.nc'
     STOP
  ENDIF

  CALL getarg (1, cdum)
  READ (cdum,*) refdepth
  CALL getarg (2, cdum)
  READ (cdum,*) sigmaref
  CALL getarg (3, cfilet)

  PRINT *, 'Potential density referenced at ', refdepth , ' meters'
  PRINT *, 'Isopycn for projection is ', sigmaref

  npiglo= getdim (cfilet,'x')
  npjglo= getdim (cfilet,'y')
  npk   = getdim (cfilet,'depth')
  npt   = getdim (cfilet,'time')

  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk
  PRINT *, 'npt   =', npt

  ALLOCATE ( prof(npk) , tim(npt) )
  ALLOCATE ( e1t(npiglo,npjglo)  , e2t(npiglo,npjglo) )

  e1t(:,:) = getvar(coordhgr, 'e1t'  ,1,npiglo,npjglo)
  e2t(:,:) = getvar(coordhgr, 'e2t'  ,1,npiglo,npjglo)

  !--------------------------------------------------------------------
  !!  Output file
  ipk(:)= 1  ! all variables are 2d
  typvar(1)%name= 'votemper_interp'
  typvar(1)%units='DegC'
  typvar(1)%missing_value=0.
  typvar(1)%valid_min= -2.
  typvar(1)%valid_max= 45.
  typvar(1)%long_name='Temperature interpolated on isopycnal layer'
  typvar(1)%short_name='votemper_interp'
  typvar(1)%online_operation='N/A'
  typvar(1)%axis='TZYX'

  typvar(2)%name= 'vosaline_interp'
  typvar(2)%units='PSU'
  typvar(2)%missing_value=0.
  typvar(2)%valid_min= 0.
  typvar(2)%valid_max= 50.
  typvar(2)%long_name='Salinity interpolated on isopycnal layer'
  typvar(2)%short_name='vosaline_interp'
  typvar(2)%online_operation='N/A'
  typvar(2)%axis='TZYX'

  typvar(3)%name= 'depth_interp'
  typvar(3)%units='meters'
  typvar(3)%missing_value=0.
  typvar(3)%valid_min= 0.0
  typvar(3)%valid_max= 8000.
  typvar(3)%long_name='Depth of the isopycnal layer'
  typvar(3)%short_name='depth_interp'
  typvar(3)%online_operation='N/A'
  typvar(3)%axis='TZYX'

  typvar(4)%name= 'soisopsi'
  typvar(4)%units=' m2s-2 (to be verified)'
  typvar(4)%missing_value=0.
  typvar(4)%valid_min= -500.
  typvar(4)%valid_max=  500.
  typvar(4)%long_name='Total streamfunction on the isopycnal layer'
  typvar(4)%short_name='soisopsi'
  typvar(4)%online_operation='N/A'
  typvar(4)%axis='TZYX'

  typvar(5)%name= 'soisopsi0'
  typvar(5)%units=' m2s-2 (to be verified)'
  typvar(5)%missing_value=0.
  typvar(5)%valid_min= -500.
  typvar(5)%valid_max=  500.
  typvar(5)%long_name='Contribution of the SSH'
  typvar(5)%short_name='soisopsi'
  typvar(5)%online_operation='N/A'
  typvar(5)%axis='TZYX'

  typvar(6)%name= 'soisopsi1'
  typvar(6)%units=' m2s-2 (to be verified)'
  typvar(6)%missing_value=0.
  typvar(6)%valid_min= -500.
  typvar(6)%valid_max=  500.
  typvar(6)%long_name='Contribution of specific volume anomaly vertical integration'
  typvar(6)%short_name='soisopsi'
  typvar(6)%online_operation='N/A'
  typvar(6)%axis='TZYX'

  typvar(7)%name= 'soisopsi2'
  typvar(7)%units=' m2s-2 (to be verified)'
  typvar(7)%missing_value=0.
  typvar(7)%valid_min= -500.
  typvar(7)%valid_max=  500.
  typvar(7)%long_name='Contribution of pressure term on the isopycnal layer'
  typvar(7)%short_name='soisopsi'
  typvar(7)%online_operation='N/A'
  typvar(7)%axis='TZYX'

  ! create output fileset

  ncout =create(cfileout, cfilet, npiglo,npjglo,npk)

  ierr= createvar(ncout ,typvar,nvars, ipk,id_varout ) 
  ierr= putheadervar(ncout, cfilet,npiglo, npjglo,npk)
  prof(:)=getvar1d(cfilet,'deptht',npk)
  tim=getvar1d(cfilet,'time_counter',npt)
  ierr=putvar1d(ncout,tim,npt,'T')

  spval=getatt(cfilet,'vosaline','missing_value')

  !---------------------------------------------------------------------------
  !! BEGIN LOOP ON TIME COUNTER
  DO jt=1,npt
     PRINT *,'time ',jt, tim(jt)/86400.,' days'

  !------------------------------------------------------------------------------
  ! 1. First we compute the potential density and store it into a 3d array
  ALLOCATE (ztemp(npiglo,npjglo), zsal(npiglo,npjglo), zmask(npiglo,npjglo))
  ALLOCATE (v3d(npiglo,npjglo,npk))

  DO jk = 1, npk
     zmask(:,:)=1.

     ztemp(:,:)= getvar(cfilet, 'votemper',  jk ,npiglo, npjglo, ktime=jt)
     zsal(:,:) = getvar(cfilet, 'vosaline',  jk ,npiglo, npjglo, ktime=jt)

     WHERE(zsal == spval ) zmask = 0

     v3d(:,:,jk) = sigmai ( ztemp,zsal,refdepth,npiglo,npjglo ) * zmask(:,:)

  END DO  ! loop to next level
  DEALLOCATE ( ztemp, zsal, zmask )

  !------------------------------------------------------------------------------
  ! 2. Projection of T,S and p on the chosen isopycnal layer (from cdfrhoproj)

  ALLOCATE ( alpha(npiglo,npjglo) )

  !! Compute coefficients
  DO ji=1,npiglo
     DO jj = 1, npjglo
        jk = 1
        !  Assume that rho (z) is increasing downward (no inversion)
        !     Caution with sigma0 at great depth !
           DO WHILE (sigmaref >=  v3d(ji,jj,jk) .AND. jk <= npk &
     &                .AND. v3d(ji,jj,jk) /=  spval )
              jk=jk+1
           END DO
           jk=jk-1
           k0=jk
           IF (jk .EQ. 0) THEN
              jk=1
              alpha(ji,jj) = 0.
           ELSE IF (v3d(ji,jj,jk+1) .EQ. spval ) THEN
              k0=0
              alpha(ji,jj) = 0.
           ELSE
           ! ... alpha is always in [0,1]. Adding k0 ( >=1 ) for saving space for k0
              alpha(ji,jj)= &
     &               (sigmaref-v3d(ji,jj,jk))/(v3d(ji,jj,jk+1)-v3d(ji,jj,jk)) + k0
           ENDIF
     END DO
  END DO

  DEALLOCATE (v3d)

  ! Working on temperature first
  ALLOCATE( ztempint(npiglo, npjglo), zint(npiglo, npjglo), pint(npiglo, npjglo) )
  ALLOCATE( ztemp3(npiglo, npjglo,npk) )

  DO jk=1,npk
     ztemp3(:,:,jk) = getvar(cfilet, 'votemper',  jk ,npiglo, npjglo, ktime=jt)
  ENDDO

  DO ji=1,npiglo
     DO jj=1,npjglo
        ! k0 is retrieved from alpha, taking the integer part.
        ! The remnant is alpha. 
        k0=INT(alpha(ji,jj))
        alpha(ji,jj) =  alpha(ji,jj) - k0
        IF (k0 /= 0) THEN
           P1=ztemp3(ji,jj,k0)
           P2=ztemp3(ji,jj,k0+1)
           IF (P1 /= spval .AND. P2 /= spval) THEN
               ztempint(ji,jj) = alpha(ji,jj)*P2  &
  &                            +(1-alpha(ji,jj))*P1
               zint(ji,jj)     = alpha(ji,jj)*prof(k0+1) &
  &                            +(1-alpha(ji,jj))*prof(k0)
           ELSE
               ztempint(ji,jj)=spval
               zint  (ji,jj)=spval
           ENDIF
        ELSE
            ztempint(ji,jj)=spval
            zint    (ji,jj)=spval
        ENDIF
        ! re-add k0 to alpha for the next computation
        alpha(ji,jj) =  alpha(ji,jj) + k0
     END DO
  END DO

  pint = zint / 10. ! pressure on the isopycnal layer = depth / 10.

  ierr = putvar(ncout, id_varout(1) ,ztempint, 1,npiglo, npjglo,ktime=jt)
  ierr = putvar(ncout, id_varout(3) ,zint, 1,npiglo, npjglo,ktime=jt)

  ! Working on salinity
  DEALLOCATE( ztemp3 )  
  ALLOCATE( zsalint(npiglo, npjglo) )
  ALLOCATE( zsal3(npiglo, npjglo,npk) )

  DO jk=1,npk
     zsal3(:,:,jk) = getvar(cfilet, 'vosaline',  jk ,npiglo, npjglo, ktime=jt)
  ENDDO

  DO ji=1,npiglo
     DO jj=1,npjglo
        ! k0 is retrieved from alpha, taking the integer part.
        ! The remnant is alpha. 
        k0=INT(alpha(ji,jj))
        alpha(ji,jj) =  alpha(ji,jj) - k0
        IF (k0 /= 0) THEN
           P1=zsal3(ji,jj,k0)
           P2=zsal3(ji,jj,k0+1)
           IF (P1 /= spval .AND. P2 /= spval) THEN
               zsalint(ji,jj) = alpha(ji,jj)*P2  &
  &                            +(1-alpha(ji,jj))*P1
           ELSE
               zsalint(ji,jj)=spval
           ENDIF
        ELSE
            zsalint(ji,jj)=spval
        ENDIF
        ! re-add k0 to alpha for the next computation
        alpha(ji,jj) =  alpha(ji,jj) + k0
     END DO
  END DO

  ierr = putvar(ncout, id_varout(2) ,zsalint, 1,npiglo, npjglo,ktime=jt)
  DEALLOCATE( zsal3 )

  ! 3. Compute means for T,S and depth on the isopycnal layer

  ALLOCATE( zmask(npiglo, npjglo) )
  zmask=1. ! define a new mask which correspond to the isopycnal layer
  WHERE( zint == 0. ) zmask = 0.

  tmean = SUM( ztempint * e1t * e2t * zmask ) / SUM( e1t * e2t * zmask )
  smean = SUM( zsalint  * e1t * e2t * zmask ) / SUM( e1t * e2t * zmask )
  hmean = SUM( zint     * e1t * e2t * zmask ) / SUM( e1t * e2t * zmask )
  pmean = SUM( pint     * e1t * e2t * zmask ) / SUM( e1t * e2t * zmask )

  DEALLOCATE ( ztempint, zsalint )

  ! 4. Compute specific volume anomaly
  ALLOCATE( sva3(npiglo,npjglo,npk) )
  ALLOCATE( zsiginsitu(npiglo,npjglo), zsig0(npiglo,npjglo) )
  ALLOCATE( ztemp(npiglo,npjglo),  zsal(npiglo,npjglo) )
  ALLOCATE( ztemp0(npiglo,npjglo), zsal0(npiglo,npjglo) )
  
  DO jk=1,npk
   
     ztemp(:,:) = getvar(cfilet, 'votemper',  jk ,npiglo, npjglo,ktime=jt)
     zsal (:,:) = getvar(cfilet, 'vosaline',  jk ,npiglo, npjglo,ktime=jt)

     ztemp0(:,:) = tmean
     zsal0 (:,:) = smean

     ! again land/sea mask
     zmask (:,:) = 1.
     WHERE( zsal == spval ) zmask = 0.

     zsiginsitu(:,:) = sigmai ( ztemp , zsal , prof(jk),npiglo,npjglo ) * zmask(:,:) ! in-situ density
     zsig0(:,:)      = sigmai ( ztemp0, zsal0, prof(jk),npiglo,npjglo ) * zmask(:,:) ! density of reference profile

     sva3(:,:,jk)    = ( 1. / zsiginsitu(:,:) ) - ( 1. / zsig0(:,:) )
     
  ENDDO

  DEALLOCATE( zsiginsitu, zsig0, ztemp0, zsal0 )

  ! 5. Integrates from surface to depth of isopycnal layer
  ALLOCATE( zdep(npiglo, npjglo), deltapsi1(npiglo, npjglo) )

  deltapsi1(:,:) = 0.

  DO jk=1, npk

    zdep(:,:) = getvar(coordzgr, 'e3t_ps', jk,npiglo,npjglo,ldiom=.true.)

    ! For each point we integrate from surface to zint(ji,jj) which is the depth
    ! of the isopycnal layer

    ! If isopycnal layer depth is below the current level
    WHERE( zint >= prof(jk) ) deltapsi1 = deltapsi1 - sva3(:,:,jk) * zdep / 10.
    ! If isopycnal layer is between current level and previous level
    WHERE( zint < prof(jk) .AND. zint > prof(jk-1) ) deltapsi1 = deltapsi1 &
                                                   & - sva3(:,:,jk) * ( zint - prof(jk-1) ) / 10.

  ENDDO

  ierr = putvar(ncout, id_varout(6) ,deltapsi1, 1,npiglo, npjglo,ktime=jt)

  DEALLOCATE( zdep )

  ! 6. Projection of the specific volume anomaly on the isopycnal layer
  ALLOCATE( sva2(npiglo,npjglo), deltapsi2(npiglo,npjglo) )

  DO ji=1,npiglo
     DO jj=1,npjglo
        ! k0 is retrieved from alpha, taking the integer part.
        ! The remnant is alpha. 
        k0=INT(alpha(ji,jj))
        alpha(ji,jj) =  alpha(ji,jj) - k0
        IF (k0 /= 0) THEN
           P1=sva3(ji,jj,k0)
           P2=sva3(ji,jj,k0+1)
           IF (P1 /= spval .AND. P2 /= spval) THEN
               sva2(ji,jj) = alpha(ji,jj)*P2  &
  &                            +(1-alpha(ji,jj))*P1
           ELSE
               sva2(ji,jj)=spval
           ENDIF
        ELSE
            sva2(ji,jj)=spval
        ENDIF
        ! re-add k0 to alpha for the next computation
        alpha(ji,jj) =  alpha(ji,jj) + k0
     END DO
  END DO

  deltapsi2 = ( pint - pmean ) * sva2

  ierr = putvar(ncout, id_varout(7) ,deltapsi2, 1,npiglo, npjglo,ktime=jt)

  DEALLOCATE ( sva3, sva2, alpha, zint, pint )

  ! 6. Finally we compute the surface streamfunction

  ALLOCATE(zssh(npiglo,npjglo) , zsigsurf(npiglo,npjglo), psi0(npiglo,npjglo) )
  
  ztemp   (:,:) = getvar(cfilet, 'votemper',  1 ,npiglo, npjglo,ktime=jt)
  zsal    (:,:) = getvar(cfilet, 'vosaline',  1 ,npiglo, npjglo,ktime=jt)
  zssh    (:,:) = getvar(cfilet, 'sossheig',  1 ,npiglo, npjglo,ktime=jt)

  ! land/sea mask at surface
  zmask (:,:) = 1.
  WHERE( zsal == spval ) zmask = 0.

  zsigsurf(:,:) = sigmai ( ztemp,zsal,prof(1),npiglo,npjglo ) * zmask(:,:)

  psi0 = zsigsurf * zssh * (9.81 / 1020. )
  ierr = putvar(ncout, id_varout(5) ,psi0, 1,npiglo, npjglo,ktime=jt)

  DEALLOCATE(zssh, zsigsurf, ztemp, zsal )

  ! 7. At least we are done with the computations
  ALLOCATE( psi(npiglo,npjglo) )

  ! final mask for output : mask the contribution of SSH where isopycn outcrops
  zmask=1.
  WHERE(deltapsi1 == spval ) zmask = 0.
  psi = ( psi0 * zmask ) + deltapsi1 + deltapsi2

  ierr = putvar(ncout, id_varout(4) ,psi, 1,npiglo, npjglo,ktime=jt)

  DEALLOCATE( psi, psi0, deltapsi1, deltapsi2, zmask )

  END DO  ! loop to next time

  istatus = closeout(ncout)

END PROGRAM cdfisopsi
