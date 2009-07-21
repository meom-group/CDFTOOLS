PROGRAM cdfmoc_gsop_x
  !!-------------------------------------------------------------------
  !!               ***  PROGRAM cdfmoc_gsop  ***
  !!
  !!  **  Purpose  :  Compute the Meridional Overturning Cell (MOC)
  !!                  Components for GSOP intercomparison
  !!                  PARTIAL STEPS
  !!  
  !!  **  Method   :  The MOC is computed from the V velocity field, integrated
  !!                  from the bottom to the surface, then zonally averaged with
  !!                  eventual masking for oceanic basins.
  !!                  In the present version the masking corresponds to the global
  !!                  configuration. MOC for Global, Atlantic, Indo-Pacific, Indian,Pacific ocean
  !!                  Results are saved on moc.nc file with variables name respectively
  !!                  zomsfglo, zomsfatl, zomsfinp, zomsfind, zomsfpac
  !!
  !!
  !! history ;
  !!  Original :  J.M. Molines (jul. 2005) 
  !!              G.C. Smith ( Sep 2007) Added MOC decomposition following :
  !!                 Lee & Marotzke (1998), Baehr, Hirschi, Beismann, &  Marotzke (2004), Cabanes, Lee, & Fu (2007),  !!                  Koehl & Stammer (2007).
  !!                 See also the powerpoint presentation by Tony Lee at the third CLIVAR-GSOP intercomparison
  !!    available at : http://www.clivar.org/organization/gsop/synthesis/mit/talks/lee_MOC_comparison.ppt
  !!-------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------

  !! * Modules used
  USE cdfio
  USE eos

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: jpbasins  ! =5 modif Alb 29/11/08 pour fonctionner avec MERA
  INTEGER, PARAMETER :: jpgsop=4
  INTEGER   :: jgsop, jbasin, jj, jk ,ji                  !: dummy loop index
  INTEGER   :: ierr                                !: working integer
  INTEGER   :: narg, iargc                         !: command line 
  INTEGER   :: npiglo,npjglo, npk                  !: size of the domain
  INTEGER   :: ncout, np
  INTEGER   :: numout=10
  INTEGER, DIMENSION(:), ALLOCATABLE ::  ipk, id_varout         ! 
  INTEGER, DIMENSION(jpgsop) ::  ipk_gsop, id_varout_gsop         !
  INTEGER, DIMENSION(2)          ::  iloc

  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  e1u, e1v, e3v, gphiv, glamv, zv !:  metrics, velocity
  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  Hdep, vbt
  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  btsf
  REAL(KIND=4), DIMENSION (:,:,:),   ALLOCATABLE ::  btsf_x
  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  ztemp,zsal,tmask,umask,vmask, vgeoz
  REAL(KIND=8), DIMENSION (:,:),     ALLOCATABLE ::  zsig0
  REAL(KIND=4), DIMENSION (:,:,:),   ALLOCATABLE ::  vgeo,vgeosh,vfull,vmaskz,tmaskz
  REAL(KIND=4)   ::  rau0, grav, f0, fcor, zmsv, zphv, rpi
  REAL(KIND=4), DIMENSION (:,:,:),   ALLOCATABLE ::  e3vz
  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  dumlon              !: dummy longitude = 0.
  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  dumlat              !: latitude for i = north pole
  REAL(KIND=4), DIMENSION (:),       ALLOCATABLE ::  deptht, gdepw       !: deptw
  REAL(KIND=4), DIMENSION (:,:,:),   ALLOCATABLE ::  zmask               !:  jpbasins x npiglo x npjglo
  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  zzmask              !:  npiglo x npjglo
  REAL(KIND=4), DIMENSION (1)                    ::  tim

  REAL(KIND=8) ,DIMENSION(:,:,:) , ALLOCATABLE   ::  zomsf               !: jpbasins x npjglo x npk
  REAL(KIND=8) ,DIMENSION(:,:,:,:) , ALLOCATABLE ::  zomsf_x             !: jpbasins x npiglo x npjglo x npk
  REAL(KIND=8) ,DIMENSION(:,:,:) , ALLOCATABLE   ::  zomsf_gsop          !: jpgsop x npjglo x npk
  REAL(KIND=8) ,DIMENSION(:,:,:,:) , ALLOCATABLE ::  zomsf_gsop_x        !: jpgsop x npiglo x npjglo x npk

  CHARACTER(LEN=256) :: cfilet, cfilev , cfileoutnc='gsopmoc.nc'
  CHARACTER(LEN=256) :: coordhgr='mesh_hgr.nc',  coordzgr='mesh_zgr.nc',cbasinmask='new_maskglo.nc'
  CHARACTER(LEN=256) ,DIMENSION(jpgsop)     :: cvarname_gsop              !: array of var name for output
  TYPE(variable), DIMENSION(jpgsop) :: typvar       !: modif Alb 26/11/08 structure for attributes
  LOGICAL    :: llglo = .false.                            !: indicator for presence of new_maskglo.nc file
  INTEGER    :: istatus

  ! constants

  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfmoc  V file Tfile'
     PRINT *,' Computes the MOC for oceanic basins as described in new_maskglo.nc'
     PRINT *,' PARTIAL CELLS VERSION'
     PRINT *,' Files mesh_hgr.nc, mesh_zgr.nc ,new_maskglo.nc ,mask.nc '
     PRINT *,'  must be in the current directory'
     PRINT *,' Output on moc.nc: '
     PRINT *,'      variables zomsfglo  : Global ocean '
     PRINT *,'      variables zomsfatl  : Atlantic Ocean '
     PRINT *,'      variables zomsfinp  : Indo Pacific '
     PRINT *,'      variables zomsfind  : Indian Ocean alone'
     PRINT *,'      variables zomsfpac  : Pacific Ocean alone'
     STOP
  ENDIF

  CALL getarg (1, cfilev)
  npiglo= getdim (cfilev,'x')
  npjglo= getdim (cfilev,'y')
  npk   = getdim (cfilev,'depth')
  CALL getarg (2, cfilet)

 ! Detects newmaskglo file modif Alb 29/11/08 pour MERA
  INQUIRE( FILE='new_maskglo.nc', EXIST=llglo )
  IF (llglo) THEN
     jpbasins = 5
  ELSE
     jpbasins = 1
  ENDIF 

 ! define new variables for output ( must update att.txt)

!  typvar(1)%name= 'zobtmsfa'
!  typvar(1)%units='Sverdrup'
!  typvar(1)%missing_value=99999.
!  typvar(1)%valid_min= -1000.
!  typvar(1)%valid_max= 1000.
!  typvar(1)%scale_factor= 1.
!  typvar(1)%add_offset= 0.
!  typvar(1)%savelog10= 0.
!  typvar(1)%long_name='Barotropic_Merid_StreamFunction'
!  typvar(1)%short_name='zobtmsfa'
!  typvar(1)%online_operation='N/A'
!  typvar(1)%axis='TZY'

!  typvar(2)%name= 'zoshmsfa'
!  typvar(2)%units='Sverdrup'
!  typvar(2)%missing_value=99999.
!  typvar(2)%valid_min= -1000.
!  typvar(2)%valid_max= 1000.
!  typvar(2)%scale_factor= 1.
!  typvar(2)%add_offset= 0.
!  typvar(2)%savelog10= 0.
!  typvar(2)%long_name='GeoShear_Merid_StreamFunction'
!  typvar(2)%short_name='zoshmsfa'
!  typvar(2)%online_operation='N/A'
!  typvar(2)%axis='TZY'

!  typvar(3)%name= 'zoagmsfa'
!  typvar(3)%units='Sverdrup'
!  typvar(3)%missing_value=99999.
!  typvar(3)%valid_min= -1000.
!  typvar(3)%valid_max= 1000.
!  typvar(3)%scale_factor= 1.
!  typvar(3)%add_offset= 0.
!  typvar(3)%savelog10= 0.
!  typvar(3)%long_name='Ageo_Merid_StreamFunction'
!  typvar(3)%short_name='zoagmsfa'
!  typvar(3)%online_operation='N/A'
!  typvar(3)%axis='TZY'

!  typvar(4)%name= 'zomsfatl'
!  typvar(4)%units='Sverdrup'
!  typvar(4)%missing_value=99999.
!  typvar(4)%valid_min= -1000.
!  typvar(4)%valid_max= 1000.
!  typvar(4)%scale_factor= 1.
!  typvar(4)%add_offset= 0.
!  typvar(4)%savelog10= 0.
!  typvar(4)%long_name='Meridional_Overt.Cell_Atlantic'
!  typvar(4)%short_name='zomsfatl'
!  typvar(4)%online_operation='N/A'
!  typvar(4)%axis='TZY'

!  ipk_gsop(1) = npk
!  ipk_gsop(2) = npk
!  ipk_gsop(3) = npk
!  ipk_gsop(4) = npk

  typvar(1)%name= 'zobtmsfa_x'
  typvar(1)%units='Sverdrup'
  typvar(1)%missing_value=99999.
  typvar(1)%valid_min= -1000.
  typvar(1)%valid_max= 1000.
  typvar(1)%scale_factor= 1.
  typvar(1)%add_offset= 0.
  typvar(1)%savelog10= 0.
  typvar(1)%long_name='Barotropic_Merid_StreamFunction_x'
  typvar(1)%short_name='zobtmsfa_x'
  typvar(1)%online_operation='N/A'
  typvar(1)%axis='TZYX'

  typvar(2)%name= 'zoshmsfa_x'
  typvar(2)%units='Sverdrup'
  typvar(2)%missing_value=99999.
  typvar(2)%valid_min= -1000.
  typvar(2)%valid_max= 1000.
  typvar(2)%scale_factor= 1.
  typvar(2)%add_offset= 0.
  typvar(2)%savelog10= 0.
  typvar(2)%long_name='GeoShear_Merid_StreamFunction_x'
  typvar(2)%short_name='zoshmsfa_x'
  typvar(2)%online_operation='N/A'
  typvar(2)%axis='TZYX'

  typvar(3)%name= 'zoagmsfa_x'
  typvar(3)%units='Sverdrup'
  typvar(3)%missing_value=99999.
  typvar(3)%valid_min= -1000.
  typvar(3)%valid_max= 1000.
  typvar(3)%scale_factor= 1.
  typvar(3)%add_offset= 0.
  typvar(3)%savelog10= 0.
  typvar(3)%long_name='Ageo_Merid_StreamFunction_x'
  typvar(3)%short_name='zoagmsfa_x'
  typvar(3)%online_operation='N/A'
  typvar(3)%axis='TZYX'

  typvar(4)%name= 'zomsfatl_x'
  typvar(4)%units='Sverdrup'
  typvar(4)%missing_value=99999.
  typvar(4)%valid_min= -1000.
  typvar(4)%valid_max= 1000.
  typvar(4)%scale_factor= 1.
  typvar(4)%add_offset= 0.
  typvar(4)%savelog10= 0.
  typvar(4)%long_name='Meridional_Overt.Cell_Atlantic_x'
  typvar(4)%short_name='zomsfatl_x'
  typvar(4)%online_operation='N/A'
  typvar(4)%axis='TZYX'

  ipk_gsop(1) = npk
  ipk_gsop(2) = npk
  ipk_gsop(3) = npk
  ipk_gsop(4) = npk

  ! Allocate arrays
  ALLOCATE ( zmask(jpbasins,npiglo,npjglo) )
  ALLOCATE ( tmask(npiglo,npjglo) )
  ALLOCATE ( umask(npiglo,npjglo) )
  ALLOCATE ( vmask(npiglo,npjglo) )
  ALLOCATE ( vmaskz(npiglo,npjglo,npk) )
  ALLOCATE ( tmaskz(npiglo,npjglo,npk) )
  ALLOCATE ( zv(npiglo,npjglo) )
  ALLOCATE ( e1u(npiglo,npjglo),e1v(npiglo,npjglo),e3v(npiglo,npjglo), gphiv(npiglo,npjglo) ,glamv(npiglo,npjglo),gdepw(npk) )
  ALLOCATE ( Hdep(npiglo,npjglo), vbt(npiglo,npjglo) )
  ALLOCATE ( ztemp(npiglo,npjglo), zsal(npiglo,npjglo), zsig0(npiglo,npjglo) )
  ALLOCATE ( btsf(npjglo,npk) )
  ALLOCATE ( btsf_x(npiglo,npjglo,npk) )
  ALLOCATE ( vgeo(npiglo,npjglo,npk) )
  ALLOCATE ( vfull(npiglo,npjglo,npk) )
  ALLOCATE ( vgeosh(npiglo,npjglo,npk) )
  ALLOCATE ( vgeoz(npiglo,npjglo) )
  ALLOCATE ( e3vz(npiglo,npjglo,npk) )
  ALLOCATE ( zomsf(jpbasins, npjglo, npk) )
  ALLOCATE ( zomsf_x(jpbasins, npiglo,npjglo, npk) )
  ALLOCATE ( zomsf_gsop(jpgsop, npjglo, npk) )
  ALLOCATE ( zomsf_gsop_x(jpgsop, npiglo,npjglo, npk) )
  ALLOCATE ( dumlon(npiglo,npjglo) , dumlat(npiglo,npjglo))
  ALLOCATE ( deptht(npk) )
  ALLOCATE ( zzmask(npiglo,npjglo) )

  e1v(:,:) = getvar(coordhgr, 'e1v', 1,npiglo,npjglo) 
  e1u(:,:) = getvar(coordhgr, 'e1u', 1,npiglo,npjglo) 
  gphiv(:,:) = getvar(coordhgr, 'gphiv', 1,npiglo,npjglo)
  glamv(:,:) = getvar(coordhgr, 'glamv', 1,npiglo,npjglo)
  deptht(:) = getvare3(coordzgr, 'gdept',npk)
  gdepw(:) = getvare3(coordzgr, 'gdepw',npk)
  gdepw(:) = -1.*  gdepw(:)

  iloc=maxloc(gphiv)
  dumlat(:,:) = gphiv(:,:)
  dumlon(:,:) = glamv(:,:) 

  ! create output fileset
!  ncout =create(cfileoutnc, cfilev,1,npjglo,npk,cdep='depthw')
  ncout =create(cfileoutnc, cfilev, npiglo,npjglo,npk,cdep='depthw')
  ierr= createvar(ncout ,typvar,jpgsop, ipk_gsop,id_varout_gsop )
!  ierr= putheadervar(ncout, cfilev,1, npjglo,npk,pnavlon=dumlon,pnavlat=dumlat,pdep=gdepw)
  ierr= putheadervar(ncout, cfilev,npiglo, npjglo,npk,pnavlon=dumlon,pnavlat=dumlat,pdep=gdepw)
  tim=getvar1d(cfilev,'time_counter',1)
  ierr=putvar1d(ncout,tim,1,'T')

  ! reading the masks
  ! 1 : global ; 2 : Atlantic ; 3 : Indo-Pacif ; 4 : Indian ; 5 : Pacif

zmask=0
zmask(1,:,:)=getvar('mask.nc','vmask',1,npiglo,npjglo)

IF (llglo) THEN

  zmask(2,:,:)=getvar(cbasinmask,'tmaskatl',1,npiglo,npjglo)
  zmask(4,:,:)=getvar(cbasinmask,'tmaskind',1,npiglo,npjglo)
  zmask(5,:,:)=getvar(cbasinmask,'tmaskpac',1,npiglo,npjglo)
  zmask(3,:,:)=zmask(5,:,:)+zmask(4,:,:)
  ! ensure that there are no overlapping on the masks
  WHERE(zmask(3,:,:) > 0 ) zmask(3,:,:) = 1

ELSE

zmask(2,:,:)=getvar('mask.nc','tmask',1,npiglo,npjglo)

ENDIF

  ! initialize moc to 0
  zomsf(:,:,:) = 0.
  zomsf_x(:,:,:,:) = 0.
  zomsf_gsop(:,:,:) = 0.
  zomsf_gsop_x(:,:,:,:) = 0.
  vbt(:,:) = 0.0
  Hdep(:,:) = 0.0
  btsf(:,:) = 0.0
  btsf_x(:,:,:) = 0.0
  vgeo(:,:,:)=0.0
  vfull(:,:,:)=0.0

  ! Constants for geostrophic calc
  rau0 = 1025.0
  grav = 9.81
  rpi = 3.14159
  f0 = 2.0*(2.0*rpi)/(24.0*3600.0)

  ! Get velocities v and e3v_ps and masks at all levels
  DO jk = 1,npk
  zv(:,:)= getvar(cfilev, 'vomecrty',  jk ,npiglo,npjglo)
  vfull(:,:,jk) = zv(:,:)
  e3v(:,:) = getvar(coordzgr, 'e3v_ps', jk,npiglo,npjglo)
  e3vz(:,:,jk) = e3v(:,:)
  vmask(:,:)=getvar('mask.nc','vmask',jk,npiglo,npjglo)
  vmaskz(:,:,jk) = vmask(:,:)
  tmask(:,:)=getvar('mask.nc','tmask',jk,npiglo,npjglo)
  tmaskz(:,:,jk) = tmask(:,:)
  ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! CALCUL OF THE TOTAL AMOC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  DO jk = 1,npk-1
     ! Integrates 'zonally' (along i-coordinate)
     DO ji=1,npiglo
       ! For all basins 
       DO jbasin = 1, jpbasins
         DO jj=1,npjglo
           zomsf(jbasin,jj,jk)=zomsf(jbasin,jj,jk) - vmaskz(ji,jj,jk)*e1v(ji,jj)*e3vz(ji,jj,jk)*zmask(jbasin,ji,jj)*vfull(ji,jj,jk)/1.e6
           zomsf_x(jbasin,ji,jj,jk) = zomsf(jbasin,jj,jk)
         ENDDO ! loop to next latitude
       END DO  ! loop to next basin
     END DO    ! loop to next longitude
  ENDDO        ! loop to next level
  ! Integrates vertically from bottom to surface
  DO jk=npk-1 , 1 , -1
     zomsf(:,:,jk) = zomsf(:,:,jk+1) + zomsf(:,:,jk)
     zomsf_x(:,:,:,jk) = zomsf_x(:,:,:,jk+1) + zomsf_x(:,:,:,jk)
  END DO  ! loop to next level
  ! Save variable in zomsf_gsop
  zomsf_gsop(4,:,:) = zomsf(2,:,:)
  zomsf_gsop_x(4,:,:,:) = zomsf_x(2,:,:,:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! CALCUL OF THE BAROTROPIC AMOC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Calcul barotropic velocity, assuming barotropic velocity is zero at the
  ! bottom
  DO jk = 1,npk-1
     DO ji=1,npiglo
       DO jj=1,npjglo
         vbt(ji,jj) = vbt(ji,jj) + e3vz(ji,jj,jk)*zmask(2,ji,jj)*vfull(ji,jj,jk)*vmaskz(ji,jj,jk)  ! hardwire to jbasin=2
         Hdep(ji,jj) = Hdep(ji,jj) + e3vz(ji,jj,jk)*zmask(2,ji,jj)*vmaskz(ji,jj,jk)
       ENDDO ! loop to next latitude
     ENDDO   ! loop to next longitude
  ENDDO      ! loop to next level

  ! Normalize barotropic velocity
  DO ji=1,npiglo
    DO jj=1,npjglo
      IF ( Hdep(ji,jj) > 0.0 ) THEN
      vbt(ji,jj) = vbt(ji,jj)/Hdep(ji,jj)
      ELSE
        IF ( vbt(ji,jj) /= 0.0 ) THEN
        print *, 'Is something wrong?, ji,jj=',ji,jj
        ENDIF
           vbt(ji,jj) = 0.0
      ENDIF
    ENDDO ! loop to next latitude
  ENDDO   ! loop to next longitude
  
  ! Integrate zonally the barotropic velocity
  DO jk=1, npk
    DO jj=1,npjglo
      DO ji=1,npiglo
        btsf(jj,jk) = btsf(jj,jk) - vmaskz(ji,jj,jk)*e1v(ji,jj)*e3vz(ji,jj,jk)*zmask(2,ji,jj)*vbt(ji,jj)/1.e6
        btsf_x(ji,jj,jk) = btsf(jj,jk)
      END DO
    ENDDO
  ENDDO
  ! Now Integrate vertically to get barotropic AMOC
  DO jk=npk-1 , 1 , -1
     btsf(:,jk) = btsf(:,jk+1) + btsf(:,jk)
     btsf_x(:,:,jk) = btsf_x(:,:,jk+1) + btsf_x(:,:,jk)
  END DO  ! loop to next level
  ! Save variable in zomsf_gsop
  zomsf_gsop(1,:,:) = btsf(:,:)
  zomsf_gsop_x(1,:,:,:) = btsf_x(:,:,:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! CALCUL OF THE GEOSTROPHIC AMOC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  DO jk = 1,npk-1 
  ! Calculate density
  ztemp(:,:)= getvar(cfilet, 'votemper',  jk ,npiglo, npjglo)
  zsal(:,:) = getvar(cfilet, 'vosaline',  jk ,npiglo, npjglo)
  zzmask=1
  WHERE(zsal(:,:)* zmask(2,:,:) == 0 ) zzmask = 0
  ! geostrophic calculation must use in situ density gradient
  zsig0(:,:) = sigmai ( ztemp,zsal,deptht(jk),npiglo,npjglo )* zzmask(:,:)

  ! Calculate Geostrophic velocity
  ! value at v points is average of values at u points
    DO ji = 2, npiglo-1
      DO jj = 2, npjglo-1
         IF ( gphiv(ji,jj) == 0.0 ) THEN 
           vgeo(ji,jj,jk) = 0.0
         ELSE
           zmsv = 1. / MAX(  tmaskz(ji,jj+1,jk)*tmaskz(ji-1,jj+1,jk) + tmaskz(ji+1,jj+1,jk)*tmaskz(ji,jj+1,jk)   &
                         + tmaskz(ji,jj,jk)*tmaskz(ji-1,jj,jk) + tmaskz(ji+1,jj,jk)*tmaskz(ji,jj,jk) , 1.  )
           zphv = ( zsig0(ji  ,jj+1) - zsig0(ji-1,jj+1) ) * tmaskz(ji  ,jj+1,jk)*tmaskz(ji-1,jj+1,jk) / e1u(ji-1,jj+1)   &
                + ( zsig0(ji+1,jj+1) - zsig0(ji  ,jj+1) ) * tmaskz(ji+1,jj+1,jk)*tmaskz(ji  ,jj+1,jk) / e1u(ji  ,jj+1)   &
                + ( zsig0(ji  ,jj  ) - zsig0(ji-1,jj  ) ) * tmaskz(ji  ,jj  ,jk)*tmaskz(ji-1,jj  ,jk) / e1u(ji-1,jj  )   &
                + ( zsig0(ji+1,jj  ) - zsig0(ji  ,jj  ) ) * tmaskz(ji+1,jj  ,jk)*tmaskz(ji  ,jj  ,jk) / e1u(ji  ,jj  )
           zphv = (1. / rau0) * zphv * zmsv * vmaskz(ji,jj,jk)
           fcor = f0*SIN(rpi*gphiv(ji,jj)/180.0)
           vgeo(ji,jj,jk) = -grav*zphv/fcor*e3vz(ji,jj,jk)*zmask(2,ji,jj)
         ENDIF
      ENDDO ! loop to next latitude
    ENDDO   ! loop to next longitude
  ENDDO     ! loop to next level

  ! Vertical shear-velocity: Remove vertical average
  vgeoz(:,:) = 0.0
  vgeosh(:,:,:)=0.0
  DO ji=1, npiglo
    DO jj = 1, npjglo
      ! Integrate vertically to get geostrophic velocity referenced  to bottom
      DO jk = npk-1,1,-1
        vgeo(ji,jj,jk) = vgeo(ji,jj,jk+1) + vgeo(ji,jj,jk)
      ENDDO
      ! Calculate vertical sum
      DO jk = 1, npk
       vgeoz(ji,jj) = vgeoz(ji,jj) + vgeo(ji,jj,jk)*zmask(2,ji,jj)*e3vz(ji,jj,jk)*vmaskz(ji,jj,jk)
      ENDDO
      ! Remove total depth to get vertical mean
      IF ( Hdep(ji,jj) > 0.0 ) THEN
        vgeoz(ji,jj) = vgeoz(ji,jj)/Hdep(ji,jj)
      ELSE
        vgeoz(ji,jj) = 0.0
      ENDIF
      ! Remove vertical mean from geostrophic velocity to get geostrophic vertical shear velocity.
      DO jk = 1, npk
        vgeosh(ji,jj,jk) = zmask(2,ji,jj)*vgeo(ji,jj,jk) - vgeoz(ji,jj)
      ENDDO
    ENDDO  ! loop to next latitude
  ENDDO    ! loop to next longitude
  ! Calculate vertical shear geostrophic AMOC - integrate over x
  DO jk=1, npk 
    DO jj=1,npjglo
      DO ji=1,npiglo
        zomsf_gsop(2,jj,jk) = zomsf_gsop(2,jj,jk) - vmaskz(ji,jj,jk)*e1v(ji,jj)*e3vz(ji,jj,jk)*zmask(2,ji,jj)*vgeosh(ji,jj,jk)/1.e6
        zomsf_gsop_x(2,ji,jj,jk) = zomsf_gsop(2,jj,jk)
      ENDDO 
    ENDDO
  ENDDO 
  ! Integrate vertically to get GEOSTROPHIC AMOC
  DO jk=npk-1 , 1 , -1 
     zomsf_gsop(2,:,jk) = zomsf_gsop(2,:,jk+1) + zomsf_gsop(2,:,jk)
     zomsf_gsop_x(2,:,:,jk) = zomsf_gsop_x(2,:,:,jk+1) + zomsf_gsop_x(2,:,:,jk)
  END DO  ! loop to next level

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALCUL OF THE AGEOSTROPHIC AMOC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Calculate Ageostrophic meridional transport as residual
  !   SFag              =  MOC                  -     SFshear           -  SFbarotropic
  zomsf_gsop(3,:,:)     = zomsf_gsop(4,:,:)     - zomsf_gsop(2,:,:)     - zomsf_gsop(1,:,:)
  zomsf_gsop_x(3,:,:,:) = zomsf_gsop_x(4,:,:,:) - zomsf_gsop_x(2,:,:,:) - zomsf_gsop_x(1,:,:,:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  jj = 190 ; jk = 26 ; ji=50
  FIND26: DO jj=1,npjglo
    IF ( dumlat(1,jj) > 26.0 ) EXIT FIND26
  ENDDO FIND26
  print *, 'MOC:zomsf(2,jj,jk) = ', zomsf(2,jj,jk)
  print *, 'BT:zomsf_gsop(1,jj,jk) = ', zomsf_gsop(1,jj,jk)
  print *, 'SH:zomsf_gsop(2,jj,jk) = ', zomsf_gsop(2,jj,jk)
  print *, 'AG:zomsf_gsop(3,jj,jk) = ', zomsf_gsop(3,jj,jk)
  print *, 'MOC:zomsf_x(2,ji,jj,jk) = ', zomsf_x(2,ji,jj,jk)
  print *, 'BT:zomsf_gsop_x(1,ji,jj,jk) = ', zomsf_gsop_x(1,ji,jj,jk)
  print *, 'SH:zomsf_gsop_x(2,ji,jj,jk) = ', zomsf_gsop_x(2,ji,jj,jk)
  print *, 'AG:zomsf_gsop_x(3,ji,jj,jk) = ', zomsf_gsop_x(3,ji,jj,jk)

  !---------------------------------
  ! netcdf output 
  !---------------------------------

  !print *, 'Writing netcdf...'
  DO jgsop = 1, jpgsop
    DO jk=1,npk
!      ierr = putvar (ncout, id_varout_gsop(jgsop),REAL(zomsf_gsop(jgsop,:,jk)), jk,1,npjglo)
      ierr = putvar (ncout, id_varout_gsop(jgsop),REAL(zomsf_gsop_x(jgsop,:,:,jk)),jk,npiglo,npjglo)
    ENDDO
  ENDDO

  ierr = closeout(ncout)
 
END PROGRAM cdfmoc_gsop_x
