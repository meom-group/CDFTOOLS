PROGRAM cdfmht_gsop
  !!-------------------------------------------------------------------
  !!               ***  PROGRAM cdfmht_gsop  ***
  !!
  !!  **  Purpose  :  Compute the Meridional Heat Transport (MHT)
  !!                  Components for GSOP intercomparison
  !!                  PARTIAL STEPS
  !!  
  !!  **  Method   :  The MHT is computed from the V velocity field and T temperature field, integrated
  !!                  from the bottom to the surface.
  !!                  The MHT is decomposed into 3 components : BT, SH, AG.
  !!                  Results are saved on gsopmht.nc file with variables name respectively
  !!                  zomhtatl, zobtmhta, zoshmhta, zoagmhta
  !!
  !!
  !! history ;
  !!  Original :  J.M. Molines (jul. 2005) 
  !!              G.C. Smith ( Sep 2007) Added MOC decomposition following :
  !!                 Lee & Marotzke (1998), Baehr, Hirschi, Beismann, &  Marotzke (2004), Cabanes, Lee, & Fu (2007),  !!                  Koehl & Stammer (2007).
  !!                 See also the powerpoint presentation by Tony Lee at the third CLIVAR-GSOP intercomparison
  !!    available at : http://www.clivar.org/organization/gsop/synthesis/mit/talks/lee_MOC_comparison.ppt
  !!
  !!              A. Lecointre (Dec 2008) Replaced by a MHT decomposition
  !!
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

  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  e1u, e1v, e3v, gphiv, zv !:  metrics, velocity
  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  Hdep, vbt
  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  btht
  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  zt,zt_v,zsal,tmask,umask,vmask, vgeoz
  REAL(KIND=8), DIMENSION (:,:),     ALLOCATABLE ::  zsig0
  REAL(KIND=4), DIMENSION (:,:,:),   ALLOCATABLE ::  vgeo,vgeosh,vageosh,vfull,tfull,vmaskz,tmaskz
  REAL(KIND=4)   ::  rau0, grav, f0, fcor, zmsv, zphv, rpi
!  REAL(KIND=4)   ::  grav, f0, fcor, zmsv, zphv, rpi
  REAL(KIND=4), DIMENSION (:,:,:),   ALLOCATABLE ::  e3vz
  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  dumlon              !: dummy longitude = 0.
  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  dumlat              !: latitude for i = north pole
  REAL(KIND=4), DIMENSION (:),       ALLOCATABLE ::  deptht, gdepw       !: deptw
  REAL(KIND=4), DIMENSION (:,:,:),   ALLOCATABLE ::  zmask               !:  jpbasins x npiglo x npjglo
  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  zzmask              !:  npiglo x npjglo
  REAL(KIND=4), DIMENSION (1)                    ::  tim

  REAL(KIND=8) ,DIMENSION(:,:) , ALLOCATABLE ::  zomht                 !: jpbasins x npjglo
  REAL(KIND=8) ,DIMENSION(:,:) , ALLOCATABLE ::  zomht_gsop            !: jpgsop x npjglo
  REAL(KIND=8) ,DIMENSION(:,:) , ALLOCATABLE ::  zomht_geos_full       !: npjglo x npk
  REAL(KIND=8) ,DIMENSION(:,:) , ALLOCATABLE ::  zomht_ageos_full      !: npjglo x npk
  REAL(KIND=8) ,DIMENSION(:,:,:) , ALLOCATABLE ::  zomhtfull           !: jpbasin x npjglo x npk


  CHARACTER(LEN=100) :: cfilet, cfilev , cfileoutnc='gsopmht.nc'
  CHARACTER(LEN=80) :: coordhgr='mesh_hgr.nc',  coordzgr='mesh_zgr.nc',cbasinmask='new_maskglo.nc'
  CHARACTER(LEN=80) ,DIMENSION(jpgsop)     :: cvarname_gsop              !: array of var name for output
  TYPE(variable), DIMENSION(jpgsop) :: typvar       !: modif Alb 26/11/08 structure for attributes
  LOGICAL    :: llglo = .false.                            !: indicator for presence of new_maskglo.nc file
  INTEGER    :: istatus

  ! constants
  REAL(KIND=4),PARAMETER   ::  rho0=1000.,   rcp=4000.   ! rau0 en kg x m-3 et rcp en m2 x s-2 x degC-1

  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfmht_gsop  V file Tfile'
     PRINT *,' Computes the MHT for atlantic basin'
     PRINT *,' PARTIAL CELLS VERSION'
     PRINT *,' Files mesh_hgr.nc, mesh_zgr.nc ,new_maskglo.nc ,mask.nc '
     PRINT *,'  must be in the current directory'
     PRINT *,' Output on gsopmht.nc: '
     PRINT *,'      variables zomhtatl  : MHT Atlantic Ocean '
     PRINT *,'      variables zobtmhta  : Barotropic component '
     PRINT *,'      variables zoshmhta  : Vertical shear geostrophic component '
     PRINT *,'      variables zoagmhta  : vertical shear ageostrophic component (Ekman + residu)'
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

 ! define new variables for output

  typvar(1)%name= 'zobtmhta'
  typvar(1)%units='PetaWatt'
  typvar(1)%missing_value=99999.
  typvar(1)%valid_min= -1000.
  typvar(1)%valid_max= 1000.
  typvar(1)%scale_factor= 1.
  typvar(1)%add_offset= 0.
  typvar(1)%savelog10= 0.
  typvar(1)%long_name='Barotropic_Merid_HeatTransport'
  typvar(1)%short_name='zobtmhta'
  typvar(1)%online_operation='N/A'
  typvar(1)%axis='TY'

  typvar(2)%name= 'zoshmhta'
  typvar(2)%units='PetaWatt'
  typvar(2)%missing_value=99999.
  typvar(2)%valid_min= -1000.
  typvar(2)%valid_max= 1000.
  typvar(2)%scale_factor= 1.
  typvar(2)%add_offset= 0.
  typvar(2)%savelog10= 0.
  typvar(2)%long_name='GeoShear_Merid_HeatTransport'
  typvar(2)%short_name='zoshmhta'
  typvar(2)%online_operation='N/A'
  typvar(2)%axis='TY'

  typvar(3)%name= 'zoagmhta'
  typvar(3)%units='PetaWatt'
  typvar(3)%missing_value=99999.
  typvar(3)%valid_min= -1000.
  typvar(3)%valid_max= 1000.
  typvar(3)%scale_factor= 1.
  typvar(3)%add_offset= 0.
  typvar(3)%savelog10= 0.
  typvar(3)%long_name='Ageo_Merid_HeatTransport'
  typvar(3)%short_name='zoagmhta'
  typvar(3)%online_operation='N/A'
  typvar(3)%axis='TY'

  typvar(4)%name= 'zomhtatl'
  typvar(4)%units='PetaWatt'
  typvar(4)%missing_value=99999.
  typvar(4)%valid_min= -1000.
  typvar(4)%valid_max= 1000.
  typvar(4)%scale_factor= 1.
  typvar(4)%add_offset= 0.
  typvar(4)%savelog10= 0.
  typvar(4)%long_name='Meridional_HeatTransport_Atlantic'
  typvar(4)%short_name='zomhtatl'
  typvar(4)%online_operation='N/A'
  typvar(4)%axis='TY'

  ipk_gsop(1) = npk
  ipk_gsop(2) = npk
  ipk_gsop(3) = npk
  ipk_gsop(4) = npk

  ! Allocate arrays
  ALLOCATE ( zmask(jpbasins,npiglo,npjglo) )
  ALLOCATE ( tmask(npiglo,npjglo) )
  ALLOCATE ( tmaskz(npiglo,npjglo,npk) )
  ALLOCATE ( umask(npiglo,npjglo) )
  ALLOCATE ( vmask(npiglo,npjglo) )
  ALLOCATE ( vmaskz(npiglo,npjglo,npk) )
  ALLOCATE ( zv(npiglo,npjglo) )
  ALLOCATE ( vfull(npiglo,npjglo,npk) )
  ALLOCATE ( zt(npiglo,npjglo), zt_v(npiglo,npjglo) ) ! temperature au point T et au point V
  ALLOCATE ( tfull(npiglo,npjglo,npk) )  ! temperature au point V
  ALLOCATE ( e1u(npiglo,npjglo),e1v(npiglo,npjglo),e3v(npiglo,npjglo), gphiv(npiglo,npjglo) ,gdepw(npk) )
  ALLOCATE ( Hdep(npiglo,npjglo), vbt(npiglo,npjglo) )
  ALLOCATE ( e3vz(npiglo,npjglo,npk) )
  ALLOCATE ( zomhtfull(jpbasins,npjglo,npk) )
  ALLOCATE ( zomht(jpbasins, npjglo) )
  ALLOCATE ( zomht_gsop(jpgsop, npjglo) )
  ALLOCATE ( btht(npjglo,npk) )
  ALLOCATE ( zsal(npiglo,npjglo), zsig0(npiglo,npjglo) )
  ALLOCATE ( deptht(npk) )
  ALLOCATE ( dumlon(1,npjglo) , dumlat(1,npjglo))
  ALLOCATE ( zzmask(npiglo,npjglo) )
  ALLOCATE ( vgeo(npiglo,npjglo,npk) )
  ALLOCATE ( vgeoz(npiglo,npjglo) )
  ALLOCATE ( vgeosh(npiglo,npjglo,npk) )
  ALLOCATE ( zomht_geos_full(npjglo,npk) )
  ALLOCATE ( vageosh(npiglo,npjglo,npk) )
  ALLOCATE ( zomht_ageos_full(npjglo,npk) )

  e1v(:,:) = getvar(coordhgr, 'e1v', 1,npiglo,npjglo) 
  e1u(:,:) = getvar(coordhgr, 'e1u', 1,npiglo,npjglo) 
  gphiv(:,:) = getvar(coordhgr, 'gphiv', 1,npiglo,npjglo)
  deptht(:) = getvare3(coordzgr, 'gdept',npk)
  gdepw(:) = getvare3(coordzgr, 'gdepw',npk)
  gdepw(:) = -1.*  gdepw(:)

  iloc=maxloc(gphiv)
  dumlat(1,:) = gphiv(iloc(1),:)
  dumlon(:,:) = 0.   ! set the dummy longitude to 0

  ! create output fileset
  ncout =create(cfileoutnc, cfilev,1,npjglo,1,cdep='depthw')
  ierr= createvar(ncout ,typvar,jpgsop, ipk_gsop,id_varout_gsop )
  ierr= putheadervar(ncout, cfilev,1, npjglo,1,pnavlon=dumlon,pnavlat=dumlat,pdep=gdepw)
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

  ! initialize mht to 0
  zomht(:,:) = 0.
  zomhtfull(:,:,:) = 0.
  zomht_gsop(:,:) = 0.
  vbt(:,:) = 0.0
  Hdep(:,:) = 0.0
  btht(:,:) = 0.0
  vgeo(:,:,:)=0.0
  vfull(:,:,:)=0.0
  tfull(:,:,:)=0.0

  ! Constants for geostrophic calc
  rau0 = 1025.0
  grav = 9.81
  rpi = 3.14159
  f0 = 2.0*(2.0*rpi)/(24.0*3600.0)

  ! Get velocities v and temperature T and e3v_ps and masks at all levels
  DO jk = 1,npk
  vmask(:,:)=getvar('mask.nc','vmask',jk,npiglo,npjglo)
  vmaskz(:,:,jk) = vmask(:,:)
  tmask(:,:)=getvar('mask.nc','tmask',jk,npiglo,npjglo)
  tmaskz(:,:,jk) = tmask(:,:)
  zv(:,:)= getvar(cfilev, 'vomecrty',  jk ,npiglo,npjglo)  ! au point V
  vfull(:,:,jk) = zv(:,:)                                  ! au point V
  zt(:,:)= getvar(cfilet, 'votemper', jk,npiglo,npjglo)    ! au point T
    DO ji = 1,npiglo                                       ! mettre la temperature au point V
      DO jj = 1,npjglo-1
        zt_v(ji,jj)= ((zt(ji,jj) + zt(ji,jj+1)) * tmask(ji,jj) * tmask(ji,jj+1))/2
      END DO
    END DO
  tfull(:,:,jk)= zt_v(:,:)                                ! au point V
  e3v(:,:) = getvar(coordzgr, 'e3v_ps', jk,npiglo,npjglo)
  e3vz(:,:,jk) = e3v(:,:)
  ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! CALCUL OF THE TOTAL MHT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  DO jk = 1,npk-1
     ! MHT totale au point V
     ! integrates 'zonally' (along i-coordinate)
     DO ji=1,npiglo
       ! For all basins 
       DO jbasin = 1, jpbasins
         DO jj=1,npjglo
           zomhtfull(jbasin,jj,jk) = zomhtfull(jbasin,jj,jk) + vmaskz(ji,jj,jk)*e1v(ji,jj)*e3vz(ji,jj,jk)*zmask(jbasin,ji,jj)*vfull(ji,jj,jk)*tfull(ji,jj,jk)*rho0*rcp/1.e15
         ENDDO ! loop to next latitude
       END DO  ! loop to next basin
     END DO    ! loop to next longitude
  ENDDO        ! loop to next level

 ! integrates vertically from bottom to surface the total MHT
  DO jk=npk , 1 , -1 
     zomht(:,:) = zomht(:,:) + zomhtfull(:,:,jk)
  END DO  ! loop to next level
  ! Save variable in zomht_gsop
  zomht_gsop(4,:) = zomht(2,:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! CALCUL OF THE BAROTROPIC MHT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! Calculate ATLANTIC Barotropic velocity au point V
  DO jk = 1,npk-1
     DO ji=1,npiglo
       DO jj=1,npjglo
         vbt(ji,jj) = vbt(ji,jj) + e3vz(ji,jj,jk)*zmask(2,ji,jj)*vfull(ji,jj,jk)*vmaskz(ji,jj,jk)  ! hardwire to jbasin=2
         Hdep(ji,jj) = Hdep(ji,jj) + e3vz(ji,jj,jk)*zmask(2,ji,jj)*vmaskz(ji,jj,jk)
       ENDDO ! loop to next latitude
     ENDDO   ! loop to next longitude
  ENDDO      ! loop to next level

  ! Normalize Barotropic velocity
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
    ENDDO  ! loop to next latitude
  ENDDO    ! loop to next longitude

  ! Integrate zonally the barotropic velocity
  DO jk=1, npk 
    DO jj=1,npjglo
      DO ji=1,npiglo
        btht(jj,jk) = btht(jj,jk) + vmaskz(ji,jj,jk)*e1v(ji,jj)*e3vz(ji,jj,jk)*zmask(2,ji,jj)*vbt(ji,jj)*tfull(ji,jj,jk)*rho0*rcp/1.e15
      ENDDO 
    ENDDO
  ENDDO 

  ! Now Integrate vertically to get Barotropic Meridional Heat Transport
  DO jk=npk , 1 , -1 
     zomht_gsop(1,:)=zomht_gsop(1,:) + btht(:,jk)
  END DO  ! loop to next level

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! CALCUL OF THE GEOSTROPHIC MHT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! reinitialiser la temperature au point T a 0
zt(:,:)=0.0
  DO jk = 1,npk-1 
     ! Calculate density !! attention, density est au point U, il faut la mettre au point V
     zsal(:,:) = getvar(cfilet, 'vosaline',  jk ,npiglo, npjglo)
     zt(:,:)= getvar(cfilet, 'votemper', jk,npiglo,npjglo)    ! au point T

     zzmask=1
     WHERE(zsal(:,:)* zmask(2,:,:) == 0 ) zzmask = 0
     ! geostrophic calculation must use in situ density gradient
     ! la il faut prendre la temperature au point T
     zsig0(:,:) = sigmai ( zt,zsal,deptht(jk),npiglo,npjglo )* zzmask(:,:)

     ! Calculate Geostrophic velocity
     ! value at v points is average of values at u points
     DO ji = 2, npiglo-1
       DO jj = 2, npjglo-1
         IF ( gphiv(ji,jj) == 0.0 ) THEN 
           vgeo(ji,jj,jk) = 0.0
         ELSE
           zmsv = 1. / MAX(  tmaskz(ji  ,jj+1,jk)*tmaskz(ji-1,jj+1,jk) + tmaskz(ji+1,jj+1,jk)*tmaskz(ji  ,jj+1,jk)   &
                         + tmaskz(ji,  jj  ,jk)*tmaskz(ji-1,jj  ,jk) + tmaskz(ji+1,jj  ,jk)*tmaskz(ji  ,jj  ,jk) , 1.  )
           zphv = ( zsig0(ji  ,jj+1) - zsig0(ji-1,jj+1) ) * tmaskz(ji  ,jj+1,jk)*tmaskz(ji-1,jj+1,jk) / e1u(ji-1,jj+1)   &
                + ( zsig0(ji+1,jj+1) - zsig0(ji  ,jj+1) ) * tmaskz(ji+1,jj+1,jk)*tmaskz(ji  ,jj+1,jk) / e1u(ji  ,jj+1)   &
                + ( zsig0(ji  ,jj  ) - zsig0(ji-1,jj  ) ) * tmaskz(ji  ,jj  ,jk)*tmaskz(ji-1,jj  ,jk) / e1u(ji-1,jj  )   &
                + ( zsig0(ji+1,jj  ) - zsig0(ji  ,jj  ) ) * tmaskz(ji+1,jj ,jk )*tmaskz(ji  ,jj  ,jk) / e1u(ji  ,jj  )
           zphv = (1. / rau0) * zphv * zmsv * vmaskz(ji,jj,jk)
           fcor = f0*SIN(rpi*gphiv(ji,jj)/180.0)
           vgeo(ji,jj,jk) = -grav*zphv/fcor*e3vz(ji,jj,jk)*zmask(2,ji,jj)
         ENDIF
       ENDDO  ! loop to next latitude
     ENDDO    ! loop to next longitude
  ENDDO       ! loop to next level

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
  ! Calculate vertical shear MHT - integrate over x
  zomht_geos_full(:,:) = 0.0
  DO jk=1, npk 
    DO jj=1,npjglo
      DO ji=1,npiglo
        zomht_geos_full(jj,jk) = zomht_geos_full(jj,jk) + &
            & vmaskz(ji,jj,jk)*e1v(ji,jj)*e3vz(ji,jj,jk)*zmask(2,ji,jj)*vgeosh(ji,jj,jk)*tfull(ji,jj,jk)*rho0*rcp/1.e15
      END DO 
    ENDDO
  ENDDO 
  ! Integrate vertically the geostrophic MHT
  DO jk=npk , 1 , -1 
     zomht_gsop(2,:) = zomht_gsop(2,:) + zomht_geos_full(:,jk)
  END DO  ! loop to next level

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALCUL OF THE AGEOSTROPHIC MHT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  vageosh(:,:,:)=0.0
  ! Calculate Vageostrophique au point V
  DO jk=1,npk
  vageosh(:,:,jk)=vfull(:,:,jk)-vgeosh(:,:,jk)-vbt(:,:)
  END DO

  ! Calculate vertical shear ageostrophique streamfunction - integrate over x
  zomht_ageos_full(:,:) = 0.0
  DO jk=1, npk 
    DO jj=1,npjglo
      DO ji=1,npiglo
        zomht_ageos_full(jj,jk) = zomht_ageos_full(jj,jk) + vmaskz(ji,jj,jk)*e1v(ji,jj)*e3vz(ji,jj,jk)*zmask(2,ji,jj)*vageosh(ji,jj,jk)*tfull(ji,jj,jk)*rho0*rcp/1.e15
      END DO 
    ENDDO
  ENDDO 

  ! Now Integrate vertically to get streamfunction AGEOSTROPHIE
  DO jk=npk , 1 , -1 
     zomht_gsop(3,:) = zomht_gsop(3,:) + zomht_ageos_full(:,jk)
  END DO  ! loop to next level
  



! ! integrates vertically from bottom to surface the total MHT
!  DO jk=npk-1 , 1 , -1 
!     zomht(:,:,jk) = zomht(:,:,jk+1) + zomht(:,:,jk)
!  END DO  ! loop to next level

!  ! Normalize Barotropic velocity
!  DO ji=1,npiglo
!    DO jj=1,npjglo
!      IF ( Hdep(ji,jj) > 0.0 ) THEN
!        vbt(ji,jj) = vbt(ji,jj)/Hdep(ji,jj)
!      ELSE
!        IF ( vbt(ji,jj) /= 0.0 ) THEN
!         print *, 'Is something wrong?, ji,jj=',ji,jj
!        ENDIF
!        vbt(ji,jj) = 0.0
!      ENDIF
!    END DO 
!  ENDDO

!  ! Calculate Barotropic Meridional Heat Transport - integrate over x
!  DO jk=1, npk 
!    DO jj=1,npjglo
!      DO ji=1,npiglo
!        btht(jj,jk) = btht(jj,jk) - vmaskz(ji,jj,jk)*e1v(ji,jj)*e3vz(ji,jj,jk)*zmask(2,ji,jj)*vbt(ji,jj)*ztfull(ji,jj,jk)/1.e15
!      END DO 
!    ENDDO
!  ENDDO 

  ! Now Integrate vertically to get Barotropic Meridional Heat Transport
!  DO jk=npk-1 , 1 , -1 
!     btht(:,jk) = btht(:,jk+1) + btht(:,jk)
!  END DO  ! loop to next level

!  ! Calculate Vageostrophique au point V
!  DO jk=1,npk
!  vageosh(:,:,jk)=vfull(:,:,jk)-vgeosh(:,:,jk)-vbt(:,:)
!  END DO

!  ! Calculate vertical shear ageostrophique streamfunction - integrate over x
!  DO jk=1, npk 
!    DO jj=1,npjglo
!      DO ji=1,npiglo
!        zomht_gsop(3,jj,jk) = zomht_gsop(3,jj,jk) - vmaskz(ji,jj,jk)*e1v(ji,jj)*e3vz(ji,jj,jk)*zmask(2,ji,jj)*vageosh(ji,jj,jk)*ztfull(ji,jj,jk)/1.e15
!      END DO 
!    ENDDO
!  ENDDO 

!  ! Now Integrate vertically to get streamfunction AGEOSTROPHIE
!  DO jk=npk-1 , 1 , -1 
!     zomht_gsop(3,:,jk) = zomht_gsop(3,:,jk+1) + zomht_gsop(3,:,jk)
!  END DO  ! loop to next level
  
!  ! Save variables in zomht_gsop
!  zomht_gsop(1,:,:) = btht(:,:)
!  zomht_gsop(4,:,:) = zomht(2,:,:)

  jj = 190
  FIND26: DO jj=1,npjglo
    IF ( dumlat(1,jj) > 26.0 ) EXIT FIND26
  ENDDO FIND26
  print *, 'MHT:zomht_gsop(4,jj) = ', zomht_gsop(4,jj)
  print *, 'BT:zomht_gsop(1,jj) = ', zomht_gsop(1,jj)
  print *, 'SH:zomht_gsop(2,jj) = ', zomht_gsop(2,jj)
  print *, 'AG:zomht_gsop(3,jj) = ', zomht_gsop(3,jj)

  !---------------------------------
  ! netcdf output 
  !---------------------------------

  !print *, 'Writing netcdf...'
  DO jgsop = 1, jpgsop
      ierr = putvar (ncout, id_varout_gsop(jgsop),REAL(zomht_gsop(jgsop,:)), 1,1,npjglo)
  ENDDO

  ierr = closeout(ncout)
 
END PROGRAM cdfmht_gsop
