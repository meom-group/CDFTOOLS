PROGRAM cdfmht_gsop
  !!======================================================================
  !!                     ***  PROGRAM  cdfmht_gsop  ***
  !!=====================================================================
  !!  **  Purpose  :  Compute the Meridional Heat Transport (MHT) components
  !!                  GSOP intercomparison project.
  !!  
  !!  **  Method   :  The MHT is computed from the V velocity field and T temperature field,
  !!                  integrated from the bottom to the surface. The velocity field is 
  !!                  decomposed into 3 components : Barotropic (BT), Shear geostrophic(SH)
  !!                  and Ageostrophic (AG). The program computes the corresponding MHT, thus
  !!                  saving 4 variables : Total, BT, SH, and AG MHT component.
  !!
  !!
  !! History : 2.0  : 07/2005 : J.M. Molines    : original MHT code
  !!                : 09/2007 : G.C. Smith      : Added MOC decomposition following :
  !!                : 12/2008 : A. Lecointre    : Replaced by a MHT decomposition
  !!         : 4.0  : 03/2017 : J.M. Molines    : rewriting at cdftools 4 standards
  !!
  !! References :  Lee & Marotzke (1998), Baehr, Hirschi, Beismann, &  Marotzke (2004),
  !!               Cabanes, Lee, & Fu (2007),  Koehl & Stammer (2007).
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!   varchk2       : check if variable is candidate for square mean
  !!   varchk3       : check if variable is candidate for cubic mean
  !!   zeromean      : substract mean value from input field
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  USE eos
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class transport
  !!-----------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4), PARAMETER                  :: jpgsop=4                   ! number of output variable
  INTEGER(KIND=4)                             :: nbasins                    ! number of basins to deal with
  INTEGER(KIND=4)                             :: jgsop, jbasin, jj, jk ,ji  ! dummy loop index
  INTEGER(KIND=4)                             :: ierr                       ! working integer
  INTEGER(KIND=4)                             :: narg, iargc, ijarg         ! command line 
  INTEGER(KIND=4)                             :: npiglo,npjglo, npk, npt    ! size of the domain
  INTEGER(KIND=4)                             :: ncout, np
  INTEGER(KIND=4)                             :: numout=10                  ! logical unit
  INTEGER(KIND=4), DIMENSION(jpgsop)          :: ipk_gsop, id_varout_gsop   ! netcdf output
  INTEGER(KIND=4), DIMENSION(2)               :: iloc

  REAL(KIND=4)                                :: rho0=1000.,   rcp=4000.    ! rau0 en kg x m-3 et rcp en m2 x s-2 x degC-1
  REAL(KIND=4)                                :: rau0, grav, f0, fcor       ! physical parameters
  REAL(KIND=4)                                :: zmsv, zphv, rpi            !
  REAL(KIND=4), DIMENSION(:),     ALLOCATABLE :: gdept, gdepw               ! deptht, deptw
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: e1u, e1v, e3v, gphiv, zv   !  metrics, velocity
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: hdep, vbt
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: btht
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: zt,zt_v,zsal
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: tmask,umask,vmask
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: vgeoz
  REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: zsig0
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: rlon                      ! dummy longitude = 0.
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: rlat                      ! latitude for i = north pole
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: zzmask                    ! npiglo x npjglo
  REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE :: zmask                     ! nbasins x npiglo x npjglo
  REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE :: vgeo,vgeosh,vageosh       ! npiglo x npjglo x npk
  REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE :: vfull,tfull               ! npiglo x npjglo x npk
  REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE :: vmaskz,tmaskz             ! npiglo x npjglo x npk
  REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE :: e3vz                      ! npiglo x npjglo x npk

  REAL(KIND=8), DIMENSION(:),     ALLOCATABLE :: dtim
  REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: dzomht                    ! nbasins x npjglo
  REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: dzomht_gsop               ! jpgsop x npjglo
  REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: dzomht_geos_full          ! npjglo x npk
  REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: dzomht_ageos_full         ! npjglo x npk
  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: dzomhtfull                ! nbasins x npjglo x npk

  CHARACTER(LEN=256)                          :: cf_tfil
  CHARACTER(LEN=256)                          :: cf_sfil
  CHARACTER(LEN=256)                          :: cf_vfil 
  CHARACTER(LEN=256)                          :: cf_out='gsopmht.nc'
  CHARACTER(LEN=256)                          :: cldum

  TYPE(variable), DIMENSION(jpgsop)           :: stypvar       

  LOGICAL                                     :: llglo = .FALSE.            ! indicator for presence of new_maskglo.nc file
  LOGICAL                                     :: lchk  = .FALSE.            ! missing files flag
  !!-----------------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()

  IF ( narg == 0 ) THEN
     PRINT *,' usage :  cdfmht_gsop -v V-file -t T-file [-s S-file] [-o OUT-file]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the meridional heat transport(MHT) for the Atlantic basin.'
     PRINT *,'       Compute 3 components of the MHT :'
     PRINT *,'          - Barotropic component'
     PRINT *,'          - Vertical shear geostrophic component '
     PRINT *,'          - Vertical shear ageostrophic component (Ekman + residual)'
     PRINT *,'      '
     PRINT *,'     REMARKS :'
     PRINT *,'       This program has been ported to CDFTOOLS4, without major changes.'
     PRINT *,'       It should work as before but is not optimized for memory (lot of'
     PRINT *,'       3D arrays declared). '
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -v V-file : name of the meridional velocity file.' 
     PRINT *,'       -t T-file : name of the temperature file.'
     PRINT *,'          If salinity not in T-file use -s option.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-s S-file]   : specify salinity file if not T-file.'
     PRINT *,'       [-o OUT-file] : output file name instead of ',TRIM(cf_out)
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       ',TRIM(cn_fhgr),', ',TRIM(cn_fzgr),' and ',TRIM(cn_fbasins) 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables :  zomhtatl  : MHT Atlantic Ocean '
     PRINT *,'                      zobtmhta  : Barotropic component '
     PRINT *,'                      zoshmhta  : Vertical shear geostrophic component '
     PRINT *,'                      zoagmhta  : vertical shear ageostrophic component'
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'       cdfmhst (compute MHT without decomposition), cdfmoc' 
     PRINT *,'      '
     STOP 
  ENDIF

  ijarg   = 1 
  cf_sfil = 'none'
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-v'   ) ; CALL getarg(ijarg, cf_vfil ) ; ijarg=ijarg+1
     CASE ( '-t'   ) ; CALL getarg(ijarg, cf_tfil ) ; ijarg=ijarg+1
        ! option
     CASE ( '-s'   ) ; CALL getarg(ijarg, cf_sfil ) ; ijarg=ijarg+1
     CASE ( '-o'   ) ; CALL getarg(ijarg, cf_out     ) ; ijarg=ijarg+1
     CASE DEFAULT    ; PRINT *, ' ERROR : ', TRIM(cldum),' : unknown option.'; STOP 99 
     END SELECT
  ENDDO
  IF ( cf_sfil == 'none' ) cf_sfil = cf_tfil

  lchk = lchk .OR. chkfile(cn_fzgr)
  lchk = lchk .OR. chkfile(cn_fhgr)
  lchk = lchk .OR. chkfile(cn_fmsk)
  lchk = lchk .OR. chkfile(cn_fbasins)
  lchk = lchk .OR. chkfile(cf_vfil )
  lchk = lchk .OR. chkfile(cf_tfil )
  lchk = lchk .OR. chkfile(cf_sfil )
  IF ( lchk ) STOP 99 

  npiglo= getdim (cf_vfil,cn_x)
  npjglo= getdim (cf_vfil,cn_y)
  npk   = getdim (cf_vfil,cn_z)
  npt   = getdim (cf_vfil,cn_t)

  ! sanity check on npt
  IF (npt > 1) THEN
     PRINT *, 'E R R O R npt > 1 not yet coded, STOP'
     STOP 99
  END IF

  ! Detects newmaskglo file modif Alb 29/11/08 pour MERA
  INQUIRE( FILE=cn_fbasins, EXIST=llglo )
  IF (llglo) THEN ; nbasins = 5
  ELSE            ; nbasins = 1
  ENDIF

  ! Allocate arrays
  ALLOCATE ( zmask(nbasins,npiglo,npjglo) )
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
  ALLOCATE ( hdep(npiglo,npjglo), vbt(npiglo,npjglo) )
  ALLOCATE ( e3vz(npiglo,npjglo,npk) )
  ALLOCATE ( dzomhtfull(nbasins,npjglo,npk) )
  ALLOCATE ( dzomht(nbasins, npjglo) )
  ALLOCATE ( dzomht_gsop(jpgsop, npjglo) )
  ALLOCATE ( btht(npjglo,npk) )
  ALLOCATE ( zsal(npiglo,npjglo), zsig0(npiglo,npjglo) )
  ALLOCATE ( gdept(npk) ,dtim(npt))
  ALLOCATE ( rlon(1,npjglo) , rlat(1,npjglo))
  ALLOCATE ( zzmask(npiglo,npjglo) )
  ALLOCATE ( vgeo(npiglo,npjglo,npk) )
  ALLOCATE ( vgeoz(npiglo,npjglo) )
  ALLOCATE ( vgeosh(npiglo,npjglo,npk) )
  ALLOCATE ( dzomht_geos_full(npjglo,npk) )
  ALLOCATE ( vageosh(npiglo,npjglo,npk) )
  ALLOCATE ( dzomht_ageos_full(npjglo,npk) )

  e1v(:,:) = getvar(cn_fhgr, cn_ve1v, 1,npiglo,npjglo) 
  e1u(:,:) = getvar(cn_fhgr, cn_ve1u, 1,npiglo,npjglo) 
  gphiv(:,:) = getvar(cn_fhgr, cn_gphiv, 1,npiglo,npjglo)
  gdept(:) = getvare3(cn_fzgr, cn_gdept,npk)
  gdepw(:) = getvare3(cn_fzgr, cn_gdepw,npk)
  gdepw(:) = -1.*  gdepw(:)

  iloc=MAXLOC(gphiv)
  rlat(1,:) = gphiv(iloc(1),:)
  rlon(:,:) = 0.   ! set the dummy longitude to 0

  CALL CreateOutput

  ! reading the masks
  ! 1 : global ; 2 : Atlantic ; 3 : Indo-Pacif ; 4 : Indian ; 5 : Pacif

  zmask=0
  zmask(1,:,:)=getvar(cn_fmsk,cn_vmask,1,npiglo,npjglo)

  IF (llglo) THEN
     zmask(2,:,:)=getvar(cn_fbasins,cn_tmaskatl,1,npiglo,npjglo)
     zmask(4,:,:)=getvar(cn_fbasins,cn_tmaskind,1,npiglo,npjglo)
     zmask(5,:,:)=getvar(cn_fbasins,cn_tmaskpac,1,npiglo,npjglo)
     zmask(3,:,:)=zmask(5,:,:)+zmask(4,:,:)
     ! ensure that there are no overlapping on the masks
     WHERE(zmask(3,:,:) > 0 ) zmask(3,:,:) = 1
  ELSE
     zmask(2,:,:)=getvar(cn_fmsk, cn_tmask,1,npiglo,npjglo)
  ENDIF

  ! initialize mht to 0
  dzomht(:,:) = 0.
  dzomhtfull(:,:,:) = 0.
  dzomht_gsop(:,:) = 0.
  vbt(:,:) = 0.0
  hdep(:,:) = 0.0
  btht(:,:) = 0.0
  vgeo(:,:,:)=0.0
  vfull(:,:,:)=0.0
  tfull(:,:,:)=0.0

  ! Constants for geostrophic calc
  rau0 = 1025.0
  grav = 9.81
  rpi = ACOS(-1.)
  f0 = 2.0*(2.0*rpi)/(24.0*3600.0)

  ! Get velocities v and temperature T and e3v_ps and masks at all levels
  DO jk = 1,npk
     vmask(:,:)=getvar(cn_fmsk, cn_vmask,    jk,npiglo,npjglo) ; vmaskz(:,:,jk) = vmask(:,:)
     tmask(:,:)=getvar(cn_fmsk, cn_tmask,    jk,npiglo,npjglo) ; tmaskz(:,:,jk) = tmask(:,:)
     zv(:,:)  = getvar(cf_vfil, cn_vomecrty, jk,npiglo,npjglo) ; vfull(:,:,jk) = zv(:,:)
     zt(:,:)  = getvar(cf_tfil, cn_votemper, jk,npiglo,npjglo)  ! au point T
     DO ji = 1,npiglo   ! mettre la temperature au point V
        DO jj = 1,npjglo-1
           zt_v(ji,jj)= ((zt(ji,jj) + zt(ji,jj+1)) * tmask(ji,jj) * tmask(ji,jj+1))/2
        END DO
     END DO
     tfull(:,:,jk)= zt_v(:,:)                                ! au point V
     e3v(:,:) = getvar(cn_fzgr, cn_ve3v, jk,npiglo,npjglo) ; e3vz(:,:,jk) = e3v(:,:)
  ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! CALCUL OF THE TOTAL MHT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  DO jk = 1,npk-1
     ! MHT totale au point V
     ! integrates 'zonally' (along i-coordinate)
     DO ji=1,npiglo
        ! For all basins 
        DO jbasin = 1, nbasins
           DO jj=1,npjglo
              dzomhtfull(jbasin,jj,jk) = dzomhtfull(jbasin,jj,jk) + vmaskz(ji,jj,jk)*e1v(ji,jj)*e3vz(ji,jj,jk)*zmask(jbasin,ji,jj)*vfull(ji,jj,jk)*tfull(ji,jj,jk)*rho0*rcp/1.d15
           ENDDO ! loop to next latitude
        END DO  ! loop to next basin
     END DO    ! loop to next longitude
  ENDDO        ! loop to next level

  ! integrates vertically from bottom to surface the total MHT
  DO jk=npk , 1 , -1 
     dzomht(:,:) = dzomht(:,:) + dzomhtfull(:,:,jk)
  END DO  ! loop to next level
  ! Save variable in dzomht_gsop
  dzomht_gsop(4,:) = dzomht(2,:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! CALCUL OF THE BAROTROPIC MHT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Calculate ATLANTIC Barotropic velocity au point V
  DO jk = 1,npk-1
     DO ji=1,npiglo
        DO jj=1,npjglo
           vbt(ji,jj) = vbt(ji,jj) + e3vz(ji,jj,jk)*zmask(2,ji,jj)*vfull(ji,jj,jk)*vmaskz(ji,jj,jk)  ! hardwire to jbasin=2
           hdep(ji,jj) = hdep(ji,jj) + e3vz(ji,jj,jk)*zmask(2,ji,jj)*vmaskz(ji,jj,jk)
        ENDDO ! loop to next latitude
     ENDDO   ! loop to next longitude
  ENDDO      ! loop to next level

  ! Normalize Barotropic velocity
  DO ji=1,npiglo
     DO jj=1,npjglo
        IF ( hdep(ji,jj) > 0.0 ) THEN
           vbt(ji,jj) = vbt(ji,jj)/hdep(ji,jj)
        ELSE
           IF ( vbt(ji,jj) /= 0.0 ) THEN
              PRINT *, 'Is something wrong?, ji,jj=',ji,jj
           ENDIF
           vbt(ji,jj) = 0.0
        ENDIF
     ENDDO  ! loop to next latitude
  ENDDO    ! loop to next longitude

  ! Integrate zonally the barotropic velocity
  DO jk=1, npk 
     DO jj=1,npjglo
        DO ji=1,npiglo
           btht(jj,jk) = btht(jj,jk) + vmaskz(ji,jj,jk)*e1v(ji,jj)*e3vz(ji,jj,jk)*zmask(2,ji,jj)*vbt(ji,jj)*tfull(ji,jj,jk)*rho0*rcp/1.d15
        ENDDO
     ENDDO
  ENDDO

  ! Now Integrate vertically to get Barotropic Meridional Heat Transport
  DO jk=npk , 1 , -1 
     dzomht_gsop(1,:)=dzomht_gsop(1,:) + btht(:,jk)
  END DO  ! loop to next level

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! CALCUL OF THE GEOSTROPHIC MHT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! reinitialiser la temperature au point T a 0
  zt(:,:)=0.0
  DO jk = 1,npk-1 
     ! Calculate density !! attention, density est au point U, il faut la mettre au point V
     zsal(:,:) = getvar(cf_sfil, cn_vosaline,  jk ,npiglo, npjglo)
     zt(:,:)= getvar(cf_tfil, cn_votemper, jk,npiglo,npjglo)    ! au point T

     zzmask=1
     WHERE(zsal(:,:)* zmask(2,:,:) == 0 ) zzmask = 0
     ! geostrophic calculation must use in situ density gradient
     ! la il faut prendre la temperature au point T
     zsig0(:,:) = sigmai ( zt,zsal,gdept(jk),npiglo,npjglo )* zzmask(:,:)

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
        IF ( hdep(ji,jj) > 0.0 ) THEN
           vgeoz(ji,jj) = vgeoz(ji,jj)/hdep(ji,jj)
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
  dzomht_geos_full(:,:) = 0.0
  DO jk=1, npk 
     DO jj=1,npjglo
        DO ji=1,npiglo
           dzomht_geos_full(jj,jk) = dzomht_geos_full(jj,jk) + &
                & vmaskz(ji,jj,jk)*e1v(ji,jj)*e3vz(ji,jj,jk)*zmask(2,ji,jj)*vgeosh(ji,jj,jk)*tfull(ji,jj,jk)*rho0*rcp/1.e15
        END DO
     ENDDO
  ENDDO
  ! Integrate vertically the geostrophic MHT
  DO jk=npk , 1 , -1 
     dzomht_gsop(2,:) = dzomht_gsop(2,:) + dzomht_geos_full(:,jk)
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
  dzomht_ageos_full(:,:) = 0.0
  DO jk=1, npk 
     DO jj=1,npjglo
        DO ji=1,npiglo
           dzomht_ageos_full(jj,jk) = dzomht_ageos_full(jj,jk) + vmaskz(ji,jj,jk)*e1v(ji,jj)*e3vz(ji,jj,jk)*zmask(2,ji,jj)*vageosh(ji,jj,jk)*tfull(ji,jj,jk)*rho0*rcp/1.e15
        END DO
     ENDDO
  ENDDO

  ! Now Integrate vertically to get streamfunction AGEOSTROPHIE
  DO jk=npk , 1 , -1 
     dzomht_gsop(3,:) = dzomht_gsop(3,:) + dzomht_ageos_full(:,jk)
  END DO  ! loop to next level

  jj = 190
  FIND26: DO jj=1,npjglo
     IF ( rlat(1,jj) > 26.0 ) EXIT FIND26
  ENDDO FIND26
  PRINT *, 'MHT:dzomht_gsop(4,jj) = ', dzomht_gsop(4,jj)
  PRINT *, 'BT:dzomht_gsop(1,jj) = ', dzomht_gsop(1,jj)
  PRINT *, 'SH:dzomht_gsop(2,jj) = ', dzomht_gsop(2,jj)
  PRINT *, 'AG:dzomht_gsop(3,jj) = ', dzomht_gsop(3,jj)

  !---------------------------------
  ! netcdf output 
  !---------------------------------

  !print *, 'Writing netcdf...'
  DO jgsop = 1, jpgsop
     ierr = putvar (ncout, id_varout_gsop(jgsop),REAL(dzomht_gsop(jgsop,:)), 1,1,npjglo)
  ENDDO

  ierr = closeout(ncout)

CONTAINS

  SUBROUTINE CreateOutput
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutput  ***
    !!
    !! ** Purpose :  Create netcdf output file(s) 
    !!
    !! ** Method  :  Use stypvar global description of variables
    !!
    !!----------------------------------------------------------------------
    ipk_gsop(:) = npk

    stypvar(1)%cname= 'zobtmhta'
    stypvar(1)%cunits='PetaWatt'
    stypvar(1)%rmissing_value=99999.
    stypvar(1)%valid_min= -1000.
    stypvar(1)%valid_max= 1000.
    stypvar(1)%scale_factor= 1.
    stypvar(1)%add_offset= 0.
    stypvar(1)%savelog10= 0.
    stypvar(1)%clong_name='Barotropic_Merid_HeatTransport'
    stypvar(1)%cshort_name='zobtmhta'
    stypvar(1)%conline_operation='N/A'
    stypvar(1)%caxis='TY'

    stypvar(2)%cname= 'zoshmhta'
    stypvar(2)%cunits='PetaWatt'
    stypvar(2)%rmissing_value=99999.
    stypvar(2)%valid_min= -1000.
    stypvar(2)%valid_max= 1000.
    stypvar(2)%scale_factor= 1.
    stypvar(2)%add_offset= 0.
    stypvar(2)%savelog10= 0.
    stypvar(2)%clong_name='GeoShear_Merid_HeatTransport'
    stypvar(2)%cshort_name='zoshmhta'
    stypvar(2)%conline_operation='N/A'
    stypvar(2)%caxis='TY'

    stypvar(3)%cname= 'zoagmhta'
    stypvar(3)%cunits='PetaWatt'
    stypvar(3)%rmissing_value=99999.
    stypvar(3)%valid_min= -1000.
    stypvar(3)%valid_max= 1000.
    stypvar(3)%scale_factor= 1.
    stypvar(3)%add_offset= 0.
    stypvar(3)%savelog10= 0.
    stypvar(3)%clong_name='Ageo_Merid_HeatTransport'
    stypvar(3)%cshort_name='zoagmhta'
    stypvar(3)%conline_operation='N/A'
    stypvar(3)%caxis='TY'

    stypvar(4)%cname= 'zomhtatl'
    stypvar(4)%cunits='PetaWatt'
    stypvar(4)%rmissing_value=99999.
    stypvar(4)%valid_min= -1000.
    stypvar(4)%valid_max= 1000.
    stypvar(4)%scale_factor= 1.
    stypvar(4)%add_offset= 0.
    stypvar(4)%savelog10= 0.
    stypvar(4)%clong_name='Meridional_HeatTransport_Atlantic'
    stypvar(4)%cshort_name='zomhtatl'
    stypvar(4)%conline_operation='N/A'
    stypvar(4)%caxis='TY'

    ! create output fileset
    ncout = create(cf_out, cf_vfil,1,npjglo,1,cdep=cn_vdepthw)
    ierr  = createvar(ncout ,stypvar,jpgsop, ipk_gsop,id_varout_gsop )
    ierr  = putheadervar(ncout, cf_vfil,1, npjglo,1,pnavlon=rlon,pnavlat=rlat,pdep=gdepw)
    dtim  = getvar1d(cf_vfil,cn_vtimec,npt)
    ierr  = putvar1d(ncout,dtim,npt,'T')

  END SUBROUTINE CreateOutput

END PROGRAM cdfmht_gsop
