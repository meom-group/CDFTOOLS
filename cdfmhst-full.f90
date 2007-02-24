PROGRAM cdfmhst_full
  !!--------------------------------------------------------------------
  !!               ***  PROGRAM  cdfmhst_full  ***
  !!
  !!  **  Purpose  : Compute Meridional Heat Salt  Transport. 
  !!                 FULL STEP version
  !!  
  !!  **  Method   : Starts from the mean VT, VS fields computed by cdfvT
  !!                 Use a basin mask file
  !!
  !!
  !! history :
  !!      Original : J.M. Molines (jan. 2005)
  !!                 J.M. Molines apr. 2005 : use modules
  !!--------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !! * Modules used
  USE cdfio

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: jj,jk                            !: dummy loop index
  INTEGER   :: narg, iargc                         !: command line 
  INTEGER   :: npiglo,npjglo, npk                  !: size of the domain
  INTEGER   :: numout = 10

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE ::  zmask, e1v, gphiv, zvt, zvs !: mask, metrics
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE ::  e3t
  REAL(KIND=4) ,DIMENSION(:) ,  ALLOCATABLE ::  gphimean_glo, gphimean_atl, gphimean_pac, & 
       &                                        gphimean_ind, gphimean_aus, gphimean_med
  REAL(KIND=4)                              ::  zpoints

  REAL(KIND=8), DIMENSION (:,:), ALLOCATABLE :: zwk , ztrp, ztrps, zwks
  REAL(KIND=8) ,DIMENSION(:) , ALLOCATABLE ::  zonal_heat_glo, zonal_heat_atl, zonal_heat_pac,&
       &                                       zonal_heat_ind, zonal_heat_aus, zonal_heat_med
  REAL(KIND=8) ,DIMENSION(:) , ALLOCATABLE ::  zonal_salt_glo, zonal_salt_atl, zonal_salt_pac,&
       &                                       zonal_salt_ind, zonal_salt_aus, zonal_salt_med

  CHARACTER(LEN=80) :: cfilet ,cfileout='zonal_heat_trp.dat', cfileouts='zonal_salt_trp.dat'
  ! to be put in namelist eventually
  CHARACTER(LEN=80) :: coordhgr='mesh_hgr.nc',  coordzgr='mesh_zgr.nc', cbasinmask='new_maskglo.nc'


  ! constants
  REAL(KIND=4),PARAMETER   ::  rau0=1000.,   rcp=4000.

  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfmhst_full  VTfile '
     PRINT *,' Files mesh_hgr.nc, mesh_zgr.nc ,mask.nc, new_maskglo.nc must be in te current directory'
     PRINT *,' Output on zonal_heat_trp.dat and zonal_salt_trp.dat'
     STOP
  ENDIF

  CALL getarg (1, cfilet)
  npiglo= getdim (cfilet,'x')
  npjglo= getdim (cfilet,'y')
  npk   = getdim (cfilet,'depth')


  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk

  ! Allocate arrays
  ALLOCATE ( zwk(npiglo,npjglo) ,zmask(npiglo,npjglo),zvt(npiglo,npjglo) )
  ALLOCATE ( zwks(npiglo,npjglo) ,zvs(npiglo,npjglo) )
  ALLOCATE ( e1v(npiglo,npjglo),e3t(npk), gphiv(npiglo,npjglo))
  ALLOCATE ( ztrp(npiglo,npjglo))
  ALLOCATE ( ztrps(npiglo,npjglo))
  ALLOCATE ( zonal_heat_glo(npjglo), zonal_heat_atl(npjglo), zonal_heat_pac(npjglo))
  ALLOCATE ( zonal_heat_ind(npjglo), zonal_heat_aus(npjglo), zonal_heat_med(npjglo) )
  ALLOCATE ( zonal_salt_glo(npjglo), zonal_salt_atl(npjglo), zonal_salt_pac(npjglo))
  ALLOCATE ( zonal_salt_ind(npjglo), zonal_salt_aus(npjglo), zonal_salt_med(npjglo) )
  ALLOCATE ( gphimean_glo(npjglo) , gphimean_atl(npjglo), gphimean_pac(npjglo))
  ALLOCATE ( gphimean_ind(npjglo),gphimean_aus(npjglo),gphimean_med(npjglo))

  ! Read metrics and latitudes
  e1v(:,:)   = getvar(coordhgr, 'e1v', 1,npiglo,npjglo)
  gphiv(:,:) = getvar(coordhgr, 'gphiv', 1,npiglo,npjglo)
  e3t(:)     = getvare3(coordzgr,'e3t',npk) 

  ztrp(:,:)= 0
  ztrps(:,:)= 0
  DO jk = 1,npk
     PRINT *,'level ',jk
     ! Get temperature and salinity at jk
     zvt(:,:)= getvar(cfilet, 'vomevt',  jk ,npiglo,npjglo)
     zvs(:,:)= getvar(cfilet, 'vomevs',  jk ,npiglo,npjglo)

     ! get e3v at level jk
     zwk(:,:)  = zvt(:,:)*e1v(:,:)*e3t(jk)
     zwks(:,:) = zvs(:,:)*e1v(:,:)*e3t(jk)

     ! integrates vertically 
     ztrp(:,:)  = ztrp(:,:)  + zwk(:,:) * rau0*rcp
     ztrps(:,:) = ztrps(:,:) + zwks(:,:)  

  END DO  ! loop to next level

  ! Zonal mean with mask
  ! Atlantic 
  zmask(:,:)=getvar(cbasinmask,'tmaskatl',1,npiglo,npjglo)
  DO jj=1,npjglo
     zonal_heat_atl(jj) = SUM( ztrp(:,jj) *zmask(:,jj))
     zonal_salt_atl(jj) = SUM( ztrps(:,jj) *zmask(:,jj))
     zpoints =  SUM(zmask(:,jj))
     IF (  zpoints == 0 ) THEN
        gphimean_atl(jj)   = 999.
     ELSE
        gphimean_atl(jj)   = SUM( gphiv(:,jj)*zmask(:,jj)) / zpoints
     END IF
  END DO

  ! Pacific
  zmask(:,:)=getvar(cbasinmask,'tmaskpac',1,npiglo,npjglo)
  DO jj=1,npjglo
     zonal_heat_pac(jj)= SUM( ztrp(:,jj)*zmask(:,jj))
     zonal_salt_pac(jj)= SUM( ztrps(:,jj)*zmask(:,jj))
     zpoints =  SUM(zmask(:,jj))
     IF (  zpoints == 0 ) THEN
        gphimean_pac(jj)   = 999.
     ELSE
        gphimean_pac(jj)   = SUM( gphiv(:,jj)*zmask(:,jj)) / zpoints
     END IF
  END DO

  ! Indian
  zmask(:,:)=getvar(cbasinmask,'tmaskind',1,npiglo,npjglo)
  DO jj=1,npjglo
     zonal_heat_ind(jj)= SUM( ztrp(:,jj)*zmask(:,jj))
     zonal_salt_ind(jj)= SUM( ztrps(:,jj)*zmask(:,jj))
     zpoints =  SUM(zmask(:,jj))
     IF (  zpoints == 0 ) THEN
        gphimean_ind(jj)   = 999.
     ELSE
        gphimean_ind(jj)   = SUM( gphiv(:,jj)*zmask(:,jj)) / zpoints
     END IF
  END DO

  ! Austral
  zmask(:,:)=getvar(cbasinmask,'tmaskant',1,npiglo,npjglo)
  DO jj=1,npjglo
     zonal_heat_aus(jj)= SUM( ztrp(:,jj)*zmask(:,jj))
     zonal_salt_aus(jj)= SUM( ztrps(:,jj)*zmask(:,jj))
     zpoints =  SUM(zmask(:,jj))
     IF (  zpoints == 0 ) THEN
        gphimean_aus(jj)   = 999.
     ELSE
        gphimean_aus(jj)   = SUM( gphiv(:,jj)*zmask(:,jj)) / zpoints
     END IF
  END DO

  ! Med
! zmask(:,:)=getvar(cbasinmask,'tmaskmed',1,npiglo,npjglo)
! DO jj=1,npjglo
!    zonal_heat_med(jj)= SUM( ztrp(:,jj)*zmask(:,jj))
!    zonal_salt_med(jj)= SUM( ztrps(:,jj)*zmask(:,jj))
!    zpoints =  SUM(zmask(:,jj))
!    IF (  zpoints == 0 ) THEN
!       gphimean_med(jj)   = 999.
!    ELSE
!       gphimean_med(jj)   = SUM( gphiv(:,jj)*zmask(:,jj)) / zpoints
!    END IF
! END DO

  ! global 
  zmask(:,:)=getvar('mask.nc','vmask',1,npiglo,npjglo)
  DO jj=1,npjglo
     zonal_heat_glo(jj)= SUM( ztrp(:,jj)*zmask(:,jj))
     zonal_salt_glo(jj)= SUM( ztrps(:,jj)*zmask(:,jj))
     zpoints =  SUM(zmask(:,jj))
     IF (  zpoints == 0 ) THEN
        gphimean_glo(jj)   = 999.
     ELSE
        gphimean_glo(jj)   = SUM( gphiv(:,jj)*zmask(:,jj)) / zpoints
     END IF
  END DO

  ! Output file
  OPEN(numout,FILE=cfileout)
  WRITE(numout,*)' FULL STEP COMPUTATION '
  WRITE(numout,*)' Zonal heat transport (integrated along I-model coordinate) (in Pw)'
  WRITE(numout,*)' J        Global          Atlantic         Pacific          Indian           Mediteranean     Austral  '
  !               '  -89.123 -5.1234  -89.123 -5.1234  -89.123 -5.1234  -89.123 -5.1234  -89.123 -5.1234  -89.123 -5.1234
  DO jj=npjglo, 1, -1
     WRITE(numout,9000) jj, &
          gphimean_glo(jj),  zonal_heat_glo(jj)/1e15 , &
          gphimean_atl(jj),  zonal_heat_atl(jj)/1e15, &
          gphimean_pac(jj),  zonal_heat_pac(jj)/1e15, &
          gphimean_ind(jj),  zonal_heat_ind(jj)/1e15, &
          gphimean_med(jj),  zonal_heat_med(jj)/1e15, &
          gphimean_aus(jj),  zonal_heat_aus(jj)/1e15
  END DO
  CLOSE(numout)

  OPEN(numout,FILE=cfileouts)
  WRITE(numout,*)' Zonal salt transport (integrated alon I-model coordinate) (in 10^6 kg/s)'
  WRITE(numout,*)' J        Global          Atlantic         Pacific          Indian           Mediteranean     Austral  '
  !               '  -89.123 -5.1234  -89.123 -5.1234  -89.123 -5.1234  -89.123 -5.1234  -89.123 -5.1234  -89.123 -5.1234
  DO jj=npjglo, 1, -1
     WRITE(numout,9001) jj, &
          gphimean_glo(jj),  zonal_salt_glo(jj)/1e6 , &
          gphimean_atl(jj),  zonal_salt_atl(jj)/1e6, &
          gphimean_pac(jj),  zonal_salt_pac(jj)/1e6, &
          gphimean_ind(jj),  zonal_salt_ind(jj)/1e6, &
          gphimean_med(jj),  zonal_salt_med(jj)/1e6, &
          gphimean_aus(jj),  zonal_salt_aus(jj)/1e6
  END DO
  CLOSE(numout)


9000 FORMAT(I4,6(f9.3,f8.4))
9001 FORMAT(I4,6(f9.2,f9.3))


END PROGRAM cdfmhst_full
