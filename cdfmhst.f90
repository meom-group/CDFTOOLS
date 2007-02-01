PROGRAM cdfmhst
  !!--------------------------------------------------------------------
  !!               ***  PROGRAM  cdfmhst  ***
  !!
  !!  **  Purpose  : Compute Meridional Heat Salt  Transport. 
  !!  
  !!  **  Method   : Starts from the mean VT, VS fields computed by cdfvT
  !!                 The program looks for the file "new_maskglo.nc". If it does not exist, 
  !!                 only the calculation over all the domain is performed (this is adequate 
  !!                 for a basin configuration like NATL4).
  !!
  !!
  !! history :
  !!      Original : J.M. Molines (jan. 2005)
  !!                 J.M. Molines apr. 2005 : use modules
  !!                 A.M. Treguier (april 2006) adaptation to NATL4 case 
  !!--------------------------------------------------------------------
  !! * Modules used
  USE cdfio

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: jj,jk                            !: dummy loop index
  INTEGER   :: narg, iargc                         !: command line 
  INTEGER   :: npiglo,npjglo, npk                  !: size of the domain
  INTEGER   :: numout = 10
  INTEGER, DIMENSION(2)          ::  iloc
  LOGICAL    :: llglo = .false.                !: indicator for presence of new_maskglo.nc file 


  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE ::  zmask, e1v, e3v ,gphiv, zvt, zvs !: mask, metrics
  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  dumlon              !: dummy longitude = 0.
  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  dumlat              !: latitude for i = north pole

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
     PRINT *,' Usage : cdfmhst  VTfile '
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
  ALLOCATE ( e1v(npiglo,npjglo),e3v(npiglo,npjglo), gphiv(npiglo,npjglo))
  ALLOCATE ( ztrp(npiglo,npjglo))
  ALLOCATE ( ztrps(npiglo,npjglo))
  ALLOCATE ( zonal_heat_glo(npjglo), zonal_heat_atl(npjglo), zonal_heat_pac(npjglo) )
  ALLOCATE ( zonal_heat_ind(npjglo), zonal_heat_aus(npjglo) , zonal_heat_med(npjglo) )
  ALLOCATE ( zonal_salt_glo(npjglo), zonal_salt_atl(npjglo), zonal_salt_pac(npjglo) )
  ALLOCATE ( zonal_salt_ind(npjglo), zonal_salt_aus(npjglo), zonal_salt_med(npjglo) )
  ALLOCATE ( dumlon(1,npjglo) , dumlat(1,npjglo))

  ! create output fileset
  e1v(:,:)   = getvar(coordhgr, 'e1v', 1,npiglo,npjglo)
  gphiv(:,:) = getvar(coordhgr, 'gphiv', 1,npiglo,npjglo)

  iloc=maxloc(gphiv)
  dumlat(1,:) = gphiv(iloc(1),:)
  dumlon(:,:) = 0.   ! set the dummy longitude to 0


  ztrp(:,:)= 0
  ztrps(:,:)= 0
  DO jk = 1,npk
     PRINT *,'level ',jk
     ! Get temperature and salinity at jk
     zvt(:,:)= getvar(cfilet, 'vomevt',  jk ,npiglo,npjglo)
     zvs(:,:)= getvar(cfilet, 'vomevs',  jk ,npiglo,npjglo)

     ! get e3v at level jk
     e3v(:,:)  = getvar(coordzgr, 'e3v_ps', jk,npiglo,npjglo)
     zwk(:,:)  = zvt(:,:)*e1v(:,:)*e3v(:,:)
     zwks(:,:) = zvs(:,:)*e1v(:,:)*e3v(:,:)

     ! integrates vertically 
     ztrp(:,:)  = ztrp(:,:)  + zwk(:,:) * rau0*rcp
     ztrps(:,:) = ztrps(:,:) + zwks(:,:)  

  END DO  ! loop to next level

  ! global 
  zmask(:,:)=getvar('mask.nc','vmask',1,npiglo,npjglo)
  DO jj=1,npjglo
     zonal_heat_glo(jj)= SUM( ztrp(2:npiglo-1,jj)*zmask(2:npiglo-1,jj))
     zonal_salt_glo(jj)= SUM( ztrps(2:npiglo-1,jj)*zmask(2:npiglo-1,jj))
     zpoints =  SUM(zmask(2:npiglo-1,jj))
  END DO

 !  Detects newmaskglo file 
  INQUIRE( FILE=cbasinmask, EXIST=llglo )
  
  IF ( llglo ) THEN
     ! Zonal mean with mask
     ! Atlantic 
     zmask(:,:)=getvar(cbasinmask,'tmaskatl',1,npiglo,npjglo)
     DO jj=1,npjglo
        zonal_heat_atl(jj) = SUM( ztrp(:,jj) *zmask(:,jj))
        zonal_salt_atl(jj) = SUM( ztrps(:,jj) *zmask(:,jj))
        zpoints =  SUM(zmask(:,jj))
     END DO

     ! Pacific
     zmask(:,:)=getvar(cbasinmask,'tmaskpac',1,npiglo,npjglo)
     DO jj=1,npjglo
        zonal_heat_pac(jj)= SUM( ztrp(:,jj)*zmask(:,jj))
        zonal_salt_pac(jj)= SUM( ztrps(:,jj)*zmask(:,jj))
        zpoints =  SUM(zmask(:,jj))
     END DO

     ! Indian
     zmask(:,:)=getvar(cbasinmask,'tmaskind',1,npiglo,npjglo)
     DO jj=1,npjglo
        zonal_heat_ind(jj)= SUM( ztrp(:,jj)*zmask(:,jj))
        zonal_salt_ind(jj)= SUM( ztrps(:,jj)*zmask(:,jj))
        zpoints =  SUM(zmask(:,jj))
     END DO

     ! Austral
     zonal_heat_aus = 0.
     zonal_salt_aus = 0.
!    zmask(:,:)=getvar(cbasinmask,'tmaskant',1,npiglo,npjglo)
!    DO jj=1,npjglo
!       zonal_heat_aus(jj)= SUM( ztrp(:,jj)*zmask(:,jj))
!       zonal_salt_aus(jj)= SUM( ztrps(:,jj)*zmask(:,jj))
!       zpoints =  SUM(zmask(:,jj))
!    END DO

!   ! Med
     zonal_heat_med = 0.
     zonal_salt_med = 0.

!    zmask(:,:)=getvar(cbasinmask,'tmaskmed',1,npiglo,npjglo)
!    DO jj=1,npjglo
!       zonal_heat_med(jj)= SUM( ztrp(:,jj)*zmask(:,jj))
!       zonal_salt_med(jj)= SUM( ztrps(:,jj)*zmask(:,jj))
!       zpoints =  SUM(zmask(:,jj))
!    END DO
  ENDIF

  ! Output file
  OPEN(numout,FILE=cfileout)
  WRITE(numout,*)'! Zonal heat transport (integrated alon I-model coordinate) (in Pw)'
  IF ( llglo ) THEN
     WRITE(numout,*)'! J        Global          Atlantic         Pacific          Indian           Mediteranean     Austral  '
     DO jj=npjglo, 1, -1
        WRITE(numout,9000) jj, &
          dumlat(1,jj),  zonal_heat_glo(jj)/1e15 , &
            zonal_heat_atl(jj)/1e15, &
            zonal_heat_pac(jj)/1e15, &
            zonal_heat_ind(jj)/1e15, &
            zonal_heat_med(jj)/1e15, &
            zonal_heat_aus(jj)/1e15
     END DO
  ELSE
     WRITE(numout,*)'! J        Global        '
     DO jj=npjglo, 1, -1
        WRITE(numout,9000) jj, &
          dumlat(1,jj),  zonal_heat_glo(jj)/1e15  
     END DO
  ENDIF       
  !               
  CLOSE(numout)

  OPEN(numout,FILE=cfileouts)
  WRITE(numout,*)' ! Zonal salt transport (integrated alon I-model coordinate) (in 10^6 kg/s)'
  IF ( llglo ) THEN
     WRITE(numout,*)' ! J        Global          Atlantic         Pacific          Indian           Mediteranean     Austral  '
  !               
     DO jj=npjglo, 1, -1
        WRITE(numout,9001) jj, &
          dumlat(1,jj),  zonal_salt_glo(jj)/1e6 , &
           zonal_salt_atl(jj)/1e6, &
           zonal_salt_pac(jj)/1e6, &
           zonal_salt_ind(jj)/1e6, &
           zonal_salt_med(jj)/1e6, &
           zonal_salt_aus(jj)/1e6
     END DO
  ELSE
     WRITE(numout,*)' J        Global  '
     DO jj=npjglo, 1, -1
        WRITE(numout,9001) jj, &
          dumlat(1,jj),  zonal_salt_glo(jj)/1e6  
     ENDDO
   ENDIF

   CLOSE(numout)


9000 FORMAT(I4,6(1x,f9.3,1x,f8.4))
9001 FORMAT(I4,6(1x,f9.2,1x,f9.3))


END PROGRAM cdfmhst
