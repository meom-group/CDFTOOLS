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
  INTEGER, DIMENSION(2)          ::  iloc
  LOGICAL    :: llglo = .false.                !: indicator for presence of new_maskglo.nc file

  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE ::  zmask, e1v, gphiv, zvt, zvs !: mask, metrics
  REAL(KIND=4), DIMENSION (:,:), ALLOCATABLE ::  dumlon              !: dummy longitude = 0.
  REAL(KIND=4), DIMENSION (:,:), ALLOCATABLE ::  dumlat              !: latitude for i = north pole

  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE ::  e3t
  REAL(KIND=4)                              ::  zpoints

  REAL(KIND=8), DIMENSION (:,:), ALLOCATABLE :: zwk , ztrp, ztrps, zwks
  REAL(KIND=8) ,DIMENSION(:) , ALLOCATABLE ::  zonal_heat_glo, zonal_heat_atl, zonal_heat_pac,&
       &                                       zonal_heat_ind, zonal_heat_aus, zonal_heat_med
  REAL(KIND=8) ,DIMENSION(:) , ALLOCATABLE ::  zonal_salt_glo, zonal_salt_atl, zonal_salt_pac,&
       &                                       zonal_salt_ind, zonal_salt_aus, zonal_salt_med, zmtrp

  CHARACTER(LEN=80) :: cfilet ,cfileout='zonal_heat_trp.dat', cfileouts='zonal_salt_trp.dat'
  ! to be put in namelist eventually
  CHARACTER(LEN=80) :: coordhgr='mesh_hgr.nc',  coordzgr='mesh_zgr.nc', cbasinmask='new_maskglo.nc'

 ! NC output
  INTEGER            :: npvar=1
  INTEGER            :: jbasins, js, jvar   !: dummy loop index
  INTEGER            :: ncout,  nbasins, ierr
  INTEGER, DIMENSION(:), ALLOCATABLE :: ipk, id_varout

  REAL(KIND=4), PARAMETER :: rpspval=9999.99
  REAL(KIND=4), DIMENSION(1) :: gdep
  REAL(KIND=4), DIMENSION (1)                    ::  tim

  CHARACTER(LEN=80) :: cfileoutnc='mhst.nc', cdum
  CHARACTER(LEN=80), DIMENSION(:), ALLOCATABLE   :: cvarname             !: array of var name for input
  CHARACTER(LEN=4),DIMENSION(5) :: cbasin=(/'_glo','_atl','_inp','_ind','_pac'/)
  TYPE(variable), DIMENSION(:), ALLOCATABLE   :: typvar                  !: structure for attributes


  ! constants
  REAL(KIND=4),PARAMETER   ::  rau0=1000.,   rcp=4000.

  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfmhst_full  VTfile [MST]'
     PRINT *,' Files mesh_hgr.nc, mesh_zgr.nc ,mask.nc, new_maskglo.nc must be in te current directory'
     PRINT *,' Output on zonal_heat_trp.dat and zonal_salt_trp.dat'
     STOP
  ENDIF

  CALL getarg (1, cfilet)
  npiglo= getdim (cfilet,'x')
  npjglo= getdim (cfilet,'y')
  npk   = getdim (cfilet,'depth')
  npvar=1
  IF ( narg == 2 ) THEN
   CALL getarg(2,cdum)
   IF ( cdum /= 'MST' ) THEN
      PRINT *,' unknown option :', TRIM(cdum)  ; STOP
   ENDIF
   npvar=2
  ENDIF

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
  ALLOCATE ( zmtrp(npjglo) )
  ALLOCATE ( dumlon(1,npjglo) , dumlat(1,npjglo))

  ! Read metrics and latitudes
  e1v(:,:)   = getvar(coordhgr, 'e1v', 1,npiglo,npjglo)
  gphiv(:,:) = getvar(coordhgr, 'gphiv', 1,npiglo,npjglo)
  e3t(:)     = getvare3(coordzgr,'e3t',npk) 
  gdep(:) = getvare3(coordzgr, 'depthv' ,1)

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
     zwk(:,:)  = zvt(:,:)*e1v(:,:)*e3t(jk)
     zwks(:,:) = zvs(:,:)*e1v(:,:)*e3t(jk)

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
 
  nbasins=1
  IF ( llglo) THEN ! 5 basins
       nbasins=5
  ENDIF
  
  ! Allocate output variables
  ALLOCATE(typvar(nbasins*npvar),cvarname(nbasins*npvar))
  ALLOCATE(ipk(nbasins*npvar),id_varout(nbasins*npvar))
  ipk(:)=1               ! all output variables have only 1 level !
  DO jbasins = 1,nbasins
   SELECT CASE ( npvar )
   CASE ( 1 )   ! only MHT is output
   cvarname(jbasins) = 'zomht'//TRIM(cbasin(jbasins))
   typvar(jbasins)%name=cvarname(jbasins)
   typvar(jbasins)%units='PW'
   typvar(jbasins)%missing_value=rpspval
   typvar(jbasins)%valid_min=-10.
   typvar(jbasins)%valid_max=20
   typvar(jbasins)%long_name='Meridional Heat Transport '//TRIM(cbasin(jbasins))
   typvar(jbasins)%short_name=cvarname(jbasins)
   typvar(jbasins)%online_operation='N/A'
   typvar(jbasins)%axis='TY'
   CASE ( 2 )   ! both MHT and MST (meridional Salt Transport )
   cvarname(jbasins) = 'zomht'//TRIM(cbasin(jbasins))
   typvar(jbasins)%name=cvarname(jbasins)
   typvar(jbasins)%units='PW'
   typvar(jbasins)%missing_value=rpspval
   typvar(jbasins)%valid_min=-10.
   typvar(jbasins)%valid_max=20
   typvar(jbasins)%long_name='Meridional Heat Transport '//TRIM(cbasin(jbasins))
   typvar(jbasins)%short_name=cvarname(jbasins)
   typvar(jbasins)%online_operation='N/A'
   typvar(jbasins)%axis='TY'
   ! MST
   cvarname(nbasins+jbasins) = 'zomst'//TRIM(cbasin(jbasins))
   typvar(nbasins+jbasins)%name=cvarname(nbasins+jbasins)
   typvar(nbasins+jbasins)%units='T/sec'
   typvar(nbasins+jbasins)%missing_value=rpspval
   typvar(nbasins+jbasins)%valid_min=-10.e9
   typvar(nbasins+jbasins)%valid_max=20.e9
   typvar(nbasins+jbasins)%long_name='Meridional Salt Transport '//TRIM(cbasin(jbasins))
   typvar(nbasins+jbasins)%short_name=cvarname(nbasins+jbasins)
   typvar(nbasins+jbasins)%online_operation='N/A'
   typvar(nbasins+jbasins)%axis='TY'
   CASE DEFAULT
      PRINT * ,'   This program is not ready for npvar > 2 ' ; STOP
   END SELECT 
  END DO

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
  ! create output fileset
  ncout = create(cfileoutnc, cfilet, 1,npjglo,1,cdep='depthv')
  ierr  = createvar(ncout ,typvar,nbasins*npvar, ipk,id_varout )
  ierr  = putheadervar(ncout, cfilet,1,npjglo,1,pnavlon=dumlon,pnavlat=dumlat,pdep=gdep)
  tim   = getvar1d(cfilet,'time_counter',1)
  ierr  = putvar1d(ncout,tim,1,'T')

  DO jvar=1,npvar   !  MHT [ and MST ]  (1 or 2 )
    IF ( jvar == 1 ) THEN
       ! MHT
       js=1
       zmtrp(:)=zonal_heat_glo(:)/1.e15                        ! GLO
       WHERE ( zmtrp == 0 ) zmtrp=rpspval
       ierr=putvar(ncout,id_varout(js),REAL(zmtrp), 1,1,npjglo)
       js=js+1
       IF ( nbasins == 5 ) THEN
         zmtrp(:)=zonal_heat_atl(:)/1.e15                      ! ATL
         WHERE ( zmtrp == 0 ) zmtrp=rpspval         
         ierr=putvar(ncout,id_varout(js),REAL(zmtrp), 1,1,npjglo)
         js=js+1
         zmtrp(:)=zonal_heat_ind(:) + zonal_heat_pac(:)/1.e15  ! INP
         WHERE ( zmtrp == 0 ) zmtrp=rpspval
         ierr=putvar(ncout,id_varout(js),REAL(zmtrp), 1,1,npjglo)
         js=js+1
         zmtrp(:)=zonal_heat_ind(:)/1.e15                      ! IND
         WHERE ( zmtrp == 0 ) zmtrp=rpspval
         ierr=putvar(ncout,id_varout(js),REAL(zmtrp), 1,1,npjglo)
         js=js+1
         zmtrp(:)=zonal_heat_pac(:)/1.e15                      ! PAC
         WHERE ( zmtrp == 0 ) zmtrp=rpspval
         ierr=putvar(ncout,id_varout(js),REAL(zmtrp), 1,1,npjglo)
         js=js+1
       ENDIF
    ELSE
       ! MST
       zmtrp(:)=zonal_salt_glo(:)/1.e6                        ! GLO
       WHERE ( zmtrp == 0 ) zmtrp=rpspval
       ierr=putvar(ncout,id_varout(js),REAL(zmtrp), 1,1,npjglo)
       js = js + 1
       IF ( nbasins == 5 ) THEN
         zmtrp(:)=zonal_salt_atl(:)/1.e6                      ! ATL
         WHERE ( zmtrp == 0 ) zmtrp=rpspval
         ierr=putvar(ncout,id_varout(js),REAL(zmtrp), 1,1,npjglo)
         js = js + 1
         zmtrp(:)=zonal_salt_ind(:) + zonal_salt_pac(:)/1.e6  ! INP
         WHERE ( zmtrp == 0 ) zmtrp=rpspval
         ierr=putvar(ncout,id_varout(js),REAL(zmtrp), 1,1,npjglo)
         js = js + 1
         zmtrp(:)=zonal_salt_ind(:)/1.e6                      ! IND
         WHERE ( zmtrp == 0 ) zmtrp=rpspval
         ierr=putvar(ncout,id_varout(js),REAL(zmtrp), 1,1,npjglo)
         js = js + 1
         zmtrp(:)=zonal_salt_pac(:)/1.e6                      ! PAC
         WHERE ( zmtrp == 0 ) zmtrp=rpspval
         ierr=putvar(ncout,id_varout(js),REAL(zmtrp), 1,1,npjglo)
       js = js + 1
       ENDIF
    ENDIF
  END DO 
  ierr=closeout(ncout)

  OPEN(numout,FILE=cfileout)
  WRITE(numout,*)'% FULL STEP COMPUTATION'
  WRITE(numout,*)'% Zonal heat transport (integrated alon I-model coordinate) (in Pw)'
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
  WRITE(numout,*)'% FULL STEP COMPUTATION'
  WRITE(numout,*)'% Zonal salt transport (integrated alon I-model coordinate) (in 10^6 kg/s)'
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

END PROGRAM cdfmhst_full
