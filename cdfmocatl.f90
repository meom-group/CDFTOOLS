PROGRAM cdfmocatl
  !!-------------------------------------------------------------------
  !!               ***  PROGRAM cdfmocatl  ***
  !!
  !!  **  Purpose  :  Compute the Meridional Overturning Cell (MOC)
  !!                  PARTIAL STEPS
  !!  
  !!  **  Method   :  The MOC is computed from the V velocity field, integrated
  !!                  from the bottom to the surface, then zonally averaged with
  !!                  eventual masking for oceanic basins.
  !!                  Results are saved on moc.nc file with variable name  zomsfatl
  !!                  This version is intended to be used with Atlantic-only configurations
  !!
  !!
  !! history ;
  !!  Original :  J.M. Molines (jul. 2005) 
  !!              Atlantic Only :  J.M. Molines (Nov. 2005) 
  !!-------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !! * Modules used
  USE cdfio

  !! * Local variables
  IMPLICIT NONE
  INTEGER, PARAMETER :: jpbasins=1                 !: atlantic only !!
  INTEGER   :: jbasin, jj, jk ,ji                  !: dummy loop index
  INTEGER   :: ierr                                !: working integer
  INTEGER   :: narg, iargc                         !: command line 
  INTEGER   :: npiglo,npjglo, npk                  !: size of the domain
  INTEGER   :: ncout, np
  INTEGER   :: numout=10
  INTEGER, DIMENSION(jpbasins) ::  ipk, id_varout         !
  INTEGER, DIMENSION(2)          ::  iloc

  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  e1v, e3v, gphiv, zv !:  metrics, velocity
  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  dumlon              !: dummy longitude = 0.
  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  dumlat              !: latitude for i = north pole
  REAL(KIND=4), DIMENSION (:),       ALLOCATABLE ::  gdepw               !: deptw
  REAL(KIND=4), DIMENSION (:,:,:),   ALLOCATABLE ::  zmask               !:  jpbasins x npiglo x npjglo
  REAL(KIND=4), DIMENSION (1)                    ::  tim

  REAL(KIND=8) ,DIMENSION(:,:,:) , ALLOCATABLE ::  zomsf                 !: jpbasins x npjglo x npk

  CHARACTER(LEN=256) :: cfilev , cfileoutnc='moc.nc'
  CHARACTER(LEN=256) :: coordhgr='mesh_hgr.nc',  coordzgr='mesh_zgr.nc'
  TYPE(variable), DIMENSION(jpbasins)      :: typvar                   !: structure for attribures


  INTEGER    :: istatus

  ! constants

  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfmoc  V file '
     PRINT *,' Computes the MOC  for a mono basin oceanic configuration'
     PRINT *,' PARTIAL CELLS VERSION'
     PRINT *,' Files mesh_hgr.nc, mesh_zgr.nc ,mask.nc '
     PRINT *,'  must be in the current directory'
     PRINT *,' Output on moc.nc: '
     PRINT *,'      variables zomsfatl  : Atlantic Ocean '
     STOP
  ENDIF

  CALL getarg (1, cfilev)
  npiglo= getdim (cfilev,'x')
  npjglo= getdim (cfilev,'y')
  npk   = getdim (cfilev,'depth')

 ! define new variables for output ( must update att.txt)
  typvar(1)%name= 'zomsfatl'
  typvar(1)%units='Sverdrup'
  typvar(1)%missing_value=99999.
  typvar(1)%valid_min= -1000.
  typvar(1)%valid_max= 1000.
  typvar(1)%long_name='Meridional_Overt.Cell_Atlantic'
  typvar(1)%short_name='zomsfatl'
  typvar(1)%online_operation='N/A'
  typvar(1)%axis='TZY'


  ipk(1) = npk  !  

  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk

  ! Allocate arrays
  ALLOCATE ( zmask(jpbasins,npiglo,npjglo) )
  ALLOCATE ( zv(npiglo,npjglo) )
  ALLOCATE ( e1v(npiglo,npjglo),e3v(npiglo,npjglo), gphiv(npiglo,npjglo) ,gdepw(npk) )
  ALLOCATE ( zomsf(jpbasins, npjglo, npk) )
  ALLOCATE ( dumlon(1,npjglo) , dumlat(1,npjglo))


  e1v(:,:) = getvar(coordhgr, 'e1v', 1,npiglo,npjglo) 
  gphiv(:,:) = getvar(coordhgr, 'gphiv', 1,npiglo,npjglo)
  gdepw(:) = getvare3(coordzgr, 'gdepw',npk)
  gdepw(:) = -1.*  gdepw(:)

  iloc=maxloc(gphiv)
  dumlat(1,:) = gphiv(iloc(1),:)
  dumlon(:,:) = 0.   ! set the dummy longitude to 0

  ! create output fileset
   ncout =create(cfileoutnc, cfilev, 1,npjglo,npk,cdep='depthw')
   ierr= createvar(ncout ,typvar ,jpbasins, ipk,id_varout )
   ierr= putheadervar(ncout, cfilev,1, npjglo,npk,pnavlon=dumlon,pnavlat=dumlat,pdep=gdepw)
   tim=getvar1d(cfilev,'time_counter',1)
   ierr=putvar1d(ncout,tim,1,'T')


  ! reading the masks
  zmask(1,:,:)=getvar('mask.nc','vmask',1,npiglo,npjglo)

  ! initialize moc to 0
  zomsf(:,:,:) = 0.
  
  DO jk = 1,npk-1
     PRINT *,'level ',jk
     ! Get velocities v at jk
     zv(:,:)= getvar(cfilev, 'vomecrty',  jk ,npiglo,npjglo)

     ! get e3v at level jk ( ps...)
     e3v(:,:) = getvar(coordzgr, 'e3v_ps', jk,npiglo,npjglo, ldiom=.true.)
     
     ! integrates 'zonally' (along i-coordinate)
     DO ji=1,npiglo
       ! For all basins 
       DO jbasin = 1, jpbasins
         DO jj=1,npjglo
            zomsf(jbasin,jj,jk)=zomsf(jbasin,jj,jk) - e1v(ji,jj)*e3v(ji,jj)* zmask(jbasin,ji,jj)*zv(ji,jj)
         ENDDO
       END DO
     END DO
  END DO

  ! integrates vertically   from bottom to surface
  DO jk=npk-1 , 1 , -1
     zomsf(:,:,jk) = zomsf(:,:,jk+1) + zomsf(:,:,jk)/1.e6
  END DO  ! loop to next level

  ! netcdf output 
  DO jbasin= 1, jpbasins
    DO jk =1, npk
    ierr = putvar (ncout, id_varout(jbasin),REAL(zomsf(jbasin,:,jk)), jk,1,npjglo)
    END DO
  END DO

  ierr = closeout(ncout)
 
   END PROGRAM cdfmocatl
