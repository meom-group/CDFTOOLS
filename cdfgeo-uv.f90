PROGRAM cdfgeo_uv
  !!---------------------------------------------------------------------------
  !!         ***  PROGRAM  cdfgeo_uv  ***
  !!
  !!  **  Purpose: Compute the ug and vg component of the geostrophic velocity
  !!               from the SSH field
  !! 
  !!  **  Method :  ug = -g/f * d(ssh)/dy
  !!                vg =  g/f * d(ssh)/dx
  !!  
  !!  **  Note : ug is located on a V grid point
  !!             vg                 U grid point 
  !!
  !!
  !! history :
  !!   Original :  J. Jouanno (Feb 2008)
  !     remark JMM : use of fmask ? use of ff ?  
  !!---------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !! * Modules used
  USE cdfio

  !! * Local variables
  IMPLICIT NONE
  INTEGER :: ji,jj                                !: dummy loop index
  INTEGER :: npiglo, npjglo , npk                 !: size of the domain
  INTEGER :: narg, iargc, ncoutu, ncoutv , ierr              !: 
  INTEGER, DIMENSION(1) ::  ipk, id_varoutu ,  id_varoutv         ! 

  REAL(kind=8), DIMENSION(:,:), ALLOCATABLE  :: e1u, e2v , ff !, e1t, e2t  !: metrics
  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE  :: glamt, gphit , glamv , gphiv , glamu , gphiu      !: longitude latitude
  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE  :: un, vn  
  REAL(KIND=4) , DIMENSION (:), ALLOCATABLE  :: dep 
  REAL(kind=4), DIMENSION(:,:), ALLOCATABLE  :: sshn , fmask 
  REAL(KIND=4) ,DIMENSION(1)                 :: tim
  REAL(KIND=4)                               :: g 
  CHARACTER(LEN=80) :: cfilt
  CHARACTER(LEN=80) :: coordhgr='mesh_hgr.nc',  coordzgr='mesh_zgr.nc'
  CHARACTER(LEN=80) :: cfiloutu='ugeo.nc' , cfileoutv='vgeo.nc'
  CHARACTER(LEN=80) :: cvart='sossheig', cvaru='vozocrtx', cvarv='vomecrty' 

  TYPE(variable), DIMENSION(1)      :: typvaru ,typvarv     !: structure for attributes

  g=9.81

  !!
  narg = iargc()
  IF ( narg < 1 ) THEN
     PRINT *,' USAGE : cdfgeo-uv fileT'
     PRINT *,'        Read sossheig on grid T'
     PRINT *,'        Produce 2 cdf file ugeo.nc and vgeo.nc with vozocrtx and vomecrty variables'
     PRINT *,'        Names of the variable have been chosen to be compatible with cdfeke, but note that Ugeo and Vgeo are now respectively on V and U grid points'
     PRINT *,'        Need mesh_hgr.nc mesh_zgr.nc'
     STOP
  ENDIF

  CALL getarg(1, cfilt)

  npiglo = getdim(cfilt,'x')
  npjglo = getdim(cfilt,'y')
  npk    = getdim(cfilt,'depth') 

  ipk(1)=1
  
  typvaru(1)%name=TRIM(cvaru)
  typvaru(1)%units='m/s'
  typvaru(1)%missing_value=0.
  typvaru(1)%valid_min= 0.
  typvaru(1)%valid_max= 20.
  typvaru(1)%long_name='Zonal_Geostrophic_Velocity'
  typvaru(1)%short_name=TRIM(cvaru)
  typvaru(1)%online_operation='N/A'
  typvaru(1)%axis='TYX'
  
  typvarv(1)%name=TRIM(cvarv)
  typvarv(1)%units='m/s'
  typvarv(1)%missing_value=0.
  typvarv(1)%valid_min= 0.
  typvarv(1)%valid_max= 20.
  typvarv(1)%long_name='Meridional_Geostrophic_Velocity'
  typvarv(1)%short_name=TRIM(cvarv)
  typvarv(1)%online_operation='N/A'
  typvarv(1)%axis='TYX'


  ! Allocate the memory
  ALLOCATE ( e1u(npiglo,npjglo) ,e2v(npiglo,npjglo) )
  ALLOCATE ( ff(npiglo,npjglo) )
  ALLOCATE ( glamt(npiglo,npjglo), gphit(npiglo,npjglo)  )
  ALLOCATE ( glamu(npiglo,npjglo), gphiu(npiglo,npjglo)  )
  ALLOCATE ( glamv(npiglo,npjglo), gphiv(npiglo,npjglo)  )
  ALLOCATE ( un(npiglo,npjglo)  , vn(npiglo,npjglo)  )
  ALLOCATE ( sshn(npiglo,npjglo) , fmask(npiglo,npjglo) )

  ! Read the metrics from the mesh_hgr file
  e2v=  getvar(coordhgr, 'e2v', 1,npiglo,npjglo)
  e1u=  getvar(coordhgr, 'e1u', 1,npiglo,npjglo)
  ff(:,:)  = getvar(coordhgr, 'ff', 1,npiglo,npjglo) 
 
  ! and the coordinates   from the mesh_hgr file
  glamt = getvar(coordhgr, 'glamt', 1,npiglo,npjglo)
  gphit = getvar(coordhgr, 'gphit', 1,npiglo,npjglo)
 
  ! create output fileset
  glamu=getvar(coordhgr,'glamu',1,npiglo,npjglo)
  gphiu=getvar(coordhgr,'gphiu',1,npiglo,npjglo)
  dep=getvare3(coordzgr,'gdept',1)

  ncoutu =create(cfiloutu,cfilt,npiglo,npjglo,0)
  ierr= createvar(ncoutu,typvaru(1),1,ipk,id_varoutu )
  ierr= putheadervar(ncoutu,cfilt,npiglo,npjglo,0,pnavlon=glamv,pnavlat=gphiv)

  tim=getvar1d(cfilt,'time_counter',1)
  ierr=putvar1d(ncoutu,tim,1,'T')

  glamv=getvar(coordhgr,'glamv',1,npiglo,npjglo)
  gphiv=getvar(coordhgr,'gphiv',1,npiglo,npjglo)
  dep=getvare3(coordzgr,'gdept',1)

  ncoutv =create(cfileoutv,cfilt,npiglo,npjglo,0)
  ierr= createvar(ncoutv,typvarv(1),1,ipk,id_varoutv )
  ierr= putheadervar(ncoutv,cfilt,npiglo,npjglo,0,pnavlon=glamu,pnavlat=gphiu)

  tim=getvar1d(cfilt,'time_counter',1)
  ierr=putvar1d(ncoutv,tim,1,'T')

  ! Read ssh
  sshn(:,:) =  getvar(cfilt, cvart,1, npiglo,npjglo)

  ! compute the masks
  fmask=0.
  DO jj = 1, npjglo - 1
     DO ji = 1, npiglo - 1
        fmask(ji,jj)=0.
        fmask(ji,jj)= sshn(ji,jj)*sshn(ji,jj+1)*sshn(ji+1,jj)
        IF (fmask(ji,jj) /= 0.) fmask(ji,jj)=1.
     ENDDO
  ENDDO

  ! Calculation of geostrophic velocity :
  un(:,:) = 0.
  vn(:,:) = 0.

  DO jj = 2, npjglo - 1
    DO ji = 2, npiglo -1 
               vn(ji,jj) =   g * fmask(ji,jj) * ( sshn(ji+1,jj  ) -sshn(ji,jj) )  / ( ff(ji,jj) * e1u(ji,jj) ) 
               un(ji,jj) = - g * fmask(ji,jj) * ( sshn(ji  ,jj+1) -sshn(ji,jj) )  / ( ff(ji,jj) * e2v(ji,jj) )
    END DO
  END DO
     
  ! write un and vn  ...
  ierr = putvar(ncoutu,id_varoutu(1),un(:,:),1,npiglo,npjglo)
  ierr = putvar(ncoutv,id_varoutv(1),vn(:,:),1,npiglo,npjglo)

  ierr = closeout(ncoutu)
  ierr = closeout(ncoutv)

END PROGRAM cdfgeo_uv

