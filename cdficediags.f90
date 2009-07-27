PROGRAM cdficediag
  !!-------------------------------------------------------------------
  !!               ***  PROGRAM cdficediag  ***
  !!
  !!  **  Purpose  :  Compute the Ice volume, area and extend for each hemisphere
  !!  
  !!  **  Method   :   Read the ice output and integrates (2D)
  !!                   determine the hemisphere by the sign of ff (coriolis)
  !!
  !! history ;
  !!  Original :  J.M. Molines (Jan. 2006)
  !!              R. Dussin (Jul. 2009) : Add netcdf output
  !!-------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !! * Modules used
  USE cdfio

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: jk, jj
  INTEGER   :: ierr                                !: working integer
  INTEGER   :: narg, iargc                         !: command line 
  INTEGER   :: npiglo,npjglo                       !: size of the domain
  INTEGER   :: nvpk                                !: vertical levels in working variable
  INTEGER   :: nperio = 4                          !: boundary condition ( periodic, north fold)
  ! added to write in netcdf
  INTEGER :: kx=1, ky=1, kz=1          ! dims of netcdf output file
  INTEGER :: nboutput=8                ! number of values to write in cdf output
  INTEGER :: ncout                     ! for netcdf output
  INTEGER, DIMENSION(:), ALLOCATABLE ::  ipk, id_varout

  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::  e1, e2            !:  metrics
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::  zmask ,ff         !:   npiglo x npjglo
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::  ricethick, riceldfra !: thickness, leadfrac (concentration)

  REAL(KIND=8)      :: zvols,  zareas, zextends,zextends2        !: volume, area extend South hemisphere
  REAL(KIND=8)      :: zvoln,  zarean, zextendn,zextendn2        !: volume, area extend North hemisphere
  ! added to write in netcdf
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE ::  dumlon, dumlat
  REAL(KIND=4), DIMENSION (1)               ::  tim ! time counter
  TYPE(variable), DIMENSION(:), ALLOCATABLE :: typvar  ! structure of output
  !  REAL(KIND=4), DIMENSION(:), ALLOCATABLE   ::   gdepw    ! depth read 
  !
  CHARACTER(LEN=256) :: cfilev , cdum
  CHARACTER(LEN=256) :: coordhgr='mesh_hgr.nc',  cmask='mask.nc'
  ! added to write in netcdf
  CHARACTER(LEN=256) :: cfileoutnc='icediags.nc' , cflagcdf
  ! added to write in netcdf
  LOGICAL :: lwrtcdf=.FALSE.

  INTEGER    :: istatus

  ! constants

  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg == 0 .OR. narg >= 3 ) THEN
     PRINT *,' Usage : cdficediag ncfile [cdfout]'
     PRINT *,' Files mesh_hgr.nc, mask.nc '
     PRINT *,'  must be in the current directory'
     PRINT *,' Output on standard output'
     PRINT *,' Optional Output in NetCDF with cdfout option'
     STOP
  ENDIF

  CALL getarg (1, cfilev)

  npiglo= getdim (cfilev,'x')
  npjglo= getdim (cfilev,'y')

  IF ( narg == 2 ) THEN
     CALL getarg (2, cflagcdf)
     IF (cflagcdf=='cdfout') THEN
        lwrtcdf=.TRUE.
     ELSE
        PRINT *, 'unknown option'
     ENDIF
  ENDIF

  ALLOCATE ( zmask(npiglo,npjglo) ,ff(npiglo,npjglo) )
  ALLOCATE ( ricethick(npiglo,npjglo) )
  ALLOCATE ( riceldfra(npiglo,npjglo) )
  ALLOCATE ( e1(npiglo,npjglo),e2(npiglo,npjglo) )

  IF(lwrtcdf) THEN

     ALLOCATE ( typvar(nboutput), ipk(nboutput), id_varout(nboutput) )
     ALLOCATE (dumlon(1,1) , dumlat(1,1) )

     dumlon(:,:)=0.
     dumlat(:,:)=0.

     DO jj=1,nboutput
        ipk(jj)=1
     ENDDO

     ! define new variables for output 
     typvar(1)%name='NVolume'
     typvar(1)%units='10^9 m3'
     typvar%scale_factor= 1.
     typvar%add_offset= 0.
     typvar%savelog10= 0.
     typvar(1)%long_name='Ice_volume_in_Northern_Hemisphere'
     typvar(1)%short_name='NVolume'
     typvar%online_operation='N/A'
     typvar%axis='T'

     typvar(2)%name='NArea'
     typvar(2)%units='10^9 m2'
     typvar(2)%long_name='Ice_area_in_Northern_Hemisphere'
     typvar(2)%short_name='NArea'

     typvar(3)%name='NExtent'
     typvar(3)%units='10^9 m2'
     typvar(3)%long_name='Ice_extent_in_Northern_Hemisphere'
     typvar(3)%short_name='NExtent'

     typvar(4)%name='NExnsidc'
     typvar(4)%units='10^9 m2'
     typvar(4)%long_name='Ice_extent_similar_to_NSIDC_in_Northern_Hemisphere'
     typvar(4)%short_name='NExnsidc'

     typvar(5)%name='SVolume'
     typvar(5)%units='10^9 m3'
     typvar(5)%long_name='Ice_volume_in_Southern_Hemisphere'
     typvar(5)%short_name='SVolume'

     typvar(6)%name='SArea'
     typvar(6)%units='10^9 m2'
     typvar(6)%long_name='Ice_area_in_Southern_Hemisphere'
     typvar(6)%short_name='SArea'

     typvar(7)%name='SExtent'
     typvar(7)%units='10^9 m2'
     typvar(7)%long_name='Ice_extent_in_Southern_Hemisphere'
     typvar(7)%short_name=''

     typvar(8)%name='SExnsidc'
     typvar(8)%units='10^9 m2'
     typvar(8)%long_name='Ice_extent_similar_to_NSIDC_in_Southern_Hemisphere'
     typvar(8)%short_name='SExnsidc'


  ENDIF

  e1(:,:) = getvar(coordhgr, 'e1t', 1,npiglo,npjglo)
  e2(:,:) = getvar(coordhgr, 'e2t', 1,npiglo,npjglo)
  ! only the sign of ff is important
  ff(:,:) = getvar(coordhgr, 'gphit' , 1,npiglo,npjglo)


  ricethick(:,:)= getvar(cfilev, 'iicethic',  1 ,npiglo,npjglo)
  riceldfra(:,:)= getvar(cfilev, 'ileadfra',  1 ,npiglo,npjglo)

  ! modify the mask for periodic and north fold condition (T pivot, F Pivot ...)
  ! in fact should be nice to use jperio as in the code ...

  zmask(:,:)=getvar(cmask,'tmask',1,npiglo,npjglo)
  SELECT CASE (nperio)
  CASE (0) ! closed boundaries
     ! nothing to do
  CASE (4) ! ORCA025 type boundary
     zmask(1:2,:)=0.
     zmask(:,npjglo)=0.
     zmask(npiglo/2+1:npiglo,npjglo-1)= 0.
  CASE (6)
     zmask(1:2,:)=0.
     zmask(:,npjglo)=0.
  CASE DEFAULT
     PRINT *,' Nperio=', nperio,' not yet coded'
     STOP
  END SELECT

  ! North : ff > 0 
  zvoln=SUM( ricethick (:,:)* e1(:,:) * e2(:,:) * riceldfra (:,:) * zmask (:,:) , (ff > 0 ) )
  zarean=SUM( e1(:,:) * e2(:,:) * riceldfra (:,:) * zmask (:,:) ,( ff > 0 ) )
  zextendn=SUM( e1(:,:) * e2(:,:) * riceldfra (:,:) * zmask (:,:), (riceldfra > 0.15 .AND. ff > 0 ) )
  ! JMM added 22/01/2007 : to compute same extent than the NSIDC
  zextendn2=SUM( e1(:,:) * e2(:,:) * zmask (:,:), (riceldfra > 0.15 .AND. ff > 0 ) )

  ! South : ff < 0
  zvols=SUM( ricethick (:,:)* e1(:,:) * e2(:,:) * riceldfra (:,:) * zmask (:,:) ,(ff < 0 ) )
  zareas=SUM( e1(:,:) * e2(:,:) * riceldfra (:,:) * zmask (:,:), ( ff < 0 ) )
  zextends=SUM( e1(:,:) * e2(:,:) * riceldfra (:,:) * zmask (:,:), (riceldfra > 0.15 .AND. ff < 0  ) )
  zextends2=SUM( e1(:,:) * e2(:,:)* zmask (:,:), (riceldfra > 0.15 .AND. ff < 0  ) )

  PRINT *,' Northern Hemisphere ' 
  PRINT *,'          NVolume (10^9 m3)  ', zvoln /1.d9
  PRINT *,'          NArea (10^9 m2)    ', zarean /1.d9
  PRINT *,'          NExtend (10^9 m2)  ', zextendn /1.d9
  PRINT *,'          NExnsidc (10^9 m2)  ', zextendn2 /1.d9
  PRINT *
  PRINT *,' Southern Hemisphere ' 
  PRINT *,'          SVolume (10^9 m3)  ', zvols /1.d9
  PRINT *,'          SArea (10^9 m2)    ', zareas /1.d9
  PRINT *,'          SExtend (10^9 m2)  ', zextends /1.d9
  PRINT *,'          SExnsidc (10^9 m2)  ', zextends2 /1.d9

  IF (lwrtcdf) THEN

     ! create output fileset
     ncout =create(cfileoutnc,'none',kx,ky,kz,cdep='depthw')
     ierr= createvar(ncout,typvar,nboutput,ipk,id_varout )
     ierr= putheadervar(ncout, cfilev,kx, &
          ky,kz,pnavlon=dumlon,pnavlat=dumlat)
     tim=getvar1d(cfilev,'time_counter',1)
     ierr=putvar1d(ncout,tim,1,'T')

     ! netcdf output 
     ierr = putvar0d(ncout,id_varout(1), REAL(zvoln /1.d9) )
     ierr = putvar0d(ncout,id_varout(2), REAL(zarean /1.d9) )
     ierr = putvar0d(ncout,id_varout(3), REAL(zextendn /1.d9) )
     ierr = putvar0d(ncout,id_varout(4), REAL(zextendn2 /1.d9) )
     ierr = putvar0d(ncout,id_varout(5), REAL(zvols /1.d9) )
     ierr = putvar0d(ncout,id_varout(6), REAL(zareas /1.d9) )
     ierr = putvar0d(ncout,id_varout(7), REAL(zextends /1.d9) )
     ierr = putvar0d(ncout,id_varout(8), REAL(zextends2 /1.d9) )

     ierr = closeout(ncout)

  ENDIF


END PROGRAM cdficediag
