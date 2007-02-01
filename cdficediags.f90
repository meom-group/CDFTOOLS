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
  !!-------------------------------------------------------------------
  !! * Modules used
  USE cdfio

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: jk
  INTEGER   :: ierr                                !: working integer
  INTEGER   :: narg, iargc                         !: command line 
  INTEGER   :: npiglo,npjglo                       !: size of the domain
  INTEGER   :: nvpk                                !: vertical levels in working variable
  INTEGER   :: nperio = 4                          !: boundary condition ( periodic, north fold)

  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::  e1, e2            !:  metrics
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::  zmask ,ff         !:   npiglo x npjglo
  REAL(KIND=4), DIMENSION (:,:),   ALLOCATABLE ::  ricethick, riceldfra !: thickness, leadfrac (concentration)

  REAL(KIND=8)      :: zvols,  zareas, zextends,zextends2        !: volume, area extend South hemisphere
  REAL(KIND=8)      :: zvoln,  zarean, zextendn,zextendn2        !: volume, area extend North hemisphere

  CHARACTER(LEN=80) :: cfilev , cdum
  CHARACTER(LEN=80) :: coordhgr='mesh_hgr.nc',  cmask='mask.nc'

  INTEGER    :: istatus

  ! constants

  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdficediag ncfile '
     PRINT *,' Files mesh_hgr.nc, mask.nc '
     PRINT *,'  must be in the current directory'
     PRINT *,' Output on standard output'
     STOP
  ENDIF

  CALL getarg (1, cfilev)

  npiglo= getdim (cfilev,'x')
  npjglo= getdim (cfilev,'y')

  ALLOCATE ( zmask(npiglo,npjglo) ,ff(npiglo,npjglo) )
  ALLOCATE ( ricethick(npiglo,npjglo) )
  ALLOCATE ( riceldfra(npiglo,npjglo) )
  ALLOCATE ( e1(npiglo,npjglo),e2(npiglo,npjglo) )

  e1(:,:) = getvar(coordhgr, 'e1t', 1,npiglo,npjglo)
  e2(:,:) = getvar(coordhgr, 'e2t', 1,npiglo,npjglo)
  ff(:,:) = getvar(coordhgr, 'ff' , 1,npiglo,npjglo)


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


   END PROGRAM cdficediag
