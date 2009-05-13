PROGRAM cdfhflx
  !!-------------------------------------------------------------------
  !!               ***  PROGRAM cdfhflx  ***
  !!
  !!  **  Purpose  :  Compute the Meridional Heat Transport from the forcing fluxes
  !!                  PARTIAL STEPS
  !!  
  !!  **  Method   :   Compute the zonaly integrated heat flux.
  !!                  The program looks for the file "new_maskglo.nc". If it does not exist, 
  !!                  only the calculation over all the domain is performed (this is adequate 
  !!                  for a basin configuration like NATL4).
  !!                  In new_maskglo.nc the masking corresponds to the global
  !!                  configuration. (Global, Atlantic, Indo-Pacific, Indian,Pacific ocean)
  !!
  !!
  !! history ;
  !!  Original :  J.M. Molines (jul. 2005) 
  !!              A.M. Treguier (april 2006) adaptation to NATL4 case 
  !!-------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !! * Modules used
  USE cdfio

  !! * Local variables
  IMPLICIT NONE
  INTEGER   :: jpbasins
  INTEGER   :: jbasin, jj, jk ,ji                  !: dummy loop index
  INTEGER   :: ierr                                !: working integer
  INTEGER   :: narg, iargc                         !: command line 
  INTEGER   :: npiglo,npjglo, npk                  !: size of the domain
  INTEGER   :: ncout, np, imean
  INTEGER   :: numout=10
  INTEGER, DIMENSION(2)          ::  iloc

  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  e1t, e2t, gphit, zflx !:  metrics, velocity
  REAL(KIND=4), DIMENSION (:,:,:),   ALLOCATABLE ::  zmask             !:  jpbasins x npiglo x npjglo
  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  dumlon              !: dummy longitude = 0.
  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  dumlat              !: latitude for i = north pole
  REAL(KIND=4) ,DIMENSION(1)                    ::  tim
  REAL(KIND=4) ,DIMENSION(:,:) , ALLOCATABLE ::  gphimean,htrp    !: jpbasins x npjglo

  REAL(KIND=8) ,DIMENSION(:,:) , ALLOCATABLE ::  zmht    !: jpbasins x npjglo 

  CHARACTER(LEN=80) :: cfilet , cfileout='hflx.out'
  CHARACTER(LEN=80) :: coordhgr='mesh_hgr.nc',cbasinmask='new_maskglo.nc'
  LOGICAL    :: llglo = .false.                          !: indicator for presence of new_maskglo.nc file 


  INTEGER    :: istatus

  ! constants

  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfhflx  T file '
     PRINT *,' Computes the MHT from heat fluxes '
     PRINT *,' Files mesh_hgr.nc, new_maskglo.nc must be in the current directory'
     PRINT *,' Output on hflx.out (ascii file )'
     STOP
  ENDIF

  CALL getarg (1, cfilet)
  npiglo= getdim (cfilet,'x')
  npjglo= getdim (cfilet,'y')
  npk   = getdim (cfilet,'depth')

  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk

 !  Detects newmaskglo file 
  INQUIRE( FILE='new_maskglo.nc', EXIST=llglo )
  IF (llglo) THEN
     jpbasins = 5
  ELSE
     jpbasins = 1
  ENDIF 

  ! Allocate arrays
  ALLOCATE ( zmask(jpbasins,npiglo,npjglo) )
  ALLOCATE ( zflx(npiglo,npjglo) )
  ALLOCATE ( e1t(npiglo,npjglo),e2t(npiglo,npjglo), gphit(npiglo,npjglo) )
  ALLOCATE ( htrp (jpbasins,npjglo) )
  ALLOCATE ( zmht(jpbasins, npjglo) )
  ALLOCATE ( dumlon(1,npjglo) , dumlat(1,npjglo))

  e1t(:,:) = getvar(coordhgr, 'e1t', 1,npiglo,npjglo) 
  e2t(:,:) = getvar(coordhgr, 'e2t', 1,npiglo,npjglo) 
  gphit(:,:) = getvar(coordhgr, 'gphit', 1,npiglo,npjglo)

  iloc=maxloc(gphit)
  dumlat(1,:) = gphit(iloc(1),:)
  dumlon(:,:) = 0.   ! set the dummy longitude to 0

  ! reading the masks
  ! 1 : global ; 2 : Atlantic ; 3 : Indo-Pacif ; 4 : Indian ; 5 : Pacif
  zmask(1,:,:)=getvar('mask.nc','vmask',1,npiglo,npjglo)

  IF (llglo) THEN
     zmask(2,:,:)=getvar(cbasinmask,'tmaskatl',1,npiglo,npjglo)
     zmask(4,:,:)=getvar(cbasinmask,'tmaskind',1,npiglo,npjglo)
     zmask(5,:,:)=getvar(cbasinmask,'tmaskpac',1,npiglo,npjglo)
     zmask(3,:,:)=zmask(5,:,:)+zmask(4,:,:)
     ! ensure that there are no overlapping on the masks
     WHERE(zmask(3,:,:) > 0 ) zmask(3,:,:) = 1
     !       change global mask for GLOBAL periodic condition
     zmask(1,1,:) = 0.
     zmask(1,npiglo,:) = 0.
  ENDIF
 
  ! initialize zmht
  zmht(:,:) = 0.
  htrp(:,:) = 0.
  

     ! Get fluxes
     zflx(:,:)= getvar(cfilet, 'sohefldo',  1 ,npiglo,npjglo)

     ! integrates 'zonally' (along i-coordinate)
     DO ji=1,npiglo
       ! For all basins 
       DO jbasin = 1, jpbasins
         DO jj=1,npjglo
            zmht(jbasin,jj)=zmht(jbasin,jj) + e1t(ji,jj)*e2t(ji,jj)* zmask(jbasin,ji,jj)*zflx(ji,jj)
         ENDDO
       END DO
     END DO

     ! cumulates transport from north to south
     DO jj=npjglo-1,1,-1
        DO jbasin=1, jpbasins
           htrp(jbasin,jj) = htrp(jbasin,jj+1) - zmht(jbasin,jj)
        END DO
     END DO

  OPEN(numout,FILE=cfileout,FORM='FORMATTED', RECL=256)  ! to avoid wrapped line with ifort
  WRITE(numout,*)'! Zonal heat transport (integrated from surface fluxes) (in Pw)'
  IF (llglo) THEN
     WRITE(numout,*)'! J        Global          Atlantic         INDO-PACIF    INDIAN  PACIF '
     DO jj=npjglo, 1, -1
        WRITE(numout,9000) jj, &
          dumlat(1,jj),  htrp(1,jj)/1e15 , &
            htrp(2,jj)/1e15, &
            htrp(3,jj)/1e15, &
            htrp(4,jj)/1e15, &
            htrp(5,jj)/1e15
     ENDDO
  ELSE
     WRITE(numout,*)'! J        Global   '
     DO jj=npjglo, 1, -1
        WRITE(numout,9000) jj, &
          dumlat(1,jj),  htrp(1,jj)/1e15  
     ENDDO
  ENDIF

  CLOSE(numout)
9000 FORMAT(I4,5(1x,f9.3,1x,f8.4))


   END PROGRAM cdfhflx
