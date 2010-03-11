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
  INTEGER   :: jpbasins
  INTEGER   :: jbasin, jj, jk ,ji                  !: dummy loop index
  INTEGER   :: ierr                                !: working integer
  INTEGER   :: narg, iargc                         !: command line 
  INTEGER   :: npiglo,npjglo, npk                  !: size of the domain
  INTEGER   :: ncout, np, imean
  INTEGER   :: numout=10
  INTEGER, DIMENSION(2)          ::  iloc
  ! added to write in netcdf
  INTEGER :: kx=1, ky=1                ! dims of netcdf output file
  INTEGER :: nboutput                  ! number of values to write in cdf output
  INTEGER, DIMENSION(:), ALLOCATABLE ::  ipk, id_varout

  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  e1t, e2t, gphit, zflx !:  metrics, velocity
  REAL(KIND=4), DIMENSION (:,:,:),   ALLOCATABLE ::  zmask             !:  jpbasins x npiglo x npjglo
  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  dumlon              !: dummy longitude = 0.
  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  dumlat              !: latitude for i = north pole
  REAL(KIND=4) ,DIMENSION(:,:) , ALLOCATABLE ::  gphimean,htrp    !: jpbasins x npjglo

  REAL(KIND=8) ,DIMENSION(:,:) , ALLOCATABLE ::  zmht    !: jpbasins x npjglo 
  ! added to write in netcdf
  REAL(KIND=4) :: threedmeanout, pmissing_value
  REAL(KIND=4), DIMENSION (1)               ::  tim ! time counter
  REAL(KIND=4), DIMENSION(:), ALLOCATABLE :: meanout

  TYPE(variable), DIMENSION(:), ALLOCATABLE :: typvar  ! structure of output

  CHARACTER(LEN=256) :: cfilet , cfileout='hflx.out'
  CHARACTER(LEN=256) :: coordhgr='mesh_hgr.nc',cbasinmask='new_maskglo.nc'
  ! added to write in netcdf
  CHARACTER(LEN=256) :: cfileoutnc='cdfhflx.nc' , cflagcdf
  CHARACTER(LEN=256) :: cdunits, cdlong_name, cdshort_name

  LOGICAL    :: llglo = .FALSE.                          !: indicator for presence of new_maskglo.nc file 
  ! added to write in netcdf
  LOGICAL :: lwrtcdf=.FALSE.

  INTEGER    :: istatus

  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfhflx  T file [cdfout]'
     PRINT *,' Computes the MHT from heat fluxes '
     PRINT *,' Files mesh_hgr.nc, new_maskglo.nc must be in the current directory'
     PRINT *,' Output on hflx.out (ascii file )'
     STOP
  ENDIF

  IF ( narg >= 1 ) THEN
     CALL getarg (1, cfilet)
     npiglo= getdim (cfilet,'x')
     npjglo= getdim (cfilet,'y')
     npk   = getdim (cfilet,'depth')

     PRINT *, 'npiglo=', npiglo
     PRINT *, 'npjglo=', npjglo
     PRINT *, 'npk   =', npk
  ENDIF

  IF ( narg == 2 ) THEN
     CALL getarg (2, cflagcdf)
     IF (cflagcdf == 'cdfout') THEN
        lwrtcdf=.TRUE.
     ELSE
        PRINT *,'Uncorrect second argument'
        PRINT *,'second argument must be "cdfout" to write in NetCDF'
     ENDIF
  ENDIF

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

  IF (lwrtcdf) THEN
     nboutput=jpbasins 
     ALLOCATE (typvar(nboutput), ipk(nboutput), id_varout(nboutput))

     DO jj=1,jpbasins
        ipk(jj)=1
     ENDDO

     ! define new variables for output 
     typvar(1)%name='hflx_glo'
     typvar%units=TRIM(cdunits)
     typvar%missing_value=99999.
     typvar%valid_min= -1000.
     typvar%valid_max= 1000.
     typvar%scale_factor= 1.
     typvar%add_offset= 0.
     typvar%savelog10= 0.
     typvar(1)%long_name='Heat_Fluxes_Global'
     typvar(1)%short_name='hflx_glo'
     typvar%online_operation='N/A'
     typvar%axis='ZT'

     IF (llglo) THEN

        typvar(1)%name='hflx_atl'
        typvar(1)%long_name='Heat_Fluxes_Atlantic'
        typvar(1)%short_name='hflx_atl'

        typvar(2)%name='hflx_indopacif'
        typvar(2)%long_name='Heat_Fluxes_Indo-Pacific'
        typvar(2)%short_name='hflx_indopacif'

        typvar(3)%name='hflx_indian'
        typvar(3)%long_name='Heat_Fluxes_Indian'
        typvar(3)%short_name='hflx_indian'

        typvar(4)%name='hflx_pacif'
        typvar(4)%long_name='Heat_Fluxes_Pacific'
        typvar(4)%short_name='hflx_pacif'

     ENDIF
  ENDIF

  e1t(:,:) = getvar(coordhgr, 'e1t', 1,npiglo,npjglo) 
  e2t(:,:) = getvar(coordhgr, 'e2t', 1,npiglo,npjglo) 
  gphit(:,:) = getvar(coordhgr, 'gphit', 1,npiglo,npjglo)

  iloc=MAXLOC(gphit)
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

  IF (lwrtcdf) THEN

     ! create output fileset
     ncout =create(cfileoutnc,'none',kx,npjglo,npk)
     ierr= createvar(ncout,typvar,nboutput,ipk,id_varout )
     ierr= putheadervar(ncout, cfilet ,kx, npjglo,npk,pnavlon=dumlon,pnavlat=dumlat)
     tim=getvar1d(cfilet,'time_counter',1)
     ierr=putvar1d(ncout,tim,1,'T')

     ! netcdf output 
     DO jj=1, jpbasins
        ierr = putvar(ncout, id_varout(jj), htrp(jj,:), ipk(jj), kx, npjglo )
     END DO

     ierr = closeout(ncout)

  ENDIF

END PROGRAM cdfhflx
