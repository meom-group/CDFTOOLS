PROGRAM cdfmocsig
  !!-------------------------------------------------------------------
  !!               ***  PROGRAM cdfmocsig  ***
  !!
  !!  **  Purpose  :  Compute the Meridional Overturning Cell (MOC)
  !!                  PARTIAL STEPS  in density coordinates.
  !!  
  !!  **  Method   :  The MOC is computed from the V velocity field, collected in density bins,
  !!                  (reference depth is given as the 3rd argument) and integrated
  !!                  throughout the density bins, then zonally averaged with
  !!                  eventual masking for oceanic basins.
  !!                  In the present version the masking corresponds to the global
  !!                  configuration. MOC for Global, Atlantic, Indo-Pacific, Indian,Pacific ocean
  !!                  Results are saved on mocsig.nc file with variables name respectively
  !!                  zomsfglo, zomsfatl, zomsfinp, zomsfind, zomsfpac.
  !!                  If no new_maskglo.nc file found, then the mask.nc file is used and
  !!                  only zomsfglo is computed.
  !!
  !! history :
  !!  Original :  J.M. Molines (jul. 2005) 
  !!     moc in density : A.M. Treguier ( Nov. 2005)
  !!     choice of reference depth, improvements: C. Dufour (March 2010)
  !!-------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------
  !! * Modules used
  USE cdfio
  USE eos

  !! * Local variables
  IMPLICIT NONE
  !  FOR sigma 1 as the density coordinate 
  CHARACTER(LEN=256) :: cdref                   !:  depth reference for density read in argument
  REAL(KIND=4)       :: pref=0.                 !:  depth reference for density used in the code
  INTEGER            :: jpbin                   !:   density  bins
  REAL(KIND=4)       :: s1min, s1scal           !:  reference for density

  INTEGER   :: jpbasins
  INTEGER   :: jbasin, jj, jk ,ji, jkk             !: dummy loop index
  INTEGER   :: ierr                                !: working integer
  INTEGER   :: narg, iargc, iarg                   !: command line 
  INTEGER   :: npiglo,npjglo, npk                  !: size of the domain
  INTEGER   :: ncout, np
  INTEGER   :: numout=10
  INTEGER, DIMENSION(:), ALLOCATABLE ::  ipk, id_varout         !
  INTEGER, DIMENSION(2)          ::  iloc

  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  e1v,  gphiv          !:  metrics, velocity
  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  zt, zs, zv, zveiv, e3v       !:  metrics, velocity
  INTEGER,      DIMENSION (:,:),     ALLOCATABLE ::  ibin   !: integer value corresponding to the density for binning
  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  dumlon               !: dummy longitude = 0.
  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  dumlat               !: latitude for i = north pole
  REAL(KIND=4), DIMENSION (:,:),     ALLOCATABLE ::  zttmp,zstmp, zmask2d !: arrays to call sigmai and mask it
  REAL(KIND=4), DIMENSION (:),       ALLOCATABLE ::  gdepw                !: deptw
  REAL(KIND=4), DIMENSION (:),       ALLOCATABLE ::  sigma                !: density coordinate
  REAL(KIND=4), DIMENSION (:,:,:),   ALLOCATABLE ::  zmask                !:  jpbasins x npiglo x npjglo
  REAL(KIND=4), DIMENSION (:,:,:),   ALLOCATABLE ::  zread                !:  jpi,1,jpk
  REAL(KIND=4), DIMENSION (1)                    ::  tim

  REAL(KIND=8) ,DIMENSION(:,:)   , ALLOCATABLE   ::  zdens                 !: density
  REAL(KIND=8) ,DIMENSION(:,:)   , ALLOCATABLE   ::  zomsftmp              !: temporary transport array
  REAL(KIND=8) ,DIMENSION(:,:,:) , ALLOCATABLE   ::  zomsf                 !: jpbasins x npjglo x npk

  CHARACTER(LEN=256) :: cfilev , cfilet, cfileoutnc='mocsig.nc'
  CHARACTER(LEN=256) :: coordhgr='mesh_hgr.nc',  coordzgr='mesh_zgr.nc',cbasinmask='new_maskglo.nc'
  CHARACTER(LEN=255) :: cglobal         !: Global attribute to trace the building command.
  TYPE(variable)    ,DIMENSION(:), ALLOCATABLE   :: typvar
  LOGICAL    :: llglo = .false.                            !: indicator for presence of new_maskglo.nc file

  REAL (KIND=4)      :: ztrans 

  INTEGER    :: istatus 
  LOGICAL    :: lprint = .false.
  LOGICAL    :: leiv = .false.

  ! constants
  lprint = .false. 
  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' Usage : cdfmocsig  V_file T_file depth_ref [-eiv]'
     PRINT *,' Computes the MOC for oceanic basins as described in new_maskglo.nc'
     PRINT *,' Last arguments is the reference depth for potential density, in m'
     PRINT *,' actually only 0 1000 or 2000 are allowed'
     PRINT *,'     Option -eiv : takes into account VEIV Meridional eddy induced velocity'
     PRINT *,'     -> To be used only if Gent and McWilliams parameterization has been used '
     PRINT *,' PARTIAL CELLS VERSION'
     PRINT *,' Files mesh_hgr.nc, mesh_zgr.nc ,new_maskglo.nc ,mask.nc '
     PRINT *,'  must be in the current directory'
     PRINT *,' Output on mocsig.nc: '
     PRINT *,'      variables zomsfglo  : Global ocean '
     PRINT *,'      variables zomsfatl  : Atlantic Ocean '
     PRINT *,'      variables zomsfinp  : Indo Pacific '
     PRINT *,'      variables zomsfind  : Indian Ocean alone'
     PRINT *,'      variables zomsfpac  : Pacific Ocean alone'
     PRINT *,' If new_maskglo.nc is not present mask file is used and only zomsfglo is'
     PRINT *,' produced'
     STOP
  ENDIF
  !!   density coordinate is sigma1. 
  !!   bins are by 0.1 from 24 to 32.6  = 87 bins 

  CALL getarg (1, cfilev)
  CALL getarg (2, cfilet)
  CALL getarg (3, cdref)
  READ(cdref,*) pref 
  IF (narg > 3 ) THEN
    iarg=4
    leiv=.TRUE.
  ENDIF

  npiglo= getdim (cfilev,'x')
  npjglo= getdim (cfilev,'y')
  npk   = getdim (cfilev,'depth')

  !setting up the building command in global attribute
  WRITE(cglobal,'(a,1x,a,1x,a,1x,a) ') 'cdfmocsig', TRIM(cfilev), TRIM(cfilet), TRIM(cdref)
  !  Detects newmaskglo file 
  INQUIRE( FILE='new_maskglo.nc', EXIST=llglo )
  IF (llglo) THEN
     jpbasins = 5
  ELSE
     jpbasins = 1
  ENDIF

  ALLOCATE ( typvar(jpbasins), ipk(jpbasins), id_varout(jpbasins) )

  ! Define density parameters
  IF ( pref==0 ) THEN
     IF (lprint) PRINT *, ' to be implemented soon... '
     jpbin = 52
     s1min = 23
     s1scal=0.1
  ELSE IF ( pref==1000 ) THEN
     jpbin = 88
     s1min = 24
     s1scal=0.1
  ELSE IF ( pref==2000 ) THEN
     jpbin = 158
     s1min = 30
     s1scal=0.05
  ELSE
     IF (lprint) PRINT *, ' this value of depth_ref is not yet implemented '
  ENDIF

  ALLOCATE ( sigma(jpbin) )  

  ! define densities at middle of bins
  DO ji=1,jpbin
     sigma(ji)  = s1min +(ji-0.5)*s1scal
  ENDDO
  IF (lprint) print *, ' min density:',sigma(1), ' max density:', sigma(jpbin)
  ! define new variables for output ( must update att.txt)

  typvar(1)%name= 'zomsfglo'
  typvar%units='Sverdrup'
  typvar%missing_value=99999.
  typvar%valid_min= -1000.
  typvar%valid_max= 1000.
  typvar%scale_factor= 1.
  typvar%add_offset= 0.
  typvar%savelog10= 0.
  typvar(1)%long_name='Meridional_Overt.Cell_Global'
  typvar(1)%short_name='zomsfglo'
  typvar%online_operation='N/A'
  typvar%axis='TZY'

  ipk(1) = npk  !  2D

  IF (llglo) THEN
     typvar(2)%name= 'zomsfatl'
     typvar(2)%long_name='Meridional_Overt.Cell_Atlantic'
     typvar(2)%short_name='zomsfatl'

     typvar(3)%name= 'zomsfinp'
     typvar(3)%long_name='Meridional_Overt.Cell_IndoPacif'
     typvar(3)%short_name='zomsfinp'

     typvar(4)%name= 'zomsfind'
     typvar(4)%long_name='Meridional_Overt.Cell_Indian'
     typvar(4)%short_name='zomsfind'

     typvar(5)%name= 'zomsfpac'
     typvar(5)%long_name='Meridional_Overt.Cell_pacif'
     typvar(5)%short_name='zomspac'

     ipk(2) = npk  !  2D
     ipk(3) = npk  !  2D
     ipk(4) = npk  !  2D
     ipk(5) = npk  !  2D
  ENDIF


  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk

  ! Allocate arrays
  ALLOCATE ( zmask(jpbasins,npiglo,npjglo) )
  ALLOCATE ( zv (npiglo,npjglo), zt(npiglo,npjglo), zs(npiglo,npjglo))
  IF ( leiv ) THEN
    ALLOCATE ( zveiv (npiglo,npjglo))
  END IF
  ALLOCATE ( e3v(npiglo,npjglo) )
  ALLOCATE ( ibin(npiglo, npjglo) )
  ALLOCATE ( e1v(npiglo,npjglo), gphiv(npiglo,npjglo) ,gdepw(npk) )
  ALLOCATE ( zomsf(jpbasins, npjglo, jpbin) )
  ALLOCATE ( zomsftmp(jpbin,npiglo) )
  ALLOCATE ( dumlon(1,npjglo) , dumlat(1,npjglo))
  ALLOCATE ( zdens(npiglo,npjglo))
  ALLOCATE ( zmask2d(npiglo,npjglo), zttmp(npiglo,npjglo))


  e1v(:,:) = getvar(coordhgr, 'e1v', 1,npiglo,npjglo) 
  gphiv(:,:) = getvar(coordhgr, 'gphiv', 1,npiglo,npjglo)
  gdepw(:)   = getvare3(coordzgr, 'gdepw',npk)
  gdepw(:) = -1.*  gdepw(:)

  iloc=maxloc(gphiv)
  dumlat(1,:) = gphiv(iloc(1),:)
  dumlon(:,:) = 0.   ! set the dummy longitude to 0

  ! create output fileset
  IF (lprint) PRINT *, ' ready to create file:',trim( cfileoutnc), ' from reference:',trim(cfilev )
  ncout =create(cfileoutnc, cfilev, 1,npjglo,jpbin,cdep='sigma')
  IF (lprint) PRINT *, ' ready to create variables:'
  ierr= createvar(ncout ,typvar ,jpbasins, ipk ,id_varout ,cdglobal=cglobal)
  IF (lprint) PRINT *, ' writing variables headers:'
  ierr= putheadervar(ncout, cfilev,1, npjglo,jpbin,pnavlon=dumlon,pnavlat=dumlat,pdep=sigma)
  IF (lprint) PRINT *, ' writing time_counter:'
  tim=getvar1d(cfilev,'time_counter',1)
  ierr=putvar1d(ncout,tim,1,'T')


  ! reading the masks
  ! 1 : global ; 2 : Atlantic ; 3 : Indo-Pacif ; 4 : Indian ; 5 : Pacif
  zmask(1,:,:)=getvar('mask.nc','vmask',1,npiglo,npjglo)
  !  zmask(1,:,:)=  1.
  IF ( llglo ) THEN
     zmask(2,:,:)=getvar(cbasinmask,'tmaskatl',1,npiglo,npjglo)
     zmask(4,:,:)=getvar(cbasinmask,'tmaskind',1,npiglo,npjglo)
     zmask(5,:,:)=getvar(cbasinmask,'tmaskpac',1,npiglo,npjglo)
     zmask(3,:,:)=zmask(5,:,:)+zmask(4,:,:)
     ! ensure that there are no overlapping on the masks
     WHERE(zmask(3,:,:) > 0 ) zmask(3,:,:) = 1
     ! change global mask for GLOBAL periodic condition
     zmask(1,1,:) = 0.
     zmask(1,npiglo,:) = 0.
  ENDIF

  ! initialize moc to 0
  zomsf(:,:,:) = 0.

  DO jk = 1,npk-1
     !               for testing purposes only loop from 2 to 400
     IF (lprint) PRINT *,' working at depth ',jk
     ! Get velocities v at jj
     zv(:,:)= getvar(cfilev, 'vomecrty', jk,npiglo,npjglo)
     IF ( leiv ) THEN
       zveiv(:,:)= getvar(cfilev, 'vomeeivv', jk,npiglo,npjglo)
     END IF
     zt(:,:)= getvar(cfilet, 'votemper', jk,npiglo,npjglo)
     zs(:,:)= getvar(cfilet, 'vosaline', jk,npiglo,npjglo)
     ! get e3v at latitude jj
     e3v(:,:) = getvar(coordzgr, 'e3v_ps', jk,npiglo,npjglo, ldiom=.true.)
     !
     !  finds density 
     ! 
     zmask2d =  1
     WHERE(zt ==0) zmask2d = 0
     zdens = sigmai(zt,zs,pref,npiglo,npjglo)
     zttmp= zdens* zmask2d !: convert to single precision 
     ibin(:,:) = ifix( (zttmp-s1min)/s1scal )
     ibin(:,:) = max( ibin(:,:) ,1)
     ibin(:,:) = min(ibin(:,:),jpbin)
     DO jj=2,npjglo-1
        zomsftmp = 0
        !  converts transport in "k" to transport in "sigma"
        !  indirect adresssing - do it once and not for each basin!
        DO ji=2,npiglo-1
           IF ( leiv ) THEN
             ztrans =   e1v(ji,jj)*e3v(ji,jj)*(zv(ji,jj)+zveiv(ji,jj))
           ELSE
             ztrans =   e1v(ji,jj)*e3v(ji,jj)*zv(ji,jj)
           END IF
           zomsftmp(ibin(ji,jj),ji)=zomsftmp(ibin(ji,jj),ji) - ztrans
        END DO
        ! integrates 'zonally' (along i-coordinate) 
        ! add to zomsf the contributions from level jk  at all densities jkk
        DO jkk =1,jpbin  
           DO ji=2,npiglo-1
              DO jbasin= 1, jpbasins
                 ! For all basins 
                 ztrans =   zomsftmp(jkk,ji) * zmask(jbasin,ji,jj) 
                 zomsf(jbasin,jj,jkk)=zomsf(jbasin,jj,jkk ) + ztrans
              ENDDO
           END DO
        END DO
        !               end of loop on latitude for filling zomsf
     END DO
     !  end of loop on depths for calculating transports     
  END DO

  ! integrates vertically   from bottom to surface
  zomsf(:,:,jpbin) = zomsf(:,:,jpbin)/1.e6
  DO jk=jpbin-1,1,-1
     zomsf(:,:,jk) = zomsf(:,:,jk+1) + zomsf(:,:,jk)/1.e6
  END DO  ! loop to next level

  ! netcdf output 
  DO jbasin= 1, jpbasins
     DO jk =1, jpbin
        ierr = putvar (ncout, id_varout(jbasin),REAL(zomsf(jbasin,:,jk)), jk,1,npjglo)
     END DO
  END DO

  ierr = closeout(ncout)

END PROGRAM cdfmocsig
   
