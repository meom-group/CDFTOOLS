PROGRAM cdfmxlhcsc
  !!---------------------------------------------------------------------
  !!               ***  PROGRAM cdfmxlhcsc  ***
  !!
  !!  **  Purpose: Compute mixed layer depth and the heat and salt contents
  !!               in the mixed layer. There is an option to limit this computation
  !!               between hmin and ml depth. For that, hmin is given as last argument (>0)
  !!               with no arguments, hmin os supposed to be 0.
  !!  
  !!  **  Method to compute MLD: 
  !!         Try to avoid 3 d arrays.
  !!            - compute surface properties
  !!            - initialize depths and model levels number
  !!            - from bottom to top compute rho and
  !!              check if rho > rho_surf +rho_c
  !!              where rho_c is a density criteria given as argument
  !! **  Method to compute HC/SC in MLD:
  !!         Compute the sum ( rho cp T * e1 * e2 * e3 * mask )
  !!                         (     cp S * e1 * e2 * e3 * mask )
  !!         for the MLD computed before, or in the water volume between hmin and MLD
  !!
  !! history :
  !!   cdfmxl.f90 : Original :  J.M. Molines (October 2005)
  !!   cdfmxlheatc.f90 : Original :  J.M. Molines (April 2006)
  !!         Merging programs : M. Juza ( April 2007 )
  !!---------------------------------------------------------------------
  !!  $Rev$
  !!  $Date$
  !!  $Id$
  !!--------------------------------------------------------------

  !! * Modules used
  USE cdfio
  USE eos

  !! * Local variables
  IMPLICIT NONE
  INTEGER :: ji,jj,jk, ik
  INTEGER :: narg, iargc                               !: command line
  INTEGER :: npiglo, npjglo, npk                       !: size of the domain
  INTEGER ,  DIMENSION(:,:), ALLOCATABLE :: mbathy     !: number of w levels in water <= npk
  INTEGER ,  DIMENSION(:,:), ALLOCATABLE :: nmln       !: last level where rho > rho + val_crit
  !:              or temp > temp + val_crit
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE ::temp, & !: temperatures
       &                                       sal,  & !: salinity
       &                                       rho,  & !: current density
       &                                   rho_surf, & !: surface density
       &                                   tem_surf, & !: surface temperature
       &                                   hmld    , & !: mixed layer depth based on criterium
       &                             zmask_surf    , & !: mixed layer depth
       &                                   zmask       !: tmask at current level
  REAL(KIND=4), DIMENSION(:), ALLOCATABLE   :: gdepw  !: depth of w levels
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: mask,zt,zs,e3
  REAL(KIND=8), PARAMETER :: rprho0=1020., rpcp=4000.
  REAL(KIND=8) :: zvol,zsum,zvol2d,zsurf, zsal3d
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: zmxlheatc,zmxlsaltc
  REAL(KIND=4) :: val, hmin = 0.

  TYPE(variable),DIMENSION(3) :: typvar  !: extension for attributes

  CHARACTER(LEN=80) :: cfilet,critere,cdum
  CHARACTER(LEN=80) :: coordzgr='mesh_zgr.nc', coordhgr='mesh_hgr.nc' , cmask='mask.nc'

  ! output stuff
  INTEGER                         :: ncout, ierr
  INTEGER,           DIMENSION(3) :: ipk, id_varout  !: only one output variable
  REAL(KIND=4),      DIMENSION(1) :: tim,dep       !: time output
  CHARACTER(LEN=80), DIMENSION(3) :: cvarname
  CHARACTER(LEN=80)               :: cfileout='mxlhcsc.nc'

  !! 0- Get started ..
  !!
  !  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg < 1  ) THEN
     PRINT *,' Usage : cdfmxlhcsc gridTfile crit val [hmin] '
     PRINT *,'           crit = ''temperature'' or ''density'' criterium '
     PRINT *,'           val  = value of the criterium '
     PRINT *,'           [hmin] = optional. If given limit depth integral from hmin to mld'
     PRINT *,'           [hmin] = 0 by defaul, ie, compute over the whole mixed layer'
     PRINT *,'Compute MLD and HC/SC in the MLD'
     PRINT *,' Files mesh_hgr.nc, mesh_zgr.nc , mask.nc must be in the current directory'
     PRINT *,' Output on mxlhcsc.nc '
     PRINT *,' Output variables : -  somxl010 (mld based on density criterium 0.01)'
     PRINT *,'      (2D)          or somxl030 (mld on density criterium 0.03)'
     PRINT *,'                    or somxlt02 (mld on temperature criterium -0.2)'
     PRINT *,'                    -  somxlheatc (heat content computed in the MLD)'
     PRINT *,'                    -  somxlsaltc (salt content computed in the MLD)'
     STOP
  ENDIF
  CALL getarg (1, cfilet)
  CALL getarg (2, critere)
  CALL getarg (3, cdum) ;  READ(cdum,*) val
  IF ( narg == 4 ) THEN  ; CALL getarg (4, cdum) ;  READ(cdum,*) hmin ; ENDIF

  ! read dimensions 
  npiglo= getdim (cfilet,'x')
  npjglo= getdim (cfilet,'y')
  npk   = getdim (cfilet,'depth')

  dep(1) = 0.
  ipk(:) = 1 

  ! Variable Mixed Layer Depth
  IF ( critere == 'temperature' .AND. val==-0.2) THEN
     typvar(1) %name = 'somxlt02'
     typvar(1) %short_name = 'somxlt02'
  ELSE IF ( critere == 'density' .AND. val==0.01) THEN
     typvar(1) %name = 'somxl010'
     typvar(1) %short_name = 'somxl010'
  ELSE IF ( critere == 'density' .AND. val==0.03) THEN
     typvar(1) %name = 'somxl030'
     typvar(1) %short_name = 'somxl030'
  ENDIF
  typvar(1) %units = 'm'
  typvar(1) %long_name = ' Mixed Layer Depth'
  ! Variable Heat Content
  typvar(2) %name = 'somxlheatc'
  typvar(2) %units = '10^9 J/m2'
  typvar(2) %long_name = 'Mixed_Layer_Heat_Content'
  typvar(2) %short_name = 'somxlheatc'
  ! Variable Salt Content
  typvar(3) %name = 'somxlsaltc'
  typvar(3) %units = '10^6 kg/m2'
  typvar(3) %long_name = 'Mixed_Layer_Salt_Content'
  typvar(3) %short_name = 'somxlsaltc'


  !PRINT *, 'npiglo=', npiglo
  !PRINT *, 'npjglo=', npjglo
  !PRINT *, 'npk   =', npk


  ! Allocate arrays
  ALLOCATE (temp(npiglo,npjglo),sal(npiglo,npjglo))
  ALLOCATE (zmask(npiglo,npjglo),zmask_surf(npiglo,npjglo))
  ALLOCATE (mbathy(npiglo,npjglo))
  ALLOCATE (nmln(npiglo,npjglo),hmld(npiglo,npjglo))
  ALLOCATE (mask(npiglo,npjglo))
  ALLOCATE (zmxlheatc(npiglo,npjglo),zmxlsaltc(npiglo,npjglo))
  ALLOCATE (zt(npiglo,npjglo),zs(npiglo,npjglo))
  ALLOCATE (e3(npiglo,npjglo))
  ALLOCATE (gdepw(npk))

  ! read mbathy and gdepw use real temp(:,:) as template (getvar is used for real only)
  temp(:,:) = getvar(coordzgr,'mbathy',1, npiglo, npjglo)
  mbathy(:,:) = temp(:,:)
  gdepw(:) = getvare3(coordzgr, 'gdepw',npk) 

  ! Initialization to the number of w ocean point mbathy
  nmln(:,:) = mbathy(:,:)

  !! 1- Get surface properties
  !!
  ! read surface T and S and deduce land-mask from salinity
  temp(:,:) = getvar(cfilet, 'votemper',  1 ,npiglo,npjglo)
  sal (:,:) = getvar(cfilet, 'vosaline',  1 ,npiglo,npjglo)
  zmask(:,:) = 1.; WHERE ( sal == 0. ) zmask = 0.
  zmask_surf(:,:) = zmask(:,:)


   SELECT CASE ( critere )
   CASE ( 'temperature','Temperature' ) !!!!! Temperature criterium
     ! temp_surf
     ALLOCATE (tem_surf(npiglo,npjglo))
     tem_surf(:,:) = temp(:,:)   

     ! Last w-level at which ABS(temp-temp_surf)>=ABS(val) (starting from jpk-1)
     ! (temp defined at t-point, thus jk-1 for w-level just above)
     DO jk = npk-1, 2, -1
        temp(:,:) = getvar(cfilet, 'votemper',  jk ,npiglo,npjglo)
        DO jj = 1, npjglo
           DO ji = 1, npiglo
              IF( ABS(temp(ji,jj) - tem_surf(ji,jj)) > ABS(val)  )   nmln(ji,jj) = jk
           ENDDO
        ENDDO
     ENDDO


   CASE ( 'density', 'Density' ) !!!!! Density criterium
     ! compute rho_surf
     ALLOCATE (rho(npiglo,npjglo))
     ALLOCATE (rho_surf(npiglo,npjglo))
     rho_surf(:,:) = sigma0 ( temp,sal,npiglo,npjglo )* zmask(:,:)

     ! Last w-level at which rhop>=rho surf+rho_c (starting from jpk-1)
     ! (rhop defined at t-point, thus jk-1 for w-level just above)
     DO jk = npk-1, 2, -1
        temp(:,:) = getvar(cfilet, 'votemper',  jk ,npiglo,npjglo)
        sal (:,:) = getvar(cfilet, 'vosaline',  jk ,npiglo,npjglo)
        zmask(:,:) = 1.  ; WHERE ( sal == 0. ) zmask = 0.
        rho(:,:)  = sigma0 ( temp,sal,npiglo,npjglo )* zmask(:,:)
        DO jj = 1, npjglo
           DO ji = 1, npiglo
              IF( rho(ji,jj)  > rho_surf(ji,jj) + val )   nmln(ji,jj) = jk
           ENDDO
        ENDDO
     ENDDO
  
    CASE DEFAULT
        PRINT *,' ERROR: Criterium on ', TRIM(critere),' not suported' ; STOP
    END SELECT

  !! 2- Determine mixed layer depth
  DO jj = 1, npjglo
     DO ji = 1, npiglo
        ik = nmln(ji,jj)
        hmld (ji,jj) = gdepw(ik) * zmask_surf(ji,jj)
     ENDDO
  ENDDO

  !! 3- Compute heat and salt contents in the mixed layer depth
  !!
  zvol=0.d0
  zsum=0.d0
  zmxlheatc(:,:)=0.d0
  zmxlsaltc(:,:)=0.d0
  DO jk = 1,npk
     ! Get temperature and salinity at jk
     zt(:,:)  = getvar(cfilet, 'votemper',jk ,npiglo,npjglo)
     zs(:,:)  = getvar(cfilet, 'vosaline',jk ,npiglo,npjglo)
     mask(:,:)= getvar(cmask,  'tmask',   jk ,npiglo,npjglo)
     ! Get e3 at level jk (ps...)
     e3(:,:)  = getvar(coordzgr,'e3t_ps',  jk ,npiglo,npjglo)
     ! e3 is used as a flag for the mixed layer; it is 0 outside the mixed layer
     e3(:,:)=MAX(0.,MIN(e3,hmld-gdepw(jk)) + MIN(e3,gdepw(jk)+e3-hmin)-e3)

     ! JMM : I think next line is useless as the masking of the product by e3...
     WHERE ( e3 == 0 ) mask = 0.
     ! Heat and salt contents
     zvol2d=SUM(e3*mask)
     zmxlheatc(:,:)=zmxlheatc(:,:) + zt*e3*mask
     zmxlsaltc(:,:)=zmxlsaltc(:,:) + zs*e3*mask
     ! We want to scan all deptht from to to bottom and as we eventually skip between surf and hmin ...
!    IF ( zvol2d /= 0 ) THEN
!       ! go on !
!    ELSE
!       ! no more layer below !
!       EXIT  ! get out of the jk loop
!    ENDIF

  ENDDO

  !! Heat and salt contents (10^9.J/m2 and 10^6.kg/m2)
  zmxlheatc = zmxlheatc*rprho0*rpcp*(10.)**(-9)
  zmxlsaltc = zmxlsaltc*rprho0*(10.)**(-6)


  !! 4- Write output file
  !!
  ncout = create(cfileout, cfilet, npiglo,npjglo,1)
  ierr = createvar(ncout ,typvar,3, ipk,id_varout )
  ierr = putheadervar(ncout, cfilet,npiglo, npjglo,1,pdep=dep)
  tim  = getvar1d(cfilet,'time_counter',1)
  ierr = putvar(ncout, id_varout(1) ,hmld, 1,npiglo, npjglo)
  ierr = putvar(ncout, id_varout(2) ,REAL(zmxlheatc), 1,npiglo, npjglo)
  ierr = putvar(ncout, id_varout(3) ,REAL(zmxlsaltc), 1,npiglo, npjglo)
  ierr = putvar1d(ncout,tim,1,'T')

  ierr = closeout(ncout)


END PROGRAM cdfmxlhcsc
