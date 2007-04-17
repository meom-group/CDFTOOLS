PROGRAM cdfmxl
  !!---------------------------------------------------------------------
  !!               ***  PROGRAM cdfmxl  ***
  !!
  !!  **  Purpose: Compute mixed layer depth 
  !!  
  !!  **  Method: Try to avoid 3 d arrays.
  !!            - compute surface properties
  !!            - initialize depths and model levels number
  !!            - from bottom to top compute rho and
  !!              check if rho > rho_surf +rho_c
  !!              where rho_c is a density criteria given as argument
  !!
  !! history :
  !!   Original :  J.M. Molines (October 2005)
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
  INTEGER :: ji,jj,jk, ik1, ik2, ikt
  INTEGER :: narg, iargc
  INTEGER :: npiglo, npjglo, npk
  INTEGER ,  DIMENSION(:,:), ALLOCATABLE :: mbathy       !: number of w levels in water <= npk
  INTEGER ,  DIMENSION(:,:), ALLOCATABLE :: nmln1  ,&    !: last level where rho > rho + rho_c1
      &                                     nmln2  ,&    !: last level where rho > rho + rho_c2
      &                                     nmlnt        !: last level where temp > temp +temp_c  (temp_c<0)

  REAL(KIND=4) :: rho_c1=0.01, rho_c2=0.03,  temp_c=-0.2
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: temp, &   !: temperatures
      &                                        sal,  &   !: salinity
      &                                        rho,  &   !: current density
      &                                    rho_surf, &   !: surface density
      &                                    tem_surf, &   !: surface temperature
      &                                    hmlp1   , &   !: mixed layer depth based on density criterium 1
      &                                    hmlp2   , &   !: mixed layer depth based on density criterium 2
      &                                    hmlt    , &   !: mixed layer depth based on Temp Criterium
      &                              zmask_surf    , &   !: mixed layer depth
      &                                    zmask         !: tmask at current level
  REAL(KIND=4), DIMENSION(:), ALLOCATABLE   :: gdepw     !: depth of w levels

  CHARACTER(LEN=80) :: cfilet,  coordzgr='mesh_zgr.nc'
  CHARACTER(LEN=80) :: cbathy='bathy_level.nc'

  ! output stuff
  INTEGER                         :: ncout, ierr
  INTEGER,           DIMENSION(3) :: ipk, id_varout  !: only one output variable
  REAL(KIND=4),      DIMENSION(1) :: tim,dep       !: time output
  CHARACTER(LEN=80)               :: cfileout='mxl.nc'

  TYPE(variable), DIMENSION(3)    :: typvar        !: structure for attributes

  LOGICAL  :: lexist                               !: flag for existence of bathy_level file

  !! 0- Get started ..
  !!
  !  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg < 1  ) THEN
     PRINT *,' Usage : cdfmxl gridTfile '
     PRINT *,' Files  mesh_zgr.nc must be in the current directory ^**'
     PRINT *,' Output on mxl.nc '
     PRINT *,'          variable somxl010 = mld on density criterium 0.01'
     PRINT *,'          variable somxl030 = mld on density criterium 0.03'
     PRINT *,'          variable somxlt02 = mld on temperature criterium -0.2'
     PRINT *,'  ^** : In case of FULL STEP run, bathy_level.nc must also be in the directory'
     STOP
  ENDIF
  CALL getarg (1, cfilet)

  ! read dimensions 
  npiglo= getdim (cfilet,'x')
  npjglo= getdim (cfilet,'y')
  npk   = getdim (cfilet,'depth')

  dep(1) = 0.
  ipk(:)= npk  ! all variables ( output are 2D)

  typvar(1)%name= 'somxl010'
  typvar(2)%name= 'somxl030'
  typvar(3)%name= 'somxlt02'
  typvar%units='m'
  typvar%missing_value=0.
  typvar%valid_min= 0.
  typvar%valid_max= 7000.
  typvar(1)%long_name='Mixed_Layer_Depth_on_0.01_rho_crit'
  typvar(2)%long_name='Mixed_Layer_Depth_on_0.03_rho_crit'
  typvar(3)%long_name='Mixed_Layer_Depth_on_-0.2_temp_crit'
  typvar(1)%short_name='somxl010'
  typvar(2)%short_name='somxl030'
  typvar(3)%short_name='somxlt02'
  typvar%online_operation='N/A'
  typvar%axis='TYX'


  PRINT *, 'npiglo=', npiglo
  PRINT *, 'npjglo=', npjglo
  PRINT *, 'npk   =', npk

  ALLOCATE (temp(npiglo,npjglo), sal(npiglo,npjglo), rho(npiglo,npjglo))
  ALLOCATE (rho_surf(npiglo,npjglo) ,tem_surf(npiglo,npjglo))
  ALLOCATE (zmask(npiglo,npjglo),zmask_surf(npiglo,npjglo)    )
  ALLOCATE (hmlp1(npiglo,npjglo),hmlp2(npiglo,npjglo), hmlt(npiglo,npjglo) )
  ALLOCATE (mbathy(npiglo,npjglo) )
  ALLOCATE (nmln1(npiglo,npjglo), nmln2(npiglo,npjglo),nmlnt(npiglo,npjglo) )
  ALLOCATE ( gdepw(npk) )
  
  ! read mbathy and gdepw use real temp(:,:) as template (getvar is used for real only)
  INQUIRE (FILE=cbathy, EXIST=lexist)
  IF ( lexist ) THEN
    temp(:,:) = getvar(cbathy,'Bathy_level',1, npiglo, npjglo)
  ELSE
    temp(:,:) = getvar(coordzgr,'mbathy',1, npiglo, npjglo)
  ENDIF
  mbathy(:,:) = temp(:,:)
  gdepw(:) = getvare3(coordzgr, 'gdepw',npk)

  !! 1- Get surface properties
  !!
  ! read surface T and S and deduce land-mask from salinity
  temp(:,:) = getvar(cfilet, 'votemper',  1 ,npiglo,npjglo)
  sal (:,:) = getvar(cfilet, 'vosaline',  1 ,npiglo,npjglo)
  zmask(:,:) = 1.; WHERE ( sal == 0. ) zmask = 0.
  zmask_surf(:,:) = zmask(:,:)
 
  ! compute rho_surf
  rho_surf(:,:) = sigma0 ( temp,sal,npiglo,npjglo )* zmask(:,:)
  tem_surf(:,:) = temp(:,:)

  ! Initialization to the number of w ocean point mbathy
  nmln1(:,:) = mbathy(:,:)
  nmln2(:,:) = mbathy(:,:)
  nmlnt(:,:) = mbathy(:,:)

 !! 2- determine mixed layer
 !!
      ! Last w-level at which rhop>=rho surf+rho_c (starting from jpk-1)
      ! (rhop defined at t-point, thus jk-1 for w-level just above)
      DO jk = npk-1, 2, -1
         temp(:,:) = getvar(cfilet, 'votemper',  jk ,npiglo,npjglo)
         sal (:,:) = getvar(cfilet, 'vosaline',  jk ,npiglo,npjglo)
         zmask(:,:) = 1.  ; WHERE ( sal == 0. ) zmask = 0.
         rho(:,:)  = sigma0 ( temp,sal,npiglo,npjglo )* zmask(:,:)

         DO jj = 1, npjglo
            DO ji = 1, npiglo
               IF( rho(ji,jj)  > rho_surf(ji,jj) + rho_c1 )   nmln1(ji,jj) = jk
               IF( rho(ji,jj)  > rho_surf(ji,jj) + rho_c2 )   nmln2(ji,jj) = jk
               IF( ABS(temp(ji,jj) - tem_surf(ji,jj)) > ABS( temp_c)  )   nmlnt(ji,jj) = jk
            END DO
         END DO
      END DO

      ! Mixed layer depth
      DO jj = 1, npjglo
         DO ji = 1, npiglo
            ik1 = nmln1(ji,jj) ; ik2 = nmln2(ji,jj) ; ikt = nmlnt(ji,jj)
            hmlp1 (ji,jj) = gdepw(ik1) * zmask_surf(ji,jj)
            hmlp2 (ji,jj) = gdepw(ik2) * zmask_surf(ji,jj)
            hmlt (ji,jj)  = gdepw(ikt) * zmask_surf(ji,jj)
         END DO
      END DO

 !! 3- Write output file
 !!
 ncout = create(cfileout, cfilet, npiglo,npjglo,1)
 ierr = createvar(ncout ,typvar ,3, ipk,id_varout )
 ierr= putheadervar(ncout, cfilet,npiglo, npjglo,1,pdep=dep)
 tim=getvar1d(cfilet,'time_counter',1)
 ierr = putvar(ncout, id_varout(1) ,hmlp1, 1,npiglo, npjglo)
 ierr = putvar(ncout, id_varout(2) ,hmlp2, 1,npiglo, npjglo)
 ierr = putvar(ncout, id_varout(3) ,hmlt , 1,npiglo, npjglo)
 ierr=putvar1d(ncout,tim,1,'T')

 ierr=closeout(ncout)


END PROGRAM cdfmxl
