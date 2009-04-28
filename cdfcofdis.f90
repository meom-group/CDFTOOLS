PROGRAM cdfcofdis
  !!--------------------------------------------------------------------------
  !!                   ***  PROGRAM  cdfcofdis   ***
  !!
  !! ** Purpose : wrap for standalone cofdis routine from OPA
  !! 
  !! ** Method : define required arrays and variables in main
  !!             call cofdis
  !!             write results using cdfio instead of ioipsl
  !!
  !!-------------------------------------------------------------------------
  USE cdfio
  IMPLICIT NONE
  ! Global variables
  INTEGER :: jpi,jpj,jpk, jpim1, jpjm1, nperio=4
  INTEGER :: narg, iargc
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: glamt, glamu,glamv, glamf
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: gphit, gphiu,gphiv, gphif
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: tmask, umask, vmask, fmask
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: pdct
  ! from phycst
  REAL(KIND=4) ::   rpi = 3.141592653589793          !: pi
  REAL(KIND=4) ::   rad = 3.141592653589793 / 180.   !: conversion from degre into radian
  REAL(KIND=4) ::   ra  = 6371229.                   !: earth radius (meter)

  CHARACTER(LEN=80) :: coordhgr='mesh_hgr.nc', cmask='mask.nc'
  CHARACTER(LEN=255) :: cfilet

  ! output stuff
  INTEGER, DIMENSION(1) :: ipk, id_varout
  TYPE(variable), DIMENSION(1) :: typvar
  REAL(KIND=4) ,DIMENSION(1)                  :: timean
  CHARACTER(LEN=80) :: cfileout='dist.coast'
  INTEGER :: ncout, ierr

  !
  narg=iargc()
  IF ( narg == 0 ) THEN
    PRINT *,' USAGE: cdfcofdis mesh_hgr.nc  mask.nc gridT.nc'
    PRINT *,'   where mesh_hgr.nc and mask.nc stand for the name of the mesh_hgr'
    PRINT *,'   and mask files respectively'
    PRINT *,'   gridT.nc is used for size and depth references'
    PRINT *,' Program will output dist.coast with variable Tcoast, representing the distance of every'
    PRINT *,' T points to the coast line '
    STOP
  ENDIF
 
  CALL getarg(1,coordhgr)
  CALL getarg(2,cmask)
  CALL getarg(3,cfilet)

  ! read domain dimensions in the mask file
  jpi=getdim(cfilet,'x')
  jpj=getdim(cfilet,'y')
  jpk=getdim(cfilet,'depth')
  IF (jpk == 0 ) THEN
    jpk=getdim(cfilet,'z')
    IF ( jpk == 0 ) THEN
      PRINT *,' ERROR in determining jpk form gridT file ....'
      STOP
    ENDIF
  ENDIF
  PRINT *, jpi,jpj,jpk
  jpim1=jpi-1 ; jpjm1=jpj-1

  ! ALLOCATION of the arrays
  ALLOCATE ( glamt(jpi,jpj), glamu(jpi,jpj), glamv(jpi,jpj), glamf(jpi,jpj) )
  ALLOCATE ( gphit(jpi,jpj), gphiu(jpi,jpj), gphiv(jpi,jpj), gphif(jpi,jpj) )
  ALLOCATE ( tmask(jpi,jpj), umask(jpi,jpj), vmask(jpi,jpj), fmask(jpi,jpj) )
  ALLOCATE ( pdct(jpi,jpj) )
  
  PRINT *, 'ALLOCATION DONE.'

  ! read latitude an longitude
  glamt(:,:) = getvar(coordhgr,'glamt',1,jpi,jpj)
  PRINT *,' READ GLAMT done.'
  glamu(:,:) = getvar(coordhgr,'glamu',1,jpi,jpj)
  PRINT *,' READ GLAMU done.'
  glamv(:,:) = getvar(coordhgr,'glamv',1,jpi,jpj)
  PRINT *,' READ GLAMV done.'
  glamf(:,:) = getvar(coordhgr,'glamf',1,jpi,jpj)
  PRINT *,' READ GLAMF done.'

  gphit(:,:) = getvar(coordhgr,'gphit',1,jpi,jpj)
  PRINT *,' READ GPHIT done.'
  gphiu(:,:) = getvar(coordhgr,'gphiu',1,jpi,jpj)
  PRINT *,' READ GPHIU done.'
  gphiv(:,:) = getvar(coordhgr,'gphiv',1,jpi,jpj)
  PRINT *,' READ GPHIV done.'
  gphif(:,:) = getvar(coordhgr,'gphif',1,jpi,jpj)
  PRINT *,' READ GPHIF done.'

  ! prepare file output
  ipk(1) = jpk
  typvar(1)%name='Tcoast'
  typvar(1)%units='m'
  typvar(1)%missing_value=0
  typvar(1)%valid_min= 0.
  typvar(1)%valid_max= 1.
  typvar(1)%long_name='Tcoast'
  typvar(1)%short_name='Tcoast'
  typvar(1)%online_operation='N/A'
  typvar(1)%axis='TZYX'
  typvar(1)%precision='r4'
  PRINT *,' CREATE ...'
  ncout=create(cfileout, cfilet,jpi,jpj,jpk)

  PRINT *,' CREATEVAR ...'
  ierr= createvar(ncout ,typvar,1, ipk,id_varout )
  PRINT *,' PUTHEADERVAR ...'
  ierr= putheadervar(ncout, cfilet, jpi,jpj,jpk)
  PRINT *, 'CALL to cofdis ...'
  CALL cofdis
  
  CONTAINS
  SUBROUTINE cofdis
    !!----------------------------------------------------------------------
    !!                 ***  ROUTINE cofdis  ***
    !!
    !! ** Purpose :   Compute the distance between ocean T-points and the
    !!      ocean model coastlines. Save the distance in a NetCDF file.
    !!
    !! ** Method  :   For each model level, the distance-to-coast is 
    !!      computed as follows : 
    !!       - The coastline is defined as the serie of U-,V-,F-points
    !!      that are at the ocean-land bound.
    !!       - For each ocean T-point, the distance-to-coast is then 
    !!      computed as the smallest distance (on the sphere) between the 
    !!      T-point and all the coastline points.
    !!       - For land T-points, the distance-to-coast is set to zero.
    !!      C A U T I O N : Computation not yet implemented in mpp case.
    !!
    !! ** Action  : - pdct, distance to the coastline (argument)
    !!              - NetCDF file 'dist.coast' 
    !!----------------------------------------------------------------------
    !!
    !!
    INTEGER ::   ji, jj, jk, jl      ! dummy loop indices
    INTEGER ::   iju, ijt            ! temporary integers
    INTEGER ::   icoast, itime
    INTEGER ::   icot         ! logical unit for file distance to the coast
    LOGICAL, DIMENSION(jpi,jpj) ::   llcotu, llcotv, llcotf   ! ???
    CHARACTER (len=32) ::   clname
    REAL(KIND=4) ::   zdate0
    REAL(KIND=4), DIMENSION(jpi,jpj)   ::   zxt, zyt, zzt, zmask   ! cartesian coordinates for T-points
    REAL(KIND=4), DIMENSION(3*jpi*jpj) ::   zxc, zyc, zzc, zdis    ! temporary workspace
    !!----------------------------------------------------------------------

    ! 0. Initialization
    ! -----------------
    PRINT *, 'COFDIS init'
    zxt(:,:) = COS( rad * gphit(:,:) ) * COS( rad * glamt(:,:) )
    zyt(:,:) = COS( rad * gphit(:,:) ) * SIN( rad * glamt(:,:) )
    zzt(:,:) = SIN( rad * gphit(:,:) )


    ! 1. Loop on vertical levels
    ! --------------------------
    !                                             ! ===============
    DO jk = 1, jpk                                ! Horizontal slab
       !                                          ! ===============
       PRINT *,'WORKING for level ', jk, nperio
       pdct(:,:) = 0.e0
       ! read the masks
       !    temp(:,:) = getvar(cbathy,'Bathy_level',1, npiglo, npjglo)

       tmask(:,:)=getvar(cmask,'tmask',jk,jpi,jpj)
       umask(:,:)=getvar(cmask,'umask',jk,jpi,jpj)
       vmask(:,:)=getvar(cmask,'vmask',jk,jpi,jpj)
       fmask(:,:)=getvar(cmask,'fmask',jk,jpi,jpj)
       PRINT *, '    READ masks done.'
       ! Define the coastline points (U, V and F)
       DO jj = 2, jpjm1
          DO ji = 2, jpim1
             zmask(ji,jj) =  ( tmask(ji,jj+1) + tmask(ji+1,jj+1) &
                  &           + tmask(ji,jj  ) + tmask(ji+1,jj  ) )
             llcotu(ji,jj) = ( tmask(ji,jj  ) + tmask(ji+1,jj  ) == 1. ) 
             llcotv(ji,jj) = ( tmask(ji,jj  ) + tmask(ji  ,jj+1) == 1. ) 
             llcotf(ji,jj) = ( zmask(ji,jj) > 0. ) .AND. ( zmask(ji,jj) < 4. )
          END DO
       END DO
       PRINT *,'  llcot? set now.'

       ! Lateral boundaries conditions
       llcotu(:, 1 ) = umask(:,  2  ) == 1
       llcotu(:,jpj) = umask(:,jpjm1) == 1
       llcotv(:, 1 ) = vmask(:,  2  ) == 1
       llcotv(:,jpj) = vmask(:,jpjm1) == 1
       llcotf(:, 1 ) = fmask(:,  2  ) == 1
       llcotf(:,jpj) = fmask(:,jpjm1) == 1

       IF( nperio == 1 .OR. nperio == 4 .OR. nperio == 6 ) THEN
          llcotu( 1 ,:) = llcotu(jpim1,:)
          llcotu(jpi,:) = llcotu(  2  ,:)
          llcotv( 1 ,:) = llcotv(jpim1,:)
          llcotv(jpi,:) = llcotv(  2  ,:)
          llcotf( 1 ,:) = llcotf(jpim1,:)
          llcotf(jpi,:) = llcotf(  2  ,:)
       ELSE
          llcotu( 1 ,:) = umask(  2  ,:) == 1
          llcotu(jpi,:) = umask(jpim1,:) == 1
          llcotv( 1 ,:) = vmask(  2  ,:) == 1
          llcotv(jpi,:) = vmask(jpim1,:) == 1
          llcotf( 1 ,:) = fmask(  2  ,:) == 1
          llcotf(jpi,:) = fmask(jpim1,:) == 1
       ENDIF
       IF( nperio == 3 .OR. nperio == 4 ) THEN
          DO ji = 1, jpim1
             iju = jpi - ji + 1
             llcotu(ji,jpj  ) = llcotu(iju,jpj-2)
             llcotf(ji,jpjm1) = llcotf(iju,jpj-2)
             llcotf(ji,jpj  ) = llcotf(iju,jpj-3)
          END DO
          DO ji = jpi/2, jpim1
             iju = jpi - ji + 1
             llcotu(ji,jpjm1) = llcotu(iju,jpjm1)
          END DO
          DO ji = 2, jpi
             ijt = jpi - ji + 2
             llcotv(ji,jpjm1) = llcotv(ijt,jpj-2)
             llcotv(ji,jpj  ) = llcotv(ijt,jpj-3)
          END DO
       ENDIF
       IF( nperio == 5 .OR. nperio == 6 ) THEN
          DO ji = 1, jpim1
             iju = jpi - ji
             llcotu(ji,jpj  ) = llcotu(iju,jpjm1)
             llcotf(ji,jpj  ) = llcotf(iju,jpj-2)
          END DO
          DO ji = jpi/2, jpim1
             iju = jpi - ji
             llcotf(ji,jpjm1) = llcotf(iju,jpjm1)
          END DO
          DO ji = 1, jpi
             ijt = jpi - ji + 1
             llcotv(ji,jpj  ) = llcotv(ijt,jpjm1)
          END DO
          DO ji = jpi/2+1, jpi
             ijt = jpi - ji + 1
             llcotv(ji,jpjm1) = llcotv(ijt,jpjm1)
          END DO
       ENDIF

       ! Compute cartesian coordinates of coastline points
       ! and the number of coastline points

       icoast = 0
       PRINT *,' START computing cartesian coord of coastlines '
       DO jj = 1, jpj
          DO ji = 1, jpi
             IF( llcotf(ji,jj) ) THEN
                icoast = icoast + 1
                zxc(icoast) = COS( rad*gphif(ji,jj) ) * COS( rad*glamf(ji,jj) )
                zyc(icoast) = COS( rad*gphif(ji,jj) ) * SIN( rad*glamf(ji,jj) )
                zzc(icoast) = SIN( rad*gphif(ji,jj) )
             ENDIF
             IF( llcotu(ji,jj) ) THEN
                icoast = icoast+1
                zxc(icoast) = COS( rad*gphiu(ji,jj) ) * COS( rad*glamu(ji,jj) )
                zyc(icoast) = COS( rad*gphiu(ji,jj) ) * SIN( rad*glamu(ji,jj) )
                zzc(icoast) = SIN( rad*gphiu(ji,jj) )
             ENDIF
             IF( llcotv(ji,jj) ) THEN
                icoast = icoast+1
                zxc(icoast) = COS( rad*gphiv(ji,jj) ) * COS( rad*glamv(ji,jj) )
                zyc(icoast) = COS( rad*gphiv(ji,jj) ) * SIN( rad*glamv(ji,jj) )
                zzc(icoast) = SIN( rad*gphiv(ji,jj) )
             ENDIF
          END DO
       END DO
       PRINT *,' END computing cartesian coord of coastlines '

       ! Distance for the T-points

       PRINT *,' START computing distance for T points', icoast
       DO jj = 1, jpj
             print *, jj
          DO ji = 1, jpi
             IF( tmask(ji,jj) == 0. ) THEN
                pdct(ji,jj) = 0.
             ELSE
                DO jl = 1, icoast
                   zdis(jl) = ( zxt(ji,jj) - zxc(jl) )**2   &
                        &     + ( zyt(ji,jj) - zyc(jl) )**2   &
                        &     + ( zzt(ji,jj) - zzc(jl) )**2
                END DO
                pdct(ji,jj) = ra * SQRT( MINVAL( zdis(1:icoast) ) )
             ENDIF
          END DO
       END DO
       PRINT *,' END computing distance for T points'
       
       ierr=putvar(ncout,id_varout(1),pdct,jk,jpi,jpj)
       !                                                ! ===============
    END DO                                              !   End of slab
    !                                                   ! ===============
    timean(:)=0.
    ierr=putvar1d(ncout,timean,1,'T')
    ierr = closeout(ncout)

  END SUBROUTINE cofdis

  END PROGRAM cdfcofdis
