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
 INTEGER :: jpi,jpj,jpk,jpkm1
 INTEGER :: narg, iargc
 REAL(kind=4), DIMENSION(:,:), ALLOCATABLE :: glamt, glamu,glamv, glamf
 REAL(kind=4), DIMENSION(:,:), ALLOCATABLE :: gphit, gphiu,gphiv, gphif
 REAL(kind=4), DIMENSION(:,:), ALLOCATABLE :: tmask, umask, vmask, fmask
 REAL(kind=4), DIMENSION(:,:), ALLOCATABLE :: pdct
   ! from phycst
   REAL(kind=4) ::   rpi = 3.141592653589793          !: pi
   REAL(kind=4) ::   rad = 3.141592653589793 / 180.   !: conversion from degre into radian
   REAL(KIND=4) ::   ra  = 6371229.                   !: earth radius (meter)

 CHARACTER(LEN=80) :: coordhgr='mesh_hgr.nc', cmask='mask.nc'

   


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
      !!              - NetCDF file 'dist.coast.nc' 
      !!----------------------------------------------------------------------
      !!
      !!
      INTEGER ::   ji, jj, jk, jl      ! dummy loop indices
      INTEGER ::   iju, ijt            ! temporary integers
      INTEGER ::   icoast, itime
      INTEGER ::   icot         ! logical unit for file distance to the coast
      LOGICAL, DIMENSION(jpi,jpj) ::   llcotu, llcotv, llcotf   ! ???
      CHARACTER (len=32) ::   clname
      REAL(wp) ::   zdate0
      REAL(wp), DIMENSION(jpi,jpj)   ::   zxt, zyt, zzt, zmask   ! cartesian coordinates for T-points
      REAL(wp), DIMENSION(3*jpi*jpj) ::   zxc, zyc, zzc, zdis    ! temporary workspace
      !!----------------------------------------------------------------------

      ! 0. Initialization
      ! -----------------
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'cofdis : compute the distance to coastline'
      IF(lwp) WRITE(numout,*) '~~~~~~'
      IF(lwp) WRITE(numout,*)
      IF( lk_mpp ) &
           & CALL ctl_stop('         Computation not yet implemented with key_mpp_...', &
           &               '         Rerun the code on another computer or ', &
           &               '         create the "dist.coast.nc" file using IDL' )

      zxt(:,:) = cos( rad * gphit(:,:) ) * cos( rad * glamt(:,:) )
      zyt(:,:) = cos( rad * gphit(:,:) ) * sin( rad * glamt(:,:) )
      zzt(:,:) = sin( rad * gphit(:,:) )


      ! 1. Loop on vertical levels
      ! --------------------------
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         pdct(:,:) = 0.e0
         ! read the masks
!    temp(:,:) = getvar(cbathy,'Bathy_level',1, npiglo, npjglo)

         tmask(:,:)=getvar(cmask,'tmask',jk,jpi,jpj)
         umask(:,:)=getvar(cmask,'umask',jk,jpi,jpj)
         vmask(:,:)=getvar(cmask,'vmask',jk,jpi,jpj)
         fmask(:,:)=getvar(cmask,'fmask',jk,jpi,jpj)
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

         ! Distance for the T-points

         DO jj = 1, jpj
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
         !                                                ! ===============
      END DO                                              !   End of slab
      !                                                   ! ===============


      ! 2. Create the  distance to the coast file in NetCDF format
      ! ----------------------------------------------------------    
      clname = 'dist.coast'
      itime = 0
      CALL ymds2ju( 0     , 1      , 1     , 0.e0 , zdate0 )
      CALL restini( 'NONE', jpi    , jpj   , glamt, gphit ,   &
         &          jpk   , gdept_0, clname, itime, zdate0,   &
         &          rdt   , icot                         )
      CALL restput( icot, 'Tcoast', jpi, jpj, jpk, 0, pdct )
      CALL restclo( icot )

   END SUBROUTINE cofdis
