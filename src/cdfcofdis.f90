PROGRAM cdfcofdis
  !!======================================================================
  !!                     ***  PROGRAM  cdfcofdis  ***
  !!=====================================================================
  !!  ** Purpose : A wrapper for NEMO routine cofdis: create a file 
  !!               with the distance to coast variable
  !!
  !!  ** Method  : Mimic some NEMO global variables to be able to use
  !!               NEMO cofdis with minimum changes. Use cdfio instead
  !!               of IOIPSL for the output file. Due to this constaint
  !!               DOCTOR norm is not fully respected (eg jpi not PARAMETER) 
  !!               pdct is not a routine argument ...
  !!
  !! History : 2.1  : 11/2009  : J.M. Molines : Original code
  !!           3.0  : 12/2010  : J.M. Molines : Doctor norm + Lic.
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   cofdis       : compute distance to coast (NEMO routine )
  !!----------------------------------------------------------------------

  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2010, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                           :: jpi, jpj, jpk, npk
  INTEGER(KIND=4)                           :: jpim1, jpjm1, nperio=4
  INTEGER(KIND=4)                           :: narg, iargc, ijarg
  INTEGER(KIND=4)                           :: ncout, ierr
  INTEGER(KIND=4), DIMENSION(1)             :: ipk, id_varout

  ! from phycst
  REAL(KIND=4)            :: rpi = 3.141592653589793          !: pi
  REAL(KIND=4)            :: rad = 3.141592653589793 / 180.   !: conv. from degre into radian
  REAL(KIND=4)            :: ra  = 6371229.                   !: earth radius (meter)

  REAL(KIND=4) ,DIMENSION(1)                :: timean
  ! to be read in mesh_hgr
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: glamt, glamu,glamv, glamf
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: gphit, gphiu,gphiv, gphif
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: tmask, umask, vmask, fmask
  ! 
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: pdct                       ! 2D only in this version
  ! It is a 3D arg in original cofdis
  CHARACTER(LEN=256)                        :: cf_out='dist.coast'
  CHARACTER(LEN=256)                        :: cf_tfil
  CHARACTER(LEN=256)                        :: cv_out='Tcoast'
  CHARACTER(LEN=256)                        :: cldum

  TYPE(variable), DIMENSION(1)              :: stypvar

  LOGICAL                                   :: lchk
  LOGICAL                                   :: lsurf = .FALSE.
  LOGICAL                                   :: lnc4  = .FALSE.     ! Use nc4 with chunking and deflation
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()
  !
  narg=iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage :  cdfcofdis -H HGR-file -M MSK-file -T gridT.nc [-jperio jperio ]...'
     PRINT *,'               ... [-surf] [-o OUT-file[ [-nc4] '
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'        Compute the distance to the coast and create a file with the ',TRIM(cv_out)
     PRINT *,'        variable, indicating the distance to the coast. This computation is '
     PRINT *,'        done for every model level, unless -surf option is used.'
     PRINT *,'        This file is used in NEMO tradmp routine for fading out restoring'
     PRINT *,'        in vicinity of the coast line.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -H HGR-file : name of the mesh_hgr file '
     PRINT *,'       -M MSK-file : name of the mask file '
     PRINT *,'       -T T-file   : netcdf file at T point ( used for looking at jpk)'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-jperio jperio ] : define the NEMO jperio variable for north fold '
     PRINT *,'           condition. Default is  4.'
     PRINT *,'       [-surf ] : only compute  distance at the surface.'
     PRINT *,'       [-o OUT-file ] : specify name of the output file instead of ', TRIM(cf_out)
     PRINT *,'       [-nc4 ]     : Use netcdf4 output with chunking and deflation level 1'
     PRINT *,'                 This option is effective only if cdftools are compiled with'
     PRINT *,'                 a netcdf library supporting chunking and deflation.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       none' 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) ,' unless -o option is used.'
     PRINT *,'         variables : ', TRIM(cv_out),' (m)'
     PRINT *,'      '
     PRINT *,'      '
     STOP
  ENDIF

  ijarg = 1
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-H'     ) ; CALL getarg(ijarg,cn_fhgr) ; ijarg=ijarg+1
     CASE ( '-M'     ) ; CALL getarg(ijarg,cn_fmsk) ; ijarg=ijarg+1
     CASE ( '-T'     ) ; CALL getarg(ijarg,cf_tfil) ; ijarg=ijarg+1
        ! options
     CASE ( '-jperio') ; CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1 ;  READ(cldum, * ) nperio
     CASE ( '-surf'  ) ; lsurf = .TRUE.
     CASE ( '-o'     ) ; CALL getarg(ijarg, cf_out) ; ijarg=ijarg+1 
     CASE ( '-nc4'   ) ; lnc4  = .TRUE.
     CASE DEFAULT      ; PRINT *,' ERROR : ', TRIM(cldum),' : unknown option.' ; STOP
     END SELECT
  ENDDO

  lchk =           chkfile ( cn_fhgr )
  lchk = lchk .OR. chkfile ( cn_fmsk )
  lchk = lchk .OR. chkfile ( cf_tfil )
  IF ( lchk ) STOP ! missing files

  ! read domain dimensions in the mask file
  jpi = getdim(cf_tfil,cn_x)
  jpj = getdim(cf_tfil,cn_y)
  jpk = getdim(cf_tfil,cn_z)

  IF (jpk == 0 ) THEN
     jpk = getdim(cf_tfil,'z')
     IF ( jpk == 0 ) THEN
        PRINT *,' ERROR in determining jpk form gridT file ....'
        STOP
     ENDIF
  ENDIF

  PRINT *, ' JPI = ', jpi
  PRINT *, ' JPJ = ', jpj
  PRINT *, ' JPK = ', jpk

  jpim1=jpi-1 ; jpjm1=jpj-1

  ! ALLOCATION of the arrays
  ALLOCATE ( glamt(jpi,jpj), glamu(jpi,jpj), glamv(jpi,jpj), glamf(jpi,jpj) )
  ALLOCATE ( gphit(jpi,jpj), gphiu(jpi,jpj), gphiv(jpi,jpj), gphif(jpi,jpj) )
  ALLOCATE ( tmask(jpi,jpj), umask(jpi,jpj), vmask(jpi,jpj), fmask(jpi,jpj) )
  ALLOCATE ( pdct(jpi,jpj) )

  PRINT *, 'ALLOCATION DONE.'

  ! read latitude an longitude
  glamt(:,:) = getvar(cn_fhgr,cn_glamt,1,jpi,jpj)
  glamu(:,:) = getvar(cn_fhgr,cn_glamu,1,jpi,jpj)
  glamv(:,:) = getvar(cn_fhgr,cn_glamv,1,jpi,jpj)
  glamf(:,:) = getvar(cn_fhgr,cn_glamf,1,jpi,jpj)

  gphit(:,:) = getvar(cn_fhgr,cn_gphit,1,jpi,jpj)
  gphiu(:,:) = getvar(cn_fhgr,cn_gphiu,1,jpi,jpj)
  gphiv(:,:) = getvar(cn_fhgr,cn_gphiv,1,jpi,jpj)
  gphif(:,:) = getvar(cn_fhgr,cn_gphif,1,jpi,jpj)

  ! prepare file output
  IF ( lsurf ) THEN ; npk                       = 1
  ELSE              ; npk = jpk
  ENDIF
  
  CALL CreateOutput

  CALL cofdis (npk)

CONTAINS

  SUBROUTINE cofdis(kpk)
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
    INTEGER(KIND=4), INTENT(in) ::   kpk                 ! number of level to deal with

    INTEGER(KIND=4) ::   ji, jj, jk, jl      ! dummy loop indices
    INTEGER(KIND=4) ::   iju, ijt            ! temporary integers
    INTEGER(KIND=4) ::   icoast, itime
    INTEGER(KIND=4) ::   icot         ! logical unit for file distance to the coast
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
    DO jk = 1, kpk                                ! Horizontal slab
       !                                          ! ===============
       PRINT *,'WORKING for level ', jk, nperio
       pdct(:,:) = 0.e0
       tmask(:,:)=getvar(cn_fmsk,cn_tmask,jk,jpi,jpj)
       DO jj = 1, jpjm1
          DO ji = 1, jpim1   ! vector loop
             umask(ji,jj) = tmask(ji,jj ) * tmask(ji+1,jj  )
             vmask(ji,jj) = tmask(ji,jj ) * tmask(ji  ,jj+1)
          END DO
          DO ji = 1, jpim1      ! NO vector opt.
             fmask(ji,jj) = tmask(ji,jj  ) * tmask(ji+1,jj  )   &
                  &         * tmask(ji,jj+1) * tmask(ji+1,jj+1)
          END DO
       END DO
       umask(jpi,:)=umask(2,:)
       vmask(jpi,:)=vmask(2,:)
       fmask(jpi,:)=fmask(2,:)


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
          PRINT *, jj
          DO ji = 1, jpi
             IF( tmask(ji,jj) == 0. ) THEN
                pdct(ji,jj) = 0.
             ELSE
                DO jl = 1, icoast
                   zdis(jl) = ( zxt(ji,jj) - zxc(jl) )**2   &
                        &     + ( zyt(ji,jj) - zyc(jl) )**2 &
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

  SUBROUTINE CreateOutput
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutput  ***
    !!
    !! ** Purpose :  Create netcdf output file(s) 
    !!
    !! ** Method  :  Use stypvar global description of variables
    !!
    !!----------------------------------------------------------------------
    ipk(1)                       = npk
    stypvar(1)%ichunk            = (/jpi,MAX(1,jpj/30),1,1 /)
    stypvar(1)%cname             = cv_out
    stypvar(1)%cunits            = 'm'
    stypvar(1)%rmissing_value    = 0
    stypvar(1)%valid_min         = 0.
    stypvar(1)%valid_max         = 1.
    stypvar(1)%clong_name        = cv_out
    stypvar(1)%cshort_name       = cv_out
    stypvar(1)%conline_operation = 'N/A'
    stypvar(1)%caxis             = 'TZYX'
    stypvar(1)%cprecision        = 'r4'

    ncout = create      (cf_out, cf_tfil, jpi, jpj, npk       , ld_nc4=lnc4 )
    ierr  = createvar   (ncout,  stypvar, 1,   ipk, id_varout , ld_nc4=lnc4 )
    ierr  = putheadervar(ncout,  cf_tfil, jpi, jpj, npk       )

  END SUBROUTINE CreateOutput

END PROGRAM cdfcofdis
