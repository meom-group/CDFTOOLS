PROGRAM cdfvita
  !!======================================================================
  !!                     ***  PROGRAM  cdfvita  ***
  !!=====================================================================
  !!  ** Purpose : Compute velocity on t grid
  !!
  !!  ** Method  : Read velocity component on input gridU and gridV file
  !!               Use gridT file for the proper location of T points
  !!               The velocity module is also output (same function than
  !!               cdfspeed) If a gridW file is given, (fifth argument)
  !!               then w is also computed on the T grid
  !!
  !! History : 2.1  : 11/2006  : J.M. Molines : Original code
  !!           3.0  : 01/2011  : J.M. Molines : Doctor norm + Lic.
  !!                : 03/2013  : J.M. Molines : add -geo option
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  USE cdfio 
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class derived_fields
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                            :: ji, jj, jk, jt, jlev    ! dummy loop index
  INTEGER(KIND=4)                            :: narg, iargc, ijarg      ! browse line
  INTEGER(KIND=4)                            :: npiglo,npjglo           ! size of the domain
  INTEGER(KIND=4)                            :: npk, npt                ! size of the domain
  INTEGER(KIND=4)                            :: nlev, ik                ! number of selected levels, current lev
  INTEGER(KIND=4)                            :: ncout                   ! ncid of output file
  INTEGER(KIND=4)                            :: ierr                    ! error status for cdfio
  INTEGER(KIND=4)                            :: ivar                    ! variable index
  INTEGER(KIND=4)                            :: ivar_vert               ! variable index of vertical velocity
  INTEGER(KIND=4)                            :: ivar_cub                ! variable index of cube of velocity module
  INTEGER(KIND=4)                            :: nvar                    ! number of variable
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: nklevel                 ! selected levels
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: ipk, id_varout          ! output stuff

  REAL(KIND=4)                               :: pi                      ! pi
  REAL(KIND=4), DIMENSION(:),    ALLOCATABLE :: gdeptall, gdept         ! depths and selected depths
  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: uc, vc                  ! velocity component on C grid
  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: ua, va, vmod, vdir      ! velocity component on A grid

  REAL(KIND=8), DIMENSION(:),    ALLOCATABLE :: dtim                    ! time counter array

  TYPE(variable), DIMENSION(:),  ALLOCATABLE :: stypvar                 ! data attributes

  CHARACTER(LEN=256)                         :: cf_ufil, cf_vfil        ! velocity files on C grid
  CHARACTER(LEN=256)                         :: cf_wfil                 ! optional W file on C grid
  CHARACTER(LEN=256)                         :: cf_tfil                 ! GridT file for T position
  CHARACTER(LEN=256)                         :: cf_out='vita.nc'        ! output file name
  CHARACTER(LEN=256)                         :: cldum                   ! dummy char variable

  LOGICAL                                    :: lvertical = .FALSE.     ! vertical velocity  flag 
  LOGICAL                                    :: lperio    = .FALSE.     ! E_W periodicity flag 
  LOGICAL                                    :: lgeo      = .FALSE.     ! input U V files are geostrophic files
  LOGICAL                                    :: lcub      = .FALSE.     ! save U*U*U on  A grid
  LOGICAL                                    :: lnc4      = .FALSE.     ! Use nc4 with chunking and deflation
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfvita -u U-file -v V-file -t T-file [-w W-file] [-geo] [-cubic]'
     PRINT *,'               ... [-o OUT-file] [-nc4] [-lev LST-level]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Creates a file with velocity components, module  and direction'
     PRINT *,'       at T points from file on C-grid. T-file is used only for' 
     PRINT *,'       getting the header of the output file. Any file on T grid'
     PRINT *,'       can be used.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -u U-file  : netcdf file with zonal component of velocity' 
     PRINT *,'       -v V-file  : netcdf file with meridional component of velocity' 
     PRINT *,'       -t T-file  : netcdf file with T points header OK.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-w W-file ]: if used, also compute vertical velocities at T points.' 
     PRINT *,'       [-geo ]     : indicate that input velocity files are produced by '
     PRINT *,'              cdfgeo-uv, hence ugeo on V-point, vgeo on U-points. '
     PRINT *,'       [-cubic ]   : Save the cube of the velocity module.'
     PRINT *,'       [-nc4 ]     : Use netcdf4 output with chunking and deflation level 1.'
     PRINT *,'              This option is effective only if cdftools are compiled with'
     PRINT *,'              a netcdf library supporting chunking and deflation.'
     PRINT *,'       [-o OUT-file ] : Specify name of output file instead of ',TRIM(cf_out)
     PRINT *,'       [-lev LST-level] : specify a blank-separated list of levels to be used.'
     PRINT *,'              (default option is to use all input levels).'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        none'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) ,' unless -o option is used'
     PRINT *,'         variables : sovitua, sovitva, sovitmod, sovitdir, [sovitmod3] and'
     PRINT *,'                     [sovitwa]. Note that the direction is relative to the'
     PRINT *,'                     model grid.'
     PRINT *,'      '
     PRINT *,'     SEE ALSO : '
     PRINT *,'        cdfspeed only computes the speed (velocity module).'
     PRINT *,'      '
     STOP 
  ENDIF

  nlev = 0
  ijarg=1
  pi=ACOS(-1.)
  DO WHILE ( ijarg <= narg )
     CALL getarg( ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-u'   ) ; CALL getarg(ijarg, cf_ufil ) ; ijarg=ijarg+1
     CASE ( '-v'   ) ; CALL getarg(ijarg, cf_vfil ) ; ijarg=ijarg+1
        ! options
     CASE ( '-lev' ) ; CALL GetLevList
     CASE ( '-t'   ) ; CALL getarg(ijarg, cf_tfil ) ; ijarg=ijarg+1
     CASE ( '-w'   ) ; CALL getarg(ijarg, cf_wfil ) ; ijarg=ijarg+1
        ;              lvertical=.TRUE.
     CASE ( '-geo' ) ; lgeo = .TRUE.
     CASE ('-cubic') ; lcub = .TRUE.
     CASE ( '-nc4' ) ; lnc4 = .TRUE.
     CASE ( '-o'   ) ; CALL getarg(ijarg, cf_out  ) ; ijarg=ijarg+1
     CASE DEFAULT    ; PRINT *,' ERROR : ', TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  ENDDO

  IF ( chkfile(cf_ufil) .OR. chkfile(cf_vfil) .OR. chkfile(cf_tfil) ) STOP 99 ! missing file
  IF ( lvertical ) THEN 
     IF ( chkfile(cf_wfil) ) STOP 99 ! missing file
  ENDIF

  npiglo = getdim (cf_ufil,cn_x)
  npjglo = getdim (cf_ufil,cn_y)
  npk    = getdim (cf_ufil,cn_z)
  npt    = getdim (cf_ufil,cn_t)

  IF (npk == 0 ) THEN
     npk  = 1
  ENDIF

  IF ( nlev == 0 ) THEN ! take all levels
     nlev = npk
     ALLOCATE (nklevel(nlev) )
     DO jlev = 1, nlev
        nklevel(jlev) = jlev
     ENDDO
  ENDIF

  ! adjust number of variable according to  option
  nvar=4
  IF ( lvertical )  nvar = nvar + 1
  IF ( lcub      )  nvar = nvar + 1

  ALLOCATE ( ipk(nvar), id_varout(nvar), stypvar(nvar) )
  ALLOCATE ( gdept(nlev) )

  PRINT *, 'npiglo =', npiglo
  PRINT *, 'npjglo =', npjglo
  PRINT *, 'npk    =', npk
  PRINT *, 'npt    =', npt
  PRINT *, 'nlev   =', nlev

  ALLOCATE( uc(npiglo,npjglo), vc(npiglo,npjglo)  )
  ALLOCATE( ua(npiglo,npjglo), va(npiglo,npjglo)  )
  ALLOCATE( vmod(npiglo,npjglo), vdir(npiglo, npjglo)  )
  ALLOCATE( dtim(npt), gdeptall(npk) )

  gdeptall(:) = getvar1d(cf_tfil,cn_vdeptht, npk)
  DO jlev = 1, nlev
     ik = nklevel(jlev)
     gdept(jlev) = gdeptall(ik)
  ENDDO

  ! check E-W periodicity using uc array as working space
  uc(:,:) = getvar(cf_tfil, cn_vlon2d, 1, npiglo, npjglo )
  IF ( uc(1,1) == uc(npiglo-1,1) )  THEN 
     lperio = .TRUE.
     PRINT *,' E-W periodicity detected.'
  ENDIF

  CALL CreateOutput

  DO jt = 1, npt
     DO jlev = 1, nlev
        ik = nklevel(jlev)
        uc(:,:) = getvar(cf_ufil, cn_vozocrtx, ik ,npiglo, npjglo, ktime=jt )
        vc(:,:) = getvar(cf_vfil, cn_vomecrty, ik ,npiglo, npjglo, ktime=jt )

        ua = 0. ; va = 0. ; ua(:,:) = 0. ; va(:,:)=0. ; vmod(:,:)=0.
        IF ( lgeo ) THEN  ! geostrophic velocities
           DO ji=2, npiglo
              DO jj=2,npjglo
                 ua(ji,jj)   = 0.5* (uc(ji,jj)+ uc(ji  ,jj-1))
                 va(ji,jj)   = 0.5* (vc(ji,jj)+ vc(ji-1,jj  ))
                 vmod(ji,jj) = SQRT( ua(ji,jj)*ua(ji,jj) + va(ji,jj)*va(ji,jj) )
                 ! the direction is relative to the model mesh, not the true North
                 vdir(ji,jj) = 90. - ATAN2(va(ji,jj),ua(ji,jj))*180./pi
                 IF ( vdir(ji,jj) < 0. ) vdir(ji,jj) = 360.+vdir(ji,jj)
              END DO
           END DO
        ELSE
           DO ji=2, npiglo
              DO jj=2,npjglo
                 ua(ji,jj)   = 0.5* (uc(ji,jj)+ uc(ji-1,jj))
                 va(ji,jj)   = 0.5* (vc(ji,jj)+ vc(ji,jj-1))
                 vmod(ji,jj) = SQRT( ua(ji,jj)*ua(ji,jj) + va(ji,jj)*va(ji,jj) )
                 ! the direction is relative to the model mesh, not the true North
                 vdir(ji,jj) = 90. - ATAN2(va(ji,jj),ua(ji,jj))*180./pi
                 IF ( vdir(ji,jj) < 0. ) vdir(ji,jj) = 360.+vdir(ji,jj)
              END DO
           END DO
        ENDIF
        IF ( lperio) THEN  ! periodic E-W boundary ...
           ua  (1,:) = ua  (npiglo-1,:)
           va  (1,:) = va  (npiglo-1,:)
           vmod(1,:) = vmod(npiglo-1,:)
           vdir(1,:) = vdir(npiglo-1,:)
        ENDIF
        ivar = 1
        ierr=putvar(ncout, id_varout(ivar), ua,   jlev ,npiglo, npjglo, ktime=jt ) ; ivar = ivar +1
        ierr=putvar(ncout, id_varout(ivar), va,   jlev ,npiglo, npjglo, ktime=jt ) ; ivar = ivar +1
        ierr=putvar(ncout, id_varout(ivar), vmod, jlev ,npiglo, npjglo, ktime=jt ) ; ivar = ivar +1
        ierr=putvar(ncout, id_varout(ivar), vdir, jlev ,npiglo, npjglo, ktime=jt ) 
        IF ( lcub ) THEN
           ierr=putvar(ncout, id_varout(ivar_cub), vmod*vmod*vmod, jlev ,npiglo, npjglo, ktime=jt ) 
        ENDIF
     END DO
  END DO

  IF ( lvertical ) THEN
     ! reuse uc an vc arrays to store Wk and Wk+1
     DO jt = 1, npt
        DO jlev=1, nlev - 1
           uc(:,:) = getvar(cf_wfil, cn_vovecrtz, nklevel(jlev),   npiglo, npjglo, ktime=jt )
           vc(:,:) = getvar(cf_wfil, cn_vovecrtz, nklevel(jlev)+1, npiglo, npjglo, ktime=jt )
           ua(:,:) = 0.5*(uc(:,:) + vc(:,:))*1000.  ! mm/sec
           ierr    = putvar(ncout, id_varout(ivar_vert), ua, jlev,      npiglo, npjglo, ktime=jt )
           uc(:,:) = vc(:,:)
        END DO
        IF ( nlev == npk ) THEN
           ua(:,:) = 0.e0  ! npk
        ELSE
           uc(:,:) = getvar(cf_wfil, cn_vovecrtz, nklevel(nlev),   npiglo, npjglo, ktime=jt )
           vc(:,:) = getvar(cf_wfil, cn_vovecrtz, nklevel(nlev)+1, npiglo, npjglo, ktime=jt )
           ua(:,:) = 0.5*(uc(:,:) + vc(:,:))*1000.  ! mm/sec
        ENDIF
        ierr = putvar(ncout, id_varout(ivar_vert), ua, nlev ,npiglo, npjglo, ktime=jt )
     ENDDO
  ENDIF

  ierr = closeout(ncout)

CONTAINS

  SUBROUTINE CreateOutput
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutput  ***
    !!
    !! ** Purpose :  Create netcdf output file(s) 
    !!
    !! ** Method  :  Use stypvar global description of variables
    !!
    !!----------------------------------------------------------------------
    ivar=0
    ! Zonal Velocity T point
    ivar                            = ivar+1
    ipk(ivar)                       = nlev
    stypvar(ivar)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(ivar)%cname             = 'sovitua'
    stypvar(ivar)%cunits            = 'm/s'
    stypvar(ivar)%rmissing_value    = 0.
    stypvar(ivar)%valid_min         = 0.
    stypvar(ivar)%valid_max         = 10000.
    stypvar(ivar)%clong_name        = 'Zonal Velocity T point'
    stypvar(ivar)%cshort_name       = 'sovitua'
    stypvar(ivar)%conline_operation = 'N/A'
    stypvar(ivar)%caxis             = 'TZYX'

    ! Meridional Velocity T point
    ivar                            = ivar+1
    ipk(ivar)                       = nlev
    stypvar(ivar)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(ivar)%cname             = 'sovitva'
    stypvar(ivar)%cunits            = 'm/s'
    stypvar(ivar)%rmissing_value    = 0.
    stypvar(ivar)%valid_min         = 0.
    stypvar(ivar)%valid_max         = 10000.
    stypvar(ivar)%clong_name        = 'Meridional Velocity T point'
    stypvar(ivar)%cshort_name       = 'sovitva'
    stypvar(ivar)%conline_operation = 'N/A'
    stypvar(ivar)%caxis             = 'TZYX'

    ! Velocity module T point
    ivar                            = ivar+1
    ipk(ivar)                       = nlev
    stypvar(ivar)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(ivar)%cname             = 'sovitmod'
    stypvar(ivar)%cunits            = 'm/s'
    stypvar(ivar)%rmissing_value    = 0.
    stypvar(ivar)%valid_min         = 0.
    stypvar(ivar)%valid_max         = 10000.
    stypvar(ivar)%clong_name        = 'Velocity module T point'
    stypvar(ivar)%cshort_name       = 'sovitmod'
    stypvar(ivar)%conline_operation = 'N/A'
    stypvar(ivar)%caxis             = 'TZYX'

    ! Velocity direction  T point
    ivar                            = ivar+1
    ipk(ivar)                       = nlev
    stypvar(ivar)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(ivar)%cname             = 'sovitdir'
    stypvar(ivar)%cunits            = 'deg N'
    stypvar(ivar)%rmissing_value    = 0.
    stypvar(ivar)%valid_min         = 0.
    stypvar(ivar)%valid_max         = 360.
    stypvar(ivar)%clong_name        = 'Velocity direction T point'
    stypvar(ivar)%cshort_name       = 'sovitdir'
    stypvar(ivar)%conline_operation = 'N/A'
    stypvar(ivar)%caxis             = 'TZYX'

    IF ( lvertical ) THEN
       ! Vertical Velocity at T point
       ivar                            = ivar+1
       ivar_vert                       = ivar
       ipk(ivar)                       = nlev
       stypvar(ivar)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
       stypvar(ivar)%cname             = 'sovitwa'
       stypvar(ivar)%cunits            = 'mm/s'
       stypvar(ivar)%rmissing_value    = 0.
       stypvar(ivar)%valid_min         = 0.
       stypvar(ivar)%valid_max         = 10000.
       stypvar(ivar)%clong_name        = 'Vertical Velocity at T point'
       stypvar(ivar)%cshort_name       = 'sovitwa'
       stypvar(ivar)%conline_operation = 'N/A'
       stypvar(ivar)%caxis             = 'TZYX'
    ENDIF

    IF ( lcub ) THEN
       ! Cube of velocity module
       ivar                            = ivar+1
       ivar_cub                        = ivar
       ipk(ivar)                       = nlev
       stypvar(ivar)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
       stypvar(ivar)%cname             = 'sovitmod3'
       stypvar(ivar)%cunits            = 'm3/s3'
       stypvar(ivar)%rmissing_value    = 0.
       stypvar(ivar)%valid_min         = 0.
       stypvar(ivar)%valid_max         = 10000.
       stypvar(ivar)%clong_name        = 'cube of velocity module'
       stypvar(ivar)%cshort_name       = 'sovitmod3'
       stypvar(ivar)%conline_operation = 'N/A'
       stypvar(ivar)%caxis             = 'TZYX'
    ENDIF

    ncout = create      (cf_out,   cf_tfil,  npiglo, npjglo, nlev     , ld_nc4=lnc4 )
    ierr  = createvar   (ncout ,   stypvar,  nvar,   ipk,    id_varout, ld_nc4=lnc4 )
    ierr  = putheadervar(ncout,    cf_tfil,  npiglo, npjglo, nlev,     pdep=gdept )

    dtim = getvar1d(cf_ufil, cn_vtimec, npt     )
    ierr = putvar1d(ncout,  dtim,       npt, 'T')

  END SUBROUTINE CreateOutput

  SUBROUTINE GetLevList
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE GetLevList  ***
    !!
    !! ** Purpose :  Set up a level list given on the command line as 
    !!               blank separated list
    !!
    !! ** Method  :  Scan the command line until a '-' is found
    !!----------------------------------------------------------------------
    INTEGER (KIND=4)  :: ji
    INTEGER (KIND=4)  :: icur
    !!----------------------------------------------------------------------
    !!
    nlev=0
    ! need to read a list of level ( number unknow ) 
    ! loop on argument till a '-' is found as first char
    icur=ijarg                          ! save current position of argument number
    DO ji = icur, narg                  ! scan arguments till - found
       CALL getarg ( ji, cldum )
       IF ( cldum(1:1) /= '-' ) THEN ; nlev = nlev+1
       ELSE                          ; EXIT
       ENDIF
    ENDDO
    ALLOCATE (nklevel(nlev) )
    DO ji = icur, icur + nlev -1
       CALL getarg(ji, cldum ) ; ijarg=ijarg+1 ; READ(cldum, * ) nklevel( ji -icur +1 )
    END DO
  END SUBROUTINE GetLevList

END PROGRAM cdfvita
