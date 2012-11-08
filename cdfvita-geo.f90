PROGRAM cdfvita_geo
  !!======================================================================
  !!                     ***  PROGRAM  cdfvita_geo  ***
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
  !!----------------------------------------------------------------------
  USE cdfio 
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2011, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                            :: ji, jj, jk, jt, jlev    ! dummy loop index
  INTEGER(KIND=4)                            :: narg, iargc, ijarg      ! browse line
  INTEGER(KIND=4)                            :: npiglo,npjglo           ! size of the domain
  INTEGER(KIND=4)                            :: npk, npt                ! size of the domain
  INTEGER(KIND=4)                            :: nlev, ik                ! number of selected levels, current lev
  INTEGER(KIND=4)                            :: ncout                   ! ncid of output file
  INTEGER(KIND=4)                            :: ierr                    ! error status for cdfio
  INTEGER(KIND=4)                            :: nvar                    ! number of variable
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: nklev                   ! selected levels
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: ipk, id_varout          ! output stuff

  REAL(KIND=4), DIMENSION(:),    ALLOCATABLE :: tim                     ! time counter array
  REAL(KIND=4), DIMENSION(:),    ALLOCATABLE :: gdeptall, gdept         ! depths and selected depths
  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: uc, vc                  ! velocity component on C grid
  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: ua, va, vmod            ! velocity component on A grid

  TYPE(variable), DIMENSION(:),  ALLOCATABLE :: stypvar                 ! data attributes

  CHARACTER(LEN=256)                         :: cf_ufil, cf_vfil        ! velocity files on C grid
  CHARACTER(LEN=256)                         :: cf_wfil                 ! optional W file on C grid
  CHARACTER(LEN=256)                         :: cf_tfil                 ! GridT file for T position
  CHARACTER(LEN=256)                         :: cf_out='vita.nc'        ! output file name
  CHARACTER(LEN=256)                         :: cldum                   ! dummy char variable

  LOGICAL                                    :: lvertical = .FALSE.     ! vertical velocity  flag 
  LOGICAL                                    :: lperio    = .FALSE.     ! E_W periodicity flag 
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfvita-geo  Ugeo-file Vgeo_file T-file [-w W-file] [-lev level_list]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Create a file with velocity components and module computed'
     PRINT *,'       at T points from file on C-grid. T-file is used only for' 
     PRINT *,'       getting the header of the output file. Any file on T grid'
     PRINT *,'       can be used.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       Ugeo-file  : netcdf file with zonal component of velocity' 
     PRINT *,'       Vgeo-file  : netcdf file with meridional component of velocity' 
     PRINT *,'       T-file  : netcdf file with T points header OK.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [ -w W-file ] : if used, also compute vertical velocities at' 
     PRINT *,'                       T points.'
     PRINT *,'       [ -lev level_list] : specify a list of level to be used '
     PRINT *,'                   (default option is to use all input levels).'
     PRINT *,'                   This option MUST be the last on the command line !!'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        none'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : sovitua, sovitva, sovitmod, [sovitwa]'
     STOP
  ENDIF

  nlev = 0
  ijarg=1
  DO WHILE ( ijarg <= narg )
     CALL getarg( ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-lev' )
        nlev= narg - ijarg + 1
        ALLOCATE (nklev(nlev) )
        DO jlev = 1, nlev
           CALL getarg( ijarg, cldum ) ; ijarg=ijarg+1 ; READ(cldum,* ) nklev(jlev)
        ENDDO
     CASE ( '-w' )
        CALL getarg( ijarg, cf_wfil ) ; ijarg=ijarg+1
        lvertical=.TRUE.
     CASE DEFAULT
        cf_ufil=cldum
        CALL getarg( ijarg, cf_vfil ) ; ijarg=ijarg+1
        CALL getarg( ijarg, cf_tfil ) ; ijarg=ijarg+1
     END SELECT
  ENDDO

  ! adjust number of variable according to -w option
  nvar=3
  IF ( lvertical )  nvar = 4

  ALLOCATE ( ipk(nvar), id_varout(nvar), stypvar(nvar) )

  IF ( chkfile(cf_ufil) .OR. chkfile(cf_vfil) .OR. chkfile(cf_tfil) ) STOP ! missing file

  IF ( lvertical ) THEN 
     IF ( chkfile(cf_wfil) ) STOP ! missing file
  ENDIF

  npiglo = getdim (cf_ufil,cn_x)
  npjglo = getdim (cf_ufil,cn_y)
  npk    = getdim (cf_ufil,cn_z)
  npt    = getdim (cf_ufil,cn_t)

  IF ( npk == 0 ) THEN ; npk = 1 ; ENDIF

  IF ( nlev == 0 ) THEN ! take all levels
     nlev = npk
     ALLOCATE (nklev(nlev) )
     DO jlev = 1, nlev
        nklev(jlev) = jlev
     ENDDO
  ENDIF

  ALLOCATE ( gdept(nlev) )

  ! Zonal Velocity T point
  ipk(1)                       = nlev
  stypvar(1)%cname             = 'sovitua'
  stypvar(1)%cunits            = 'm/s'
  stypvar(1)%rmissing_value    = 0.
  stypvar(1)%valid_min         = 0.
  stypvar(1)%valid_max         = 10000.
  stypvar(1)%clong_name        = 'Zonal Velocity T point'
  stypvar(1)%cshort_name       = 'sovitua'
  stypvar(1)%conline_operation = 'N/A'
  stypvar(1)%caxis             = 'TZYX'

  ! Meridional Velocity T point
  ipk(2)                       = nlev
  stypvar(2)%cname             = 'sovitva'
  stypvar(2)%cunits            = 'm/s'
  stypvar(2)%rmissing_value    = 0.
  stypvar(2)%valid_min         = 0.
  stypvar(2)%valid_max         = 10000.
  stypvar(2)%clong_name        = 'Meridional Velocity T point'
  stypvar(2)%cshort_name       = 'sovitva'
  stypvar(2)%conline_operation = 'N/A'
  stypvar(2)%caxis             = 'TZYX'

  ! Velocity module T point
  ipk(3)                       = nlev
  stypvar(3)%cname             = 'sovitmod'
  stypvar(3)%cunits            = 'm/s'
  stypvar(3)%rmissing_value    = 0.
  stypvar(3)%valid_min         = 0.
  stypvar(3)%valid_max         = 10000.
  stypvar(3)%clong_name        = 'Velocity module T point'
  stypvar(3)%cshort_name       = 'sovitmod'
  stypvar(3)%conline_operation = 'N/A'
  stypvar(3)%caxis             = 'TZYX'

  IF ( lvertical ) THEN
     ! Vertical Velocity at T point
     ipk(nvar)                       = nlev
     stypvar(nvar)%cname             = 'sovitwa'
     stypvar(nvar)%cunits            = 'mm/s'
     stypvar(nvar)%rmissing_value    = 0.
     stypvar(nvar)%valid_min         = 0.
     stypvar(nvar)%valid_max         = 10000.
     stypvar(nvar)%clong_name        = 'Vertical Velocity at T point'
     stypvar(nvar)%cshort_name       = 'sovitwa'
     stypvar(nvar)%conline_operation = 'N/A'
     stypvar(nvar)%caxis             = 'TZYX'
  ENDIF

  PRINT *, 'npiglo =', npiglo
  PRINT *, 'npjglo =', npjglo
  PRINT *, 'npk    =', npk
  PRINT *, 'npt    =', npt
  PRINT *, 'nlev   =', nlev

  ALLOCATE( uc(npiglo,npjglo), vc(npiglo,npjglo)  )
  ALLOCATE( ua(npiglo,npjglo), va(npiglo,npjglo), vmod(npiglo,npjglo) )
  ALLOCATE( tim(npt), gdeptall(npk) )

  gdeptall(:) = getvar1d(cf_tfil,cn_vdeptht, npk)
  DO jlev = 1, nlev
     ik = nklev(jlev)
     gdept(jlev) = gdeptall(ik)
  ENDDO

  ! check E-W periodicity using uc array as working space
  uc(:,:) = getvar(cf_tfil, cn_vlon2d, 1, npiglo, npjglo )
  IF ( uc(1,1) == uc(npiglo-1,1) )  THEN 
    lperio = .TRUE.
    PRINT *,' E-W periodicity detected.'
  ENDIF

  ncout = create      (cf_out,   cf_tfil,  npiglo, npjglo, nlev                 )
  ierr  = createvar   (ncout ,   stypvar,  nvar,   ipk,    id_varout            )
  ierr  = putheadervar(ncout,    cf_tfil,  npiglo, npjglo, nlev,     pdep=gdept )

  DO jt = 1, npt
     DO jlev = 1, nlev
        ik = nklev(jlev)
        uc(:,:) = getvar(cf_ufil, cn_vozocrtx, ik ,npiglo, npjglo, ktime=jt )
        vc(:,:) = getvar(cf_vfil, cn_vomecrty, ik ,npiglo, npjglo, ktime=jt )

        ua = 0. ; va = 0. ; ua(:,:) = 0. ; va(:,:)=0. ; vmod(:,:)=0.
        DO ji=2, npiglo
           DO jj=2,npjglo
              ua(ji,jj)   = 0.5* (uc(ji,jj  )+ uc(ji,jj-1))
              va(ji,jj)   = 0.5* (vc(ji-1,jj)+ vc(ji,jj  ))
              vmod(ji,jj) = SQRT( ua(ji,jj)*ua(ji,jj) + va(ji,jj)*va(ji,jj) )
           END DO
        END DO
        IF ( lperio) THEN  ! periodic E-W boundary ...
          ua  (1,:) = ua  (npiglo-1,:)
          va  (1,:) = va  (npiglo-1,:)
          vmod(1,:) = vmod(npiglo-1,:)
        ENDIF

        ierr=putvar(ncout, id_varout(1), ua,   jlev ,npiglo, npjglo, ktime=jt )
        ierr=putvar(ncout, id_varout(2), va,   jlev ,npiglo, npjglo, ktime=jt )
        ierr=putvar(ncout, id_varout(3), vmod, jlev ,npiglo, npjglo, ktime=jt )
     END DO
  END DO

  IF ( lvertical ) THEN
     ! reuse uc an vc arrays to store Wk and Wk+1
     DO jt = 1, npt
        DO jlev=1, nlev - 1
           uc(:,:) = getvar(cf_wfil, cn_vovecrtz, nklev(jlev),   npiglo, npjglo, ktime=jt )
           vc(:,:) = getvar(cf_wfil, cn_vovecrtz, nklev(jlev)+1, npiglo, npjglo, ktime=jt )
           ua(:,:) = 0.5*(uc(:,:) + vc(:,:))*1000.  ! mm/sec
           ierr    = putvar(ncout, id_varout(4), ua, jlev,      npiglo, npjglo, ktime=jt )
           uc(:,:) = vc(:,:)
        END DO
        IF ( nlev == npk ) THEN
           ua(:,:) = 0.e0  ! npk
        ELSE
           uc(:,:) = getvar(cf_wfil, cn_vovecrtz, nklev(nlev),   npiglo, npjglo, ktime=jt )
           vc(:,:) = getvar(cf_wfil, cn_vovecrtz, nklev(nlev)+1, npiglo, npjglo, ktime=jt )
           ua(:,:) = 0.5*(uc(:,:) + vc(:,:))*1000.  ! mm/sec
     ENDIF
        ierr = putvar(ncout, id_varout(4), ua, nlev ,npiglo, npjglo, ktime=jt )
     ENDDO
  ENDIF

  tim  = getvar1d(cf_ufil, cn_vtimec, npt     )
  ierr = putvar1d(ncout,  tim,       npt, 'T')
  ierr = closeout(ncout)

END PROGRAM cdfvita_geo
