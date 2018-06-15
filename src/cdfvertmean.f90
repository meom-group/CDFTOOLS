PROGRAM cdfvertmean
  !!======================================================================
  !!                     ***  PROGRAM  cdfvertmean  ***
  !!=====================================================================
  !!  ** Purpose : Compute the vertical average of a scalar quantity
  !!               between 2 z layers. Can handle full step configuration
  !!               using the -full option.
  !!
  !!  ** Method  : compute the sum ( V  * e1 *e2 * e3 *mask )
  !!
  !! History : 2.1  : 11/2008  : J.M. Molines : Original code
  !!           3.0  : 01/2011  : J.M. Molines : Doctor norm + Lic.
  !!           3.0  : 06/2013  : J.M. Molines : multi variables capability
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  USE modutils
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class integration
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                               :: jk, jvar, jvarin, jt! dummy loop index
  INTEGER(KIND=4)                               :: it                  ! time index for vvl
  INTEGER(KIND=4)                               :: ik1, ik2            ! vertical limit of integration
  INTEGER(KIND=4)                               :: narg, iargc         ! command line 
  INTEGER(KIND=4)                               :: ijarg               ! command line 
  INTEGER(KIND=4)                               :: iimin=0, iimax=0    ! horizontal window
  INTEGER(KIND=4)                               :: ijmin=0, ijmax=0    ! horizontal window
  INTEGER(KIND=4)                               :: npiglo, npjglo      ! size of the domain
  INTEGER(KIND=4)                               :: npk, npt            ! size of the domain,
  INTEGER(KIND=4)                               :: nvars, ivar         ! variables in input
  INTEGER(KIND=4)                               :: nvaro=1             ! variables for output
  INTEGER(KIND=4)                               :: ncout, ierr         ! ncid and error status
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: ipk, id_varout      ! levels and varid's of output vars

  REAL(KIND=4)                                  :: rdep_up             ! upper level of integration
  REAL(KIND=4)                                  :: rdep_down           ! lower level of integration
  REAL(KIND=4), DIMENSION(1)                    :: rdep                ! dummy depth output
  REAL(KIND=4), DIMENSION(:),       ALLOCATABLE :: gdep                ! vertical levels
  REAL(KIND=4), DIMENSION(:),       ALLOCATABLE :: e31d                ! vertical metric full
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: e3                  ! metrics
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: zv                  ! working variable
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: hdep                ! depth of the levels
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: zmask               ! mask
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: rlon                ! output longitude
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: rlat                ! output latitude

  REAL(KIND=8)                                  :: dvol                ! total volume
  REAL(KIND=8), DIMENSION(:),       ALLOCATABLE :: dtim                ! time counter
  REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: dvol2d              ! layer volume
  REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: dvertmean           ! value of integral

  CHARACTER(LEN=256)                            :: cf_in               ! input file name
  CHARACTER(LEN=256)                            :: cf_out='vertmean.nc'! output file
  CHARACTER(LEN=256)                            :: cf_e3               ! vertical metric file
  CHARACTER(LEN=256)                            :: cv_cur              ! variable name
  CHARACTER(LEN=256)                            :: cv_dep              ! depth name
  CHARACTER(LEN=256)                            :: cv_e3               ! vertical metric name (partial)
  CHARACTER(LEN=256)                            :: cv_e31d             ! vertical metric name (full)
  CHARACTER(LEN=256)                            :: cv_msk              ! mask variable name
  CHARACTER(LEN=256)                            :: ctype='T'           ! position of the variable
  CHARACTER(LEN=256)                            :: cglobal             ! global attribute
  CHARACTER(LEN=256)                            :: cldum               ! dummy string
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_names            ! name of input variables
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_in               ! name of output variables

  TYPE(variable), DIMENSION(:), ALLOCATABLE     :: stypvarin           ! stucture for attributes (input)
  TYPE(variable), DIMENSION(:), ALLOCATABLE     :: stypvar             ! stucture for attributes (output)

  LOGICAL                                       :: lfull =.FALSE.      ! full step flag
  LOGICAL                                       :: ldebug=.FALSE.      ! debug flag for extra output
  LOGICAL                                       :: lnc4  =.FALSE.      ! netcdf4 flag
  LOGICAL                                       :: lchk                ! file existence flag (true if missing)
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()
  

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage :  cdfvertmean -f IN-file -l LST-var -p C-type -zlim dep1 dep2'
     PRINT *,'              ... [-w imin imax jmin jmax ] '
     PRINT *,'              ... [-full] [-o OUT-file] [-nc4] [-vvl] [-debug]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the vertical mean between dep1 and dep2 given in m, for the '
     PRINT *,'       list of variables LST-var, belonging to the input file.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f IN-file  : netcdf input file.' 
     PRINT *,'       -l LST-var : Comma separated list of input variables to process.'
     PRINT *,'       -p C-type  : one of T U V W indicating position of variable on C-grid'
     PRINT *,'       -zlim dep1 dep2 : depths limit for vertical integration (in meters), '
     PRINT *,'           from top to bottom, positive depths.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [ -w imin imax jmin jmax ] : Define a subarea of the input domain where'
     PRINT *,'            the vertical mean will be applied.'
     PRINT *,'       [-full  ] : for full step configurations. Default is partial step.'
     PRINT *,'       [-debug ] : print some extra informations.'
     PRINT *,'       [-vvl ] : use time-varying vertical metrics.'
     PRINT *,'       [-o OUT-file ] : specify output file instead of ',TRIM(cf_out)
     PRINT *,'       [-nc4 ]: Use netcdf4 output with chunking and deflation level 1.'
     PRINT *,'           This option is effective only if cdftools are compiled with'
     PRINT *,'           a netcdf library supporting chunking and deflation.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       ', TRIM(cn_fzgr),' and ',TRIM(cn_fmsk)
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : varin_vertmean (same units as input variable)'
     PRINT *,'      '
     STOP 
  ENDIF

  ijarg = 1 
  DO WHILE ( ijarg <= narg ) 
     CALL getarg (ijarg, cldum ) ; ijarg = ijarg + 1
     SELECT CASE ( cldum )
     CASE ( '-f'     ) ; CALL getarg (ijarg, cf_in  ) ; ijarg=ijarg+1
     CASE ( '-l'     ) ; CALL getarg (ijarg, cldum  ) ; ijarg=ijarg+1 ; CALL ParseVars(cldum)
     CASE ( '-p'     ) ; CALL getarg (ijarg, ctype  ) ; ijarg=ijarg+1 
     CASE ( '-zlim'  ) ; CALL getarg (ijarg, cldum  ) ; ijarg=ijarg+1 ; READ(cldum,*) rdep_up
        ;                CALL getarg (ijarg, cldum  ) ; ijarg=ijarg+1 ; READ(cldum,*) rdep_down
        ! options
     CASE ('-w'      ) ; CALL getarg(ijarg, cldum   ) ; ijarg = ijarg + 1 ;  READ(cldum,*) iimin
        ;                CALL getarg(ijarg, cldum   ) ; ijarg = ijarg + 1 ;  READ(cldum,*) iimax
        ;                CALL getarg(ijarg, cldum   ) ; ijarg = ijarg + 1 ;  READ(cldum,*) ijmin
        ;                CALL getarg(ijarg, cldum   ) ; ijarg = ijarg + 1 ;  READ(cldum,*) ijmax
     CASE ( '-full'  ) ; lfull  = .TRUE.
     CASE ( '-debug' ) ; ldebug = .TRUE.
     CASE ( '-vvl'   ) ; lg_vvl = .TRUE.
     CASE ( '-nc4'   ) ; lnc4   = .TRUE.
     CASE ( '-o'     ) ; CALL getarg (ijarg, cf_out ) ; ijarg=ijarg+1
     CASE DEFAULT      ; PRINT *,' ERROR : ',TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  ENDDO

  lchk = chkfile (cn_fzgr)
  lchk = chkfile (cn_fmsk) .OR. lchk
  lchk = chkfile (cf_in  ) .OR. lchk
  IF ( lchk  ) STOP 99 ! missing files

  CALL SetGlobalAtt (cglobal)

  IF (rdep_down < rdep_up ) THEN
     PRINT *,'Give depth limits in increasing order !'
     STOP 99
  ENDIF

  npiglo = getdim (cf_in,cn_x)
  npjglo = getdim (cf_in,cn_y)
  npk    = getdim (cf_in,cn_z)
  npt    = getdim (cf_in,cn_t)

  IF ( iimin /= 0 ) THEN ; npiglo = iimax -iimin + 1;  ELSE ; iimin=1 ;
  ENDIF
  IF ( ijmin /= 0 ) THEN ; npjglo = ijmax -ijmin + 1;  ELSE ; ijmin=1 ;
  ENDIF

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  nvars       = getnvar(cf_in)
  ALLOCATE( cv_names(nvars), stypvarin(nvars) )
  ALLOCATE( stypvar(nvaro), ipk(nvaro), id_varout(nvaro) )
  cv_names(:) = getvarname(cf_in, nvars, stypvarin)

  ! just chck if var exist in file 
  DO jvar = 1, nvaro
     IF ( chkvar( cf_in, cv_in(jvar)) ) STOP 99  ! message is written in cdfio.chkvar
  ENDDO

  ! Allocate arrays
  ALLOCATE ( zmask(npiglo,npjglo), dvertmean(npiglo, npjglo) )
  ALLOCATE ( rlon(npiglo,npjglo), rlat(npiglo, npjglo)       )
  ALLOCATE ( zv(npiglo,npjglo), hdep(npiglo,npjglo)          )
  ALLOCATE ( e3(npiglo,npjglo), dvol2d(npiglo,npjglo)        )
  ALLOCATE ( gdep(npk), dtim(npt)                            )
  IF ( lfull ) ALLOCATE ( e31d(npk) )

  rlon =getvar(cf_in, cn_vlon2d, 1, npiglo, npjglo, kimin=iimin, kjmin=ijmin)  ! nav_lon
  rlat =getvar(cf_in, cn_vlat2d, 1, npiglo, npjglo, kimin=iimin, kjmin=ijmin)  ! nav_lat


  CALL CreateOutput

  SELECT CASE ( ctype)
  CASE( 'T','t');  cv_dep=cn_gdepw ; cv_e3 = cn_ve3t ; cv_e31d=cn_ve3t1d ;  cv_msk = cn_tmask ; cf_e3 = cn_fe3t
                   IF (lg_vvl ) cv_e3 = cn_ve3tvvl
  CASE( 'U','u');  cv_dep=cn_gdepw ; cv_e3 = cn_ve3u ; cv_e31d=cn_ve3t1d ;  cv_msk = cn_umask ; cf_e3 = cn_fe3u
                   IF (lg_vvl ) cv_e3 = cn_ve3uvvl
  CASE( 'V','v');  cv_dep=cn_gdepw ; cv_e3 = cn_ve3v ; cv_e31d=cn_ve3t1d ;  cv_msk = cn_vmask ; cf_e3 = cn_fe3v
                   IF (lg_vvl ) cv_e3 = cn_ve3vvvl
  CASE( 'W','w');  cv_dep=cn_gdept ; cv_e3 = cn_ve3w ; cv_e31d=cn_ve3w1d ;  cv_msk = cn_tmask ; cf_e3 = cn_fe3w
                   IF (lg_vvl ) cv_e3 = cn_ve3wvvl
  CASE DEFAULT ; PRINT *,'Point type ', TRIM(ctype),' not known! ' ; STOP 99
  END SELECT

  IF ( lg_vvl) cf_e3 = cf_in

               gdep(:) = getvare3(cn_fzgr, cv_dep,  npk)
  IF ( lfull ) e31d(:) = getvare3(cn_fzgr, cv_e31d, npk)

  ! Look for ik1 and ik2 as nearest level of rdep_up and rdep_down
  ik1 = 1; ik2 = npk
  DO jk=1,npk
     IF ( gdep(jk) <= rdep_up   ) ik1 = jk
     IF ( gdep(jk) <= rdep_down ) ik2 = jk
  ENDDO

  PRINT '(a,2f10.3)', 'depth limit of integration : ', rdep_up, rdep_down
  PRINT '(a,2i10  )', 'nearest level found        : ', ik1,     ik2
  PRINT '(a,2f10.3)', 'corresponding depth        : ', gdep(ik1), gdep(MIN(ik2+1,npk))

  DO jt=1,npt
     IF ( lg_vvl ) THEN ; it = jt
     ELSE               ; it = 1
     ENDIF

     DO jvar = 1,nvaro
        cv_cur = cv_in(jvar)
        dvol           = 0.d0
        dvol2d(:,:)    = 0.d0
        dvertmean(:,:) = 0.d0

        DO jk = ik1, ik2
           ! Get values at jk
           zv(   :,:) = getvar(cf_in,   cv_cur, jk, npiglo, npjglo, kimin=iimin, kjmin=ijmin, ktime=jt )
           zmask(:,:) = getvar(cn_fmsk, cv_msk, jk, npiglo, npjglo, kimin=iimin, kjmin=ijmin           )

           ! get e3 at level jk ( ps...)
           IF ( lfull ) THEN ; e3(:,:) = e31d(jk)
           ELSE              ; e3(:,:) = getvar(cf_e3, cv_e3, jk, npiglo, npjglo, kimin=iimin, kjmin=ijmin,  &
                             &                                                    ktime=it, ldiom=.NOT.lg_vvl )
           ENDIF

           IF ( jk == ik1 ) THEN
              hdep(:,:) = gdep(jk) + e3(:,:)
              e3(  :,:) = MIN(e3, hdep -rdep_up          )
           ENDIF

           IF ( jk == ik2 ) THEN
              e3(  :,:) = MIN(e3, (rdep_down) - gdep(jk) )
           ENDIF

           dvol       = SUM( DBLE(e3 * zmask) )
           dvol2d     =      e3 * zmask * 1.d0  + dvol2d
           dvertmean  = zv * e3 * zmask * 1.d0  + dvertmean

           IF (dvol == 0 )THEN
              ! no more layer below !    
              EXIT   ! get out of the jk loop
           ENDIF
        END DO

        ! Output to netcdf file 
        WHERE ( dvol2d /= 0 )  dvertmean = dvertmean/dvol2d
        ierr = putvar(ncout, id_varout(jvar), REAL(dvertmean), 1, npiglo, npjglo, ktime=jt)
     END DO
  END DO  ! loop on time

  ierr = closeout(ncout)
CONTAINS
  SUBROUTINE ParseVars (cdum)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE ParseVars  ***
    !!
    !! ** Purpose :  Decode variable name  option from command line
    !!
    !! ** Method  :  look for , in the argument string and set the number of
    !!         variable (nvaro), allocate cv_in array and fill it with the
    !!         decoded  names.
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(in) :: cdum

    CHARACTER(LEN=80), DIMENSION(100) :: cl_dum  ! 100 is arbitrary
    INTEGER  :: ji
    INTEGER  :: inchar,  i1=1
    !!----------------------------------------------------------------------
    inchar= LEN(TRIM(cdum))
    ! scan the input string and look for ',' as separator
    DO ji=1,inchar
       IF ( cdum(ji:ji) == ',' ) THEN
          cl_dum(nvaro) = cdum(i1:ji-1)
          i1=ji+1
          nvaro=nvaro+1
       ENDIF
    ENDDO

    ! last name of the list does not have a ','
    cl_dum(nvaro) = cdum(i1:inchar)

    ALLOCATE ( cv_in(nvaro) )
    IF ( ldebug) PRINT *,' SELECTED VARIABLES :'
    DO ji=1, nvaro
       cv_in(ji) = cl_dum(ji)
       IF ( ldebug) PRINT *, "    ",TRIM(cv_in(ji))
    ENDDO
  END SUBROUTINE ParseVars

  SUBROUTINE CreateOutput
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutput  ***
    !!
    !! ** Purpose :  Create netcdf output file(s) 
    !!
    !! ** Method  :  Use stypvar global description of variables
    !!
    !!----------------------------------------------------------------------

    rdep(1)                      = 0.
    ipk(:)                       = 1
    DO jvar=1, nvaro
       DO jvarin=1,nvars
          IF ( cv_in(jvar) == stypvarin(jvarin)%cname ) EXIT  ! cv_in match cv_varin.
       END DO
       stypvar(jvar)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
       stypvar(jvar)%cname             = TRIM(cv_in(jvar))//'_vert_mean'
       stypvar(jvar)%cunits            = stypvarin(jvarin)%cunits
       stypvar(jvar)%rmissing_value    = stypvarin(jvarin)%rmissing_value
       stypvar(jvar)%valid_min         = stypvarin(jvarin)%valid_min
       stypvar(jvar)%valid_max         = stypvarin(jvarin)%valid_max
       stypvar(jvar)%clong_name        = 'vertical average of '//TRIM(stypvarin(jvarin)%clong_name)
       stypvar(jvar)%cshort_name       = TRIM(cv_in(jvar))//'_vert_mean'
       stypvar(jvar)%conline_operation = 'N/A'
       stypvar(jvar)%caxis             = 'TYX'
    END DO

    ! Initialize output file
    ncout = create      (cf_out, cf_in,   npiglo, npjglo, 1                          , ld_nc4=lnc4)
    ierr  = createvar   (ncout,  stypvar, nvaro,  ipk,    id_varout, cdglobal=cglobal, ld_nc4=lnc4)
    ierr  = putheadervar(ncout,  cf_in,   npiglo, npjglo, 1, pnavlon=rlon, pnavlat=rlat, pdep=rdep)

    dtim = getvar1d(cf_in, cn_vtimec, npt     )
    ierr = putvar1d(ncout, dtim,      npt, 'T')

  END SUBROUTINE CreateOutput

END PROGRAM cdfvertmean
