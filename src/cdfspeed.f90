PROGRAM cdfspeed
  !!======================================================================
  !!                     ***  PROGRAM  cdfspeed  ***
  !!=====================================================================
  !!  ** Purpose : combine u and v to obtains the wind speed
  !!
  !!  ** Method  : speed=sqrt(u**2 + v**2)
  !!
  !! History : 2.1  : 11/2007  : P. Mathiot   : Original code
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

  INTEGER(KIND=4)                            :: ji, jj, jk, jt, jlev ! dummy loop index
  INTEGER(KIND=4)                            :: narg, iargc, ijarg   ! command line 
  INTEGER(KIND=4)                            :: ncout, ierr          ! output file stuff
  INTEGER(KIND=4)                            :: npiglo, npjglo       ! size of the domain
  INTEGER(KIND=4)                            :: npk, npt, nlev       ! size of the domain
  INTEGER(KIND=4)                            :: nvpk                 ! vertical levels in working variable
  INTEGER(KIND=4)                            :: ik                   ! level counter
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: nklevel              ! requested levels
  INTEGER(KIND=4), DIMENSION(1)              :: ipk, id_varout       ! output variable vertical level, varid

  REAL(KIND=4), DIMENSION(:),    ALLOCATABLE :: tim                  ! time counter array
  REAL(KIND=4), DIMENSION(:),    ALLOCATABLE :: gdept, gdeptall      ! deptht values for requested/all levels
  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: zu, zv, zspeed       ! working arrays, speed

  CHARACTER(LEN=256)                         :: cf_vfil, cf_ufil     ! file for u and v components
  CHARACTER(LEN=256)                         :: cf_tfil='none'       ! file for T point position
  CHARACTER(LEN=256)                         :: cv_u, cv_v           ! name of u and v variable
  CHARACTER(LEN=256)                         :: cf_out='speed.nc'    ! output file name
  CHARACTER(LEN=256)                         :: cldum                ! dummy char variable

  TYPE (variable), DIMENSION(1)              :: stypvar              ! structure for attibutes

  LOGICAL                                    :: lforcing             ! forcing flag
  LOGICAL                                    :: lnc4 = .false.       ! flag for netcdf4 chunking and deflation
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfspeed  U-file V-file U-var V-var [-t T-file] ...'
     PRINT *,'            ... [-nc4] [-o OUT-file ] [-lev level_list]' 
     PRINT *,'    PURPOSE :'
     PRINT *,'       Computes the speed of ocean currents or wind speed'
     PRINT *,'       '
     PRINT *,'       If the input files are 3D, the input is assumed to be '
     PRINT *,'       a model output on native C-grid. Speed is computed on the A-grid.'
     PRINT *,'       '
     PRINT *,'       If the input file is 2D and then we assume that this is '
     PRINT *,'       a forcing file already on the A-grid.'
     PRINT *,'    '
     PRINT *,'    ARGUMENTS :'
     PRINT *,'       U-file : netcdf file for U component'
     PRINT *,'       V-file : netcdf file for V component'
     PRINT *,'       U-var  : netcdf variable name for U component'
     PRINT *,'       V-var  : netcdf variable name for V component'
     PRINT *,'    '
     PRINT *,'    OPTIONS :'
     PRINT *,'       -t T-file  : indicate any file on gridT for correct header'
     PRINT *,'                 of the output file (usefull for 3D files)'
     PRINT *,'       -lev level_list  : indicate a list of levels to be processed'
     PRINT *,'                 If not used, all levels are processed.'
     PRINT *,'                 This option should be the last on the command line'
     PRINT *,'       -nc4 : use netcdf4 output with chunking and deflation'
     PRINT *,'       -o OUT-file : use specified output file instead of ',TRIM(cf_out)
     PRINT *,'    '
     PRINT *,'    OUTPUT :'
     PRINT *,'       Output on ',TRIM(cf_out),'  variable U '
     STOP
  ENDIF

  nlev =0
  ijarg=1
  DO WHILE ( ijarg <= narg ) 
     CALL getarg(ijarg, cldum ) ; ijarg = ijarg + 1
     SELECT CASE ( cldum )
     CASE ( '-lev' ) 
       nlev = narg -ijarg + 1
       ALLOCATE ( nklevel(nlev) )
       DO jlev = 1, nlev 
          CALL getarg(ijarg, cldum ) ; ijarg = ijarg + 1 ; READ( cldum,*) nklevel(jlev)
       END DO
     CASE ( '-t' ) 
       CALL getarg(ijarg, cf_tfil ) ; ijarg = ijarg + 1
       IF ( chkfile (cf_tfil) ) STOP ! missing file
     CASE ( '-nc4' ) 
       lnc4=.true.
     CASE ( '-o' ) 
       CALL getarg(ijarg, cf_out ) ; ijarg = ijarg + 1
     CASE DEFAULT
       cf_ufil = cldum
       CALL getarg(ijarg, cf_vfil ) ; ijarg = ijarg + 1
       IF ( chkfile(cf_ufil) .OR. chkfile(cf_vfil) ) STOP ! missing file
       CALL getarg(ijarg, cv_u ) ; ijarg = ijarg + 1
       CALL getarg(ijarg, cv_v ) ; ijarg = ijarg + 1
     END SELECT 
  ENDDO

  npiglo = getdim (cf_vfil,cn_x)
  npjglo = getdim (cf_vfil,cn_y)
  npk    = getdim (cf_vfil,cn_z)
  nvpk   = getvdim(cf_vfil,cv_v)
  npt    = getdim (cf_vfil,cn_t)

  IF ( (npk == 0) ) THEN
     lforcing=.TRUE.
     npk=1
     PRINT *, 'W A R N I N G : you used a forcing field'
  ELSE
     lforcing=.FALSE.
     IF ( TRIM(cf_tfil) == 'none' ) THEN
       PRINT *,'  ERROR: you must specify a griT file as fifth argument '
       PRINT *,'     This is for the proper header of output file '
       STOP
     ENDIF
  END IF

  IF ( nlev == 0 ) THEN 
    nlev = npk 
    ALLOCATE ( nklevel(nlev) )
    DO jlev =1, nlev
      nklevel(jlev) = jlev
    ENDDO
  ENDIF

  IF (nvpk == 2 ) nvpk = 1
  IF (nvpk == 3 ) nvpk = nlev

  PRINT *, 'npiglo =', npiglo
  PRINT *, 'npjglo =', npjglo
  PRINT *, 'npk    =', npk
  PRINT *, 'nvpk   =', nvpk
  PRINT *, 'nlev   =', nlev
  PRINT *, 'npt    =', npt

  IF ( nlev >  nvpk ) THEN
     PRINT *, 'W A R N I N G : nlev larger than nvpk, we assume nlev=nvpk'
     nlev = nvpk
  END IF

  ! choose chunk size for output ... not easy not used if lnc4=.false. but
  ! anyway ..
  stypvar(1)%ichunk=(/npiglo,MAX(1,npjglo/30),1,1 /)

  ! define new variables for output
  stypvar(1)%cname             = 'U'
  stypvar(1)%cunits            = 'm.s-1'
  stypvar(1)%rmissing_value    = 0.
  stypvar(1)%valid_min         = -1000.
  stypvar(1)%valid_max         = 1000.
  stypvar(1)%clong_name        = 'Current or wind speed'
  stypvar(1)%cshort_name       = 'U'
  stypvar(1)%conline_operation = 'N/A'
  stypvar(1)%caxis             = 'TZYX'

  ! create output fileset
  IF (lforcing ) THEN
     ipk(1) = 1 !  2D  no dep variable
     ncout  = create      (cf_out, cf_vfil, npiglo, npjglo, 0         , ld_nc4=lnc4 )
     ierr   = createvar   (ncout,  stypvar, 1,      ipk,    id_varout , ld_nc4=lnc4 )
     ierr   = putheadervar(ncout,  cf_vfil, npiglo, npjglo, 0                       )
     nlev=1 ; nklevel(nlev) = 1
  ELSE
     ALLOCATE ( gdept(nlev), gdeptall(npk) )
     gdeptall = getvar1d ( cf_tfil, cn_vdeptht, npk )
     DO jlev = 1, nlev
       gdept(jlev) = gdeptall( nklevel(jlev) )
     END DO
     ipk(1) = nlev
     ncout  = create      (cf_out, cf_tfil, npiglo, npjglo, nlev      , ld_nc4=lnc4   )
     ierr   = createvar   (ncout,  stypvar, 1,      ipk,    id_varout , ld_nc4=lnc4   )
     ierr   = putheadervar(ncout,  cf_tfil, npiglo, npjglo, nlev,  pdep=gdept         )
  END IF

  ! Allocate arrays
  ALLOCATE ( zv(npiglo,npjglo), zu(npiglo,npjglo), zspeed(npiglo,npjglo), tim(npt))

  DO jt=1,npt
     tim(jt)=jt
  END DO

  ierr=putvar1d(ncout, tim, npt, 'T')

  DO jt = 1,npt
     DO jlev = 1, nlev
           ik = nklevel(jlev)
        ! Get velocities v at jk
           zu(:,:) = getvar(cf_ufil, cv_u, ik, npiglo, npjglo, ktime=jt)
           zv(:,:) = getvar(cf_vfil, cv_v, ik, npiglo, npjglo, ktime=jt)
        IF ( lforcing ) THEN
           ! u and v are already on the T grid points
        ELSE
           ! in this case we are on the C-grid and the speed must be computed 
           ! on the A-grid. We use reverse loop in order to use only one array
           DO ji=npiglo,2,-1
             DO jj=1,npjglo
               zu(ji,jj) = 0.5*(zu(ji-1,jj)+zu(ji,jj))
             ENDDO
           ENDDO

           DO ji=1,npiglo 
             DO jj=npjglo,2,-1
               zv(ji,jj) = 0.5*(zv(ji,jj-1)+zv(ji,jj))
             ENDDO
           ENDDO
        END IF
        zspeed    = SQRT(zv*zv+zu*zu)
        ierr = putvar(ncout, id_varout(1), zspeed, jlev ,npiglo, npjglo, ktime=jt)
     END DO
  END DO
  ierr = closeout(ncout)

END PROGRAM cdfspeed
