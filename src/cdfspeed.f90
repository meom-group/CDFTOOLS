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

  INTEGER(KIND=4)                            :: ji, jj, jk, jt, jlev ! dummy loop index
  INTEGER(KIND=4)                            :: narg, iargc, ijarg   ! command line 
  INTEGER(KIND=4)                            :: ncout, ierr          ! output file stuff
  INTEGER(KIND=4)                            :: npiglo, npjglo       ! size of the domain
  INTEGER(KIND=4)                            :: npk, npt, nlev       ! size of the domain
  INTEGER(KIND=4)                            :: nvpk                 ! vertical levels in working variable
  INTEGER(KIND=4)                            :: ik                   ! level counter
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: nklevel              ! requested levels
  INTEGER(KIND=4), DIMENSION(1)              :: ipk, id_varout       ! output variable vertical level, varid

  REAL(KIND=4), DIMENSION(:),    ALLOCATABLE :: gdept, gdeptall      ! deptht values for requested/all levels
  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: zu, zv, zspeed       ! working arrays, speed

  REAL(KIND=8), DIMENSION(:),    ALLOCATABLE :: dtim                 ! time counter array

  CHARACTER(LEN=256)                         :: cf_vfil, cf_ufil     ! file for u and v components
  CHARACTER(LEN=256)                         :: cf_tfil='none'       ! file for T point position
  CHARACTER(LEN=256)                         :: cv_u, cv_v           ! name of u and v variable
  CHARACTER(LEN=256)                         :: cf_out='speed.nc'    ! output file name
  CHARACTER(LEN=256)                         :: cldum                ! dummy char variable

  TYPE (variable), DIMENSION(1)              :: stypvar              ! structure for attibutes

  LOGICAL                                    :: lforcing             ! forcing flag
  LOGICAL                                    :: lchk   = .FALSE.     ! flag for missing files
  LOGICAL                                    :: lnc4   = .FALSE.     ! flag for netcdf4 chunking and deflation
  LOGICAL                                    :: lcgrid = .FALSE.     ! flag for C-grid 2D data
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfspeed  -u U-file U-var -v V-file V-var [-t T-file] ...'
     PRINT *,'            ... [-o OUT-file] [-nc4] [-lev LST-level] [-C]' 
     PRINT *,'       '
     PRINT *,'    PURPOSE :'
     PRINT *,'       Compute the speed of ocean currents or wind speed.'
     PRINT *,'       '
     PRINT *,'       If the input files are 3D, the input is assumed to be a model'
     PRINT *,'       output on native C-grid. Speed is computed on the A-grid.'
     PRINT *,'       '
     PRINT *,'       If the input file is 2D then we assume that this is a forcing'
     PRINT *,'       file already on the A-grid, unless -C option is used.'
     PRINT *,'    '
     PRINT *,'    ARGUMENTS :'
     PRINT *,'       -u U-file U-var : netcdf file for U component and variable name.'
     PRINT *,'       -v V-file V-var : netcdf file for V componentt and variable name.'
     PRINT *,'    '
     PRINT *,'    OPTIONS :'
     PRINT *,'       [-t T-file] : indicate any file on gridT for correct header of the'
     PRINT *,'             output file (needed for 3D files or if -C option is used).'
     PRINT *,'       [-lev LST-level] : indicate a list of levels to be processed.'
     PRINT *,'             If not used, all levels are processed.'
     PRINT *,'       [-C] : indicates that data are on a C-grid even if input files are 2D.'
     PRINT *,'       [-o OUT-file] : use specified output file instead of ',TRIM(cf_out)
     PRINT *,'       [-nc4] : use netcdf4 output with chunking and deflation.'
     PRINT *,'    '
     PRINT *,'    OUTPUT :'
     PRINT *,'       Output on ',TRIM(cf_out),' unless -o option is used.'
     PRINT *,'         netcdf variable : U '
     PRINT *,'    '
     PRINT *,'    SEE ALSO :'
     PRINT *,'       cdfvita also computes the speed.'
     PRINT *,'    '
     STOP 
  ENDIF

  nlev =0
  ijarg=1
  DO WHILE ( ijarg <= narg ) 
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-u'   ) ; CALL getarg(ijarg, cf_ufil ) ; ijarg=ijarg+1
        ;              CALL getarg(ijarg, cv_u    ) ; ijarg=ijarg+1
     CASE ( '-v'   ) ; CALL getarg(ijarg, cf_vfil ) ; ijarg=ijarg+1
        ;              CALL getarg(ijarg, cv_v    ) ; ijarg=ijarg+1
        ! options
     CASE ( '-lev' ) ; CALL GetLevList
     CASE ( '-t'   ) ; CALL getarg(ijarg, cf_tfil ) ; ijarg=ijarg+1
     CASE ( '-o'   ) ; CALL getarg(ijarg, cf_out  ) ; ijarg=ijarg+1
     CASE ( '-nc4' ) ; lnc4   = .TRUE.
     CASE ( '-C'   ) ; lcgrid = .TRUE.
     CASE DEFAULT    ; PRINT *,' ERROR : ',TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  ENDDO

  IF ( cf_tfil /= 'none' ) lchk = lchk .OR. chkfile (cf_tfil)
  lchk = lchk .OR. chkfile (cf_ufil)
  lchk = lchk .OR. chkfile (cf_vfil)
  IF ( lchk ) STOP 99 ! missing file

  npiglo = getdim (cf_vfil,cn_x)
  npjglo = getdim (cf_vfil,cn_y)
  npk    = getdim (cf_vfil,cn_z)
  nvpk   = getvdim(cf_vfil,cv_v)
  npt    = getdim (cf_vfil,cn_t)

  IF ( (npk == 0) ) THEN
     IF ( lcgrid ) THEN ; lforcing = .FALSE.
     ELSE               ; lforcing = .TRUE.
     ENDIF
     npk=1
     PRINT *, 'W A R N I N G : you used a forcing field'
  ELSE
     lforcing=.FALSE.
     IF ( TRIM(cf_tfil) == 'none' ) THEN
        PRINT *,'  ERROR: you must specify a griT file as fifth argument '
        PRINT *,'     This is for the proper header of output file '
        STOP 99
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

  ! Allocate arrays
  ALLOCATE ( zv(npiglo,npjglo), zu(npiglo,npjglo), zspeed(npiglo,npjglo), dtim(npt))

  CALL CreateOutput

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
    ! define new variables for output
    stypvar(1)%ichunk=(/npiglo,MAX(1,npjglo/30),1,1 /)
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

    dtim = getvar1d(cf_vfil, cn_vtimec, npt     )
    ierr = putvar1d(ncout, dtim,        npt, 'T')

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

END PROGRAM cdfspeed
