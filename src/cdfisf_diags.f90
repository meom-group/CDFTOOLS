PROGRAM cdfisf_diags
  !!======================================================================
  !!                     ***  PROGRAM  cdfisf_diags  ***
  !!=====================================================================
  !!  ** Purpose : Compute the integrated ice shelf melt for each ice shelves 
  !!             specified in a list 
  !!
  !!  ** Method  : total_melt = sum ( melt_rate * e1 * e2 * mask )
  !!
  !! History : 3.0  : 08/2015  : P. Mathiot   : Original code (from cdfsum)
  !!           4.0  : 08/2017  : J.M. Molines : port to CDFTOOLS 4.0
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  USE modutils
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class ice_shelf_processing
  !!-----------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                           :: jt, jisf            ! dummy loop index
  INTEGER(KIND=4)                           :: iimin=0, iimax=0    ! domain limitation for computation
  INTEGER(KIND=4)                           :: ijmin=0, ijmax=0    ! domain limitation for computation
  INTEGER(KIND=4)                           :: ikx=1, iky=1        ! dims of netcdf output file
  INTEGER(KIND=4)                           :: idum                ! dummy integer for txt file reading
  INTEGER(KIND=4)                           :: id_isf              ! ice shelf id
  INTEGER(KIND=4)                           :: iunit=10            ! ice shelf text input file id
  INTEGER(KIND=4)                           :: ierr                ! working integer
  INTEGER(KIND=4)                           :: narg, iargc, ijarg  ! command line 
  INTEGER(KIND=4)                           :: npiglo, npjglo      ! size of the domain
  INTEGER(KIND=4)                           :: npk, npt, nisf      ! size of the domain
  INTEGER(KIND=4)                           :: nvpk                ! vertical levels in working variable
  INTEGER(KIND=4)                           :: ncout               ! for netcdf output
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: ipk, id_varout

  REAL(KIND=4)                              :: rdummy              ! dummy 2d variable for result
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE :: gdep                ! depth 
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: e1, e2, zmlt        ! metrics and ice shelf melt
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: zmask, zmaskisf, mask ! npiglo x npjglo
  REAL(KIND=4), DIMENSION(1,1)              :: rdumlon, rdumlat    ! dummy latitude and longitude
  REAL(KIND=4), DIMENSION(1,1)              :: rdum                ! dummy variable for txt file reading

  REAL(KIND=8)                              :: dfwfobs             ! ice shelf melt observation
  REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: dtim                ! time

  CHARACTER(LEN=256)                        :: cldum, cdum         ! dummy string
  CHARACTER(LEN=3)                          :: cid_isf             ! isf id print in netcdf output long name
  CHARACTER(LEN=256)                        :: cf_isfmsk           ! input file names
  CHARACTER(LEN=256)                        :: cf_isflst           ! input file names
  CHARACTER(LEN=256)                        :: cf_in               ! file name 
  CHARACTER(LEN=256)                        :: cf_out='cdfisf.nc'  ! output file name 
  CHARACTER(LEN=256)                        :: cv_in               ! variable name
  CHARACTER(LEN=20)                         :: cv_msk              ! name of mask variable
  CHARACTER(LEN=20)                         :: cv_isfmsk           ! name of mask variable
  CHARACTER(LEN=256)                        :: cglobal             ! netcdf global attribute


  TYPE(variable), DIMENSION(:), ALLOCATABLE :: stypvar             ! structure of output
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfisf_diags -f IN-file -v IN-var -ff ISF-fill_file -fv ISF-fill_var '
     PRINT *,'                -l ISF-list [-w imin imax jmin jmax kmin kmax] [-o OUT-file] '
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Computes the integrated ice shelf melting for all ice shelves'
     PRINT *,'       defined in the input list.' 
     PRINT *,'       This diags can be optionally limited to a sub-area.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f IN-file : netcdf input file for ice shelf melting.' 
     PRINT *,'       -v IN-var  : netcdf ice shelf melt variable to work with (kg/m2/y).'
     PRINT *,'                 convertion from kg/m2/s to Gt/y is made with 365d.'
     PRINT *,'       -ff ISF-fill_file : file built by cdfisf_fill (all the ice shelves are'
     PRINT *,'                       tagged with an id)'
     PRINT *,'       -fv ISF-fill_var  : name of fill variable to use in ISF-fill_file'
     PRINT *,'       -l ISF-list : ice shelf list outputed by cdfisf_fill'
     PRINT *,'                  (id in the list must match the list in ISF-fill_file)'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       -w imin imax jmin jmax  : limit of the sub area to work with.' 
     PRINT *,'              if imin=0 all i are taken'
     PRINT *,'              if jmin=0 all j are taken'
     PRINT *,'       -o OUT-file : specify the name of the output file instead of ',TRIM(cf_out)
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'      ', TRIM(cn_fhgr),' and ',TRIM(cn_fmsk)
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       Standard output.'
     PRINT *,'       netcdf file : ',TRIM(cf_out),' with 1 variable for each ice shelves'
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'       cdfisf_fill, cdfisf_forcing, cdfisf_poolchk, cdfisf_rnf'
     PRINT *,'      '

     STOP
  ENDIF

  ! global attribute
  ! setting up the building command in global attribute
  CALL SetGlobalAtt (cglobal)  ! append command name to global attribute

  ijarg = 1  ; iw = 0
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg = ijarg + 1
     SELECT CASE (cldum )
     CASE ( '-f'  ) ; CALL getarg (ijarg,cf_in    )  ; ijarg=ijarg+1
     CASE ( '-v'  ) ; CALL getarg (ijarg,cv_in    )  ; ijarg=ijarg+1
     CASE ( '-ff' ) ; CALL getarg (ijarg,cf_isfmsk)  ; ijarg=ijarg+1
     CASE ( '-fv' ) ; CALL getarg (ijarg,cv_isfmsk)  ; ijarg=ijarg+1
     CASE ( '-l'  ) ; CALL getarg (ijarg,cf_isflst)  ; ijarg=ijarg+1
        !
     CASE ( '-w'  ) ; CALL getarg (ijarg,cldum    )  ; ijarg=ijarg+1 ; READ(*,cldum) iimin 
        ;           ; CALL getarg (ijarg,cldum    )  ; ijarg=ijarg+1 ; READ(*,cldum) iimax 
        ;           ; CALL getarg (ijarg,cldum    )  ; ijarg=ijarg+1 ; READ(*,cldum) ijmin 
        ;           ; CALL getarg (ijarg,cldum    )  ; ijarg=ijarg+1 ; READ(*,cldum) ijmax
     CASE ( '-o'  ) ; CALL getarg (ijarg,cf_out   )  ; ijarg=ijarg+1
     CASE DEFAULT ; PRINT *,' ERROR : ', TRIM(cldum) ,' : unknown option.'; STOP 99
     END SELECT
  ENDDO

  ! check input files and variables
  CALL CheckInput

  ! get dimension
  npiglo = getdim (cf_in,cn_x)
  npjglo = getdim (cf_in,cn_y)
  npk    = getdim (cf_in,cn_z)
  npt    = getdim (cf_in,cn_t)
  nvpk   = 1

  IF (iimin /= 0 ) THEN ; npiglo = iimax - iimin + 1;
  ELSE ; iimin = 1 ;
  ENDIF
  IF (ijmin /= 0 ) THEN ; npjglo = ijmax - ijmin + 1;
  ELSE ; ijmin = 1 ;
  ENDIF

  PRINT *, 'Size of the extracted area :'
  PRINT *, '  npiglo = ', npiglo
  PRINT *, '  npjglo = ', npjglo
  PRINT *, '  npk    = ', npk
  PRINT *, '  npt    = ', npt

  IF ( (npk == 0) ) THEN
     npk      = 1
  END IF

  ! Allocate arrays
  ALLOCATE ( zmask(npiglo,npjglo), zmaskisf(npiglo,npjglo), mask(npiglo,npjglo) )
  ALLOCATE ( zmlt(npiglo,npjglo) )
  ALLOCATE ( e1   (npiglo,npjglo), e2(npiglo,npjglo) )
  ALLOCATE ( gdep (npk), dtim(npt) )

  e1(:,:) = getvar  (cn_fhgr, cn_ve1t, 1, npiglo, npjglo, kimin=iimin, kjmin=ijmin)
  e2(:,:) = getvar  (cn_fhgr, cn_ve2t, 1, npiglo, npjglo, kimin=iimin, kjmin=ijmin)

  ! open text file for isf
  OPEN(unit=iunit, file=cf_isflst, form='FORMATTED', status='OLD')

  ! get number of isf
  nisf = 0
  cdum='XXX'
  DO WHILE ( TRIM(cdum) /= 'EOF')
     READ(iunit,*) cdum
     nisf=nisf+1
  END DO
  REWIND(iunit)
  nisf = nisf - 1
  PRINT *, nisf,' ice shelf detected in the input text file'

  ! define new variables for output 
  ALLOCATE( stypvar(nisf), ipk(nisf), id_varout(nisf) )

  ! create netcdf output
  CALL CreateOutput

  ! compute number of isf
  DO jt = 1,npt
     PRINT *, jt,'/',npt

     ! Get velocities v at ik
     zmlt(:,:) = getvar(cf_in, cv_in, 1, npiglo, npjglo, ktime=jt, kimin=iimin, kjmin=ijmin)

     ! get isf mask
     zmask   (:,:) = getvar(cn_fmsk  , cn_tmaskutil, 1, npiglo, npjglo, kimin=iimin, kjmin=ijmin)
     zmaskisf(:,:) = getvar(cf_isfmsk, cv_isfmsk   , 1, npiglo, npjglo, kimin=iimin, kjmin=ijmin)

     ! compute melt for each ice shelf
     REWIND(iunit)
     DO jisf=1,nisf

        ! read text file
        READ(iunit,*) id_isf, cdum, rdum, rdum, idum, idum, idum, idum, dfwfobs 

        ! set mask
        mask(:,:) = 0.0
        WHERE (zmaskisf == -id_isf)
           mask(:,:) = zmask
        END WHERE

        ! compute total melt for each ice shelf
        rdummy = SUM(DBLE(zmlt * e1 * e2 * mask)) *86400.d0 * 365.d0 / 1.e12

        ! sanity check
        IF (ABS(rdummy) > stypvar(jisf)%valid_max) THEN
           PRINT *, 'total melt of ',TRIM(cdum),' is unexpected'
           PRINT *, 'please check that input file is kg/m2/s'
           STOP
        END IF

        ! standard output
        WRITE (6,'(i3,a5,f6.1,f12.1)') jt, TRIM(cdum), dfwfobs, rdummy

        ! netcdf output
        ierr = putvar0d(ncout, id_varout(jisf), rdummy, ktime=jt )

     END DO
  END DO

  ierr=closeout(ncout)

CONTAINS

  SUBROUTINE CreateOutput
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutput  ***
    !!
    !! ** Purpose :  Create the output file. This is done outside the main
    !!               in order to increase readability of the code. 
    !!
    !! ** Method  :  Use global variables, defined in main 
    !!----------------------------------------------------------------------

    REWIND(iunit)
    ! loop on all isf
    DO jisf=1,nisf
       ! read text file
       READ(iunit,*) cid_isf, cdum, rdum, rdum, idum, idum, idum, idum, dfwfobs

       stypvar(jisf)%cunits            = 'Gt/y'
       stypvar(jisf)%rmissing_value    = -1.e4
       stypvar(jisf)%valid_min         = -1.e4
       stypvar(jisf)%valid_max         = 1.e4
       stypvar(jisf)%scale_factor      = 1.
       stypvar(jisf)%add_offset        = 0.
       stypvar(jisf)%savelog10         = 0.
       stypvar(jisf)%conline_operation = 'N/A'

       stypvar(jisf)%cname             = 'isfmelt_'//TRIM(cdum)
       stypvar(jisf)%clong_name        = 'isf melt for '//TRIM(cdum)//' with id: '//TRIM(cid_isf)
       stypvar(jisf)%cshort_name       = 'isfmelt_'//TRIM(cdum)
       stypvar(jisf)%caxis             = 'ZT'
       stypvar(jisf)%cprecision        = 'r8'

       ipk(jisf) = 1
    END DO

    ! create output fileset
    ncout = create      (cf_out, 'none',  ikx,   iky,   nvpk, cdep=cn_gdept)
    ierr  = createvar   (ncout, stypvar, nisf, ipk,   id_varout, cdglobal=TRIM(cglobal) )
    ierr  = putheadervar(ncout, cf_in,  ikx, iky, npk, pnavlon=rdumlon, pnavlat=rdumlat, pdep=gdep(1:nvpk), cdep=cn_gdept)

    ! read/write time variable
    dtim  = getvar1d(cf_in, cn_vtimec, npt     )
    ierr  = putvar1d(ncout, dtim,      npt, 'T')

  END SUBROUTINE CreateOutput

  SUBROUTINE CheckInput
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CheckInput  ***
    !!
    !! ** Purpose :  Check the presence of input file and variable needed
    !!               by the tool. This is done outside the main in order 
    !!               to increase readability of the code. 
    !!
    !! ** Method  :  Use global variables, defined in main 
    !!----------------------------------------------------------------------
    LOGICAL :: lchk                ! flag for missing files and variables

    ! check file
    lchk = chkfile(cn_fhgr  )
    lchk = chkfile(cn_fmsk  ) .OR. lchk
    lchk = chkfile(cf_in    ) .OR. lchk
    lchk = chkfile(cf_isfmsk) .OR. lchk
    lchk = chkfile(cf_isflst) .OR. lchk
    lchk = chkfile(cf_in    ) .OR. lchk
    IF ( lchk ) STOP 96 ! missing file

    ! check var
    lchk = chkvar(cf_isfmsk, cv_isfmsk   )
    lchk = chkvar(cf_in    , cv_in       ) .OR. lchk
    lchk = chkvar(cn_fhgr  , cn_ve1t     ) .OR. lchk
    lchk = chkvar(cn_fhgr  , cn_ve2t     ) .OR. lchk
    lchk = chkvar(cn_fmsk  , cn_tmaskutil) .OR. lchk
    IF ( lchk ) STOP 96 ! missing variable

  END SUBROUTINE CheckInput

END PROGRAM cdfisf_diags
