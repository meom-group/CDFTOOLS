PROGRAM cdfmoy_weighted
  !!======================================================================
  !!                     ***  PROGRAM  cdfmoy_weighted  ***
  !!=====================================================================
  !!  ** Purpose : Compute weighted mean values from already processed
  !!               mean files (by cdfmoy)
  !!
  !!               when computing the time average. 
  !!
  !! History : 2.1  : 11/2009  : J.M. Molines : Original code
  !!           3.0  : 04/2011  : J.M. Molines : Doctor norm + Lic.
  !!----------------------------------------------------------------------
  !!      function       : comments
  !!  setweight   : return weight for given variable and file
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

  INTEGER(KIND=4)                               :: jk, jt, jvar        ! dummy loop index
  INTEGER(KIND=4)                               :: ierr                ! working integer
  INTEGER(KIND=4)                               :: narg, iargc, ijarg  ! command line
  INTEGER(KIND=4)                               :: npiglo, npjglo, npk ! size of the domain
  INTEGER(KIND=4)                               :: nvars               ! number of variables in a file
  INTEGER(KIND=4)                               :: ixtra               ! number of tags to process
  INTEGER(KIND=4)                               :: iweight             ! variable weight
  INTEGER(KIND=4)                               :: ncout               ! ncid of output file
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: id_var              ! array of input var id's
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: ipk                 ! array of output var levels
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: id_varout           ! array of output var id's

  REAL(KIND=4), DIMENSION (:,:),    ALLOCATABLE :: v2d                 ! array to read a layer of data
  REAL(KIND=4), DIMENSION(1)                    :: timean, tim         ! time counter

  REAL(KIND=8), DIMENSION (:,:),    ALLOCATABLE :: dtab                ! array for cumulated values
  REAL(KIND=8)                                  :: dtotal_time, dsumw  ! cumulated times and weights

  CHARACTER(LEN=256)                            :: cf_in               ! current input file name
  CHARACTER(LEN=256)                            :: cf_out='cdfmoy_weighted.nc' ! output file name
  CHARACTER(LEN=256)                            :: cv_dep              ! name of depth variable
  CHARACTER(LEN=256)                            :: cldum               ! dummy character variable
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_names            ! array of var name
  
  TYPE (variable), DIMENSION(:),    ALLOCATABLE :: stypvar             ! structure for output var attributes

  LOGICAL                                       :: lold5d=.false.      ! flag for old5d output
  LOGICAL                                       :: lmonth=.false.      ! flag for true month output
  LOGICAL                                       :: lleap=.false.       ! flag for leap years
  LOGICAL                                       :: lnc4=.false.        ! flag for netcdf4 output with chunking and deflation
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfmoy_weighted list of files [-old5d ] [-month] [-leap] ...'
     PRINT *,'                [-o output file]'
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute weight average of files. The weight for each file is'
     PRINT *,'       read from the iweight attribute. In particular, this attribute'
     PRINT *,'       is set to the number of elements used when computing a time'
     PRINT *,'       average (cdfmoy program). A primary application is thus for'
     PRINT *,'       computing annual mean from monthly means.' 
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       The list of files to be averaged, which are supposed to be of' 
     PRINT *,'       the same type and to contain the same variables. This list MUST'
     PRINT *,'       be given before any options'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-old5d ] : This option is used to mimic/replace the cdfmoy_annual'
     PRINT *,'                   which is no longer available. With this option, 12 monthly'
     PRINT *,'                   files must be given, and it is assumed that the monthly'
     PRINT *,'                   means were computed from 5d output of a simulation using'
     PRINT *,'                   a noleap calendar ( weights are fixed, predetermined)'
     PRINT *,'       [-month ] : This option is used to build annual mean from true month'
     PRINT *,'                   output (1mo) in XIOS output for instance.'
     PRINT *,'       [-leap ] : This option has only effect together with the -month option.'
     PRINT *,'                  When used set 29 days in february'
     PRINT *,'       [ -nc4 ] : Use netcdf4 chunking and deflation in output file.'
     PRINT *,'       [-o output file ] : Specify the name for output file instead of the'
     PRINT *,'                 default name ', TRIM(cf_out)
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       none'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'       variables : same as in the input files'
     STOP
  ENDIF

  ! scan command line and check if files exist
  ijarg = 1 ; ixtra=0
  DO WHILE ( ijarg <= narg ) 
    CALL getarg ( ijarg, cldum ) ; ijarg = ijarg +1
    SELECT CASE ( cldum )
    CASE ( '-old5d' )  ; lold5d = .TRUE.
    CASE ( '-month' )  ; lmonth = .TRUE.
    CASE ( '-leap'  )  ; lleap  = .TRUE.
    CASE ( '-nc4'   )  ; lnc4   = .TRUE.
    CASE ( '-o'     )  ; CALL getarg ( ijarg, cf_out ) ; ijarg = ijarg +1
    CASE DEFAULT
        ixtra = ixtra + 1
        cf_in = cldum
        IF ( chkfile (cldum ) ) STOP ! missing file
    END SELECT
  ENDDO

  ! additional check in case of old_5d averaged files
  IF ( lold5d .OR. lmonth ) THEN
    IF ( ixtra /= 12 ) THEN 
      PRINT *,' +++ ERROR : exactly 12 monthly files are required for -old5d/-month options.'
      STOP
    ENDIF
  ENDIF

  npiglo = getdim (cf_in, cn_x                              )
  npjglo = getdim (cf_in, cn_y                              )
  npk    = getdim (cf_in, cn_z, cdtrue=cv_dep, kstatus=ierr )

  IF (ierr /= 0 ) THEN
     npk   = getdim (cf_in,'z',cdtrue=cv_dep, kstatus=ierr   )
     IF (ierr /= 0 ) THEN
       npk   = getdim (cf_in,'sigma',cdtrue=cv_dep,kstatus=ierr)
        IF ( ierr /= 0 ) THEN
          npk = getdim (cf_in,'nav_lev',cdtrue=cv_dep,kstatus=ierr)
            IF ( ierr /= 0 ) THEN
              npk = getdim (cf_in,'levels',cdtrue=cv_dep,kstatus=ierr)
              IF ( ierr /= 0 ) THEN
                PRINT *,' assume file with no depth'
                npk=0
              ENDIF
            ENDIF
        ENDIF
     ENDIF
  ENDIF

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk

  ALLOCATE( dtab(npiglo,npjglo), v2d(npiglo,npjglo) )

  nvars = getnvar(cf_in)
  PRINT *,' nvars =', nvars

  ALLOCATE (cv_names(nvars) )
  ALLOCATE (stypvar(nvars)  )
  ALLOCATE (id_var(nvars), ipk(nvars), id_varout(nvars) )

  ! get list of variable names and collect attributes in stypvar (optional)
  cv_names(:) = getvarname(cf_in, nvars, stypvar)

  id_var(:)   = (/(jvar, jvar=1,nvars)/)

  ! ipk gives the number of level or 0 if not a T[Z]YX  variable
  ipk(:)     = getipk(cf_in, nvars, cdep=cv_dep)
  WHERE( ipk == 0 ) cv_names='none'
  stypvar(:)%cname = cv_names

  DO jk = 1, nvars
    stypvar(jk)%ichunk  = (/ npiglo, MAX(1,npjglo/30), 1, 1 /)
  ENDDO

  ! create output file taking the sizes in cf_in
  ncout = create      (cf_out, cf_in,   npiglo, npjglo, npk,      cdep=cv_dep , ld_nc4=lnc4 )
  ierr  = createvar   (ncout , stypvar, nvars,  ipk,    id_varout             , ld_nc4=lnc4 )
  ierr  = putheadervar(ncout , cf_in,   npiglo, npjglo, npk,      cdep=cv_dep               )

  DO jvar = 1,nvars
     IF ( cv_names(jvar) == cn_vlon2d .OR. &
          cv_names(jvar) == cn_vlat2d        ) THEN
        ! skip these variable
     ELSE
        PRINT *,' Working with ', TRIM(cv_names(jvar)), ipk(jvar)
        DO jk = 1, ipk(jvar)
           PRINT *,'Level ',jk
           dtab(:,:) = 0.d0 ; dtotal_time = 0.d0 ; dsumw=0.d0

           DO jt = 1, ixtra
              CALL getarg   (jt, cf_in)

              iweight   = setweight(cf_in, jt, cv_names(jvar)) 
              dsumw     = dsumw + iweight
              v2d(:,:)  = getvar(cf_in, cv_names(jvar), jk ,npiglo, npjglo )
              dtab(:,:) = dtab(:,:) + iweight * v2d(:,:)

              IF (jk == 1 .AND. jvar == nvars )  THEN
                 tim         = getvar1d(cf_in, cn_vtimec, 1 )
                 dtotal_time = dtotal_time + iweight * tim(1)
              END IF
           END DO

           ! finish with level jk ; compute mean (assume spval is 0 )
           ! store variable on outputfile
           ierr = putvar(ncout, id_varout(jvar), SNGL(dtab(:,:)/dsumw), jk, npiglo, npjglo, kwght=INT(dsumw) )
           IF (jk == 1 .AND. jvar == nvars )  THEN
              timean(1) = dtotal_time/dsumw
              ierr      = putvar1d(ncout, timean, 1, 'T')
           END IF
        END DO  ! loop to next level
     END IF
  END DO ! loop to next var in file

  ierr = closeout(ncout)

  CONTAINS

  INTEGER(KIND=4) FUNCTION setweight( cdfile, kt, cdvar )
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION setweight  ***
    !!
    !! ** Purpose : Return the weight of cdvar in cdfile
    !!
    !! ** Method  : Get attribute iweight from cdfvar in cdfile.
    !!              If lold5d is true, assume weight for 5d build monthly 
    !!              means. If iweight not found 1 is return.
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*),   INTENT(in) :: cdfile
    INTEGER(KIND=4),    INTENT(in) :: kt
    CHARACTER(LEN=*),   INTENT(in) :: cdvar

    INTEGER(KIND=4), DIMENSION(12) :: iweight5d=(/6,5,7,6,6,6,6,6,6,6,6,7/)
    INTEGER(KIND=4), DIMENSION(12) :: iweightmo=(/31,28,31,30,31,30,31,31,30,31,30,31/)
    INTEGER(KIND=4), DIMENSION(12) :: iweightleap=(/31,29,31,30,31,30,31,31,30,31,30,31/)
    !!----------------------------------------------------------------------
    IF ( lold5d ) THEN 
      setweight = iweight5d(kt)
    ELSE IF ( lmonth ) THEN
      IF ( lleap ) THEN
        setweight = iweightleap(kt)
      ELSE
        setweight = iweightmo(kt)
      ENDIF
    ELSE
      setweight = getatt( cdfile, cdvar, 'iweight') 
      IF ( setweight == 0 ) setweight = 1
    ENDIF

    END FUNCTION setweight

END PROGRAM cdfmoy_weighted
