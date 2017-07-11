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
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  !!      function       : comments
  !!  setweight   : return weight for given variable and file
  !!----------------------------------------------------------------------
  USE cdfio 
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class time_averaging
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                               :: jk, jt, jvar        ! dummy loop index
  INTEGER(KIND=4)                               :: ierr                ! working integer
  INTEGER(KIND=4)                               :: idep, idep_max      ! possible depth inde
  INTEGER(KIND=4)                               :: narg, iargc, ijarg  ! command line
  INTEGER(KIND=4)                               :: npiglo, npjglo, npk ! size of the domain
  INTEGER(KIND=4)                               :: nvars               ! number of variables in a file
  INTEGER(KIND=4)                               :: nfiles              ! number of tags to process
  INTEGER(KIND=4)                               :: iweight             ! variable weight
  INTEGER(KIND=4)                               :: ncout               ! ncid of output file
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: id_var              ! array of input var id's
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: ipk                 ! array of output var levels
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: id_varout           ! array of output var id's

  REAL(KIND=4), DIMENSION (:,:),    ALLOCATABLE :: v2d                 ! array to read a layer of data
  REAL(KIND=4), DIMENSION (:,:),    ALLOCATABLE :: e3                  ! array to read vertical metrics
  REAL(KIND=4), DIMENSION (:)  ,    ALLOCATABLE :: v1d                 ! array to read column of data
  REAL(KIND=4), DIMENSION(1)                    :: timean, tim         ! time counter

  REAL(KIND=8), DIMENSION (:,:),    ALLOCATABLE :: dtab                ! array for cumulated values
  REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: de3s                ! arrays for cumulated e3 (vvl)
  REAL(KIND=8), DIMENSION (:)  ,    ALLOCATABLE :: dtab1d              ! array for cumulated values
  REAL(KIND=8)                                  :: dtotal_time, dsumw  ! cumulated times and weights

  CHARACTER(LEN=256)                            :: cf_in               ! current input file name
  CHARACTER(LEN=256)                            :: cf_out='cdfmoy_weighted.nc' ! output file name
  CHARACTER(LEN=256)                            :: cf_e3               ! file name for reading vertical metrics (vvl)
  CHARACTER(LEN=256)                            :: cv_dep              ! name of depth variable
  CHARACTER(LEN=256)                            :: cv_e3               ! name of e3t variable for vvl
  CHARACTER(LEN=256)                            :: cv_skip             ! name of  variable to skip
  CHARACTER(LEN=256)                            :: cldum               ! dummy character variable
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cf_lst              ! file name list of input files
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_names            ! array of var name
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: clv_dep             ! array of possible depth name or 3rd dim.

  TYPE (variable), DIMENSION(:),    ALLOCATABLE :: stypvar             ! structure for output var attributes

  LOGICAL                                       :: lold5d = .FALSE.      ! flag for old5d output
  LOGICAL                                       :: lmonth = .FALSE.      ! flag for true month output
  LOGICAL                                       :: lleap  = .FALSE.      ! flag for leap years
  LOGICAL                                       :: lnc4   = .FALSE.      ! flag for netcdf4 output with chunking and deflation
  LOGICAL                                       :: ll_vvl                ! working flag for vvl AND 3D fields
  LOGICAL                                       :: lchk                  ! flag for checking file existence
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfmoy_weighted -l LST-files [-old5d ] [-month] [-leap] ...'
     PRINT *,'       ... [-skip variable] [-vvl] [-o OUT-file] [-nc4]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute weighted average of files. The weight for each file is read from'
     PRINT *,'       the iweight variable attribute. In particular, this attribute is set to'
     PRINT *,'       the number of elements used when computing a time average (''cdfmoy'').'
     PRINT *,'       A primary application is thus for computing annual mean from monthly '
     PRINT *,'       means.' 
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -l LST-files : The list of files to be averaged, which are supposed to'
     PRINT *,'             be of the same type and to contain the same variables.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-old5d ] : This option is used to mimic/replace the cdfmoy_annual which'
     PRINT *,'             is no longer available. With this option, 12 monthly files must be'
     PRINT *,'             given, and it is assumed that the monthly means were computed from'
     PRINT *,'             5d output of a simulation using a noleap  calendar (weights are '
     PRINT *,'             fixed, predetermined).'
     PRINT *,'       [-month ] : This option is used to build annual mean from true month'
     PRINT *,'             output (1mo) in XIOS output for instance.'
     PRINT *,'       [-leap ] : This option has only effect together with the -month option.'
     PRINT *,'             When used set 29 days in february.'
     PRINT *,'       [-skip variable ] : name of variable to skip, in the input file. '
     PRINT *,'       [-vvl ] : Use time-varying vertical metrics for weighted averages.'
     PRINT *,'       [-o OUT-file] : Specify the name for output file instead of'
     PRINT *,'           ', TRIM(cf_out)
     PRINT *,'       [-nc4 ] : Use netcdf4 chunking and deflation in output file.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       none'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'      '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'       variables : same as in the input files'
     PRINT *,'      '
     PRINT *,'     SEE ALSO : '
     PRINT *,'       cdfmoy, cdfmoyt, cdfmoy_freq'
     PRINT *,'      '
     STOP 
  ENDIF

  ! scan command line and check if files exist
  ijarg = 1 
  DO WHILE ( ijarg <= narg ) 
     CALL getarg ( ijarg, cldum ) ; ijarg = ijarg +1
     SELECT CASE ( cldum )
     CASE ( '-l'     ) ; CALL GetFileList
        ! options
     CASE ( '-old5d' ) ; lold5d = .TRUE.
     CASE ( '-month' ) ; lmonth = .TRUE.
     CASE ( '-leap'  ) ; lleap  = .TRUE.
     CASE ( '-vvl'   ) ; lg_vvl = .TRUE.
     CASE ( '-nc4'   ) ; lnc4   = .TRUE.
     CASE ( '-o'     ) ; CALL getarg ( ijarg, cf_out ) ; ijarg = ijarg +1
     CASE ( '-skip'  ) ; CALL getarg ( ijarg, cv_skip) ; ijarg = ijarg +1
     CASE DEFAULT      ; PRINT *,' ERROR : ', TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  ENDDO
  
  lchk = .FALSE.
  DO jt = 1,nfiles
     lchk = lchk .OR. chkfile ( cf_lst(jt) )
  ENDDO
  IF ( lchk ) STOP 99 ! missing files  
 
  ! work with 1rst file for dimension lookup
  cf_in=cf_lst(1)

  ! additional check in case of old_5d averaged files
  IF ( lold5d .OR. lmonth ) THEN
     IF ( nfiles /= 12 ) THEN 
        PRINT *,' +++ ERROR : exactly 12 monthly files are required for -old5d/-month options.'
        STOP 99
     ENDIF
  ENDIF

  npiglo = getdim (cf_in, cn_x )
  npjglo = getdim (cf_in, cn_y )

  ! looking for npk among various possible name
  idep_max=8
  ALLOCATE ( clv_dep(idep_max) )
  clv_dep(:) = (/cn_z,'z','sigma','nav_lev','levels','ncatice','icbcla','icbsect'/)
  idep=1  ; ierr=1000
  DO WHILE ( ierr /= 0 .AND. idep <= idep_max )
     npk  = getdim (cf_in, clv_dep(idep), cdtrue=cv_dep, kstatus=ierr)
     idep = idep + 1
  ENDDO

  IF ( ierr /= 0 ) THEN  ! none of the dim name was found
      PRINT *,' assume file with no depth'
      npk=0
  ENDIF

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk

  nvars = getnvar(cf_in)
  PRINT *,' nvars =', nvars

  ALLOCATE(               dtab(npiglo,npjglo), v2d(npiglo,npjglo) )
  IF ( lg_vvl ) ALLOCATE( de3s(npiglo,npjglo),  e3(npiglo,npjglo) )
  ALLOCATE (cv_names(nvars) )
  ALLOCATE (stypvar(nvars)  )
  ALLOCATE (id_var(nvars), ipk(nvars), id_varout(nvars) )

  ! Prepare output files
  CALL CreateOutput
  ! for vvl, look for e3x variable in the file
  cv_e3 = 'none'
  IF ( lg_vvl ) THEN
    DO jvar = 1, nvars
      IF ( INDEX(cv_names(jvar), 'e3' ) /= 0 ) THEN
         cv_e3=cv_names(jvar)
         EXIT
      ENDIF
    ENDDO
  ENDIF

  DO jvar = 1,nvars
     ll_vvl=lg_vvl .AND.  (ipk(jvar) > 1) ! JMM : assume 2D var are not weigthed averaged !!!
     IF ( cv_names(jvar) == cn_vlon2d .OR. &
          cv_names(jvar) == cn_vlat2d .OR. &
          cv_names(jvar) == 'none'    .OR. &
          cv_names(jvar) == cv_skip        ) THEN
        ! skip these variable
     ELSE
        PRINT *,' Working with ', TRIM(cv_names(jvar)), ipk(jvar)
        IF ( npiglo == 1 .AND. npjglo == 1 ) THEN 
         !
           ALLOCATE (v1d( ipk(jvar)), dtab1d(ipk(jvar)) )
           dtab1d(:) = 0.d0 ; dtotal_time=0.d0 ; dsumw=0.d0
           DO jt=1, nfiles
              cf_in = cf_lst(jt)
              iweight   = setweight(cf_in, jt, cv_names(jvar)) 
              dsumw     = dsumw + iweight
              tim = getvar1d(cf_in, cn_vtimec, 1 )
              dtotal_time = dtotal_time + iweight * tim(1)
              v1d=getvare3(cf_in, cv_names(jvar), ipk(jvar) )
              dtab1d(:)=dtab1d(:) + iweight * v1d(:)
           ENDDO
           timean(1) = dtotal_time/dsumw
           ierr      = putvar1d(ncout, timean, 1, 'T')
           ierr      = putvar(ncout, id_varout(jvar), SNGL(dtab1d(:)/dsumw), ipk(jvar), 'vert', ktime=1 , kwght=INT(dsumw) )  ! module interface to putvare3
           DEALLOCATE (v1d, dtab1d)
        ELSE
        DO jk = 1, ipk(jvar)
           PRINT *,'Level ',jk
           dtab(:,:) = 0.d0 ; dtotal_time = 0.d0 ; dsumw=0.d0
           IF ( ll_vvl ) THEN  ; de3s(:,:)  = 0.d0    ; ENDIF

           DO jt = 1, nfiles
              cf_in = cf_lst(jt)
              iweight   = setweight(cf_in, jt, cv_names(jvar)) 
              dsumw     = dsumw + iweight
              v2d(:,:)  = getvar(cf_in, cv_names(jvar), jk ,npiglo, npjglo )

              IF ( ll_vvl ) THEN
                 cf_e3     = cf_in
                 e3(:,:)   = getvar (cf_e3, cv_e3, jk ,npiglo, npjglo, ktime=jt )
                 de3s(:,:) = de3s(:,:) + iweight * e3(:,:)  ! cumulate e3
                 dtab(:,:) = dtab(:,:) + iweight * e3(:,:) * v2d(:,:)
              ELSE
                 dtab(:,:) = dtab(:,:) + iweight * v2d(:,:)
              ENDIF

              IF (jk == 1 .AND. jvar == nvars )  THEN
                 tim         = getvar1d(cf_in, cn_vtimec, 1 )
                 dtotal_time = dtotal_time + iweight * tim(1)
              END IF
           END DO

           ! finish with level jk ; compute mean (assume spval is 0 )
           ! store variable on outputfile
           IF ( ll_vvl ) THEN
             ierr = putvar(ncout, id_varout(jvar), SNGL(dtab(:,:)/de3s ), jk, npiglo, npjglo, kwght=INT(dsumw) )
           ELSE
             ierr = putvar(ncout, id_varout(jvar), SNGL(dtab(:,:)/dsumw), jk, npiglo, npjglo, kwght=INT(dsumw) )
           ENDIF
           IF (jk == 1 .AND. jvar == nvars )  THEN
              timean(1) = dtotal_time/dsumw
              ierr      = putvar1d(ncout, timean, 1, 'T')
           END IF
        END DO  ! loop to next level
        END IF ! 1D variable
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
  
  SUBROUTINE CreateOutput
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutput  ***
    !!
    !! ** Purpose :  Create netcdf output file 
    !!
    !! ** Method  :  Create file, define dims and variables
    !!
    !!----------------------------------------------------------------------
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
  END SUBROUTINE CreateOutput

  SUBROUTINE GetFileList
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE GetFileList  ***
    !!
    !! ** Purpose :  Set up a file list given on the command line as 
    !!               blank separated list
    !!
    !! ** Method  :  Scan the command line until a '-' is found
    !!----------------------------------------------------------------------
    INTEGER (KIND=4)  :: ji
    INTEGER (KIND=4)  :: icur
    !!----------------------------------------------------------------------
    !!
    nfiles=0
    ! need to read a list of file ( number unknow ) 
    ! loop on argument till a '-' is found as first char
    icur=ijarg                          ! save current position of argument number
    DO ji = icur, narg                  ! scan arguments till - found
       CALL getarg ( ji, cldum )
       IF ( cldum(1:1) /= '-' ) THEN ; nfiles = nfiles+1
       ELSE                          ; EXIT
       ENDIF
    ENDDO
    ALLOCATE (cf_lst(nfiles) )
    DO ji = icur, icur + nfiles -1
       CALL getarg(ji, cf_lst( ji -icur +1 ) ) ; ijarg=ijarg+1
    END DO
  END SUBROUTINE GetFileList

END PROGRAM cdfmoy_weighted
