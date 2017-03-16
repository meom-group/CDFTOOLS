PROGRAM cdfenstat
  !!======================================================================
  !!                     ***  PROGRAM  cdfenstat  ***
  !!=====================================================================
  !!  ** Purpose : Compute mean values standard dev for all the variables 
  !!               in a bunch of cdf files given as arguments (ensemble)
  !!               Store the results on a 'similar' cdf file (same number
  !!               of time frame) with same variables + stdev_variable
  !!
  !!  ** Method  : Use optimize algorithm for better accuracy
  !!
  !!  ** Reference : Accuracy of floating point arithmetic.
  !!
  !! History  3.0  : 07/2015  : J.M. Molines : from cdfmoyt
  !!        : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!----------------------------------------------------------------------
  USE cdfio 
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id: cdfenstat.f90 550 2011-09-19 16:07:45Z molines $
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class ensemble_processing
  !!-----------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                               :: jk, jfil, jrec, jv! dummy loop index
  INTEGER(KIND=4)                               :: ierr               ! working integer
  INTEGER(KIND=4)                               :: idep, idep_max     ! possible depth index, maximum
  INTEGER(KIND=4)                               :: inpt               ! working integer
  INTEGER(KIND=4)                               :: narg, iargc, ijarg ! browsing command line
  INTEGER(KIND=4)                               :: nfiles             ! number of files to average
  INTEGER(KIND=4)                               :: npiglo, npjglo     ! size of the domain
  INTEGER(KIND=4)                               :: npk, npt           ! size of the domain
  INTEGER(KIND=4)                               :: nvars              ! number of variables in a file
  INTEGER(KIND=4)                               :: ncout              ! ncid of output files
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: id_var             ! arrays of var id's
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: ipk                ! arrays of vertical level for each var
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: id_varout          ! varid's of average vars

  REAL(KIND=4), DIMENSION(:,:,:),   ALLOCATABLE :: v3d                ! array to read a layer of data for all files
  REAL(KIND=4), DIMENSION(:,:,:,:), ALLOCATABLE :: v4d                ! array to read a layer of data for all files
  REAL(KIND=4), DIMENSION(:),       ALLOCATABLE :: tim                ! time counter
  REAL(KIND=4), DIMENSION(:),       ALLOCATABLE :: zspval_in          ! input missing value

  REAL(KIND=8), DIMENSION(:,:,:),   ALLOCATABLE :: dtabn, dtab2n      ! arrays for cumulated values
  REAL(KIND=8), DIMENSION(:,:,:),   ALLOCATABLE :: dtabb, dtab2b      ! arrays for cumulated values
  REAL(KIND=8), DIMENSION(:,:,:),   ALLOCATABLE :: dtmp               ! temporary array
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: d4tabn, d4tab2n    ! arrays for cumulated values
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: d4tabb, d4tab2b    ! arrays for cumulated values
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: d4tmp              ! temporary array

  CHARACTER(LEN=256)                            :: cf_in              ! input file names
  CHARACTER(LEN=256)                            :: cf_out  = 'cdfmoy.nc'  ! output file for average
  CHARACTER(LEN=256)                            :: cv_dep             ! depth dimension name
  CHARACTER(LEN=256)                            :: cldum              ! dummy string argument
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cf_lst             ! list of input files
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_nam             ! array of var name
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: clv_dep            ! array of possible depth name (or 3rd dimension)

  TYPE (variable), DIMENSION(:),    ALLOCATABLE :: stypvar            ! attributes for average values

  LOGICAL                                       :: lspval0 = .FALSE.  ! cdfmoy_chsp flag
  LOGICAL                                       :: lnc4    = .FALSE.  ! flag for netcdf4 chinking and deflation
  LOGICAL                                       :: lv4d    = .FALSE.  ! flag for netcdf4 chinking and deflation
  !!----------------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfenstat -l LIST-mbr-files [-spval0] [-v4d] [-o OUT-file] [-nc4]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute an ensemble mean and standard deviation for a set of files'
     PRINT *,'       corresponding to different members from an ensemble run.'
     PRINT *,'       This program assumes that each of the member files have the same '
     PRINT *,'       variables, and the same number of time frames.'
     PRINT *,'       This program uses optimal algorithm for computing the mean and std dev,'
     PRINT *,'       in order to reduce truncation errors.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -l LIST-mbr-files : A list of members  model output files. '
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [ -spval0 ] :  set missing_value attribute to 0 for all output'
     PRINT *,'               variables and take care of the input missing_value.'
     PRINT *,'               This option is usefull if missing_values differ from files '
     PRINT *,'               to files; it was formely done by cdfmoy_chsp).'
     PRINT *,'       [ -v4d ] : uses 4D arrays for improved performance (use more memory !)'
     PRINT *,'       [ -o OUT-file ] : specify a name for output file instead of ',TRIM(cf_out)
     PRINT *,'       [ -nc4 ] : output file will be in netcdf4, with chunking and deflation'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       none '
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out), 'unless -o option in use'
     PRINT *,'       variables : are the same than in the input files. Standard Dev are '
     PRINT *,'        named  stdev_<variable>'
     STOP
  ENDIF

  ALLOCATE ( cf_lst(narg) )
  ! look for -spval0 option and set up cf_lst, nfiles 
  ijarg = 1 
  nfiles = 0
  DO WHILE ( ijarg <= narg )
     CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1
     SELECT CASE ( cldum )
     CASE ( '-l'      ) ;  CALL GetFileList
        ! options
     CASE ( '-spval0' ) ; lspval0 = .TRUE.
     CASE ( '-nc4   ' ) ; lnc4    = .TRUE.
     CASE ( '-v4d  '  ) ; lv4d    = .TRUE.
     CASE ( '-o   '   ) ; CALL getarg (ijarg, cf_out) ; ijarg = ijarg + 1
     CASE DEFAULT       ; PRINT *,' ERROR : ', TRIM(cldum),' : unknown option.' ; STOP
     END SELECT
  END DO

  IF ( nfiles == 0 ) THEN 
     PRINT *, ' ERROR : You must give a list of member files with -l option' ; STOP
  ENDIF

  ! Initialisation from  1rst file (all file are assume to have the same geometry)
  IF ( chkfile (cf_lst(1)) ) STOP ! missing file

  cf_in  = cf_lst(1)
  npiglo = getdim (cf_in,cn_x)
  npjglo = getdim (cf_in,cn_y)
  npt    = getdim (cf_in,cn_t)

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

  ! check that all files have the same number of time frames
  ierr = 0
  DO jfil = 1, nfiles
     IF (  chkfile (cf_lst(jfil)      ) ) STOP ! missing file
     inpt = getdim (cf_lst(jfil), cn_t)
     IF ( inpt /= npt ) THEN
        PRINT *, 'File ',TRIM(cf_lst(jfil) ),' has ',inpt,' time frames instead of ', npt
        ierr = ierr + 1
     ENDIF
  ENDDO
  IF ( ierr /= 0 ) STOP ! frame numbers mismatch

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  IF ( .NOT. lv4d ) THEN 
     ALLOCATE( v3d(npiglo,npjglo,npt) )
     ALLOCATE( dtabn(npiglo,npjglo,npt), dtab2n(npiglo,npjglo,npt)) 
     ALLOCATE( dtabb(npiglo,npjglo,npt), dtab2b(npiglo,npjglo,npt), dtmp(npiglo,npjglo,npt))
  ENDIF

  nvars = getnvar(cf_in)
  PRINT *,' nvars = ', nvars

  ALLOCATE (cv_nam(2*nvars) )
  ALLOCATE (stypvar(2*nvars))
  ALLOCATE (id_var(2*nvars), ipk(2*nvars), id_varout(2*nvars)  )
  ALLOCATE( tim(npt) )

  CALL CreateOutput

  DO jv = 1,nvars
     IF ( cv_nam(jv) == cn_vlon2d .OR. &     ! nav_lon
          & cv_nam(jv) == cn_vlat2d .OR. &
          & cv_nam(jv) == 'none'    ) THEN     ! nav_lat
        ! skip these variable
     ELSE
        PRINT *,' Working with ', TRIM(cv_nam(jv)), ipk(jv)
        IF ( lv4d) THEN
           ALLOCATE( v4d(npiglo,npjglo,ipk(jv),npt) )
           ALLOCATE( d4tabn(npiglo,npjglo,ipk(jv),npt), d4tab2n(npiglo,npjglo,ipk(jv),npt) ) 
           ALLOCATE( d4tabb(npiglo,npjglo,ipk(jv),npt), d4tab2b(npiglo,npjglo,ipk(jv),npt), d4tmp(npiglo,npjglo,ipk(jv),npt))

           DO jfil = 1, nfiles
              cf_in     = cf_lst(jfil)
              v4d(:,:,:,:)  = getvar4d(cf_in, cv_nam(jv), npiglo, npjglo, ipk(jv), npt )
              IF ( jfil == 1 ) THEN
                 d4tabb(:,:,:,:)=v4d(:,:,:,:) ; d4tab2b(:,:,:,:)=0.d0
                 CYCLE
              ENDIF
              IF ( lspval0 )  WHERE ( v4d == zspval_in(jv) )  v4d = 0.  ! change missing values to 0
              d4tmp(:,:,:,:)   = v4d(:,:,:,:) - d4tabb(:,:,:,:)
              d4tabn(:,:,:,:)  = d4tabb(:,:,:,:)  + d4tmp(:,:,:,:) / jfil
              d4tab2n(:,:,:,:) = d4tab2b(:,:,:,:) + d4tmp(:,:,:,:) * ( v4d(:,:,:,:) - d4tabn(:,:,:,:) )
              ! swap tabs
              d4tabb(:,:,:,:)  = d4tabn(:,:,:,:)
              d4tab2b(:,:,:,:) = d4tab2n(:,:,:,:)
           END DO

           ! store variable on outputfile
           DO jrec = 1, npt
              DO jk=1, ipk(jv)
                 ierr = putvar(ncout, id_varout(jv),       SNGL(d4tabn(:,:,jk,jrec)),                   jk, npiglo, npjglo, kwght=nfiles, ktime = jrec )
                 ierr = putvar(ncout, id_varout(jv+nvars), SQRT(SNGL(d4tab2n(:,:,jk,jrec))/(nfiles-1)), jk, npiglo, npjglo, kwght=nfiles, ktime = jrec )
              ENDDO
           ENDDO

           ! Deallocate array for this variable
           DEALLOCATE( v4d, d4tabn, d4tab2n, d4tabb, d4tab2b, d4tmp )
        ELSE  ! work with 3D arrays
           DO jk = 1, ipk(jv)
              PRINT *,'level ',jk
              DO jfil = 1, nfiles
                 cf_in     = cf_lst(jfil)
                 v3d(:,:,:)  = getvar3dt(cf_in, cv_nam(jv), jk, npiglo, npjglo, npt )
                 IF ( jfil == 1 ) THEN
                    dtabb(:,:,:)=v3d(:,:,:) ; dtab2b(:,:,:)=0.d0
                    CYCLE
                 ENDIF
                 IF ( lspval0 )  WHERE ( v3d == zspval_in(jv) )  v3d = 0.  ! change missing values to 0
                 dtmp(:,:,:)   = v3d(:,:,:) - dtabb(:,:,:)
                 dtabn(:,:,:)  = dtabb(:,:,:)  + dtmp(:,:,:) / jfil
                 dtab2n(:,:,:) = dtab2b(:,:,:) + dtmp(:,:,:) * ( v3d(:,:,:) - dtabn(:,:,:) )
                 ! swap tabs
                 dtabb(:,:,:)  = dtabn(:,:,:)
                 dtab2b(:,:,:) = dtab2n(:,:,:)
              END DO

              ! store variable on outputfile
              DO jrec = 1, npt
                 ierr = putvar(ncout, id_varout(jv),       SNGL(dtabn(:,:,jrec)),                   jk, npiglo, npjglo, kwght=nfiles, ktime = jrec )
                 ierr = putvar(ncout, id_varout(jv+nvars), SQRT(SNGL(dtab2n(:,:,jrec))/(nfiles-1)), jk, npiglo, npjglo, kwght=nfiles, ktime = jrec )
              ENDDO
           END DO  ! loop to next level
        ENDIF

     END IF
  END DO ! loop to next var in file

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
    ! get list of variable names and collect attributes in stypvar (optional)
    cv_nam(1:nvars) = getvarname(cf_in,nvars, stypvar)
    DO jv =1, nvars
       cv_nam(jv+nvars) ='stdev_'//TRIM(cv_nam(jv))
    ENDDO
    ! choose chunk size for output ... not easy not used if lnc4=.false. but anyway ..
    DO jv = 1, 2*nvars
       stypvar(jv)%ichunk=(/npiglo,MAX(1,npjglo/30),1,1 /)
    ENDDO

    IF ( lspval0 ) THEN 
       ALLOCATE ( zspval_in(nvars) )
       zspval_in(:) = stypvar(1:nvars)%rmissing_value
       stypvar(:)%rmissing_value = 0.
    ENDIF

    DO jv = 1, nvars
       ! variables that will not be computed or stored are named 'none'
       stypvar(jv+nvars)%cname             = cv_nam(jv+nvars)
       stypvar(jv+nvars)%cunits            = TRIM(stypvar(jv)%cunits)             ! unit
       stypvar(jv+nvars)%rmissing_value    = stypvar(jv)%rmissing_value           ! missing_value
       stypvar(jv+nvars)%valid_min         = 0.                                   ! valid_min = zero
       stypvar(jv+nvars)%valid_max         = stypvar(jv)%valid_max                ! valid_max *valid_max
       stypvar(jv+nvars)%scale_factor      = 1.
       stypvar(jv+nvars)%add_offset        = 0.
       stypvar(jv+nvars)%savelog10         = 0.
       stypvar(jv+nvars)%clong_name        = TRIM(stypvar(jv)%clong_name)//'_Std_Dev'   ! 
       stypvar(jv+nvars)%cshort_name       = cv_nam(jv+nvars)
       stypvar(jv+nvars)%conline_operation = TRIM(stypvar(jv)%conline_operation) 
       stypvar(jv+nvars)%caxis             = TRIM(stypvar(jv)%caxis) 
    END DO

    id_var(:)  = (/(jv, jv=1,2*nvars)/)
    ! ipk gives the number of level or 0 if not a T[Z]YX  variable
    ipk(1:nvars)  = getipk (cf_in,nvars,cdep=cv_dep)
    ipk(nvars+1:2*nvars) = ipk(1:nvars)
    WHERE( ipk == 0 ) cv_nam='none'
    stypvar( :)%cname = cv_nam

    ! create output file taking the sizes in cf_in
    ncout  = create      (cf_out,  cf_in,     npiglo, npjglo, npk, cdep=cv_dep, ld_nc4=lnc4)
    ierr   = createvar   (ncout ,  stypvar,  2*nvars,  ipk,   id_varout       , ld_nc4=lnc4)
    ierr   = putheadervar(ncout,   cf_in,    npiglo, npjglo, npk, cdep=cv_dep)

    ! all files shoud have the same time , take the first of the list !
    tim(:) = getvar1d(cf_in, cn_vtimec, npt)
    ierr   = putvar1d(ncout,  tim, npt, 'T')

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

END PROGRAM cdfenstat
