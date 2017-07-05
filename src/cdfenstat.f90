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
  !!               variables belonging to cn_sqdvar(:), than can be changed 
  !!               in the nam_cdf_names namelist if wished.
  !!
  !!  ** Reference : Accuracy of floating point arithmetic.
  !!
  !! History  3.0  : 07/2015  : J.M. Molines : from cdfmoyt
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!----------------------------------------------------------------------
  USE cdfio 
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id: cdfenstat.f90 550 2011-09-19 16:07:45Z molines $
  !! Copyright (c) 2010, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!-----------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                               :: jk, jfil, jrec, jv! dummy loop index
  INTEGER(KIND=4)                               :: ierr               ! working integer
  INTEGER(KIND=4)                               :: inpt               ! working integer
  INTEGER(KIND=4)                               :: narg, iargc, ijarg ! browsing command line
  INTEGER(KIND=4)                               :: nfil               ! number of files to average
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

  REAL(KIND=8), DIMENSION(:,:,:),     ALLOCATABLE :: dtabn, dtab2n      ! arrays for cumulated values
  REAL(KIND=8), DIMENSION(:,:,:),     ALLOCATABLE :: dtabb, dtab2b      ! arrays for cumulated values
  REAL(KIND=8), DIMENSION(:,:,:),     ALLOCATABLE :: dtmp               ! temporary array
  REAL(KIND=8), DIMENSION(:,:,:,:),   ALLOCATABLE :: d4tabn, d4tab2n    ! arrays for cumulated values
  REAL(KIND=8), DIMENSION(:,:,:,:),   ALLOCATABLE :: d4tabb, d4tab2b    ! arrays for cumulated values
  REAL(KIND=8), DIMENSION(:,:,:,:),   ALLOCATABLE :: d4tmp               ! temporary array
  REAL(KIND=8), DIMENSION(:),       ALLOCATABLE :: dtotal_time        ! to compute mean time

  CHARACTER(LEN=256)                            :: cf_in              ! input file names
  CHARACTER(LEN=256)                            :: cf_out  = 'cdfmoy.nc'  ! output file for average
  CHARACTER(LEN=256)                            :: cv_dep             ! depth dimension name
  CHARACTER(LEN=256)                            :: cldum              ! dummy string argument
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cf_list            ! list of input files
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_nam             ! array of var name

  TYPE (variable), DIMENSION(:),    ALLOCATABLE :: stypvar            ! attributes for average values

  LOGICAL                                       :: lspval0 = .FALSE.  ! cdfmoy_chsp flag
  LOGICAL                                       :: lnc4    = .FALSE.  ! flag for netcdf4 chinking and deflation
  LOGICAL                                       :: lv4d    = .FALSE.  ! flag for netcdf4 chinking and deflation
  !!----------------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfenstat list_of_model_files [-spval0] [-nc4] [-v4d] -o OUT-file]'
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the time average of a list of files given as arguments.' 
     PRINT *,'       This program handle multi time-frame files is such a way that'
     PRINT *,'       the output files are also multi time-frame, each frame being'
     PRINT *,'       the average across the files given in the list.'
     PRINT *,'       '
     PRINT *,'       The program assume that all files in the list are of same'
     PRINT *,'       type (shape, variables , and number of time frames ). '
     PRINT *,'       For some variables, the program also compute the time average '
     PRINT *,'       of the squared variables, which is used in other cdftools '
     PRINT *,'       (cdfeke, cdfrmsssh, cdfstdevw, cdfstddevts ... The actual variables'
     PRINT *,'       selected for squared average are :'
     PRINT '(10x,"- ",a)' , (TRIM(cn_sqdvar(jv)), jv=1, nn_sqdvar)
     PRINT *,'       This selection can be adapted with the nam_cdf_namelist process.'
     PRINT *,'       (See cdfnamelist -i for details).'
     PRINT *,'       If you want to compute the average of already averaged files,'
     PRINT *,'       consider using cdfmoy_weighted instead, in order to take into'
     PRINT *,'       account a particular weight for each file in the list.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       A list of similar model output files. '
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [ -spval0 ] :  set missing_value attribute to 0 for all output'
     PRINT *,'               variables and take care of the input missing_value.'
     PRINT *,'               This option is usefull if missing_values differ from files '
     PRINT *,'               to files; it was formely done by cdfmoy_chsp).'
     PRINT *,'       [ -nc4 ] : output file will be in netcdf4, with chunking and deflation'
     PRINT *,'       [ -v4d ] : uses 4D arrays for improved performance (use more memory !)'
     PRINT *,'       [ -o OUT-file ] : specify a name for output file instead of '//TRIM(cf_out)
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

  ALLOCATE ( cf_list(narg) )
  ! look for -spval0 option and set up cf_list, nfil 
  ijarg = 1 
  nfil = 0
  DO WHILE ( ijarg <= narg )
     CALL getarg (ijarg, cldum) ; ijarg = ijarg + 1
     SELECT CASE ( cldum )
     CASE ( '-spval0' )   ! option to reset spval to 0 in the output files
        lspval0 = .TRUE.
     CASE ( '-nc4   ' )   ! option to reset spval to 0 in the output files
        lnc4 = .TRUE.
     CASE ( '-v4d  '  )   ! option to reset spval to 0 in the output files
        lv4d = .TRUE.
     CASE ( '-o   ' )   ! option to reset spval to 0 in the output files
        CALL getarg (ijarg, cf_out) ; ijarg = ijarg + 1
     CASE DEFAULT         ! then the argument is a file
        nfil          = nfil + 1
        cf_list(nfil) = TRIM(cldum)
     END SELECT
  END DO
  ! Initialisation from  1rst file (all file are assume to have the same geometry)
  ! time counter can be different for each file in the list. It is read in the
  ! loop for files
  IF ( chkfile (cf_list(1)) ) STOP 99 ! missing file

  cf_in  = cf_list(1)
  npiglo = getdim (cf_in,cn_x)
  npjglo = getdim (cf_in,cn_y)
  npk    = getdim (cf_in,cn_z, cdtrue=cv_dep, kstatus=ierr)
  npt    = getdim (cf_in,cn_t)

  IF (ierr /= 0 ) THEN
     npk   = getdim (cf_in,'z',cdtrue=cv_dep,kstatus=ierr)
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

  ! check that all files have the same number of time frames
  ierr = 0
  DO jfil = 1, nfil
     IF (  chkfile (cf_list(jfil)      ) ) STOP 99 ! missing file
     inpt = getdim (cf_list(jfil), cn_t)
     IF ( inpt /= npt ) THEN
        PRINT *, 'File ',TRIM(cf_list(jfil) ),' has ',inpt,' time frames instead of ', npt
        ierr = ierr + 1
     ENDIF
  ENDDO
  IF ( ierr /= 0 ) STOP 99 ! frame numbers mismatch

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  IF ( .NOT. lv4d ) THEN 
  ALLOCATE( v3d(npiglo,npjglo,npt) )
  ALLOCATE( dtabn(npiglo,npjglo,npt), dtab2n(npiglo,npjglo,npt)) 
  ALLOCATE( dtabb(npiglo,npjglo,npt), dtab2b(npiglo,npjglo,npt), dtmp(npiglo,npjglo,npt))
  ENDIF

  ALLOCATE( dtotal_time(npt), tim(npt) )

  nvars = getnvar(cf_in)
  PRINT *,' nvars = ', nvars

  ALLOCATE (cv_nam(2*nvars) )
  ALLOCATE (stypvar(2*nvars))
  ALLOCATE (id_var(2*nvars), ipk(2*nvars), id_varout(2*nvars)  )

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

  ! Compute the mean time for each mean frame
  dtotal_time(:) = 0.d0
  DO jfil = 1, nfil 
     cf_in          = cf_list(jfil)
     tim(:)         = getvar1d(cf_in, cn_vtimec, npt)
     dtotal_time(:) = dtotal_time(:) + tim (:)
  END DO
  tim(:) = dtotal_time(:)/ nfil
  ierr   = putvar1d(ncout,  tim, npt, 'T')


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

              DO jfil = 1, nfil
                 cf_in     = cf_list(jfil)
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
                ierr = putvar(ncout, id_varout(jv),       SNGL(d4tabn(:,:,jk,jrec)),                 jk, npiglo, npjglo, kwght=nfil, ktime = jrec )
                ierr = putvar(ncout, id_varout(jv+nvars), SQRT(SNGL(d4tab2n(:,:,jk,jrec))/(nfil-1)), jk, npiglo, npjglo, kwght=nfil, ktime = jrec )
                ENDDO
              ENDDO

             ! Deqllocate array for this variable
             DEALLOCATE( v4d, d4tabn, d4tab2n, d4tabb, d4tab2b, d4tmp )
           ELSE
           DO jk = 1, ipk(jv)
              PRINT *,'level ',jk
              DO jfil = 1, nfil
                 cf_in     = cf_list(jfil)
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
                ierr = putvar(ncout, id_varout(jv),       SNGL(dtabn(:,:,jrec)),                 jk, npiglo, npjglo, kwght=nfil, ktime = jrec )
                ierr = putvar(ncout, id_varout(jv+nvars), SQRT(SNGL(dtab2n(:,:,jrec))/(nfil-1)), jk, npiglo, npjglo, kwght=nfil, ktime = jrec )
              ENDDO
           END DO  ! loop to next level
           ENDIF 

        END IF
     END DO ! loop to next var in file

  ierr = closeout(ncout)

END PROGRAM cdfenstat
