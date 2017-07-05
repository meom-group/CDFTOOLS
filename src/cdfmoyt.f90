PROGRAM cdfmoyt
  !!======================================================================
  !!                     ***  PROGRAM  cdfmoyt  ***
  !!=====================================================================
  !!  ** Purpose : Compute mean values for all the variables in a bunch
  !!               of cdf files given as arguments.
  !!               Store the results on a 'similar' cdf file. This version
  !!               differ from cdfmoy, because if the input files have many
  !!               time frames in it, the output file will have the same 
  !!               number of time frames, each being the average accross the
  !!               input files.
  !!
  !!  ** Method  : Also store the mean squared values for the nn_sqdvar
  !!               variables belonging to cn_sqdvar(:), than can be changed 
  !!               in the nam_cdf_names namelist if wished.
  !!
  !! History : 2.0  : 11/2004  : J.M. Molines : Original code
  !!         : 2.1  : 06/2007  : P. Mathiot   : Modif for forcing fields
  !!           3.0  : 12/2010  : J.M. Molines : Doctor norm + Lic.
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!   varchk2       : check if variable is candidate for square mean
  !!----------------------------------------------------------------------
  USE cdfio 
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2010, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!-----------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                               :: jk, jfil, jrec     ! dummy loop index
  INTEGER(KIND=4)                               :: jvar, jv, jt       ! dummy loop index
  INTEGER(KIND=4)                               :: ierr               ! working integer
  INTEGER(KIND=4)                               :: inpt               ! working integer
  INTEGER(KIND=4)                               :: narg, iargc, ijarg ! browsing command line
  INTEGER(KIND=4)                               :: nfil               ! number of files to average
  INTEGER(KIND=4)                               :: npiglo, npjglo     ! size of the domain
  INTEGER(KIND=4)                               :: npk, npt           ! size of the domain
  INTEGER(KIND=4)                               :: nvars              ! number of variables in a file
  INTEGER(KIND=4)                               :: ntframe            ! cumul of time frame
  INTEGER(KIND=4)                               :: ncout, ncout2      ! ncid of output files
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: id_var             ! arrays of var id's
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: ipk                ! arrays of vertical level for each var
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: id_varout          ! varid's of average vars
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: id_varout2         ! varid's of sqd average vars

  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: v2d                ! array to read a layer of data
  REAL(KIND=4), DIMENSION(:),       ALLOCATABLE :: tim                ! time counter
  REAL(KIND=4), DIMENSION(:),       ALLOCATABLE :: zspval_in          ! input missing value

  REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: dtab, dtab2        ! arrays for cumulated values
  REAL(KIND=8), DIMENSION(:),       ALLOCATABLE :: dtotal_time        ! to compute mean time

  CHARACTER(LEN=256)                            :: cf_in              ! input file names
  CHARACTER(LEN=256)                            :: cf_out  = 'cdfmoy.nc'  ! output file for average
  CHARACTER(LEN=256)                            :: cf_out2 = 'cdfmoy2.nc' ! output file for squared average
  CHARACTER(LEN=256)                            :: cv_dep             ! depth dimension name
  CHARACTER(LEN=256)                            :: cldum              ! dummy string argument
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cf_list            ! list of input files
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_nam             ! array of var name
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_nam2            ! array of var2 name for output

  TYPE (variable), DIMENSION(:),    ALLOCATABLE :: stypvar            ! attributes for average values
  TYPE (variable), DIMENSION(:),    ALLOCATABLE :: stypvar2           ! attributes for square averaged values

  LOGICAL                                       :: lspval0 = .FALSE.  ! cdfmoy_chsp flag
  !!----------------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfmoyt list_of_model_files [-spval0]  '
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
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       none '
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out),' and ',TRIM(cf_out2)
     PRINT *,'       variables : are the same than in the input files. For squared averages' 
     PRINT *,'       _sqd is append to the original variable name.'
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

  ALLOCATE( dtab(npiglo,npjglo), dtab2(npiglo,npjglo), v2d(npiglo,npjglo) )
  ALLOCATE( dtotal_time(npt), tim(npt) )

  nvars = getnvar(cf_in)
  PRINT *,' nvars = ', nvars

  ALLOCATE (cv_nam(nvars), cv_nam2(nvars) )
  ALLOCATE (stypvar(nvars), stypvar2(nvars) )
  ALLOCATE (id_var(nvars), ipk(nvars), id_varout(nvars), id_varout2(nvars)  )

  ! get list of variable names and collect attributes in stypvar (optional)
  cv_nam(:) = getvarname(cf_in,nvars, stypvar)

  IF ( lspval0 ) THEN 
     ALLOCATE ( zspval_in(nvars) )
     zspval_in(:) = stypvar(:)%rmissing_value
     stypvar(:)%rmissing_value = 0.
  ENDIF

  DO jvar = 1, nvars
     ! variables that will not be computed or stored are named 'none'
     IF ( varchk2 ( cv_nam(jvar) ) ) THEN 
        cv_nam2(jvar)                    = TRIM(cv_nam(jvar))//'_sqd'
        stypvar2(jvar)%cname             = TRIM(stypvar(jvar)%cname)//'_sqd'         ! name
        stypvar2(jvar)%cunits            = '('//TRIM(stypvar(jvar)%cunits)//')^2'    ! unit
        stypvar2(jvar)%rmissing_value    = stypvar(jvar)%rmissing_value              ! missing_value
        stypvar2(jvar)%valid_min         = 0.                                        ! valid_min = zero
        stypvar2(jvar)%valid_max         = stypvar(jvar)%valid_max**2                ! valid_max *valid_max
        stypvar2(jvar)%scale_factor      = 1.
        stypvar2(jvar)%add_offset        = 0.
        stypvar2(jvar)%savelog10         = 0.
        stypvar2(jvar)%clong_name        = TRIM(stypvar(jvar)%clong_name)//'_Squared'   ! 
        stypvar2(jvar)%cshort_name       = TRIM(stypvar(jvar)%cshort_name)//'_sqd'     !
        stypvar2(jvar)%conline_operation = TRIM(stypvar(jvar)%conline_operation) 
        stypvar2(jvar)%caxis             = TRIM(stypvar(jvar)%caxis) 
     ELSE
        cv_nam2(jvar) = 'none'
     END IF
  END DO

  id_var(:)  = (/(jv, jv=1,nvars)/)
  ! ipk gives the number of level or 0 if not a T[Z]YX  variable
  ipk(:)     = getipk (cf_in,nvars,cdep=cv_dep)
  WHERE( ipk == 0 ) cv_nam='none'
  stypvar( :)%cname = cv_nam
  stypvar2(:)%cname = cv_nam2

  ! create output file taking the sizes in cf_in
  ncout  = create      (cf_out,  cf_in,    npiglo, npjglo, npk, cdep=cv_dep)
  ierr   = createvar   (ncout ,  stypvar,  nvars,  ipk,    id_varout       )
  ierr   = putheadervar(ncout,   cf_in,    npiglo, npjglo, npk, cdep=cv_dep)

  ncout2 = create      (cf_out2, cf_in,    npiglo, npjglo, npk, cdep=cv_dep)
  ierr   = createvar   (ncout2,  stypvar2, nvars,  ipk,    id_varout2      )
  ierr   = putheadervar(ncout2,  cf_in,    npiglo, npjglo, npk, cdep=cv_dep)

  ! Compute the mean time for each mean frame
  dtotal_time(:) = 0.d0
  DO jfil = 1, nfil 
     cf_in          = cf_list(jfil)
     tim(:)         = getvar1d(cf_in, cn_vtimec, npt)
     dtotal_time(:) = dtotal_time(:) + tim (:)
  END DO
  tim(:) = dtotal_time(:)/ nfil
  ierr   = putvar1d(ncout,  tim, npt, 'T')
  ierr   = putvar1d(ncout2, tim, npt, 'T')

  DO jrec = 1, npt

     DO jvar = 1,nvars
        IF ( cv_nam(jvar) == cn_vlon2d .OR. &     ! nav_lon
             cv_nam(jvar) == cn_vlat2d ) THEN     ! nav_lat
           ! skip these variable
        ELSE
           PRINT *,' Working with ', TRIM(cv_nam(jvar)), ipk(jvar)
           DO jk = 1, ipk(jvar)
              PRINT *,'level ',jk
              dtab(:,:) = 0.d0 ; dtab2(:,:) = 0.d0 
              ntframe = 0
              DO jfil = 1, nfil
                 cf_in     = cf_list(jfil)
                 v2d(:,:)  = getvar(cf_in, cv_nam(jvar), jk, npiglo, npjglo, ktime=jrec )
                 IF ( lspval0 )  WHERE (v2d == zspval_in(jvar))  v2d = 0.  ! change missing values to 0
                 dtab(:,:) = dtab(:,:) + v2d(:,:)
                 IF (cv_nam2(jvar) /= 'none' ) dtab2(:,:) = dtab2(:,:) + v2d(:,:)*v2d(:,:)
              END DO

              ! store variable on outputfile
              ierr = putvar(ncout, id_varout(jvar), SNGL(dtab(:,:)/nfil), jk, npiglo, npjglo, kwght=nfil, ktime = jrec )
              IF (cv_nam2(jvar) /= 'none' )  THEN 
                 ierr = putvar(ncout2, id_varout2(jvar), SNGL(dtab2(:,:)/nfil), jk, npiglo, npjglo, kwght=nfil, ktime=jrec)
              ENDIF

           END DO  ! loop to next level
        END IF
     END DO ! loop to next var in file
  END DO ! loop to next record in input file

  ierr = closeout(ncout)
  ierr = closeout(ncout2)

CONTAINS 

  LOGICAL FUNCTION varchk2 ( cd_var ) 
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION varchk2  ***
    !!
    !! ** Purpose : Return true if cd_var is candidate for mean squared value  
    !!
    !! ** Method  : List of candidate is established in modcdfnames, and
    !!              can be changed via the nam_cdf_names namelist   
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(in) :: cd_var

    INTEGER(KIND=4)              :: jv
    !!----------------------------------------------------------------------
    varchk2 = .FALSE.
    DO jv = 1, nn_sqdvar 
       IF ( cd_var == cn_sqdvar(jv) ) THEN
          varchk2 = .TRUE.
          EXIT
       ENDIF
    ENDDO

  END FUNCTION varchk2

END PROGRAM cdfmoyt
