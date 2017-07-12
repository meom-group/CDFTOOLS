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
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!   varchk2       : check if variable is candidate for square mean
  !!----------------------------------------------------------------------
  USE cdfio 
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class time_averaging
  !!-----------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                               :: jk, jfil, jrec     ! dummy loop index
  INTEGER(KIND=4)                               :: jvar, jv, jt       ! dummy loop index
  INTEGER(KIND=4)                               :: ierr               ! working integer
  INTEGER(KIND=4)                               :: idep, idep_max     ! possible depth index, maximum
  INTEGER(KIND=4)                               :: inpt               ! working integer
  INTEGER(KIND=4)                               :: narg, iargc, ijarg ! browsing command line
  INTEGER(KIND=4)                               :: nfiles             ! number of files to average
  INTEGER(KIND=4)                               :: npiglo, npjglo     ! size of the domain
  INTEGER(KIND=4)                               :: npk, npt           ! size of the domain
  INTEGER(KIND=4)                               :: nvars              ! number of variables in a file
  INTEGER(KIND=4)                               :: ntframe            ! cumul of time frame
  INTEGER(KIND=4)                               :: ncout, ncout2      ! ncid of output files
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: id_var             ! arrays of var id's
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: ipk                ! arrays of vertical level for each var
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: id_varout          ! varid's of average vars
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: id_varout2         ! varid's of sqd average vars

  REAL(KIND=4), DIMENSION(:),       ALLOCATABLE :: zspval_in          ! input missing value
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: v2d                ! array to read a layer of data

  REAL(KIND=8), DIMENSION(:),       ALLOCATABLE :: dtim               ! time counter
  REAL(KIND=8), DIMENSION(:),       ALLOCATABLE :: dtotal_time        ! to compute mean time
  REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: dtab, dtab2        ! arrays for cumulated values

  CHARACTER(LEN=256)                            :: cf_in              ! input file names
  CHARACTER(LEN=256)                            :: cf_root = 'cdfmoyt'! root of the output file name
  CHARACTER(LEN=256)                            :: cf_out             ! output file for average
  CHARACTER(LEN=256)                            :: cf_out2            ! output file for squared average
  CHARACTER(LEN=256)                            :: cv_dep             ! depth dimension name
  CHARACTER(LEN=256)                            :: cldum              ! dummy string argument
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cf_lst             ! list of input files
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_nam             ! array of var name
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cv_nam2            ! array of var2 name for output
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: clv_dep            ! array of possible depth name (or 3rd dimension)

  TYPE (variable), DIMENSION(:),    ALLOCATABLE :: stypvar            ! attributes for average values
  TYPE (variable), DIMENSION(:),    ALLOCATABLE :: stypvar2           ! attributes for square averaged values

  LOGICAL                                       :: lspval0 = .FALSE.  ! cdfmoy_chsp flag
  LOGICAL                                       :: lnc4    = .FALSE.  ! Use nc4 with chunking and deflation
  LOGICAL                                       :: lchk               ! flag for checking presence of files
  !!----------------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfmoyt -l LST-files [-spval0] [-vvl] [-o OUT-rootname] [-nc4]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute the ''file average'' of the files listed on the command line.'
     PRINT *,'       The ''file average'' will have the same number of time frame than any'
     PRINT *,'       individual file in the list, the average being done frame by frame.'
     PRINT *,'       '
     PRINT *,'       The main use of this program is the calculation of a climatological'
     PRINT *,'       average, for instance. It can also be used for the calculation of an'
     PRINT *,'       ensemble mean, although cdfenstat is more appropriate. '
     PRINT *,'       '
     PRINT *,'       The program assumes that all files in the list are of same type (shape,'
     PRINT *,'       variables, and number of time frames).'
     PRINT *,'       '
     PRINT *,'       For some variables, the program also computes the ''file average'' of the'
     PRINT *,'       squared variables, which is used in other cdftools (cdfeke, cdfrmsssh,'
     PRINT *,'       cdfstdevw, cdfstddevts...). The actual variables selected for squared'
     PRINT *,'       average are :'
     PRINT '(10x,"- ",a)' , (TRIM(cn_sqdvar(jv)), jv=1, nn_sqdvar)
     PRINT *,'       This selection can be adapted with the nam_cdf_namelist process.'
     PRINT *,'       (See cdfnamelist -i for details).'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -l LST-files: List of files whose average will be computed.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-spval0 ] : set missing_value attribute to 0 for all variables and'
     PRINT *,'              take care of the input missing_value. This option is usefull'
     PRINT *,'              if missing_values differ from files to files.'
     PRINT *,'       [-vvl] : Use time-varying vertical metrics.'
     PRINT *,'       [-o OUT-rootname] : Define output root-name instead of ', TRIM(cf_root)
     PRINT *,'       [-nc4 ]: Use netcdf4 output with chunking and deflation level 1..'
     PRINT *,'              This option is effective only if cdftools are compiled with'
     PRINT *,'              a netcdf library supporting chunking and deflation.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       none '
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_root),'.nc and ',TRIM(cf_root),'2.nc'
     PRINT *,'       Variables name are the same than in the input files. '
     PRINT *,'       For squared averages ''_sqd'' is appended to the original variable name.'
     STOP 
  ENDIF

  ijarg = 1 
  DO WHILE ( ijarg <= narg )
     CALL getarg (ijarg, cldum) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-l'      ) ; CALL GetFileList
        ! options
     CASE ( '-spval0' ) ; lspval0 = .TRUE.
     CASE ( '-o'      ) ; CALL getarg(ijarg, cf_root) ; ijarg=ijarg+1
     CASE ( '-nc4'    ) ; lnc4    = .TRUE.
     CASE ( '-vvl'    ) ; lg_vvl  = .TRUE. ; PRINT *, ' vvl not yet supported' ; STOP 99
     CASE DEFAULT       ; PRINT *,' ERROR : ', TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  END DO

  cf_out  = TRIM(cf_root)//'.nc'
  cf_out2 = TRIM(cf_root)//'2.nc'

  lchk = .FALSE.
  DO jfil=1,nfiles
     lchk = lchk .OR. chkfile (cf_lst(jfil) )
  ENDDO
  IF ( lchk ) STOP 99 ! missing file

  ! Initialisation from  1rst file (all file are assumed to have the same geometry)
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
     IF (  chkfile (cf_lst(jfil)      ) ) STOP 99 ! missing file
     inpt = getdim (cf_lst(jfil), cn_t)
     IF ( inpt /= npt ) THEN
        PRINT *, 'File ',TRIM(cf_lst(jfil) ),' has ',inpt,' time frames instead of ', npt
        ierr = ierr + 1
     ENDIF
  ENDDO
  IF ( ierr /= 0 ) STOP 99 ! frame numbers mismatch

  PRINT *, 'npiglo = ', npiglo
  PRINT *, 'npjglo = ', npjglo
  PRINT *, 'npk    = ', npk
  PRINT *, 'npt    = ', npt

  nvars = getnvar(cf_in)
  PRINT *,' nvars = ', nvars

  ALLOCATE( dtab(npiglo,npjglo), dtab2(npiglo,npjglo), v2d(npiglo,npjglo) )
  ALLOCATE( dtotal_time(npt), dtim(npt) )
  ALLOCATE (cv_nam(nvars), cv_nam2(nvars) )
  ALLOCATE (stypvar(nvars), stypvar2(nvars) )
  ALLOCATE (id_var(nvars), ipk(nvars), id_varout(nvars), id_varout2(nvars)  )

  CALL CreateOutput

  ! Compute the mean time for each mean frame
  dtotal_time(:) = 0.d0
  DO jfil = 1, nfiles 
     cf_in          = cf_lst(jfil)
     dtim(:)        = getvar1d(cf_in, cn_vtimec, npt)
     dtotal_time(:) = dtotal_time(:) + dtim (:)
  END DO
  dtim(:)= dtotal_time(:)/ nfiles
  ierr   = putvar1d(ncout,  dtim, npt, 'T')
  ierr   = putvar1d(ncout2, dtim, npt, 'T')

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
              DO jfil = 1, nfiles
                 cf_in     = cf_lst(jfil)
                 v2d(:,:)  = getvar(cf_in, cv_nam(jvar), jk, npiglo, npjglo, ktime=jrec )
                 IF ( lspval0 )  WHERE (v2d == zspval_in(jvar))  v2d = 0.  ! change missing values to 0
                 dtab(:,:) = dtab(:,:) + v2d(:,:)
                 IF (cv_nam2(jvar) /= 'none' ) dtab2(:,:) = dtab2(:,:) + v2d(:,:)*v2d(:,:)
              END DO

              ! store variable on outputfile
              ierr = putvar(ncout, id_varout(jvar), SNGL(dtab(:,:)/nfiles), jk, npiglo, npjglo, kwght=nfiles, ktime = jrec )
              IF (cv_nam2(jvar) /= 'none' )  THEN 
                 ierr = putvar(ncout2, id_varout2(jvar), SNGL(dtab2(:,:)/nfiles), jk, npiglo, npjglo, kwght=nfiles, ktime=jrec)
              ENDIF

           END DO  ! loop to next level
        END IF
     END DO ! loop to next var in file
  END DO ! loop to next record in input file

  ierr = closeout(ncout )
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
    cv_nam(:) = getvarname(cf_in,nvars, stypvar)
    DO jvar = 1,nvars
       stypvar(jvar)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    ENDDO

    IF ( lspval0 ) THEN 
       ALLOCATE ( zspval_in(nvars) )
       zspval_in(:) = stypvar(:)%rmissing_value
       stypvar(:)%rmissing_value = 0.
    ENDIF

    DO jvar = 1, nvars
       ! variables that will not be computed or stored are named 'none'
       IF ( varchk2 ( cv_nam(jvar) ) ) THEN 
          cv_nam2(jvar)                    = TRIM(cv_nam(jvar))//'_sqd'
          stypvar2(jvar)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
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
    ncout  = create      (cf_out,  cf_in,    npiglo, npjglo, npk, cdep=cv_dep, ld_nc4=lnc4 )
    ierr   = createvar   (ncout ,  stypvar,  nvars,  ipk,    id_varout       , ld_nc4=lnc4 )
    ierr   = putheadervar(ncout,   cf_in,    npiglo, npjglo, npk, cdep=cv_dep)

    ncout2 = create      (cf_out2, cf_in,    npiglo, npjglo, npk, cdep=cv_dep, ld_nc4=lnc4 )
    ierr   = createvar   (ncout2,  stypvar2, nvars,  ipk,    id_varout2      , ld_nc4=lnc4 )
    ierr   = putheadervar(ncout2,  cf_in,    npiglo, npjglo, npk, cdep=cv_dep)

  END SUBROUTINE CreateOutput

END PROGRAM cdfmoyt
