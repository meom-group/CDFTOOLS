PROGRAM cdfisf_forcing
  !!======================================================================
  !!                     ***  PROGRAM  cdfisf_forcing  ***
  !!=====================================================================
  !!  ** Purpose : spread a specified total ice shelf melting over a specific
  !!               ice shelf melting pattern
  !!
  !!  ** Method  : read an ice shelf mask file produce by cdffill, read the
  !!               integrate melting for each ice shelf and the wanted pattern,
  !!               then compute the final melting for each ice shelf.
  !!
  !! History : 3.0  : 04/2014  : Pierre Mathiot 
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  USE cdfio 
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class ice_shelf_processing
  !!-----------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                               :: jisf               ! dummy loop counter
  INTEGER(KIND=4)                               :: ierr               ! working integer
  INTEGER(KIND=4)                               :: narg, iargc, ijarg ! browsing command line
  INTEGER(KIND=4)                               :: npiglo, npjglo     ! size of the domain
  INTEGER(KIND=4)                               :: npk,      nisf     ! size of the domain
  INTEGER(KIND=4)                               :: iunit=10           ! id file
  INTEGER(KIND=4)                               :: ncout              ! ncid of output files
  INTEGER(KIND=4)                               :: iiseed, ijseed
  INTEGER(KIND=4)                               :: ifill
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: ipk                ! arrays of vertical level for each var
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: id_varout          ! varid's of average vars

  INTEGER(KIND=2), DIMENSION(:,:),  ALLOCATABLE :: ipoolmsk           ! mask for closed pool
  INTEGER(KIND=2), DIMENSION(:,:),  ALLOCATABLE :: isfindex           ! index of each ISF ( negative integer)
  INTEGER(KIND=2), DIMENSION(:,:),  ALLOCATABLE :: isfindex_wk        ! index of each ISF 'working variable)

  REAL(KIND=4)                                  :: rdraftmax, rdraftmin ! dummy information in input file
  REAL(KIND=4)                                  :: rlon, rlat         ! dummy information in input file

  REAL(KIND=8)                                  :: dl_fwf, dsumcoef
  REAL(KIND=8)                                  :: dfwf
  REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: dfwfisf2d
  REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: de12t
  REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: dl_fwfisf2d, dl_fwfispat

  CHARACTER(LEN=256)                            :: cf_fill            ! input file names
  CHARACTER(LEN=256)                            :: cf_isflist         ! input file names
  CHARACTER(LEN=256)                            :: cf_out='isfforcing.nc' ! output file for average
  CHARACTER(LEN=256)                            :: cf_pat='isfpattern.nc' ! pattern file
  CHARACTER(LEN=256)                            :: cf_pool='isfpool.nc'   ! pools mask file
  CHARACTER(LEN=256)                            :: cv_dep             ! depth dimension name
  CHARACTER(LEN=256)                            :: cv_pat='sowflisf'  ! pattern variable name
  CHARACTER(LEN=256)                            :: cv_pool='isfpoolmask'! pattern variable name
  CHARACTER(LEN=256)                            :: cv_fill            ! isf index variable in cf_fill
  CHARACTER(LEN=256)                            :: cldum              ! dummy string argument

  TYPE (variable), DIMENSION(:),    ALLOCATABLE :: stypvar            ! attributes for average values

  LOGICAL                                       :: lnc4 = .FALSE.     ! flag for netcdf4 chunking and deflation
  LOGICAL                                       :: lchk = .FALSE.     ! flag for missing values

  !!----------------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfisf_forcing -f ISF-fill_file  -v ISF-fill_var -l ISF-listfile '
     PRINT *,'             -m ISF-poolmask [-vm ISF-poolmask_variable] [-p PATTERN-file] '
     PRINT *,'            [-vp PATTERN-variable] [-nc4] [-o OUT-file ]'
     PRINT *,'      '
     PRINT *,'     PURPOSE : '
     PRINT *,'         Build basal melting rate file used in NEMO ISF when nn_isf=4 '
     PRINT *,'      '
     PRINT *,'     ARGUMENTS : '
     PRINT *,'          -f ISF-fill_file : file built by cdfisf_fill (all the ice shelves '
     PRINT *,'                             are tagged with an id)'
     PRINT *,'          -v ISF-fill_var  : name of fill variable to use in ISF-fill_file'
     PRINT *,'          -l ISF-listfile : text file used to build the ISF-fill_file. '
     PRINT *,'                            Only the last variable on each line is used (GT/y)'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'          -p PATTERN-file : specify the file use for patterns. '
     PRINT *,'                            [ default : ',TRIM(cf_pat),' ]'
     PRINT *,'          -vp PATTERN-variable : specify the name of the pattern variable. '
     PRINT *,'                            [ default : ',TRIM(cv_pat),' ]'
     PRINT *,'          -vm ISF-poolmask_variable : specify the name of the variable used '
     PRINT *,'                 for masking the pools. [ default : ',TRIM(cv_pool),' ]'
     PRINT *,'          -nc4 : use netcdf4 chunking and deflation'
     PRINT *,'          -o OUT-file : specify output filename. [ default : ',TRIM(cf_out),' ]'
     PRINT *,'              '
     PRINT *,'     REQUIRED FILES : '
     PRINT *,'           mesh_zgr.nc mesh_hgr.nc,'
     PRINT *,'           isfpattern.nc (ie reference file used to define the isf melting '
     PRINT *,'                 pattern), unless -p option is used to give different name.'
     PRINT *,'      '
     PRINT *,'     OUTPUT :'
     PRINT *,'         netcdf file : ', TRIM(cf_out),' unless specified with -o option'
     PRINT *,'         variable : sofwfisf '
     PRINT *,'      '
     PRINT *,'     SEE ALSO : cdfisf_fill, cdfisf_rnf, cdfisf_poolchk'
     PRINT *,'      '
     STOP 
  ENDIF

  ijarg = 1
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg = ijarg + 1
     SELECT CASE ( cldum)
     CASE ( '-f' ) ; CALL getarg(ijarg, cf_fill    ) ; ijarg = ijarg + 1
     CASE ( '-v' ) ; CALL getarg(ijarg, cv_fill    ) ; ijarg = ijarg + 1
     CASE ( '-l' ) ; CALL getarg(ijarg, cf_isflist ) ; ijarg = ijarg + 1
     CASE ( '-o' ) ; CALL getarg(ijarg, cf_out     ) ; ijarg = ijarg + 1
     CASE ( '-p' ) ; CALL getarg(ijarg, cf_pat     ) ; ijarg = ijarg + 1
     CASE ( '-vp') ; CALL getarg(ijarg, cv_pat     ) ; ijarg = ijarg + 1
     CASE ( '-m' ) ; CALL getarg(ijarg, cf_pool    ) ; ijarg = ijarg + 1
     CASE ( '-vm') ; CALL getarg(ijarg, cv_pool    ) ; ijarg = ijarg + 1
     CASE ('-nc4') ; lnc4 = .TRUE.
     CASE DEFAULT  ; PRINT *,' ERROR : ', TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  ENDDO

  lchk = lchk .OR. chkfile (cf_fill   )
  lchk = lchk .OR. chkfile (cf_isflist)
  lchk = lchk .OR. chkfile (cf_pat    )
  lchk = lchk .OR. chkfile (cf_pool   )
  lchk = lchk .OR. chkfile (cn_fzgr   )
  lchk = lchk .OR. chkfile (cn_fhgr   )
  IF ( lchk  ) STOP 99 ! missing file

  npiglo = getdim (cf_fill, cn_x)
  npjglo = getdim (cf_fill, cn_y)
  npk    = 0  ! working with 2D files only

  PRINT *, 'NPIGLO = ', npiglo
  PRINT *, 'NPJGLO = ', npjglo
  PRINT *, 'NPK    = ', npk

  ALLOCATE(de12t(npiglo,npjglo))
  ALLOCATE(ipoolmsk(npiglo, npjglo), dfwfisf2d(npiglo, npjglo)    )
  ALLOCATE(dl_fwfisf2d(npiglo,npjglo), dl_fwfispat(npiglo,npjglo))
  ALLOCATE(isfindex(npiglo, npjglo), isfindex_wk(npiglo, npjglo) )

  ALLOCATE (stypvar(1))
  ALLOCATE (ipk(1),id_varout(1))
  ! initialisation of final fwf
  dfwfisf2d(:,:) = 0.0d0

  CALL CreateOutput

  ! define variable
  ! read ice shelf draft data
  de12t(:,:)    = getvar(cn_fhgr, cn_ve1t,  1, npiglo, npjglo ) *  getvar(cn_fhgr, cn_ve2t,  1, npiglo, npjglo )

  ! open isf file
  OPEN(unit=iunit, file=cf_isflist, form='formatted', status='old')
  ! get number of isf
  nisf = 0
  cldum='XXX'
  DO WHILE ( TRIM(cldum) /= 'EOF')
     READ(iunit,*) cldum
     nisf=nisf+1
  END DO
  REWIND(iunit)
  nisf = nisf - 1

  PRINT *, '   Number of ISF found in file list : ', nisf

  ! Read the Basal melting pattern, once for all
  dl_fwfispat(:,:) = getvar(cf_pat , cv_pat, 1 ,npiglo, npjglo )
  isfindex(:,:)    = getvar(cf_fill, cv_fill,1 ,npiglo, npjglo )  
  ipoolmsk(:,:)    = getvar(cf_pool, cv_pool,1 ,npiglo, npjglo )

  ! loop over all the ice shelf
  DO jisf=1,nisf
     ! initialize working pattern with the fixed one
     dl_fwfisf2d = dl_fwfispat

     ! reset working isf index to its initial value
     isfindex_wk(:,:) = isfindex(:,:)

     ! eliminate closed pools from isfindex, using ISF-pool file
     isfindex_wk = isfindex_wk * ipoolmsk

     ! read ice shelf data
     READ(iunit,*) ifill,cldum,rlon, rlat, iiseed, ijseed ,rdraftmin, rdraftmax, dfwf

     ! convertion of total ice shelf melting from Gt/y -> kg/s
     dl_fwf = dfwf * 1.d9 * 1.d3 / 86400.d0 / 365.d0

     ! initialisation of cumulative variable
     dsumcoef = 0.0d0

     ! isolate the ice shelf data we want
     WHERE (isfindex_wk /= -ifill)
        isfindex_wk(:,:) = 0      ! eliminate all ISF not current (-ifill)
        dl_fwfisf2d(:,:) = 0.0d0
     END WHERE

     ! set the halo to 0 (to avoid double counting)  ( E-W periodicity !)
     dl_fwfisf2d(1,:)=0.0d0 ; dl_fwfisf2d(npiglo,:)=0.0d0 ;

     dsumcoef = SUM(dl_fwfisf2d * de12t)
     dl_fwfisf2d(:,:) = dl_fwfisf2d(:,:) / dsumcoef * dl_fwf

     ! Value read from the text file has the wrong sign for melting.
     ! As the shelves are disjoint, cumulate is OK !
     dfwfisf2d(:,:) = dfwfisf2d(:,:) - dl_fwfisf2d(:,:)
  END DO

  ! print isf forcing file
  ierr = putvar(ncout, id_varout(1), dfwfisf2d, 1, npiglo, npjglo)

  ! close file
  ierr = closeout(ncout)

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
    REAL(KIND=8), DIMENSION(1) :: dl_tim
    !!----------------------------------------------------------------------
    ! define new variables for output
    ipk(1) = 1  !  2D
    stypvar(1)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(1)%cname             = 'sofwfisf'
    stypvar(1)%cunits            = 'kg/s'
    stypvar(1)%rmissing_value    =  -99.d0
    stypvar(1)%valid_min         =  0.d0
    stypvar(1)%valid_max         =  2000.d0
    stypvar(1)%clong_name        = 'Ice Shelf Fresh Water Flux '
    stypvar(1)%cshort_name       = 'sofwfisf'
    stypvar(1)%conline_operation = 'N/A'
    stypvar(1)%caxis             = 'TYX'
    stypvar(1)%cprecision        = 'r8'

    ! create output file taking the sizes in cf_fill
    ncout  = create      (cf_out, cf_fill, npiglo, npjglo, npk,   ld_nc4=lnc4 )
    ierr   = createvar   (ncout , stypvar, 1,  ipk,  id_varout,   ld_nc4=lnc4 )
    ierr   = putheadervar(ncout,  cn_fzgr, npiglo, npjglo, npk                )
   
    dl_tim(1) = 0.d0
    ierr = putvar1d( ncout, dl_tim,1 , 'T')

  END SUBROUTINE CreateOutput

END PROGRAM cdfisf_forcing
