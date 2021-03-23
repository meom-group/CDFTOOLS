PROGRAM cdficb_clv
  !!======================================================================
  !!                     ***  PROGRAM  cdficb_clv  ***
  !!=====================================================================
  !!  ** Purpose : Prepare netcdf iceberg calving file
  !!
  !!  ** Method  : generate random disctribution along ice shelf front 
  !!               then scale it
  !!
  !! History : 4.0  : 04/2020  : P.Mathiot 
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  USE cdfio 
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id: cdficb_rnf.f90 668 2013-05-30 12:54:00Z molines $
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class ice_shelf_processing
  !!-----------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                               :: jicb, ji, jj       ! dummy loop index
  INTEGER(KIND=4)                               :: ierr               ! working integer
  INTEGER(KIND=4)                               :: narg, iargc, ijarg ! browsing command line
  INTEGER(KIND=4)                               :: npiglo, npjglo, npk! size of the domain
  INTEGER(KIND=4)                               :: nicb               ! number of ice shelves
  INTEGER(KIND=4)                               :: iunit=10           ! id file
  INTEGER(KIND=4)                               :: ncout              ! ncid of output files
  INTEGER(KIND=4)                               :: ifill
  INTEGER(KIND=4)                               :: iiseed, ijseed     ! position of a point within ice shelf (not used here)
  INTEGER(KIND=4)                               :: ijmin, ijmax       ! j-limits of a particular ice shelf
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: ipk                ! arrays of vertical level for each var
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: id_varout          ! varid's of average vars
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: isum               ! zonal sum of icb index for optim
  INTEGER(KIND=2), DIMENSION(:,:),  ALLOCATABLE :: icbindex, icbmask  ! 
  INTEGER(KIND=2), DIMENSION(:,:),  ALLOCATABLE :: icbindex_wk        !
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: iseed
  INTEGER(KIND=4)                               :: nseed              ! size of seed array
  INTEGER(KIND=4)                               :: irdsf              ! seed scale factor

  REAL(KIND=4)                                  :: rdraftmax, rdraftmin ! dummy information in input file
  REAL(KIND=4)                                  :: rlon, rlat         ! dummy information in input file
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: bathy, zdrft       ! bathymetry and ice shelf draft

  REAL(KIND=8)                                  :: dfwf, dclv, dclv_tot
  REAL(KIND=8)                                  :: dsumcoef
  REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: dl_icbclv2d
  REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: dicbclv2d

  !                                                FILES
  CHARACTER(LEN=256)                            :: cf_fill            ! input file names
  CHARACTER(LEN=256)                            :: cf_icblist         ! input file names
  CHARACTER(LEN=256)                            :: cf_bathy='bathy.nc'! bathymetry file name
  CHARACTER(LEN=256)                            :: cf_isfdr='isf_draft.nc'! ice_draft file name
  CHARACTER(LEN=256)                            :: cf_out='rnficb.nc' ! output file for average
  !                                                VARIABLES
  CHARACTER(LEN=256)                            :: cv_fill            ! fill var name
  CHARACTER(LEN=256)                            :: cv_bathy='Bathymetry' ! bathymetry name
  CHARACTER(LEN=256)                            :: cv_isfdr='isf_draft'  ! ice shelf draft name
  CHARACTER(LEN=256)                            :: cldum              ! dummy string argument

  TYPE (variable), DIMENSION(:),    ALLOCATABLE :: stypvar            ! attributes for average values

  LOGICAL                                       :: lchk    = .FALSE.  ! flag for missing files
  LOGICAL                                       :: lnc4    = .FALSE.  ! flag for netcdf4 chunking and deflation
  LOGICAL                                       :: lperio  = .FALSE.  !
  LOGICAL                                       :: ltot    = .FALSE.  !
  !!----------------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()

  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdficb_clv -f ISF-fill-file -v ISF-fill-var -l ISF-listfile '
     PRINT *,'    [-b BATHY-file] [-vb BATHY-var] [-i ISFDRAFT-file] [-vi ISFDRAFT-variable]'
     PRINT *,'    [-s irdsf] [-st] [-nc4] [-o OUT-file ]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'        Build a netcdf calving file, for use in NEMO ICB module.'
     PRINT *,'      '
     PRINT *,'        This file hold the 2D  calving rate field (km3/year) at points located'
     PRINT *,'        along the ice shelves edge. In this program, all points along a given'
     PRINT *,'        iceshelf contribute to the calving.  The total calving rate for a given' 
     PRINT *,'        iceshelf is randomly spread over the calving points.'
     PRINT *,'      '
     PRINT *,'        Optionally, [-st], the overall calving rate can be scaled to match'
     PRINT *,'        observational estimates.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'        -f ISF-fill-file : file built by cdffill (all the ice shelves are'
     PRINT *,'                           tagged with an id)'
     PRINT *,'        -v ISF-fill_var  : name of fill variable to use in ISF-fill_file'
     PRINT *,'        -l ISF-list : Text file with the calving rate (GT/y) given for'
     PRINT *,'                      each ice shelf.' 
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'        -s irdsf : seed random scale factor (seed = irdsf * ifill) '
     PRINT *,'               This is a way to change the seed for the random number generator'
     PRINT *,'               [ default : 1 ]'
     PRINT *,'        -b BATHY-file : give name of bathy file.'
     PRINT *,'                    [ default : ',TRIM(cf_bathy),' ]'
     PRINT *,'        -vb BATHY-var : give name of bathy variable.'
     PRINT *,'                    [ default : ',TRIM(cv_bathy),' ]'
     PRINT *,'        -i ISFDRAFT-file : give name of isf_draft file.'
     PRINT *,'                      [ default : ',TRIM(cf_isfdr),' ]'
     PRINT *,'        -vi ISFDRAFT-var : give name of isf_draft variable.'
     PRINT *,'                      [ default : ',TRIM(cv_isfdr),' ]'
     PRINT *,'        -st : scale to total iceberg calving all at the end'
     PRINT *,'                to include calving from front unresolved ice shelves'
     PRINT *,'        -nc4 : Use this option to have netcdf4 output file, with chunking'
     PRINT *,'               and deflation.'
     PRINT *,'        -o OUT-file : Specify the name of the output file instead of '
     PRINT *,'               the default name ',TRIM(cf_out)
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       ',TRIM(cn_fhgr),' and all files specified on the command line' 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) ,' unless -o option used'
     PRINT *,'         variables : soicbclv (km3/y)'
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'       cdfisb_fill'
     PRINT *,'      '
     STOP 
  ENDIF

  ijarg=1
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg,cldum) ; ijarg=ijarg+1
     SELECT CASE (cldum)
     CASE ('-f' ) ; CALL getarg(ijarg,cf_fill   ) ; ijarg=ijarg+1
     CASE ('-v' ) ; CALL getarg(ijarg,cv_fill   ) ; ijarg=ijarg+1
     CASE ('-l' ) ; CALL getarg(ijarg,cf_icblist) ; ijarg=ijarg+1
     CASE ('-s' ) ; CALL getarg(ijarg,cldum     ) ; ijarg=ijarg+1 ; READ(cldum,*) irdsf
     CASE ('-b' ) ; CALL getarg(ijarg,cf_bathy  ) ; ijarg=ijarg+1
     CASE ('-vb') ; CALL getarg(ijarg,cv_bathy  ) ; ijarg=ijarg+1
     CASE ('-i' ) ; CALL getarg(ijarg,cf_isfdr  ) ; ijarg=ijarg+1
     CASE ('-vi') ; CALL getarg(ijarg,cv_isfdr  ) ; ijarg=ijarg+1
     CASE ('-nc4'); lnc4   = .TRUE.
     CASE ('-ew' ); lperio = .TRUE.
     CASE ('-st' ); ltot   = .TRUE.
     CASE ('-o' ) ; CALL getarg(ijarg,cf_out    ) ; ijarg=ijarg+1
     CASE DEFAULT ; PRINT *,' ERROR : ',TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  ENDDO

  lchk = lchk .OR. chkfile (cf_fill   )
  lchk = lchk .OR. chkfile (cf_icblist)
  lchk = lchk .OR. chkfile (cf_bathy  )
  lchk = lchk .OR. chkfile (cf_isfdr  )
  lchk = lchk .OR. chkfile (cn_fhgr   )

  IF ( lchk ) STOP 99 ! missing file

  npiglo = getdim (cf_fill, cn_x)
  npjglo = getdim (cf_fill, cn_y)
  npk    = 0                  ! bathy file has no dep dimension

  PRINT *, 'NPIGLO = ', npiglo
  PRINT *, 'NPJGLO = ', npjglo
  PRINT *, 'NPK    = ', npk

  ALLOCATE(icbindex (npiglo, npjglo), icbmask    (npiglo, npjglo))
  ALLOCATE(icbindex_wk(npiglo, npjglo), isum(npjglo))

  ALLOCATE(bathy    (npiglo, npjglo), zdrft      (npiglo, npjglo))

  ALLOCATE(dicbclv2d(npiglo, npjglo), dl_icbclv2d(npiglo, npjglo))

  ALLOCATE (stypvar(1))
  ALLOCATE (ipk(1),id_varout(1))

  CALL CreateOutput

  ! define variable
  ! read ice shelf draft data
  icbindex(:,:) = getvar(cf_fill,  cv_fill,  1 ,npiglo, npjglo )  ! fill index
  bathy(:,:)    = getvar(cf_bathy, cv_bathy, 1, npiglo, npjglo )  ! ocean bathy
  zdrft(:,:)    = getvar(cf_isfdr, cv_isfdr, 1, npiglo, npjglo )  ! ice shelf draft

!  ! sanity check
!  zdrft = MIN(zdrft,bathy)

  ! open icb file
  OPEN(unit=iunit, file=cf_icblist, form='formatted', status='old')
  ! get number of icb
  nicb = 0
  cldum='XXX'
  DO WHILE ( TRIM(cldum) /= 'EOF')
     READ(iunit,*) cldum
     nicb=nicb+1
  END DO
  REWIND(iunit)
  nicb = nicb - 1

  PRINT *, '   Number of ISF found in file list : ', nicb

  CALL RANDOM_SEED(SIZE=nseed)
  ALLOCATE(iseed(nseed))
  dicbclv2d = 0.0
  DO jicb=1,nicb-1
     ! reset working icb index to its initial value
     icbindex_wk(:,:) = icbindex(:,:)

     ! read ice shelf data for jsf
     READ(iunit,*) ifill,cldum,rlon, rlat, iiseed, ijseed ,rdraftmin, rdraftmax, dfwf, dclv

     icbmask  (:,:) = 0
     dl_icbclv2d(:,:) = 0.0d0
     dsumcoef       = 0.0d0

     ! only deal with current ice shelf
     WHERE (icbindex_wk /=  -ifill ) icbindex_wk = 0

     ! find the j-limits for the current ice shelf 
     isum(:)=SUM(icbindex_wk,dim=1)
     ijmin = npjglo ; ijmax=2
     DO jj=npjglo, 3, -1
        IF ( isum(jj) /= 0 ) THEN
           ijmin=jj
        ENDIF
     ENDDO
     DO jj=1, npjglo-2
        IF ( isum(jj) /= 0 ) THEN
           ijmax=jj
        ENDIF
     ENDDO

     iseed = ifill * irdsf
     CALL RANDOM_SEED(PUT=iseed)
     DO ji=2,npiglo-1
        DO jj = ijmin-1, ijmax+1
           IF ( zdrft(ji,jj) == 0 .AND.  &                                          ! not under ice_shelf
                &    bathy(ji,jj) /= 0 .AND.  &                                     ! but in the ocean
                &    MINVAL(icbindex_wk(ji-1:ji+1 , jj-1:jj+1)) == -ifill  .AND. &  ! adjacent to the correct icb
                &    icbindex_wk(ji,jj) == 0                               .AND. &  ! 
                &    icbmask(ji,jj)     == 0       ) THEN                           ! not filled yet
                ! case ice shelf in pultiple part in txt list
                ! WHY, should not be needed dicbclv2d(ji,jj) = 0.0
                CALL RANDOM_NUMBER(dl_icbclv2d(ji,jj))
                dsumcoef = dsumcoef + dl_icbclv2d(ji,jj)
                icbmask(ji,jj) = 1
           END IF
        END DO
     END DO

     IF ( SUM(icbmask) == 0 .AND. ifill < 99) THEN
        PRINT *, 'E R R O R : this ice shelf ',ifill, ' do not have sea access and so cannot calve'
        STOP
     END IF

     WHERE (icbmask == 1)
        dl_icbclv2d = dl_icbclv2d / dsumcoef * dclv
     END WHERE

     PRINT *, TRIM(cldum),' calving is : ',SUM(dl_icbclv2d)

     dicbclv2d = dicbclv2d + dl_icbclv2d

  END DO

  ! scale to the total ice berg calving
  IF ( ltot ) THEN
     READ(iunit,*) ifill,cldum,rlon, rlat, iiseed, ijseed ,rdraftmin, rdraftmax, dfwf, dclv_tot
     IF (ifill /= 99) THEN
        PRINT *, 'last iceberg line must be -99 TOTA ....'
        STOP 42
     END IF
     dicbclv2d = dicbclv2d * dclv_tot / SUM(dicbclv2d)     
  END IF

  ! Diagnose total amount of fwf
  PRINT *, 'diag : '
  PRINT *, 'total sum of icb (Gt/y) = ', SUM(dicbclv2d(2:npiglo-1,2:npjglo-1))

  ! convertion of Gt/y in km3/y (icb density 850 (NEMO default))
  dicbclv2d = dicbclv2d*1.d+3/850.d0
  IF ( lperio ) THEN 
     dicbclv2d(npiglo,:) = dicbclv2d(2       ,:)
     dicbclv2d(1     ,:) = dicbclv2d(npiglo-1,:)
  ENDIF

  ierr = putvar(ncout, id_varout(1), dicbclv2d, 1, npiglo, npjglo)

  ierr = closeout(ncout)

CONTAINS

  SUBROUTINE CreateOutput 
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutput  ***
    !!
    !! ** Purpose :  Create output netcdf output file using cdfio 
    !!              We use a routine just to increase readability
    !! ** Method  :  Use global variables to know about the file to be created
    !!           
    !!----------------------------------------------------------------------
    REAL(KIND=8), DIMENSION(1) :: dl_tim
    !!----------------------------------------------------------------------

    ! define new variables for output
    stypvar(1)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
    stypvar(1)%cname             = 'soicbclv'
    stypvar(1)%rmissing_value    =  -99.d0
    stypvar(1)%valid_min         =  0.
    stypvar(1)%valid_max         =  2000.
    stypvar(1)%clong_name        = 'icberg calving'
    stypvar(1)%cshort_name       = 'soicbclv'
    stypvar(1)%conline_operation = 'N/A'
    stypvar(1)%caxis             = 'TYX'
    stypvar(1)%cprecision        = 'r8'
    stypvar(1)%cunits             = 'km3/y'
    ipk(1) = 1  !  2D

    ! create output file taking the sizes in cf_fill
    ncout  = create      (cf_out, cf_fill,   npiglo, npjglo, npk,  ld_nc4=lnc4)
    ierr   = createvar   (ncout,  stypvar, 1,   ipk, id_varout  ,  ld_nc4=lnc4)
    ierr   = putheadervar(ncout,  cf_fill,   npiglo, npjglo, npk              )

    dl_tim(1) = 0.d0
    ierr = putvar1d(ncout, dl_tim, 1, 'T') 

  END SUBROUTINE CreateOutput

END PROGRAM cdficb_clv
