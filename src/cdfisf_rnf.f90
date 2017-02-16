PROGRAM cdfisf_rnf
  !!======================================================================
  !!                     ***  PROGRAM  cdfisf_rnf  ***
  !!=====================================================================
  !!  ** Purpose : Prepare netcdf file with Ice Shelf melting parametrized
  !!               as runoff (spread on the vertical)
  !!
  !!  ** Method  : 
  !!
  !! History : 3.0  : 04/2014  : P.Mathiot 
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  USE cdfio 
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2016
  !! $Id: cdfisf_rnf.f90 668 2013-05-30 12:54:00Z molines $
  !! Copyright (c) 2016, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!-----------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                               :: jisf, ji, jj, jw   ! dummy loop index
  INTEGER(KIND=4)                               :: ierr               ! working integer
  INTEGER(KIND=4)                               :: narg, iargc, ijarg ! browsing command line
  INTEGER(KIND=4)                               :: npiglo, npjglo, npk! size of the domain
  INTEGER(KIND=4)                               :: nisf, nwidth       ! number of ice shelves
  INTEGER(KIND=4)                               :: iunit=10           ! id file
  INTEGER(KIND=4)                               :: ncout              ! ncid of output files
  INTEGER(KIND=4)                               :: ifill
  INTEGER(KIND=4)                               :: iiseed, ijseed     ! position of a point within ice shelf
  INTEGER(KIND=4)                               :: ijmin, ijmax       ! j-limits of a particular ice shelf
  INTEGER(KIND=4)                               :: iwscale=75         ! horizontal scale of the fading function (points)
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: ipk                ! arrays of vertical level for each var
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: id_varout          ! varid's of average vars
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: isum               ! zonal sum of isf index for optim
  INTEGER(KIND=2), DIMENSION(:,:),  ALLOCATABLE :: isfindex, isfmask  ! 
  INTEGER(KIND=2), DIMENSION(:,:),  ALLOCATABLE :: isfindex_wk        !

  REAL(KIND=4)                                  :: rdraftmax, rdraftmin ! dummy information in input file
  REAL(KIND=4)                                  :: rlon, rlat         ! dummy information in input file
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: bathy, zdrft       ! bathymetry and ice shelf draft
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE :: zmax2d, zmin2d     ! output arrays

  REAL(KIND=8)                                  :: dfwf
  REAL(KIND=8)                                  :: dl_fwf, dsumcoef
  REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: de12t
  REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: dl_fwfisf2d
  REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: dfwfisf2d

  !                                                FILES
  CHARACTER(LEN=2048)                            :: cf_fill            ! input file names
  CHARACTER(LEN=2048)                            :: cf_isflist         ! input file names
  CHARACTER(LEN=2048)                            :: cf_bathy='bathy.nc'! bathymetry file name
  CHARACTER(LEN=2048)                            :: cf_isfdr='isf_draft.nc'! ice_draft file name
  CHARACTER(LEN=2048)                            :: cf_out='rnfisf.nc' ! output file for average
  !                                                VARIABLES
  CHARACTER(LEN=2048)                            :: cv_fill            ! fill var name
  CHARACTER(LEN=2048)                            :: cv_bathy='Bathymetry' ! bathymetry name
  CHARACTER(LEN=2048)                            :: cv_isfdr='isf_draft'  ! ice shelf draft name
  CHARACTER(LEN=2048)                            :: cdum               ! dummy string argument
  
  TYPE (variable), DIMENSION(:),    ALLOCATABLE :: stypvar            ! attributes for average values

  LOGICAL                                       :: lchk    = .false.  ! flag for missing files
  LOGICAL                                       :: lnc4    = .false.  ! flag for netcdf4 chunking and deflation
  !!----------------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()

  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfisf_rnf -f ISF-fill-file -v ISF-fill_var -l ISF-listfile -w width '
     PRINT *,'     [-b BATHY-file] [-vb BATHY-var] [-i ISFDRAFT-file] [-vi ISFDRAFT-variable]'
     PRINT *,'     [-nc4] [-o OUT-file ]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'        Build a netcdf file runoff file using the basal melting of the '
     PRINT *,'        ice-shelves. This netcdf file is intented to be used with NEMO when'
     PRINT *,'        nn_isf namelist parameter is set to 3.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'          -f ISF-fill_file : file built by cdffill (all the ice shelves are'
     PRINT *,'                             tagged with an id)'
     PRINT *,'          -v ISF-fill_var  : name of fill variable to use in ISF-fill_file'
     PRINT *,'          -l ISF-list : Text file with the melting rate (GT/y) given for'
     PRINT *,'               each ice shelf.' 
     PRINT *,'          -w width : specify the width (in grid points) on which the run-off'
     PRINT *,'               will be applied.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'          -b BATHY-file : give name of bathy file.'
     PRINT *,'                      [ default : ',TRIM(cf_bathy),' ]'
     PRINT *,'          -vp BATHY-var : give name of bathy variable.'
     PRINT *,'                      [ default : ',TRIM(cv_bathy),' ]'
     PRINT *,'          -i ISFDRAFT-file : give name of isf_draft file.'
     PRINT *,'                      [ default : ',TRIM(cf_isfdr),' ]'
     PRINT *,'          -vi ISFDRAFT-var : give name of isf_draft variable.'
     PRINT *,'                      [ default : ',TRIM(cv_isfdr),' ]'
     PRINT *,'          -nc4 : Use this option to have netcdf4 output file, with chunking'
     PRINT *,'               and deflation.'
     PRINT *,'          -o OUT-file : Specify the name of the output file instead of '
     PRINT *,'               the default name ',TRIM(cf_out)
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       mesh_hgr.nc and all files specified on the command line' 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) ,' unless -o option used'
     PRINT *,'         variables : sozisfmax (m), sozisfmin(m), sofwfisf (kg/m2/s)'
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'       cdfisf_fill, cdfisf_forcing, cdfisf_poolchk'
     PRINT *,'      '
     STOP
  ENDIF

  ijarg=1
  DO WHILE ( ijarg <= narg )
    CALL getarg(ijarg,cdum) ; ijarg=ijarg+1
    SELECT CASE (cdum)
    CASE ('-f' ) ; CALL getarg(ijarg,cf_fill   ) ; ijarg=ijarg+1
    CASE ('-v' ) ; CALL getarg(ijarg,cv_fill   ) ; ijarg=ijarg+1
    CASE ('-l' ) ; CALL getarg(ijarg,cf_isflist) ; ijarg=ijarg+1
    CASE ('-w' ) ; CALL getarg(ijarg,cdum      ) ; ijarg=ijarg+1
                   READ(cdum,*) nwidth
    CASE ('-b' ) ; CALL getarg(ijarg,cf_bathy  ) ; ijarg=ijarg+1
    CASE ('-vb') ; CALL getarg(ijarg,cv_bathy  ) ; ijarg=ijarg+1
    CASE ('-i' ) ; CALL getarg(ijarg,cf_isfdr  ) ; ijarg=ijarg+1
    CASE ('-vi') ; CALL getarg(ijarg,cv_isfdr  ) ; ijarg=ijarg+1

    CASE ('-nc4'); lnc4 = .TRUE.
    CASE ('-o' ) ; CALL getarg(ijarg,cf_out    ) ; ijarg=ijarg+1
    CASE DEFAULT 
        PRINT *,' Option ',TRIM(cdum),' not understood'
        STOP
    END SELECT
  ENDDO

  lchk = lchk .OR. chkfile (cf_fill   )
  lchk = lchk .OR. chkfile (cf_isflist)
  lchk = lchk .OR. chkfile (cf_bathy  )
  lchk = lchk .OR. chkfile (cf_isfdr  )
  lchk = lchk .OR. chkfile (cn_fhgr   )

  IF ( lchk ) STOP ! missing file

  npiglo = getdim (cf_fill, cn_x)
  npjglo = getdim (cf_fill, cn_y)
  npk    = 0                  ! bathy file has no dep dimension

  PRINT *, 'NPIGLO = ', npiglo
  PRINT *, 'NPJGLO = ', npjglo
  PRINT *, 'NPK    = ', npk

  ALLOCATE(isfindex (npiglo, npjglo), isfmask    (npiglo, npjglo))
  ALLOCATE(isfindex_wk(npiglo, npjglo), isum(npjglo))

  ALLOCATE(bathy    (npiglo, npjglo), zdrft      (npiglo, npjglo))
  ALLOCATE(zmax2d   (npiglo, npjglo), zmin2d     (npiglo, npjglo))

  ALLOCATE(de12t    (npiglo, npjglo) )
  ALLOCATE(dfwfisf2d(npiglo, npjglo), dl_fwfisf2d(npiglo, npjglo))

  ALLOCATE (stypvar(3))
  ALLOCATE (ipk(3),id_varout(3))

  CALL CreateOutput
  
  ! define variable
  ! read ice shelf draft data
  isfindex(:,:) = getvar(cf_fill,  cv_fill,  1 ,npiglo, npjglo )  ! fill index
  bathy(:,:)    = getvar(cf_bathy, cv_bathy, 1, npiglo, npjglo )  ! ocean bathy
  zdrft(:,:)    = getvar(cf_isfdr, cv_isfdr, 1, npiglo, npjglo )  ! ice shelf draft
  de12t(:,:)    = getvar(cn_fhgr,  cn_ve1t,  1, npiglo, npjglo ) * getvar(cn_fhgr, cn_ve2t, 1, npiglo, npjglo )

  ! open isf file
  OPEN(unit=iunit, file=cf_isflist, form='formatted', status='old')
  ! get number of isf
  nisf = 0
  cdum='XXX'
  DO WHILE ( TRIM(cdum) /= 'EOF')
     READ(iunit,*) cdum
     nisf=nisf+1
  END DO
  REWIND(iunit)
  nisf = nisf - 1

  PRINT *, '   Number of ISF found in file list : ', nisf

  DO jisf=1,nisf
     ! reset working isf index to its initial value
     isfindex_wk(:,:) = isfindex(:,:)

     ! read ice shelf data for jsf
     READ(iunit,*) ifill,cdum,rlon, rlat, iiseed, ijseed ,rdraftmin, rdraftmax, dfwf

     dl_fwf = dfwf * 1.d9 * 1.d3 / 86400.d0 / 365.d0  ! convert GT/yr to kg/m2/s
     isfmask    (:,:) = 0
     dl_fwfisf2d(:,:) = 0.0d0
     dsumcoef         = 0.0d0
     zmax2d(:,:)      = rdraftmax
     zmin2d(:,:)      = rdraftmin

     ! only deal with current ice shelf
     WHERE (isfindex_wk /=  -ifill ) isfindex_wk = 0

     ! find the j-limits for the current ice shelf 
     isum(:)=SUM(isfindex_wk,dim=1)
     ijmin = npjglo ; ijmax=2
     DO jj=ijseed, 2, -1
        IF ( isum(jj) /= 0 ) THEN
          ijmin=jj
        ENDIF
     ENDDO
     DO jj=ijseed, npjglo-1
        IF ( isum(jj) /= 0 ) THEN
          ijmax=jj
        ENDIF
     ENDDO

     DO jw = 1,nwidth
        DO ji=2,npiglo-1
           DO jj = ijmin, ijmax
              IF ( zdrft(ji,jj) == 0 .AND.  &    ! not under ice_shelf
              &    bathy(ji,jj) /= 0 .AND.  &    ! but in the ocean
              &    MINVAL(isfindex_wk(ji-1:ji+1 , jj-1:jj+1)) == -ifill  .AND. &  ! 
              &    isfindex_wk(ji,jj) == 0 ) THEN
                 ! compute dfwf in mm/s  ???
                 isfmask(ji,jj)  = isfmask(ji,jj) + jw 
                 ! use an empirical ocean ward fading function, in order to distribute the
                 ! runoff on nwidth points along the coast. iwscale can be adapted for
                 ! sharper (decrease) or flatter (increase) fading function. Default = 75
                 dl_fwfisf2d(ji,jj) = exp(-((isfmask(ji,jj)-1.)/iwscale)**2)  
              END IF
           END DO
        END DO
        WHERE (isfmask >= 1) isfindex_wk = -ifill ! capture already treated points under the shelf
     END DO

     dsumcoef=SUM(dl_fwfisf2d)
     PRINT *, SUM(dl_fwfisf2d), dsumcoef, SUM(dl_fwfisf2d / dsumcoef), dl_fwf, dfwf
     WHERE (isfmask >= 1)
       dl_fwfisf2d = dl_fwfisf2d / dsumcoef * dl_fwf / de12t
     END WHERE
     dfwfisf2d = dfwfisf2d + dl_fwfisf2d
  END DO
  ierr = putvar(ncout, id_varout(1), zmax2d   , 1, npiglo, npjglo)
  ierr = putvar(ncout, id_varout(2), zmin2d   , 1, npiglo, npjglo)
  ierr = putvar(ncout, id_varout(3), dfwfisf2d, 1, npiglo, npjglo)

  ierr = closeout(ncout)

  ! Diagnose total amount of fwf
  dl_fwf = 0.0d0
  DO ji=2,npiglo-1
     DO jj=2,npjglo-1
        dl_fwf = dl_fwf + dfwfisf2d(ji,jj) * de12t(ji,jj) * 86400.d0 * 365.d0
     END DO
  END DO
  PRINT *, 'total sum of isf = ', dl_fwf

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
  ! define new variables for output
  stypvar(1)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
  stypvar(1)%cname             = 'sozisfmax'
  stypvar(1)%rmissing_value    =  -99.
  stypvar(1)%valid_min         =  0.
  stypvar(1)%valid_max         =  2000.
  stypvar(1)%clong_name        = 'max depth of isf'
  stypvar(1)%cshort_name       = 'sozisfmax'
  stypvar(1)%conline_operation = 'N/A'
  stypvar(1)%caxis             = 'TYX'
  ipk(1) = 1  !  2D
  ! define new variables for output
  stypvar(2)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
  stypvar(2)%cname             = 'sozisfmin'
  stypvar(2)%rmissing_value    =  -99.
  stypvar(2)%valid_min         =  0.
  stypvar(2)%valid_max         =  2000.
  stypvar(2)%clong_name        = 'min depth of isf'
  stypvar(2)%cshort_name       = 'sozisfmin'
  stypvar(2)%conline_operation = 'N/A'
  stypvar(2)%caxis             = 'TYX'
  ipk(2) = 1  !  2D
  ! define new variables for output
  stypvar(3)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
  stypvar(3)%cname             = 'sofwfisf'
  stypvar(3)%rmissing_value    =  -99.d0
  stypvar(3)%valid_min         =  0.
  stypvar(3)%valid_max         =  2000.
  stypvar(3)%clong_name        = 'fwfisf'
  stypvar(3)%cshort_name       = 'sofwfisf'
  stypvar(3)%conline_operation = 'N/A'
  stypvar(3)%caxis             = 'TYX'
  stypvar(3)%cprecision        = 'r8'
  ipk(3) = 1  !  2D

  ! create output file taking the sizes in cf_fill
  ncout  = create      (cf_out, cf_fill,   npiglo, npjglo, npk,  ld_nc4=lnc4)
  ierr   = createvar   (ncout,  stypvar, 3,   ipk, id_varout  ,  ld_nc4=lnc4)
  ierr   = putheadervar(ncout,  cf_fill,   npiglo, npjglo, npk              )

  END SUBROUTINE CreateOutput

END PROGRAM cdfisf_rnf
