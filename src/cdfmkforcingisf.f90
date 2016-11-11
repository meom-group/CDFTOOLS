PROGRAM cdfmkforcingisf
  !!======================================================================
  !!                     ***  PROGRAM  cdfmkforcingisf  ***
  !!=====================================================================
  !!  ** Purpose : spread a specified total ice shelf melting over a specific
  !!               ice shelf melting pattern
  !!
  !!  ** Method  : read an ice shelf mask file produce by cdffill, read the
  !!               integrate melting for each ice shelf and the wanted pattern,
  !!               then compute the final melting for each ice shelf.
  !!
  !! History : 3.0  : 04/2014  : Pierre Mathiot 
  !!----------------------------------------------------------------------
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

  INTEGER(KIND=4)                               :: ierr               ! working integer
  INTEGER(KIND=4)                               :: narg, iargc, ijarg ! browsing command line
  INTEGER(KIND=4)                               :: npiglo, npjglo         ! size of the domain
  INTEGER(KIND=4)                               :: npk, npt, nisf, nwidth ! size of the domain
  INTEGER(KIND=4)                               :: jisf               ! loop counter
  INTEGER(KIND=4)                               :: iunit=10           ! id file
  INTEGER(KIND=4)                               :: ncout              ! ncid of output files
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: iiseed, ijseed 
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: ipk                ! arrays of vertical level for each var
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: id_varout          ! varid's of average vars
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: ifill
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: izmax, izmin

  REAL(KIND=8)                                  :: dl_fwf, dsumcoef
  REAL(KIND=8), DIMENSION(:),       ALLOCATABLE :: dfwf
  REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: dicedep, de1t, de2t, dmask2d, dfwfisf2d, disfmask, dl_fwfisf2d

  CHARACTER(LEN=256)                            :: cf_in              ! input file names
  CHARACTER(LEN=256)                            :: cf_isflist         ! input file names
  CHARACTER(LEN=256)                            :: cf_out='isfforcing.nc' ! output file for average
  CHARACTER(LEN=256)                            :: cv_dep             ! depth dimension name
  CHARACTER(LEN=256)                            :: cv_in              ! depth dimension name
  CHARACTER(LEN=256)                            :: cdum               ! dummy string argument

  TYPE (variable), DIMENSION(:),    ALLOCATABLE :: stypvar            ! attributes for average values

  !!----------------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfmkforcingisf ISF-maskfile  ISF-maskvariable ISF-listfile '
     PRINT *,'      '
     PRINT *,'     PURPOSE : '
     PRINT *,'         Build basal melting rate file used in NEMO ISF when nn_isf=4 '
     PRINT *,'      '
     PRINT *,'     ARGUMENTS : '
     PRINT *,'          ISF-maskfile  : mask build by cdffill (all the ice shelves are tagged with an id)'
     PRINT *,'          ISF-maskvariable : name of mask variable to use in ISF-maskfile'
     PRINT *,'          ISF-listfile : test file used to build isfmask.nc. Use last variable only (GT/y)'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'              '
     PRINT *,'     REQUIRED FILES : '
     PRINT *,'           mesh_zgr.nc mesh_hgr.nc,'
     PRINT *,'           isfforcing.nc (ie reference file used to define the isf melting pattern)'
     PRINT *,'      '
     PRINT *,'     OUTPUT :'
     PRINT *,'         netcdf file : isfforcing.nc'
     PRINT *,'         variable : sofwfisf '
     PRINT *,'      '
     PRINT *,'     SEE ALSO : cdffill and cdfmkrnfisf'
     PRINT *,'      '
     STOP
  ENDIF

  CALL getarg (1, cf_in) ; 
  CALL getarg (2, cv_in) ; 
  CALL getarg (3, cf_isflist) ; 
  IF ( chkfile (cf_in) .OR. chkfile (cf_isflist) ) STOP ! missing file

  npiglo = getdim (cf_in, cn_x)
  npjglo = getdim (cf_in, cn_y)
  npk    = getdim (cf_in, cn_z, cdtrue=cv_dep, kstatus=ierr)

  IF (ierr /= 0 ) THEN
     npk   = getdim (cf_in, 'z',cdtrue=cv_dep,kstatus=ierr)
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

  ALLOCATE(de1t(npiglo, npjglo), de2t(npiglo, npjglo))
  ALLOCATE(dmask2d(npiglo, npjglo), dfwfisf2d(npiglo, npjglo), dl_fwfisf2d(npiglo,npjglo))
  ALLOCATE(dicedep(npiglo, npjglo), disfmask(npiglo, npjglo) )

  ALLOCATE (stypvar(1))
  ALLOCATE (ipk(1),id_varout(1))
  ! initialisation of final fwf
  dfwfisf2d(:,:) = 0.0d0

  ! define new variables for output
  stypvar(1)%cname             = 'sofwfisf'
  stypvar(1)%cunits            = 'kg/s'
  stypvar(1)%rmissing_value    =  -99.
  stypvar(1)%valid_min         =  0.
  stypvar(1)%valid_max         =  2000.
  stypvar(1)%clong_name        = 'Ice Shelf Fresh Water Flux '
  stypvar(1)%cshort_name       = 'sofwfisf'
  stypvar(1)%conline_operation = 'N/A'
  stypvar(1)%caxis             = 'TYX'
  ipk(1) = 1  !  2D

  ! create output file taking the sizes in cf_in
  ncout  = create      (cf_out,  cf_in,    npiglo, npjglo, npk, cdep=cv_dep)
  ierr   = createvar   (ncout ,  stypvar,  1,  ipk,    id_varout           )
  ierr   = putheadervar(ncout,   cn_fzgr,  npiglo, npjglo, npk, cdep=cv_dep)

  ! define variable
  ! read ice shelf draft data
  dicedep  = getvar(cn_fzgr, 'misf', 1, npiglo, npjglo )
  de1t     = getvar(cn_fhgr, cn_ve1t,  1, npiglo, npjglo )
  de2t     = getvar(cn_fhgr, cn_ve2t,  1, npiglo, npjglo )

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

  PRINT *, nisf,' ISF'
  ! allocate variable
  ALLOCATE(iiseed(nisf), ijseed(nisf), ifill(nisf), izmax(nisf), izmin(nisf), dfwf(nisf))

  ! loop over all the ice shelf
  DO jisf=1,nisf
     ! read the ice shelf melting used a pattern
     dl_fwfisf2d = getvar('isfforcing.nc'   , 'sowflisf', 1 ,npiglo, npjglo )
     ! read the sossheig use to mask all the close pool in the model (ssh > 0.0
     ! WARNING it is not an universal test, have to find something better.
     dmask2d   = getvar('isfforcing.nc'   , 'sossheig', 1 ,npiglo, npjglo )
     WHERE (dmask2d >=  0.0d0) 
        dmask2d = 0.0d0
     END WHERE
     WHERE (dmask2d <  0.0d0) 
        dmask2d = 1.0d0
     END WHERE
     ! get ice shelf mask
     disfmask = getvar(cf_in , cv_in, 1 ,npiglo, npjglo )
     ! update isf mask with your input pattern mask (WARNING: issue if your
     ! pattern mask is smaller than your isf mask, you should drown your pattern
     ! or take care during the build of the pattern file)
     disfmask = disfmask * dmask2d
     ! read ice shelf data
     READ(iunit,'(i3,a4,4i5,f7.1)') ifill(jisf),cdum,iiseed(jisf),ijseed(jisf),izmin(jisf),izmax(jisf),dfwf(jisf)
     PRINT *, 'filling isf ',TRIM(cdum), ' in progress ... (',ifill(jisf),TRIM(cdum),iiseed(jisf),ijseed(jisf),izmin(jisf),izmax(jisf),dfwf(jisf),')'
     ! convertion of total ice shelf melting from Gt/y -> kg/s
     dl_fwf = dfwf(jisf) * 1.d9 * 1.d3 / 86400.d0 / 365.d0
     ! initialisation of variable
     dsumcoef = 0.0d0
     ! isolate the ice shelf data we want
     WHERE (disfmask .NE. -ifill(jisf))
        disfmask(:,:) = 0.d0
        dl_fwfisf2d(:,:) = 0.0d0
     END WHERE

     ! set the halo to 0 (to avoid double counting) 
     dl_fwfisf2d(1,:)=0.0d0 ; dl_fwfisf2d(npiglo,:)=0.0d0 ;

     dsumcoef = SUM(dl_fwfisf2d * de1t * de2t)
     ! similar that zoldfwf * e12t / sum(zoldfwf*e12t) * newfwf / e12t
     dl_fwfisf2d = dl_fwfisf2d(:,:) / dsumcoef * dl_fwf

     ! Value read from the text file has the wrong sign for melting.
     dfwfisf2d(:,:) = dfwfisf2d(:,:) - dl_fwfisf2d(:,:)
  END DO

  ! print isf forcing file
  ierr = putvar(ncout, id_varout(1), dfwfisf2d, 1, npiglo, npjglo)

  ! close file
  ierr = closeout(ncout)

END PROGRAM cdfmkforcingisf
