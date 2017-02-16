PROGRAM cdfisf_fill
  !!======================================================================
  !!                     ***  PROGRAM  cdfisf_fill  ***
  !!=====================================================================
  !!  ** Purpose : Build a file containing one value for each closed pools
  !!               seeded by a list of points.
  !!
  !!  ** Method  : flood filling algorithm
  !!               
  !! History : 3.0  : 04/2014  : Pierre Mathiot 
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  USE cdfio 
  USE modcdfnames
  USE modutils
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2014
  !! $Id$
  !! Copyright (c) 2010, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!-----------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                               :: jisf               ! dummy loop integer 
  INTEGER(KIND=4)                               :: ierr, ipos         ! working integer
  INTEGER(KIND=4)                               :: narg, iargc, ijarg ! browsing command line
  INTEGER(KIND=4)                               :: npiglo, npjglo     ! size of the domain
  INTEGER(KIND=4)                               :: npk, npt, nisf     ! size of the domain
  INTEGER(KIND=4)                               :: iunit=10           ! file unit for txt input file
  INTEGER(KIND=4)                               :: iunitu=11          ! file unit for txt output file
  INTEGER(KIND=4)                               :: ncout              ! ncid of output files
  INTEGER(KIND=4)                               :: iiseed, ijseed
  INTEGER(KIND=4)                               :: ifill
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: ipk                ! arrays of vertical level for each var
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: id_varout          ! varid's of average vars
  INTEGER(KIND=2), DIMENSION(:,:),  ALLOCATABLE :: itab

  REAL(KIND=4)                                  :: rlon, rlat         ! longitude and latitude of one point in ISF
  REAL(KIND=4)                                  :: rdraftmin, rdraftmax
  REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: dtab

  CHARACTER(LEN=2048)                            :: cf_in              ! input file name
  CHARACTER(LEN=2048)                            :: cf_isflist         ! input file name (txt)
  CHARACTER(LEN=2048)                            :: cf_isflistup       ! output file name (update of input, with draftmin/max
  CHARACTER(LEN=2048)                            :: cf_out='fill.nc'   ! output file for average
  CHARACTER(LEN=2048)                            :: cv_dep             ! depth dimension name
  CHARACTER(LEN=2048)                            :: cv_in              ! depth dimension name
  CHARACTER(LEN=2048)                            :: cdum               ! dummy string argument

  TYPE (variable), DIMENSION(:),    ALLOCATABLE :: stypvar            ! attributes for average values
  LOGICAL                                       :: lnc4 = .FALSE.     ! flag for netcdf4 chunk and deflation

  !!----------------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfisf_fill  -f ISF-file -v ISF-var -l ISF-list [-nc4 ] [-o OUT-file]'
     PRINT *,'      '
     PRINT *,'     PURPOSE : Build a nc file with a single value for each pool around a list'
     PRINT *,'               of given point. A warning is given when neighbouring ice-shelves'
     PRINT *,'               cannot be discriminated (no gap in between). In this case, hand'
     PRINT *,'               edit on the ISF-file is required.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS : '
     PRINT *,'         -f ISF-file : netcdf file  which contains the ice shelf draft variable'
     PRINT *,'                     (mesh_zgr is OK). It is used as a mask, only.'
     PRINT *,'         -v ISF-var  : variable name corresponding to the ice shelf draft or '
     PRINT *,'                      ice shelf level'
     PRINT *,'         -l ISF-list : text file containing at least the following information: '
     PRINT *,'                 1  NAME    LON  LAT I  J '
     PRINT *,'                 ...             '
     PRINT *,'                 i  NAMEi   LON  LAT I  J '
     PRINT *,'                 ...             '
     PRINT *,'                 EOF             '
     PRINT *,'                 No NAME  X    Y   I  J '
     PRINT *,'      '
     PRINT *,'     OPTIONS : '
     PRINT *,'          -nc4 : use NetCDF4 chunking and deflation for the output'
     PRINT *,'          -o OUT-file : specify the name of the output file instead of ',TRIM(cf_out)
     PRINT *,'                 This file will be one of the input file for cdfmkforcingisf '
     PRINT *,'                 as the ISF-fill_file '
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'              netcdf file : fill.nc '
     PRINT *,'              variable : sofillvar contains for all points in ice shelf NAME '
     PRINT *,'                         the value -i (negative value)'
     PRINT *,'              text file : <ISF-list>_zmin_zmax.txt '
     PRINT *,'                        this output file is similar to <ISF-list> but updated'
     PRINT *,'                        with the minimum and maximul value of ice-draft for '
     PRINT *,'                        each shelf.'
     PRINT *,'      '
     PRINT *,'     SEE ALSO : '
     PRINT *,'           cdfisf_forcing,  cdfisf_rnf '
     PRINT *,'      '
     STOP
  ENDIF

  ijarg = 1
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cdum ) ; ijarg = ijarg + 1 
     SELECT CASE ( cdum)
     CASE ( '-f' ) ; CALL getarg(ijarg, cf_in      ) ; ijarg = ijarg + 1
     CASE ( '-v' ) ; CALL getarg(ijarg, cv_in      ) ; ijarg = ijarg + 1
     CASE ( '-l' ) ; CALL getarg(ijarg, cf_isflist ) ; ijarg = ijarg + 1
     CASE ( '-o' ) ; CALL getarg(ijarg, cf_out     ) ; ijarg = ijarg + 1
     CASE ('-nc4') ; lnc4=.TRUE.
     CASE DEFAULT
        PRINT *, ' Option ', TRIM(cdum),' not understood'
        STOP
     END SELECT
  ENDDO

  IF ( chkfile (cf_in) .OR. chkfile (cf_isflist)  ) STOP ! missing file

  ipos = INDEX(cf_isflist,'.')
  cdum=cf_isflist(ipos+1:)
  cf_isflistup=cf_isflist(1:ipos-1)//'_zmin_zmax.'//TRIM(cdum)

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

  PRINT *, 'NPIGLO = ', npiglo
  PRINT *, 'NPJGLO = ', npjglo
  PRINT *, 'NPK    = ', npk

  ALLOCATE(dtab(npiglo, npjglo))
  ALLOCATE(itab(npiglo, npjglo))

  ALLOCATE (stypvar(1))
  ALLOCATE (ipk(1),id_varout(1))

  CALL CreateOutput 

  ! initialize variable
  dtab(:,:) = 0.d0 
  ! read ice shelf draft data
  dtab = getvar(cf_in, cv_in, 1 ,npiglo, npjglo )
  itab=1
  WHERE ( dtab <=0 ) itab=0
  PRINT *, 'Maximum of ISF-draft : ', MAXVAL(dtab),' m'

  ! open isf-list file
  OPEN(unit=iunit,  file=cf_isflist  , form='formatted', status='old')
  OPEN(unit=iunitu, file=cf_isflistup, form='formatted'              )
  ! get total number of isf
  nisf = 0
  cdum='XXX'
  DO WHILE ( TRIM(cdum) /= 'EOF')
     READ(iunit,*) cdum
     nisf=nisf+1
  END DO
  REWIND(iunit)

  nisf = nisf - 1
  PRINT *, '   Number of ISF found in file list : ', nisf

  ! loop over each ice shelf
  DO jisf=1,nisf
     ! get iiseed, ijseed, ice shelf number ifill
     READ(iunit,*) ifill, cdum, rlon, rlat, iiseed, ijseed
     IF (dtab(iiseed, ijseed) < 0 ) THEN
        PRINT *,'  ==> WARNING: Likely a problem with ',TRIM(cdum)
        PRINT *,'               check separation with neighbours'
     ENDIF
     CALL FillPool2D(iiseed, ijseed,itab, -ifill)

     rdraftmax=MAXVAL(dtab, (itab == -ifill) )
     rdraftmin=MINVAL(dtab, (itab == -ifill) )

     PRINT *,'Iceshelf : ', TRIM(cdum)
     PRINT *,'  index  : ', ifill
     PRINT *,'  code   : ', INT(dtab(iiseed, ijseed ) )
     PRINT *,'  depmin : ', rdraftmin
     PRINT *,'  depmax : ', rdraftmax
     PRINT *,'   '
     WRITE(iunitu,'(i4,1x,a20,2f9.4,2i5,2f8.1)') jisf,ADJUSTL(cdum),rlon, rlat, iiseed, ijseed,rdraftmin,rdraftmax
  END DO
  WRITE(iunitu,'(a)') 'EOF  '
  WRITE(iunitu,'(a5,a20,2a9,2a5,2a8,a)' ) 'No ','NAME                           ',' X',' Y',' I ',' J ',' Zmin',' Zmax',' FWF'

  CLOSE(iunitu)
  CLOSE(iunit)

  ! set to 0 all unwanted point (why not .GE. 0.0, I don't know)
  WHERE (dtab >= 1.d0)
     dtab = 0.0d0
  END WHERE

  ierr = putvar(ncout, id_varout(1), itab, 1, npiglo, npjglo)

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
    !! ** Method  :  Use global variables, defined in mail 
    !!----------------------------------------------------------------------
  ! define new variables for output
  stypvar(1)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
  stypvar(1)%cname             = 'sofillvar'
  stypvar(1)%cunits            = 'N/A'
  stypvar(1)%rmissing_value    = 0.
  stypvar(1)%valid_min         = -1000.
  stypvar(1)%valid_max         =  1000.
  stypvar(1)%clong_name        = 'Fill var'
  stypvar(1)%cshort_name       = 'sofillvar'
  stypvar(1)%conline_operation = 'N/A'
  stypvar(1)%caxis             = 'TYX'
  stypvar(1)%cprecision        = 'i2'
  ipk(1) = 1  !  2D

  ! create output file taking the sizes in cf_in
  ncout  = create      (cf_out,  cf_in,    npiglo, npjglo, npk, cdep=cv_dep, ld_nc4=lnc4)
  ierr   = createvar   (ncout ,  stypvar,  1,  ipk,    id_varout           , ld_nc4=lnc4)
  ierr   = putheadervar(ncout,   cf_in,    npiglo, npjglo, npk, cdep=cv_dep             )

  END SUBROUTINE CreateOutput

END PROGRAM cdfisf_fill
