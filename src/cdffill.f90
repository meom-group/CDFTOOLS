PROGRAM cdffill
  !!======================================================================
  !!                     ***  PROGRAM  cdffill  ***
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
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2014
  !! $Id$
  !! Copyright (c) 2010, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!-----------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                               :: ierr               ! working integer
  INTEGER(KIND=4)                               :: narg, iargc, ijarg ! browsing command line
  INTEGER(KIND=4)                               :: npiglo, npjglo     ! size of the domain
  INTEGER(KIND=4)                               :: npk, npt, nisf     ! size of the domain
  INTEGER(KIND=4)                               :: jisf               ! loop integer 
  INTEGER(KIND=4)                               :: nid=10             ! id file
  INTEGER(KIND=4)                               :: ncout              ! ncid of output files
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: iseed, jseed 
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: ipk                ! arrays of vertical level for each var
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: id_varout          ! varid's of average vars
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: ifill

  REAL(KIND=4)                                  :: isddraftmin, isddraftmax
  REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: dtab

  CHARACTER(LEN=256)                            :: cf_in              ! input file names
  CHARACTER(LEN=256)                            :: cf_isflist         ! input file names
  CHARACTER(LEN=256)                            :: cf_out='fill.nc'   ! output file for average
  CHARACTER(LEN=256)                            :: cv_dep             ! depth dimension name
  CHARACTER(LEN=256)                            :: cv_in              ! depth dimension name
  CHARACTER(LEN=256)                            :: cdum               ! dummy string argument
  
  TYPE (variable), DIMENSION(:),    ALLOCATABLE :: stypvar            ! attributes for average values
  LOGICAL                                       :: lnc4 = .FALSE.     ! flag for netcdf4 chunk and deflation

  !!----------------------------------------------------------------------------
  CALL ReadCdfNames()

  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdffill  -f ISF-file -v ISF-var -l ISF-list [-nc4 ] [-o OUT-file]'
     PRINT *,'      '
     PRINT *,'     PURPOSE : build a nc file with a single value for each pool'
     PRINT *,'               around a list of given point'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS : '
     PRINT *,'         -f ISF-file : netcdf file  which contains the ice shelf draft variable'
     PRINT *,'                     (mesh_zgr is OK)'
     PRINT *,'         -v ISF-var  : variable name corresponding to the ice shelf draft or '
     PRINT *,'                      ice shelf level'
     PRINT *,'         -l ISF-list : text file containing at least the following information: '
     PRINT *,'                 1  NAME    I  J '
     PRINT *,'                 ...             '
     PRINT *,'                 i  NAMEi   I  J '
     PRINT *,'                 ...             '
     PRINT *,'                 EOF             '
     PRINT *,'                 No NAME  X    Y '
     PRINT *,'      '
     PRINT *,'     OPTIONS : '
     PRINT *,'          -nc4 : use NetCDF4 chunking and deflation for the output'
     PRINT *,'          -o OUT-file : specify the name of the output file instead of ',TRIM(cf_out)
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'              netcdf file : fill.nc '
     PRINT *,'              variable : sofillvar contains for all points in ice shelf NAMEi '
     PRINT *,'                         the value -i (negative value)'
     PRINT *,'      '
     PRINT *,'     SEE ALSO : '
     PRINT *,'           cdfmkforcingisf.f90 cdfmkrnfisf.f90 '
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
    CASE ('-nc4') ; lnc4=.true.
    CASE DEFAULT
       PRINT *, ' Option ', TRIM(cdum),' not understood'
       STOP
    END SELECT
  ENDDO

  IF ( chkfile (cf_in) .OR. chkfile (cf_isflist)  ) STOP ! missing file

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

  ALLOCATE (stypvar(1))
  ALLOCATE (ipk(1),id_varout(1))

  ! define new variables for output
  stypvar(1)%ichunk            = (/npiglo,MAX(1,npjglo/30),1,1 /)
  stypvar(1)%cname             = 'sofillvar'
  stypvar(1)%rmissing_value    = 0.
  stypvar(1)%valid_min         = -1000.
  stypvar(1)%valid_max         =  1000.
  stypvar(1)%clong_name        = 'Fill var'
  stypvar(1)%cshort_name       = 'sofillvar'
  stypvar(1)%conline_operation = 'N/A'
  stypvar(1)%caxis             = 'TYX'

  ipk(1) = 1  !  2D

  ! create output file taking the sizes in cf_in
  ncout  = create      (cf_out,  cf_in,    npiglo, npjglo, npk, cdep=cv_dep, ld_nc4=lnc4)
  ierr   = createvar   (ncout ,  stypvar,  1,  ipk,    id_varout           , ld_nc4=lnc4)
  ierr   = putheadervar(ncout,   cf_in,    npiglo, npjglo, npk, cdep=cv_dep             )

  ! initialize variable
  dtab(:,:) = 0.d0 
  ! read ice shelf draft data
  dtab = getvar(cf_in, cv_in, 1 ,npiglo, npjglo )
  PRINT *, 'Maximum of ISF-draft : ', MAXVAL(dtab),' m'

  ! open isf-list file
  OPEN(unit=nid, file=cf_isflist, form='formatted', status='old')
  ! get total number of isf
  nisf = 0
  cdum='XXX'
  DO WHILE ( TRIM(cdum) /= 'EOF')
     READ(nid,*) cdum
     nisf=nisf+1
  END DO
  REWIND(nid)

  nisf = nisf - 1
  PRINT *, '   Number of ISF found in file list : ', nisf

  ! allocate variable
  ALLOCATE(iseed(nisf), jseed(nisf), ifill(nisf))
  ! loop over each ice shelf
  DO jisf=1,nisf
     ! get iseed, jseed, ice shelf number ifill
     READ(nid,*) ifill(jisf), cdum, iseed(jisf), jseed(jisf)
!    READ(nid,'(i3,a4,2i5)') ifill(jisf), cdum, iseed(jisf), jseed(jisf)
     PRINT *, 'filling isf ',TRIM(cdum), ' in progress ... (',ifill(jisf), TRIM(cdum), iseed(jisf), jseed(jisf),')'
     CALL fillpool(iseed(jisf), jseed(jisf), dtab, -ifill(jisf), isddraftmax, isddraftmin)
     PRINT *,TRIM(cdum), ' depmax = ', isddraftmax
     PRINT *,TRIM(cdum), ' depmin = ', isddraftmin
  END DO

  ! set to 0 all unwanted point (why not .GE. 0.0, I don't know)
  WHERE (dtab >= 1.d0)
     dtab = 0.0d0
  END WHERE

  ! print sofillvar
  ierr = putvar(ncout, id_varout(1), REAL(dtab), 1, npiglo, npjglo)

  ! close file
  ierr = closeout(ncout)

CONTAINS

  SUBROUTINE fillpool(kiseed, kjseed, ddta, kifill, pdepmax, pdepmin)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE fillpool  ***
    !!
    !! ** Purpose :  Replace all area surrounding by mask value by mask value
    !!
    !! ** Method  :  flood fill algorithm
    !!
    !!----------------------------------------------------------------------
    INTEGER(KIND=4),                 INTENT(in)    :: kiseed, kjseed 
    INTEGER(KIND=4),                 INTENT(in)    :: kifill   ! new bathymetry
    REAL(KIND=8), DIMENSION(:,:),    INTENT(inout) :: ddta     ! new bathymetry
    REAL(KIND=4),                    INTENT(out)   :: pdepmin, pdepmax

    INTEGER :: ik                       ! number of point change
    INTEGER :: ip                       ! size of the pile
    INTEGER :: ji, jj                   ! loop index
    INTEGER :: iip1, iim1, ii, ij       ! working integer
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: ipile    ! pile variable

    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dl_data   ! new bathymetry
    !!----------------------------------------------------------------------
    PRINT *, 'WARNING North fold case not coded'
    ! allocate variable
    ! Why 240004, I don't remember. Will it be enough for eORCA12, I don't know ?
    ! to be sure, replace 240004 by npiglo*npjglo but it sould be much lower.
    ALLOCATE(ipile(240004,2))
    ALLOCATE(dl_data(npiglo,npjglo))

   ! initialise variables
   dl_data=ddta
   ipile(:,:)=0
   ipile(1,:)=[kiseed,kjseed]
   ip=1; ik=0
   pdepmax=0.0
   pdepmin=99999.9

   ! loop until the pile size is 0 or if the pool is larger than the critical size
   DO WHILE ( ip /= 0 .AND. ik < 600000);
      ik=ik+1
      ii=ipile(ip,1); ij=ipile(ip,2)

      ! update bathy and update pile size
      IF (ddta(ii,ij) <= pdepmin) pdepmin=ddta(ii,ij)
      IF (ddta(ii,ij) >= pdepmax) pdepmax=ddta(ii,ij)
      dl_data(ii,ij) =kifill
      ipile(ip,:)  =[0,0]; ip=ip-1

      ! check neighbour cells and update pile
      iip1=ii+1; IF ( iip1 == npiglo+1 ) iip1=2
      iim1=ii-1; IF ( iim1 == 0        ) iim1=npiglo-1
      IF (dl_data(ii, ij+1) > 1.0) THEN
          ip=ip+1; ipile(ip,:)=[ii  ,ij+1]
      END IF
      IF (dl_data(ii, ij-1) > 1.0) THEN
          ip=ip+1; ipile(ip,:)=[ii  ,ij-1]
      END IF
      IF (dl_data(iip1, ij) > 1.0) THEN
          ip=ip+1; ipile(ip,:)=[iip1,ij  ]
      END IF
      IF (dl_data(iim1, ij) > 1.0) THEN
          ip=ip+1; ipile(ip,:)=[iim1,ij  ]
      END IF
   END DO
   IF (ik < 600000) ddta=dl_data;

   DEALLOCATE(ipile); DEALLOCATE(dl_data)

  END SUBROUTINE fillpool

END PROGRAM cdffill
