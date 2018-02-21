PROGRAM cdficetrp
  !!======================================================================
  !!                     ***  PROGRAM  cdficetrp  ***
  !!=====================================================================
  !!  ** Purpose : compute ice transport across a section (either zonal
  !!               or meridional.
  !!
  !!  ** Method  : read horizontal metrics, ice velocities, ice thickness
  !!               and ice fraction. Then compute the following transport :
  !!               icertrp = sum( ice_frac * ice_thickness * e1/2 * ice_vel )
  !!
  !! History :  4.0  : 02/2018  : J.M. Molines : Original code
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2018
  !! $Id$
  !! Copyright (c) 2018, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER(KIND=4)                               :: narg, iargc, ijarg
  INTEGER(KIND=4)                               :: nsection             ! number of sections (overall)

  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: iimina, iimaxa       ! sections limits
  INTEGER(KIND=4), DIMENSION(:),    ALLOCATABLE :: ijmina, ijmaxa       ! sections limits



  CHARACTER(LEN=255) :: cf_ifil  ! input ice model
  CHARACTER(LEN=255) :: cf_sfil='ice_section.dat'  ! input section file (txt)
  CHARACTER(LEN=255) :: cldum    ! dummy character variable
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: csection             ! section name
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: cvarname             ! output variable name (root)
  CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: clongname            ! output long name (root)


  LOGICAL :: lchk = .FALSE.
  LOGICAL :: lim3 = .FALSE.

  !!----------------------------------------------------------------------
  CALL ReadCdfNames()
  narg = iargc()

  IF ( narg == 0 ) THEN
     PRINT *,' usage :  cdficetrp -i ICE-file [-s SECTION-file]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Compute ice transport through sections described into a section file,'
     PRINT *,'       whose default name is ',TRIM(cf_sfil),'. This file has the '
     PRINT *,'       folowing format:  It is a text file with pairs of lines foreach '
     PRINT *,'       section giving : (1) section name and (2) section location.'
     PRINT *,'     '
     PRINT *,'       First line with section name may also have 2 additional strings holding'
     PRINT *,'       a prefix for variable output, and a long name to be used as attribute in'
     PRINT *,'       the output file. '
     PRINT *,'       Second line gives the location of the section with specification of four'
     PRINT *,'       integer values (imin imax jmin jmax), relative to the model grid.'
     PRINT *,'       Only  zonal or meridional sections are allowed.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'        -i ICE-file : specify ICE-file containing ice velocity, and'
     PRINT *,'            ice concentration.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'        [-s SECTION-file] : give the name of the section file, instead'
     PRINT *,'            of default ',TRIM(cf_sfil)
     PRINT *,'        '
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        mesh_hgr.nc and mask.nc'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       Netcdf file : There is 1 netcdf file per section. File name is build'
     PRINT *,'         from section name : <SECTION>_icetrp.nc'
     PRINT *,'         variables :  icetrp '

!    PRINT *,'       netcdf file : ', TRIM(cf_out) 
!    PRINT *,'         variables : ', TRIM(cv_out),' (Sv)'
     PRINT *,'      '
     PRINT *,'     SEE ALSO : '
     PRINT *,'       cdftransport, cdficediags'
     PRINT *,'      '
     STOP
  ENDIF

  ! Parse  command line
  ijarg = 1 
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-i'   ) ; CALL getarg(ijarg, cf_ifil ) ; ijarg=ijarg+1
! Option
     CASE ( '-s'   ) ; CALL getarg(ijarg, cf_sfil ) ; ijarg=ijarg+1
     CASE DEFAULT    ; PRINT *, ' ERROR : ', TRIM(cldum),' : unknown option.'; STOP 1
     END SELECT
  ENDDO

! Check existence of files :
  lchk = chkfile (cn_fhgr) .OR. lchk
  lchk = chkfile (cn_fmsk) .OR. lchk
  lchk = chkfile (cf_ifil) .OR. lchk
  lchk = chkfile (cf_sfil) .OR. lchk
  IF ( lchk ) STOP 99 ! missing file

  ! check if cf_ifile corresponds to lim3 or lim2
  IF ( chkvar ( cf_ifil, cn_ileadfra ,.FALSE. ) ) THEN
     lim3 = .FALSE.
  ELSE
     lim3 = .TRUE.
     cn_ileadfra = cn_ileadfra3
     cn_iicethic = cn_iicethic3
     cn_iicevelu = cn_iicevelu3
     cn_iicevelv = cn_iicevelv3
  ENDIF

! Initialise section
  nsection=0
  CALL section_init(cf_sfil, csection,cvarname,clongname,iimina, iimaxa, ijmina, ijmaxa, nsection)

CONTAINS
  SUBROUTINE section_init(cdfile, cdsection, cdvarname, cdlongname, kimin, kimax, kjmin, kjmax, knumber)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE section_init  ***
    !!
    !! ** Purpose : Read input ASCII file that defines section names and limit of
    !!              sections.
    !!
    !! ** Method  : At fisrt call only return the number of sections for further
    !!              allocation.  
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*),                       INTENT(in   ) :: cdfile
    CHARACTER(LEN=256), DIMENSION(knumber), INTENT(out  ) :: cdsection
    CHARACTER(LEN=256), DIMENSION(knumber), INTENT(out  ) :: cdvarname
    CHARACTER(LEN=256), DIMENSION(knumber), INTENT(out  ) :: cdlongname
    INTEGER(KIND=4),                        INTENT(inout) :: knumber
    INTEGER(KIND=4), DIMENSION(knumber),    INTENT(out  ) :: kimin, kimax, kjmin, kjmax

    ! Local variables
    INTEGER(KIND=4)                                       :: jsec
    INTEGER(KIND=4)                                       :: ii, inum=10
    INTEGER(KIND=4)                                       :: ipos
    CHARACTER(LEN=256)                                    :: cline
    CHARACTER(LEN=80), DIMENSION(3)                       :: cldum
    LOGICAL                                               :: llfirst
    !!----------------------------------------------------------------------
    llfirst=.FALSE.
    IF ( knumber == 0 ) llfirst=.TRUE.

    OPEN(inum, FILE=cdfile)
    REWIND(inum)
    ii = 0

    ! read the file just to count the number of sections
    DO
       READ(inum,'(a)') cline
       IF (INDEX(cline,'EOF') == 0 ) THEN
          READ(inum,*)    ! skip one line
          ii = ii + 1
       ELSE
          EXIT
       ENDIF
    END DO

    knumber=ii
    IF ( llfirst ) RETURN

    REWIND(inum)
    DO jsec=1,knumber
       READ(inum,'(a)') cline
       ii = 0
       cldum(:) = 'none'
       ipos = INDEX(cline,' ')
       DO WHILE ( ipos > 1 )
          ii = ii + 1
          cldum(ii) = cline(1:ipos - 1 )
          cline = TRIM ( cline(ipos+1:) )
          ipos  = INDEX( cline,' ' )
          IF ( ii >= 3 ) EXIT
       END DO
       cdsection(jsec) = TRIM(cldum(1) )
       cdvarname(jsec) = TRIM(cldum(2) )
       cdlongname(jsec) = TRIM(cldum(3) )
       READ(inum,*    ) kimin(jsec), kimax(jsec), kjmin(jsec), kjmax(jsec)
    END DO

    CLOSE(inum)

  END SUBROUTINE section_init



END PROGRAM cdficetrp
