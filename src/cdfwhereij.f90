PROGRAM cdfwhereij
  !!======================================================================
  !!                     ***  PROGRAM  cdfwhereij  ***
  !!=====================================================================
  !!  ** Purpose : Give the values of longitude latitude for a given i, j
  !!
  !!  ** Method  : Read the coordinate/mesh_hgr file and look for the glam,
  !!               gphi variables. The point type ( T U V F ) is specified 
  !!               on the command line.
  !!
  !! History : 2.1  : 05/2005  : J.M. Molines : Original code
  !!           3.0  : 01/2011  : J.M. Molines : Doctor norm + Lic.
  !!         : 4.0  : 03/2017  : J.M. Molines  
  !!----------------------------------------------------------------------
  USE cdfio
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017 
  !! $Id$
  !! Copyright (c) 2017, J.-M. Molines 
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !! @class file_informations
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                           :: narg, iargc    ! browse line
  INTEGER(KIND=4)                           :: ijarg          ! browse line
  INTEGER(KIND=4)                           :: iimin, iimax   ! i-zoom limit
  INTEGER(KIND=4)                           :: ijmin, ijmax   ! j-zoom limit
  INTEGER(KIND=4)                           :: npiglo, npjglo ! global size
  INTEGER(KIND=4)                           :: npts, ipts     ! number of pts, loop index
  INTEGER(KIND=4)                           :: numpts, numout ! file id

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: glam, gphi     ! longitude, latitude

  CHARACTER(LEN=256)                        :: cv_lam         ! longitude name
  CHARACTER(LEN=256)                        :: cv_phi         ! latitude name
  CHARACTER(LEN=256)                        :: ctype='T'      ! type of point on C-grid
  CHARACTER(LEN=256)                        :: cldum          ! dummmy string
  CHARACTER(LEN=256)                        :: cfile, cfout   ! list of point
  CHARACTER(LEN=256)                        :: cname          ! name of the list of points
  CHARACTER(LEN=256)                        :: clcoo          ! dummy character variable

  LOGICAL :: llist = .FALSE. , lw = .FALSE. , lout = .FALSE. ! option flag 
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()
  clcoo = cn_fcoo


  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfwhereij  -w imin imax jmin jmax | -l LST-file]' 
     PRINT *,'                ...  [-c COOR-file] [-p C-type] [-o TXT-file]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Return the geographical coordinates of a model sub-area specified' 
     PRINT *,'       by a rectangular window in (i,j) space.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       One of the two following arguments is mandatory, according to your '
     PRINT *,'       wishes :'
     PRINT *,'       -w          : imin imax jmin jmax : (i,j) space window coordinates.' 
     PRINT *,'          or '
     PRINT *,'       -l LST-file : LST-file contains a list of points whse longitude and'
     PRINT *,'                     latitude will be looked up. It uses the same format than'
     PRINT *,'                     section files in other tools, with a 2-line header:'
     PRINT *,'                         NAME'
     PRINT *,'                         NPTS'
     PRINT *,'                         ii jj'
     PRINT *,'                         .....'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-c COOR_file ] : specify a coordinates file instead of ', TRIM(cn_fcoo)
     PRINT *,'       [-p C-type ]    : specify a point type on the C-grid (T U V F),  '
     PRINT *,'               default is ', TRIM(ctype),'.'
     PRINT *,'       [-o TXT-file ]  : output a text file that can be used directly in  '
     PRINT *,'                        cdf_xtract_broken tool. Only available if -l option '
     PRINT *,'                        activated'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       ',TRIM(cn_fcoo),' or COOR-file given in the -c option'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       Standard output'
     PRINT *,'      '
     PRINT *,'     SEE ALSO : '
     PRINT *,'      cdffindij, cdf_xtract_broken'
     PRINT *,'      '
     STOP 
  ENDIF

  ijarg = 1 
  DO WHILE ( ijarg <= narg) 
     CALL getarg( ijarg, cldum ) ; ijarg= ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-w' ) ; lw = .TRUE.
        ;            CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1 ;  READ(cldum,*) iimin
        ;            CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1 ;  READ(cldum,*) iimax
        ;            CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1 ;  READ(cldum,*) ijmin
        ;            CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1 ;  READ(cldum,*) ijmax
     CASE ( '-l' ) ; CALL getarg(ijarg, cfile ) ; ijarg=ijarg+1 ; llist = .TRUE. 
        ! options
     CASE ( '-c' ) ; CALL getarg(ijarg, clcoo ) ; ijarg=ijarg+1
     CASE ( '-p' ) ; CALL getarg(ijarg, ctype ) ; ijarg=ijarg+1
     CASE ( '-o' ) ; CALL getarg(ijarg, cfout ) ; ijarg=ijarg+1 ; lout  = .TRUE.
     CASE DEFAULT  ; PRINT *,' ERROR : ', TRIM(cldum),' : unknown option.' ; STOP 99
     END SELECT
  END DO

  IF ( chkfile(clcoo) ) STOP 99 ! missing file

  npiglo = getdim (clcoo, cn_x)
  npjglo = getdim (clcoo, cn_y)

  ALLOCATE (glam(npiglo,npjglo), gphi(npiglo,npjglo) )

  SELECT CASE ( ctype )
  CASE ('T' , 't' ) ; cv_lam = cn_glamt ; cv_phi = cn_gphit
  CASE ('U' , 'u' ) ; cv_lam = cn_glamu ; cv_phi = cn_gphiu
  CASE ('V' , 'v' ) ; cv_lam = cn_glamv ; cv_phi = cn_gphiv
  CASE ('F' , 'f' ) ; cv_lam = cn_glamf ; cv_phi = cn_gphif
  CASE DEFAULT
     PRINT *,' ERROR : type of point not known: ', TRIM(ctype)
  END SELECT

  glam(:,:) = getvar(clcoo, cv_lam, 1, npiglo, npjglo)
  gphi(:,:) = getvar(clcoo, cv_phi, 1, npiglo, npjglo)

  IF ( lw ) THEN
     ! check index
     IF ( chkindex(iimin, ijmin) ) STOP 99
     IF ( chkindex(iimax, ijmax) ) STOP 99
     ! print output
     PRINT '(2a)'     ,' Type of point   : ', TRIM(ctype)
     PRINT '(a,4i6)'  ,'   I J zoom      : ', iimin, iimax, ijmin, ijmax
     PRINT '(a,4f9.3)','   LON LAT zoom  : ', glam(iimin,ijmin), glam(iimax,ijmax), gphi(iimin,ijmin), gphi(iimax,ijmax)

  ELSE IF (llist) THEN
     ! open section file
     OPEN(numpts,FILE=cfile)
     READ(numpts,*) cname 
     READ(numpts,*) npts 
     ! print std output
     PRINT '(2a)'     ,' Type of point   : ', TRIM(ctype)
     PRINT '(2a)'     ,'   NAME  : ', cname
     PRINT '(a,1i6)'  ,'   NUMBER of points : ', npts
     PRINT '(a)'      ,'---------------------------'
     PRINT *,''
     ! print xtract_broken output
     IF (lout) THEN
        numout = 42
        OPEN(numout,FILE=cfout)
        WRITE(numout,'(a)')   TRIM(cname) 
        WRITE(numout,'(4i6)') npts 
     END IF

     DO ipts = 1,npts
        ! read point
        READ(numpts,*) iimin,ijmin
        IF ( chkindex(iimin, ijmin) ) STOP 99
        ! print std output
        PRINT '(a,1i6)'  ,'   POINT         : ', ipts
        PRINT '(a,2i6)'  ,'   I J zoom      : ', iimin, ijmin
        PRINT '(a,2f9.3)','   LON LAT zoom  : ', glam(iimin,ijmin), gphi(iimin,ijmin)
        PRINT *,''
        ! print xtract_broken output
        IF (lout) THEN
           WRITE(numout,'(2f9.3)') glam(iimin,ijmin), gphi(iimin,ijmin)
        END IF
     END DO
  ELSE
     PRINT  *, ' *** ERROR: you must use either -w or -l arguments.'
     STOP 99
  END IF

CONTAINS

  LOGICAL FUNCTION chkindex(kii,kjj)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION chkindex ***
    !!
    !! ** Purpose :  check if a the point specified by kii, kjj is within
    !!               the model boundaries. 
    !!
    !! ** Method  :  chkindex is set to TRUE if either kii or kjj is greater
    !!               than respectively npiglo or npjglo. In this case an error
    !!               message is printed. 
    !!               If the point is OK, chkindex is false, which allows a 
    !!               simple test like :
    !!               IF ( chkindex (ii,ij) ) STOP ! point outside the domain.
    !!
    !!----------------------------------------------------------------------
     INTEGER(KIND=4), INTENT( in ) :: kii, kjj 

     chkindex = .FALSE.
     IF ( kii > npiglo ) THEN
        PRINT *,' ERROR : imax is greater than the maximum size ', kii, npiglo
        chkindex = .TRUE.
     ENDIF

     IF ( kjj > npjglo ) THEN
        PRINT *,' ERROR : jmax is greater than the maximum size ', kjj, npjglo
        chkindex = .TRUE.
     END IF

  END FUNCTION chkindex

END PROGRAM cdfwhereij
