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

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: glam, gphi     ! longitude, latitude

  CHARACTER(LEN=256)                        :: cv_lam         ! longitude name
  CHARACTER(LEN=256)                        :: cv_phi         ! latitude name
  CHARACTER(LEN=256)                        :: ctype='T'      ! type of point on C-grid
  CHARACTER(LEN=256)                        :: cldum          ! dummmy string
  CHARACTER(LEN=256)                        :: clcoo          ! dummy character variable
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()
  clcoo = cn_fcoo


  narg= iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : cdfwhereij  -w imin imax jmin jmax [-c COOR-file] [-p C-point]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Returns the geographical coordinates of a model sub-area specified' 
     PRINT *,'       by a rectangular window in (i,j) space.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -w imin imax jmin jmax : (i,j) space window coordinates.' 
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-c COOR_file ] : specify a coordinates file instead of ', TRIM(cn_fcoo)
     PRINT *,'       [-p C-point   ] : specify a point type on the C-grid (T U V F) '
     PRINT *,'               default is ', TRIM(ctype),'.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       ',TRIM(cn_fcoo),' or COOR-file given in the -c option'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       Standard output'
     STOP
  ENDIF

  ijarg = 1 
  DO WHILE ( ijarg <= narg) 
     CALL getarg( ijarg, cldum ) ; ijarg= ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-w' ) ; CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1 ;  READ(cldum,*) iimin
        ;            CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1 ;  READ(cldum,*) iimax
        ;            CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1 ;  READ(cldum,*) ijmin
        ;            CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1 ;  READ(cldum,*) ijmax
        ! options
     CASE ( '-c' ) ; CALL getarg(ijarg, clcoo ) ; ijarg=ijarg+1
     CASE ( '-p' ) ; CALL getarg(ijarg, ctype ) ; ijarg=ijarg+1
     CASE DEFAULT  ; PRINT *,' ERROR : ', TRIM(cldum),' : unknown option.' ; STOP 1
     END SELECT
  END DO

  IF ( chkfile(clcoo) ) STOP ! missing file

  npiglo = getdim (clcoo, cn_x)
  npjglo = getdim (clcoo, cn_y)

  IF ( iimax > npiglo ) THEN
     PRINT *,' ERROR : imax is greater than the maximum size ', iimax, npiglo
     STOP
  ENDIF

  IF ( ijmax > npjglo ) THEN
     PRINT *,' ERROR : jmax is greater than the maximum size ', ijmax, npjglo
     STOP
  END IF

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

  PRINT '(2a)'     ,' Type of point   : ', TRIM(ctype)
  PRINT '(a,4i6)'  ,'   I J zoom      : ', iimin, iimax, ijmin, ijmax
  PRINT '(a,4f9.3)','   LON LAT zoom  : ', glam(iimin,ijmin), glam(iimax,ijmax), gphi(iimin,ijmin), gphi(iimax,ijmax)

END PROGRAM cdfwhereij
