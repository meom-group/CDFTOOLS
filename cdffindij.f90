PROGRAM cdffindij
  !!======================================================================
  !!                     ***  PROGRAM  cdffindij  ***
  !!=====================================================================
  !!  ** Purpose : Return the window index (imin imax jmin jmax )
  !!               for the geographical windows given on input 
  !!               (longmin longmax latmin matmax)
  !!
  !!  ** Method  : Read the coordinate/mesh_hgr file and look for the glam,
  !!               gphi variables.
  !!               Then use a search algorithm to find the corresponding I J
  !!               The point type ( T U V F ) is specified on the command 
  !!               line as well as the name of the coordinate/mesh hgr file.
  !!
  !! History : 2.1  : 11/2005  : J.M. Molines : Original code
  !!           3.0  : 12/2010  : J.M. Molines : Doctor norm + Lic.
  !!----------------------------------------------------------------------
  USE cdfio
  USE cdftools
  USE modcdfnames
  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id$
  !! Copyright (c) 2010, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)    :: narg, iargc                ! command line
  INTEGER(KIND=4)    :: ijarg, ireq                ! command line
  INTEGER(KIND=4)    :: iimin, iimax, ijmin, ijmax ! model grid window

  REAL(KIND=4)       :: xmin, xmax, ymin, ymax     ! geographical window

  CHARACTER(LEN=256) :: cltype='F'                 ! point type to search for
  CHARACTER(LEN=256) :: cldum                      ! dummy character variable
  CHARACTER(LEN=256) :: clcoo                      ! dummy character variable
  !!----------------------------------------------------------------------
  CALL ReadCdfNames()
  clcoo = cn_fcoo

  !!  Read command line and output usage message if not compliant.
  narg= iargc()
  IF ( narg < 4 ) THEN
     PRINT *,' usage :   cdffindij  xmin xmax ymin ymax  [-c COOR-file] [-p point_type]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Return the model limit (i,j space) of the geographical window ' 
     PRINT *,'       given on the input line.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       xmin xmax ymin ymax : geographical limits of the window, in lon/lat' 
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       [-c COOR-file ] : specify a particular coordinate file' 
     PRINT *,'                     default is ',TRIM(cn_fcoo)
     PRINT *,'       [-p point type] : specify the point on the C-grid (T U V F)'
     PRINT *,'                     default is ',TRIM(cltype)
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       ', TRIM(cn_fcoo),' or the specified coordinates file.' 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       Output is done on standard output.'
     STOP
  ENDIF

  ijarg = 1 ; ireq = 0
  DO WHILE ( ijarg <= narg ) 
    CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
    SELECT CASE ( cldum )
    CASE ( '-c' ) ; CALL getarg(ijarg, clcoo  ) ; ijarg=ijarg+1
    CASE ( '-p' ) ; CALL getarg(ijarg, cltype ) ; ijarg=ijarg+1
    CASE DEFAULT
       ireq=ireq+1
       SELECT CASE (ireq)
       CASE ( 1 ) ; READ(cldum,*) xmin
       CASE ( 2 ) ; READ(cldum,*) xmax
       CASE ( 3 ) ; READ(cldum,*) ymin
       CASE ( 4 ) ; READ(cldum,*) ymax
       CASE DEFAULT 
         PRINT *,' Too many arguments !' ; STOP
       END SELECT
    END SELECT
  END DO

  CALL cdf_findij ( xmin, xmax, ymin, ymax, iimin, iimax, ijmin, ijmax, cd_coord=clcoo, cd_point=cltype)

END PROGRAM cdffindij
